/*
 * Copyright 2025 NWChemEx-Project
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "libxc.hpp"
#include "xc.h"

namespace scf::xc::libxc {

void load_modules(pluginplay::ModuleManager& mm) {
    mm.add_module<LibXCEnergy>("LibXC Energy");
    mm.add_module<LibXCPotential>("LibXC Potential");
}

void set_defaults(pluginplay::ModuleManager& mm) {
    mm.change_submod("LibXC Energy", "Density on a grid", "Density2Grid");
    mm.change_submod("LibXC Potential", "Density on a grid", "Density2Grid");
    mm.change_submod("LibXC Potential", "AOs on a grid", "AOs on a Grid");
}

std::pair<int, int> to_libxc_codes(chemist::qm_operator::xc_functional func) {
    using namespace chemist::qm_operator;
    switch(func) {
        case xc_functional::SVWN3: return std::pair(XC_LDA_X, XC_LDA_C_VWN_3);
        case xc_functional::SVWN5: return std::pair(XC_LDA_X, XC_LDA_C_VWN);
        default: throw std::runtime_error("Functional not supported");
    }
}

simde::type::molecule aos2molecule(const simde::type::ao_basis_set& aos) {
    simde::type::molecule mol;
    using atom_type = typename chemist::Molecule::atom_type;
    for(const auto& atomic_bs : aos) {
        const auto& center = atomic_bs.center();
        auto Z             = atomic_bs.atomic_number().value();
        atom_type a("X", Z, 1.0, center.x(), center.y(), center.z());
        mol.push_back(std::move(a));
    }
    return mol;
}

simde::type::tensor libxc_lda_energy_density(
  chemist::qm_operator::xc_functional func,
  const simde::type::tensor& rho_on_grid) {
    const auto [id_x, id_c] = to_libxc_codes(func);
    xc_func_type func_x, func_c;

    if(xc_func_init(&func_x, id_x, XC_UNPOLARIZED) != 0)
        throw std::runtime_error("Failed to initialize libxc functional ");
    if(xc_func_init(&func_c, id_c, XC_UNPOLARIZED) != 0)
        throw std::runtime_error("Failed to initialize libxc functional ");

    const auto& rho_buffer = rho_on_grid.buffer();
    auto& rv               = rho_buffer.allocator().runtime();
    tensorwrapper::allocator::Eigen<double> allocator(rv);
    const auto* prho = allocator.rebind(rho_buffer).get_immutable_data();

    auto n_grid = rho_on_grid.logical_layout().shape().as_smooth().extent(0);

    // Buffers for exchange and correlation energy densities per particle
    tensorwrapper::shape::Smooth shape{n_grid};
    tensorwrapper::layout::Physical layout(shape);

    auto pbuffer_x = allocator.allocate(layout);
    auto pbuffer_c = allocator.allocate(layout);

    auto* pepsilon_x = pbuffer_x->get_mutable_data();
    auto* pepsilon_c = pbuffer_c->get_mutable_data();
    xc_lda_exc(&func_c, n_grid, prho, pepsilon_x);
    xc_lda_exc(&func_x, n_grid, prho, pepsilon_c);

    simde::type::tensor epsilon_x(shape, std::move(pbuffer_x));
    simde::type::tensor epsilon_c(shape, std::move(pbuffer_c));
    simde::type::tensor epsilon_xc, exc_xc;

    epsilon_xc("i") = epsilon_x("i") + epsilon_c("i");
    exc_xc("i")     = epsilon_xc("i") * rho_on_grid("i");

    xc_func_end(&func_x);
    xc_func_end(&func_c);

    return exc_xc;
}

simde::type::tensor libxc_lda_energy_density_derivative(
  chemist::qm_operator::xc_functional func,
  const simde::type::tensor& rho_on_grid) {
    const auto [id_x, id_c] = to_libxc_codes(func);
    xc_func_type func_x, func_c;

    if(xc_func_init(&func_x, id_x, XC_UNPOLARIZED) != 0)
        throw std::runtime_error("Failed to initialize libxc functional ");
    if(xc_func_init(&func_c, id_c, XC_UNPOLARIZED) != 0)
        throw std::runtime_error("Failed to initialize libxc functional ");

    const auto& rho_buffer = rho_on_grid.buffer();
    auto& rv               = rho_buffer.allocator().runtime();
    tensorwrapper::allocator::Eigen<double> allocator(rv);
    const auto* prho = allocator.rebind(rho_buffer).get_immutable_data();

    auto n_grid = rho_on_grid.logical_layout().shape().as_smooth().extent(0);

    // Buffers for derivative of exchange and correlation energy densities
    tensorwrapper::shape::Smooth shape{n_grid};
    tensorwrapper::layout::Physical layout(shape);

    auto pbuffer_x = allocator.allocate(layout);
    auto pbuffer_c = allocator.allocate(layout);

    auto* pdepsilon_x = pbuffer_x->get_mutable_data();
    auto* pdepsilon_c = pbuffer_c->get_mutable_data();
    xc_lda_vxc(&func_c, n_grid, prho, pdepsilon_x);
    xc_lda_vxc(&func_x, n_grid, prho, pdepsilon_c);

    simde::type::tensor depsilon_x(shape, std::move(pbuffer_x));
    simde::type::tensor depsilon_c(shape, std::move(pbuffer_c));
    simde::type::tensor depsilon_xc, dexc_xc;

    depsilon_xc("i") = depsilon_x("i") + depsilon_c("i");
    dexc_xc("i")     = depsilon_xc("i");

    xc_func_end(&func_x);
    xc_func_end(&func_c);

    return dexc_xc;
}

simde::type::tensor tensorify_weights(const chemist::Grid& grid,
                                      parallelzone::runtime::RuntimeView rv) {
    auto n_grid = grid.size();
    tensorwrapper::allocator::Eigen<double> allocator(rv);
    tensorwrapper::shape::Smooth shape{n_grid};
    tensorwrapper::layout::Physical layout(shape);
    auto weight_buffer = allocator.allocate(layout);
    for(std::size_t i = 0; i < n_grid; ++i) {
        weight_buffer->set_elem({i}, grid.at(i).weight());
    }
    return simde::type::tensor(shape, std::move(weight_buffer));
}

namespace {

// A(m,i) = w(i) * B(m,i);
struct WeightMatrixKernel {
    using buffer_base = tensorwrapper::buffer::BufferBase;

    template<typename FloatType>
    auto run(const buffer_base& w, const buffer_base& b,
             parallelzone::runtime::RuntimeView& rv) {
        tensorwrapper::allocator::Eigen<FloatType> allocator(rv);

        const auto& eigen_w = allocator.rebind(w);
        const auto& eigen_b = allocator.rebind(b);

        const auto* pw = eigen_w.get_immutable_data();
        const auto* pb = eigen_b.get_immutable_data();

        const auto& shape_b = eigen_b.layout().shape().as_smooth();
        auto n_aos          = shape_b.extent(0);
        auto n_grid         = shape_b.extent(1);

        tensorwrapper::shape::Smooth rv_shape{n_aos, n_grid};
        tensorwrapper::layout::Physical rv_layout(rv_shape);
        auto rv_buffer = allocator.allocate(rv_layout);

        // AOs on rows, grid points on columns
        for(std::size_t grid_i = 0; grid_i < n_grid; ++grid_i) {
            for(std::size_t ao_i = 0; ao_i < n_aos; ++ao_i) {
                auto value = pw[grid_i] * pb[ao_i * n_grid + grid_i];
                rv_buffer->set_elem(std::vector{ao_i, grid_i}, value);
            }
        }
        return simde::type::tensor(rv_shape, std::move(rv_buffer));
    }
};

// Does A(m,i) = 1/sqrt(w(m)) * B(m,i);
struct NormalizeKernel {
    using buffer_base = tensorwrapper::buffer::BufferBase;

    template<typename FloatType>
    auto run(const buffer_base& w, const buffer_base& b,
             parallelzone::runtime::RuntimeView& rv) {
        tensorwrapper::allocator::Eigen<FloatType> allocator(rv);

        const auto& eigen_w = allocator.rebind(w);
        const auto& eigen_b = allocator.rebind(b);

        const auto* pw = eigen_w.get_immutable_data();
        const auto* pb = eigen_b.get_immutable_data();

        const auto& shape_b = eigen_b.layout().shape().as_smooth();
        auto n_aos          = shape_b.extent(0);
        auto n_grid         = shape_b.extent(1);

        tensorwrapper::shape::Smooth rv_shape{n_aos, n_grid};
        tensorwrapper::layout::Physical rv_layout(rv_shape);
        auto rv_buffer = allocator.allocate(rv_layout);

        // AOs on rows, grid points on columns
        for(std::size_t grid_i = 0; grid_i < n_grid; ++grid_i) {
            for(std::size_t ao_i = 0; ao_i < n_aos; ++ao_i) {
                if constexpr(std::is_same_v<FloatType, double>) {
                    auto new_weight = 1.0 / std::sqrt(pw[ao_i]);
                    auto value      = new_weight * pb[ao_i * n_grid + grid_i];
                    rv_buffer->set_elem(std::vector{ao_i, grid_i}, value);
                } else {
                    throw std::runtime_error(
                      "Normalization not implemented for uncertain floats");
                }
            }
        }
        return simde::type::tensor(rv_shape, std::move(rv_buffer));
    }
};

// Does rho(i) = X(n,i) * aos(n,i) or rh(n) = X(n,i) * aos(n,i)
struct BatchedDotKernel {
    using buffer_base = tensorwrapper::buffer::BufferBase;

    bool m_sum_row = true;

    BatchedDotKernel() = default;
    explicit BatchedDotKernel(bool sum_row) : m_sum_row(sum_row) {}

    template<typename FloatType>
    auto run(const buffer_base& aos_on_grid, const buffer_base& X,
             parallelzone::runtime::RuntimeView& rv) {
        tensorwrapper::allocator::Eigen<FloatType> allocator(rv);

        const auto& eigen_aos_on_grid = allocator.rebind(aos_on_grid);
        const auto& eigen_X           = allocator.rebind(X);

        const auto* paos_on_grid = eigen_aos_on_grid.get_immutable_data();
        const auto* pX           = eigen_X.get_immutable_data();

        const auto& shape_X = eigen_X.layout().shape().as_smooth();
        auto n_aos          = shape_X.extent(0);
        auto n_grid         = shape_X.extent(1);

        auto out_size = m_sum_row ? n_grid : n_aos;
        tensorwrapper::shape::Smooth rv_shape{out_size};
        tensorwrapper::layout::Physical rv_layout(rv_shape);
        auto rv_buffer = allocator.allocate(rv_layout);

        // AOs on rows, grid points on columns
        if(m_sum_row) {
            for(std::size_t grid_i = 0; grid_i < n_grid; ++grid_i) {
                FloatType sum = 0.0;
                for(std::size_t ao_i = 0; ao_i < n_aos; ++ao_i) {
                    const auto idx = ao_i * n_grid + grid_i;
                    sum += paos_on_grid[idx] * pX[idx];
                }
                rv_buffer->set_elem(std::vector{grid_i}, sum);
            }
        } else {
            for(std::size_t ao_i = 0; ao_i < n_aos; ++ao_i) {
                FloatType sum = 0.0;
                for(std::size_t grid_i = 0; grid_i < n_grid; ++grid_i) {
                    const auto idx = ao_i * n_grid + grid_i;
                    sum += paos_on_grid[idx] * pX[idx];
                }
                rv_buffer->set_elem(std::vector{ao_i}, sum);
            }
        }
        return simde::type::tensor(rv_shape, std::move(rv_buffer));
    }
};

} // namespace

simde::type::tensor weight_a_matrix(const simde::type::tensor& w,
                                    const simde::type::tensor& b) {
    using tensorwrapper::utilities::floating_point_dispatch;
    WeightMatrixKernel k;
    auto runtime = b.buffer().allocator().runtime();
    return floating_point_dispatch(k, w.buffer(), b.buffer(), runtime);
}

simde::type::tensor normalize_row(const simde::type::tensor& w,
                                  const simde::type::tensor& b) {
    using tensorwrapper::utilities::floating_point_dispatch;
    NormalizeKernel k;
    auto runtime = b.buffer().allocator().runtime();
    return floating_point_dispatch(k, w.buffer(), b.buffer(), runtime);
}

simde::type::tensor batched_dot(const simde::type::tensor& aos_on_grid,
                                const simde::type::tensor& X, bool sum_row) {
    using tensorwrapper::utilities::floating_point_dispatch;
    BatchedDotKernel k(sum_row);
    auto rv = X.buffer().allocator().runtime();
    return floating_point_dispatch(std::move(k), aos_on_grid.buffer(),
                                   X.buffer(), rv);
}
} // namespace scf::xc::libxc
