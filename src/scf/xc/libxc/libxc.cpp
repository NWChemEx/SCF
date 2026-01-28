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
#ifdef BUILD_LIBXC
#include "xc.h"
#include "xc_funcs.h"
#endif

namespace scf::xc::libxc {
namespace {
template<typename FloatType>
auto* get_pointer(tensorwrapper::buffer::BufferBase& buffer) {
    using tensorwrapper::buffer::make_contiguous;
    auto& contig_buffer = make_contiguous(buffer);
    auto wtf_buffer     = contig_buffer.get_mutable_data();
    return wtf::buffer::contiguous_buffer_cast<FloatType>(wtf_buffer).data();
}

template<typename FloatType>
const auto* get_const_pointer(const tensorwrapper::buffer::BufferBase& buffer) {
    using tensorwrapper::buffer::make_contiguous;
    const auto& contig_buffer = make_contiguous(buffer);
    auto wtf_buffer           = contig_buffer.get_immutable_data();
    return wtf::buffer::contiguous_buffer_cast<const FloatType>(wtf_buffer)
      .data();
}
} // namespace

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
#ifdef BUILD_LIBXC
    switch(func) {
        case xc_functional::SVWN3: return std::pair(XC_LDA_X, XC_LDA_C_VWN_3);
        case xc_functional::SVWN5: return std::pair(XC_LDA_X, XC_LDA_C_VWN);
        default: throw std::runtime_error("Functional not supported");
    }
#else
    throw std::runtime_error("LibXC support not compiled in");
#endif
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
#ifdef BUILD_LIBXC
    const auto [id_x, id_c] = to_libxc_codes(func);
    xc_func_type func_x, func_c;

    if(xc_func_init(&func_x, id_x, XC_UNPOLARIZED) != 0)
        throw std::runtime_error("Failed to initialize libxc functional ");
    if(xc_func_init(&func_c, id_c, XC_UNPOLARIZED) != 0)
        throw std::runtime_error("Failed to initialize libxc functional ");

    const auto& rho_buffer = rho_on_grid.buffer();
    const auto* prho       = get_const_pointer<double>(rho_buffer);

    auto n_grid = rho_on_grid.logical_layout().shape().as_smooth().extent(0);

    // Buffers for exchange and correlation energy densities per particle
    tensorwrapper::shape::Smooth shape{n_grid};

    auto pbuffer_x = tensorwrapper::buffer::make_contiguous(rho_buffer, shape);
    auto pbuffer_c = tensorwrapper::buffer::make_contiguous(rho_buffer, shape);

    auto* pepsilon_x = get_pointer<double>(pbuffer_x);
    auto* pepsilon_c = get_pointer<double>(pbuffer_c);
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
#else
    throw std::runtime_error("LibXC support not compiled in");
#endif
}

simde::type::tensor libxc_lda_energy_density_derivative(
  chemist::qm_operator::xc_functional func,
  const simde::type::tensor& rho_on_grid) {
#ifdef BUILD_LIBXC
    const auto [id_x, id_c] = to_libxc_codes(func);
    xc_func_type func_x, func_c;

    if(xc_func_init(&func_x, id_x, XC_UNPOLARIZED) != 0)
        throw std::runtime_error("Failed to initialize libxc functional ");
    if(xc_func_init(&func_c, id_c, XC_UNPOLARIZED) != 0)
        throw std::runtime_error("Failed to initialize libxc functional ");

    const auto& rho_buffer = make_contiguous(rho_on_grid.buffer());
    const auto* prho       = get_const_pointer<double>(rho_buffer);

    auto n_grid = rho_on_grid.logical_layout().shape().as_smooth().extent(0);

    // Buffers for derivative of exchange and correlation energy densities
    tensorwrapper::shape::Smooth shape{n_grid};

    auto pbuffer_x = make_contiguous(rho_buffer, shape);
    auto pbuffer_c = make_contiguous(rho_buffer, shape);

    auto* pdepsilon_x = get_pointer<double>(pbuffer_x);
    auto* pdepsilon_c = get_pointer<double>(pbuffer_c);
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
#else
    throw std::runtime_error("LibXC support not compiled in");
#endif
}

simde::type::tensor tensorify_weights(const chemist::Grid& grid,
                                      parallelzone::runtime::RuntimeView) {
    auto n_grid = grid.size();
    tensorwrapper::shape::Smooth shape{n_grid};
    std::vector<double> weights(n_grid);
    for(std::size_t i = 0; i < n_grid; ++i) {
        weights[i] = grid.at(i).weight();
    }
    using namespace tensorwrapper::buffer;
    auto weight_buffer = Contiguous(std::move(weights), shape);
    return simde::type::tensor(shape, std::move(weight_buffer));
}

namespace {

// A(m,i) = w(i) * B(m,i);
struct WeightMatrixKernel {
    using buffer_base = tensorwrapper::buffer::BufferBase;
    using tensor_type = simde::type::tensor;

    WeightMatrixKernel(std::size_t n_aos, std::size_t n_grid) :
      m_n_aos(n_aos), m_n_grid(n_grid) {}

    template<typename FloatType0, typename FloatType1>
    tensor_type operator()(const std::span<FloatType0>& w,
                           const std::span<FloatType1>& b) {
        throw std::runtime_error(
          "WeightMatrixKernel: Mixed float types not supported");
    }

    template<typename FloatType>
    tensor_type operator()(const std::span<FloatType>& w,
                           const std::span<FloatType>& b) {
        tensorwrapper::shape::Smooth rv_shape{m_n_aos, m_n_grid};
        using clean_type = std::remove_cv_t<FloatType>;
        std::vector<clean_type> rv_buffer(m_n_aos * m_n_grid);

        // AOs on rows, grid points on columns
        for(std::size_t grid_i = 0; grid_i < m_n_grid; ++grid_i) {
            for(std::size_t ao_i = 0; ao_i < m_n_aos; ++ao_i) {
                auto value = w[grid_i] * b[ao_i * m_n_grid + grid_i];
                rv_buffer[ao_i * m_n_grid + grid_i] = value;
            }
        }
        auto contig_buffer =
          tensorwrapper::buffer::Contiguous(std::move(rv_buffer), rv_shape);
        return simde::type::tensor(rv_shape, std::move(contig_buffer));
    }

    std::size_t m_n_aos;
    std::size_t m_n_grid;
};

// Does A(m,i) = 1/sqrt(w(m)) * B(m,i);
struct NormalizeKernel {
    using buffer_base = tensorwrapper::buffer::BufferBase;

    using tensor_type = simde::type::tensor;

    NormalizeKernel(std::size_t n_aos, std::size_t n_grid) :
      m_n_aos(n_aos), m_n_grid(n_grid) {}

    template<typename FloatType0, typename FloatType1>
    tensor_type operator()(const std::span<FloatType0>& w,
                           const std::span<FloatType1>& b) {
        throw std::runtime_error(
          "NormalizeKernel: Mixed float types not supported");
    }

    template<typename FloatType>
    tensor_type operator()(const std::span<FloatType>& w,
                           const std::span<FloatType>& b) {
        tensorwrapper::shape::Smooth rv_shape{m_n_aos, m_n_grid};
        using clean_type = std::remove_cv_t<FloatType>;
        std::vector<clean_type> rv_buffer(m_n_aos * m_n_grid);

        // AOs on rows, grid points on columns
        for(std::size_t grid_i = 0; grid_i < m_n_grid; ++grid_i) {
            for(std::size_t ao_i = 0; ao_i < m_n_aos; ++ao_i) {
                if constexpr(std::is_same_v<FloatType, double>) {
                    auto new_weight = 1.0 / std::sqrt(w[ao_i]);
                    auto value      = new_weight * b[ao_i * m_n_grid + grid_i];
                    rv_buffer[ao_i * m_n_grid + grid_i] = value;
                } else {
                    throw std::runtime_error(
                      "Normalization not implemented for uncertain floats");
                }
            }
        }
        auto cbuffer =
          tensorwrapper::buffer::Contiguous(std::move(rv_buffer), rv_shape);
        return simde::type::tensor(rv_shape, std::move(cbuffer));
    }

    std::size_t m_n_aos;
    std::size_t m_n_grid;
};

// Does rho(i) = X(n,i) * aos(n,i) or rh(n) = X(n,i) * aos(n,i)
struct BatchedDotKernel {
    using buffer_base = tensorwrapper::buffer::BufferBase;
    using tensor_type = simde::type::tensor;

    bool m_sum_row = true;
    std::size_t m_n_grid;
    std::size_t m_n_aos;

    BatchedDotKernel(std::size_t n_aos, std::size_t n_grid,
                     bool sum_row = true) :
      m_sum_row(true), m_n_grid(n_grid), m_n_aos(n_aos) {}

    template<typename FloatType0, typename FloatType1>
    tensor_type operator()(const std::span<FloatType0>& w,
                           const std::span<FloatType1>& b) {
        throw std::runtime_error(
          "BatchedDotKernel: Mixed float types not supported");
    }

    template<typename FloatType>
    tensor_type operator()(const std::span<FloatType>& aos_on_grid,
                           const std::span<FloatType>& X) {
        auto out_size = m_sum_row ? m_n_grid : m_n_aos;
        tensorwrapper::shape::Smooth rv_shape{out_size};
        using clean_type = std::remove_cv_t<FloatType>;
        std::vector<clean_type> rv_buffer(out_size);

        // AOs on rows, grid points on columns
        if(m_sum_row) {
            for(std::size_t grid_i = 0; grid_i < m_n_grid; ++grid_i) {
                clean_type sum = 0.0;
                for(std::size_t ao_i = 0; ao_i < m_n_aos; ++ao_i) {
                    const auto idx = ao_i * m_n_grid + grid_i;
                    sum += aos_on_grid[idx] * X[idx];
                }
                rv_buffer[grid_i] = sum;
            }
        } else {
            for(std::size_t ao_i = 0; ao_i < m_n_aos; ++ao_i) {
                clean_type sum = 0.0;
                for(std::size_t grid_i = 0; grid_i < m_n_grid; ++grid_i) {
                    const auto idx = ao_i * m_n_grid + grid_i;
                    sum += aos_on_grid[idx] * X[idx];
                }
                rv_buffer[ao_i] = sum;
            }
        }
        tensorwrapper::buffer::Contiguous cbuffer(std::move(rv_buffer),
                                                  rv_shape);
        return simde::type::tensor(rv_shape, std::move(cbuffer));
    }
};

} // namespace

simde::type::tensor weight_a_matrix(const simde::type::tensor& w,
                                    const simde::type::tensor& b) {
    using tensorwrapper::buffer::make_contiguous;
    using tensorwrapper::buffer::visit_contiguous_buffer;
    auto shape  = b.logical_layout().shape().as_smooth();
    auto n_aos  = shape.extent(0);
    auto n_grid = shape.extent(1);
    WeightMatrixKernel k(n_aos, n_grid);
    const auto& b_buffer = make_contiguous(b.buffer());
    const auto& w_buffer = make_contiguous(w.buffer());
    return visit_contiguous_buffer(k, w_buffer, b_buffer);
}

simde::type::tensor normalize_row(const simde::type::tensor& w,
                                  const simde::type::tensor& b) {
    using tensorwrapper::buffer::make_contiguous;
    using tensorwrapper::buffer::visit_contiguous_buffer;
    auto shape  = b.logical_layout().shape().as_smooth();
    auto n_aos  = shape.extent(0);
    auto n_grid = shape.extent(1);
    NormalizeKernel k(n_aos, n_grid);
    const auto& b_buffer = make_contiguous(b.buffer());
    const auto& w_buffer = make_contiguous(w.buffer());
    return visit_contiguous_buffer(k, w_buffer, b_buffer);
}

simde::type::tensor batched_dot(const simde::type::tensor& aos_on_grid,
                                const simde::type::tensor& X, bool sum_row) {
    using tensorwrapper::buffer::make_contiguous;
    using tensorwrapper::buffer::visit_contiguous_buffer;
    auto shape  = X.logical_layout().shape().as_smooth();
    auto n_aos  = shape.extent(0);
    auto n_grid = shape.extent(1);
    BatchedDotKernel k(n_aos, n_grid, sum_row);
    const auto& aos_on_grid_buffer = make_contiguous(aos_on_grid.buffer());
    const auto& X_buffer           = make_contiguous(X.buffer());
    return visit_contiguous_buffer(std::move(k), aos_on_grid_buffer, X_buffer);
}
} // namespace scf::xc::libxc
