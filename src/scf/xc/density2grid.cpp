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

#include "xc.hpp"
#include <simde/integration_grids/collocation_matrix.hpp>

namespace scf::xc {
namespace {
const auto desc = R"(

DensityCollocationMatrix
-----------------
)";

struct Kernel {
    using buffer_base = tensorwrapper::buffer::BufferBase;

    template<typename FloatType>
    auto run(const buffer_base& aos_on_grid, const buffer_base& X,
             parallelzone::runtime::RuntimeView& rv) {
        tensorwrapper::allocator::Eigen<FloatType> allocator(rv);

        const auto& eigen_aos_on_grid = allocator.rebind(aos_on_grid);
        const auto* paos_on_grid      = eigen_aos_on_grid.get_immutable_data();
        const auto& eigen_X           = allocator.rebind(X);
        const auto* pX                = eigen_X.get_immutable_data();
        const auto& shape_X           = eigen_X.layout().shape().as_smooth();
        auto n_aos                    = shape_X.extent(0);
        auto n_grid                   = shape_X.extent(1);

        tensorwrapper::shape::Smooth rv_shape{n_grid};
        tensorwrapper::layout::Physical rv_layout(rv_shape);
        auto rv_buffer = allocator.allocate(rv_layout);

        // AOs on rows, grid points on columns
        for(std::size_t grid_i = 0; grid_i < n_grid; ++grid_i) {
            FloatType sum = 0;
            for(std::size_t ao_i = 0; ao_i < n_aos; ++ao_i) {
                const auto idx = ao_i * n_grid + grid_i;
                sum += paos_on_grid[idx] * pX[idx];
            }
            rv_buffer->set_elem(std::vector{grid_i}, sum);
        }
        return simde::type::tensor(rv_shape, std::move(rv_buffer));
    }
};
} // namespace

using pt         = simde::EDensityCollocationMatrix;
using ao2grid_pt = simde::AOCollocationMatrix;

MODULE_CTOR(Density2Grid) {
    satisfies_property_type<pt>();
    description(desc);

    add_submodule<ao2grid_pt>("AOs on a grid");
}

MODULE_RUN(Density2Grid) {
    const auto& [grid, density] = pt::unwrap_inputs(inputs);

    const auto& rho = density.value();
    const auto& aos = density.basis_set().ao_basis_set();

    auto& ao2grid_mod = submods.at("AOs on a grid");
    auto aos_on_grid  = ao2grid_mod.run_as<ao2grid_pt>(grid, aos);

    simde::type::tensor X;
    X("m,i") = rho("m,n") * aos_on_grid("n,i");

    using tensorwrapper::utilities::floating_point_dispatch;
    Kernel k;
    auto& runtime          = get_runtime();
    const auto& aos_buffer = aos_on_grid.buffer();
    const auto& X_buffer   = X.buffer();
    auto rho_on_grid =
      floating_point_dispatch(k, aos_buffer, X_buffer, runtime);

    auto rv = results();
    return pt::wrap_results(rv, std::move(rho_on_grid));
}

} // namespace scf::xc
