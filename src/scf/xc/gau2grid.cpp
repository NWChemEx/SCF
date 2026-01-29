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
#include <gau2grid/gau2grid.h>
#include <simde/integration_grids/collocation_matrix.hpp>

namespace scf::xc {
namespace {
const auto desc = R"(

CollocationMatrix
-----------------

Warning! Gau2Grid assumes that the primitive normalization has been folded into
the contraction coefficients!!!
)";

template<typename T>
std::vector<T> flatten_grid(const chemist::Grid& grid) {
    std::vector<T> flattened_grid;
    flattened_grid.reserve(grid.size() * 3);
    for(const auto& point : grid) {
        flattened_grid.push_back(static_cast<T>(point.point().x()));
        flattened_grid.push_back(static_cast<T>(point.point().y()));
        flattened_grid.push_back(static_cast<T>(point.point().z()));
    }
    return flattened_grid;
}

} // namespace

using pt = simde::AOCollocationMatrix;

MODULE_CTOR(Gau2Grid) {
    satisfies_property_type<pt>();

    description(desc);
}

MODULE_RUN(Gau2Grid) {
    const auto& [grid, ao_basis] = pt::unwrap_inputs(inputs);
    auto n_points                = grid.size();

    using float_type    = double;
    auto flattened_grid = flatten_grid<float_type>(grid);

    tensorwrapper::shape::Smooth matrix_shape{ao_basis.n_aos(), n_points};
    std::vector<float_type> matrix_data(matrix_shape.size());

    std::size_t ao_i = 0;
    for(const auto& atomic_basis : ao_basis) {
        for(const auto& shell_i : atomic_basis) {
            const auto& cg          = shell_i.contracted_gaussian();
            const auto L            = shell_i.l();
            const auto n_primitives = cg.size();
            const auto n_aos        = shell_i.size();

            // TODO: Expose exponent_data/coefficient_data methods for Shells
            std::vector<double> exponents(n_primitives);
            std::vector<double> coefficients(n_primitives);
            std::vector<double> center(3);

            for(std::size_t i = 0; i < n_primitives; ++i) {
                exponents[i]    = cg.at(i).exponent();
                coefficients[i] = cg.at(i).coefficient();
            }
            center[0] = shell_i.center().x();
            center[1] = shell_i.center().y();
            center[2] = shell_i.center().z();

            auto is_pure = shell_i.pure() == chemist::ShellType::pure;
            auto order   = is_pure ? GG_SPHERICAL_CCA : GG_CARTESIAN_CCA;

            auto offset       = ao_i * n_points;
            auto shell_i_data = matrix_data.data() + offset;
            gg_collocation(static_cast<int>(L), static_cast<int>(n_points),
                           flattened_grid.data(), 3,
                           static_cast<int>(n_primitives), coefficients.data(),
                           exponents.data(), center.data(), order,
                           shell_i_data);

            ao_i += n_aos;
        }
    }
    tensorwrapper::buffer::Contiguous matrix_buffer(std::move(matrix_data),
                                                    matrix_shape);
    simde::type::tensor collocation_matrix(matrix_shape,
                                           std::move(matrix_buffer));

    auto rv = results();
    return pt::wrap_results(rv, std::move(collocation_matrix));
}

} // namespace scf::xc
