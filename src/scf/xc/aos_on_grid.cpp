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

AO CollocationMatrix
-----------------
)";

} // namespace

using pt = simde::AOCollocationMatrix;

MODULE_CTOR(AOsOnGrid) {
    satisfies_property_type<pt>();
    description(desc);
}

MODULE_RUN(AOsOnGrid) {
    const auto& [grid, ao_basis] = pt::unwrap_inputs(inputs);
    auto n_points                = grid.size();
    auto n_aos                   = ao_basis.n_aos();
    using float_type             = double;

    tensorwrapper::shape::Smooth matrix_shape{n_aos, n_points};
    auto buffer =
      tensorwrapper::buffer::make_contiguous<float_type>(matrix_shape);

    std::vector<std::size_t> idx{0, 0};
    for(const auto& atomic_basis : ao_basis) {
        for(const auto& shell_i : atomic_basis) {
            assert(shell_i.l() == 0); // only s is supported for now
            const auto& cg = shell_i.contracted_gaussian();
            for(; idx[1] < n_points; ++idx[1]) {
                float_type ao_value = 0.0;
                for(const auto& prim : cg) {
                    // TODO: update when eval accounts for normalization
                    const auto val = prim.evaluate(grid.at(idx[1]).point());
                    const auto exponent = prim.exponent();
                    auto norm = std::sqrt(std::pow(2.0 * exponent / M_PI, 1.5));
                    ao_value += norm * val;
                }
                buffer.set_elem(idx, ao_value);
            }
            idx[0] += shell_i.size();
            idx[1] = 0;
        }
    }
    simde::type::tensor collocation_matrix(matrix_shape, std::move(buffer));
    auto rv = results();
    return pt::wrap_results(rv, std::move(collocation_matrix));
};

} // namespace scf::xc
