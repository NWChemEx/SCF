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

#include "libxc/libxc.hpp"
#include "xc.hpp"
#include <simde/integration_grids/collocation_matrix.hpp>

namespace scf::xc {
namespace {
const auto desc = R"(

DensityCollocationMatrix
-----------------
)";

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

    simde::type::tensor X, rho2;
    rho2("m,n") = rho("m,n") * 2.0; // Assumes restricted orbitals
    X("m,i")    = rho2("m,n") * aos_on_grid("n,i");

    auto rho_on_grid = libxc::batched_dot(aos_on_grid, X);

    auto rv = results();
    return pt::wrap_results(rv, std::move(rho_on_grid));
}

} // namespace scf::xc
