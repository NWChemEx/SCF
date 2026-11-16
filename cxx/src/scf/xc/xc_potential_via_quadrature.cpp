/*
 * Copyright 2026 NWChemEx-Project
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
#include <scf/property_types/quadrature.hpp>
#include <simde/simde.hpp>

namespace scf::xc {
namespace {

const auto desc = R"(
Exchange-Correlation Potential via Quadrature
=============================================

Computes the AO matrix of the exchange-correlation potential by quadrature.
This module builds the integration grid for the molecule the AO basis set is
centered on, and then passes the grid, along with its inputs, to the
"Quadrature" submodule, which evaluates the quadrature on that grid. The choice
of functional, and the checks on the inputs, are made by the "Quadrature"
submodule.
)";

} // namespace

using pt      = simde::aos_xc_e_aos;
using grid_pt = simde::MolecularGrid;
using quad_pt = scf::property_types::xc_potential_quadrature;

MODULE_CTOR(XCPotentialViaQuadrature) {
    satisfies_property_type<pt>();
    description(desc);
    add_submodule<grid_pt>("Integration grid");
    add_submodule<quad_pt>("Quadrature");
}

MODULE_RUN(XCPotentialViaQuadrature) {
    const auto& [braket] = pt::unwrap_inputs(inputs);

    const auto& aos = braket.bra().ao_basis_set();

    auto& grid_mod   = submods.at("Integration grid");
    const auto& grid = grid_mod.run_as<grid_pt>(libxc::aos2molecule(aos));

    auto& quad_mod = submods.at("Quadrature");
    const auto vxc = quad_mod.run_as<quad_pt>(braket, grid);

    auto rv = results();
    return pt::wrap_results(rv, vxc);
}

} // namespace scf::xc
