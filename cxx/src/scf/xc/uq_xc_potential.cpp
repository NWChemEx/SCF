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
#include <scf/property_types/xc_uncertain.hpp>
#include <simde/simde.hpp>

namespace scf::xc {
namespace {

const auto desc = R"(
Exchange-Correlation Potential with Uncertainty
===============================================

Computes the AO matrix of the exchange-correlation potential by quadrature and
attaches the uncertainty from the error of the grid. This module builds the
integration grid for the molecule the AO basis set is centered on, and passes
the grid, along with its inputs, to the "XC Potential" submodule, which
evaluates the quadrature on that grid. The inputs, the grid, and the result of
the quadrature are then passed to the "Incompleteness Error Model" submodule,
which attaches the uncertainty from the incompleteness of the grid. Since the
quadrature and the error model are given the same grid, they see the same grid.

Uncertainty in the inputs, *e.g.*, the density, is not attached here. If the
inputs are uncertain, it propagates through the quadrature.
)";

} // namespace

using pt      = simde::aos_xc_e_aos;
using grid_pt = simde::MolecularGrid;
using quad_pt = scf::property_types::xc_potential_quadrature;
using uq_pt   = scf::property_types::xc_potential_uncertainty;

MODULE_CTOR(UQXCPotential) {
    satisfies_property_type<pt>();
    description(desc);
    add_submodule<grid_pt>("Integration grid");
    add_submodule<quad_pt>("XC Potential");
    add_submodule<uq_pt>("Incompleteness Error Model");
}

MODULE_RUN(UQXCPotential) {
    const auto& [braket] = pt::unwrap_inputs(inputs);

    const auto& aos = braket.bra().ao_basis_set();

    auto& grid_mod   = submods.at("Integration grid");
    const auto& grid = grid_mod.run_as<grid_pt>(libxc::aos2molecule(aos));

    auto& quad_mod = submods.at("XC Potential");
    const auto vxc = quad_mod.run_as<quad_pt>(braket, grid);

    auto& error_mod   = submods.at("Incompleteness Error Model");
    const auto uq_vxc = error_mod.run_as<uq_pt>(braket, grid, vxc);

    auto rv = results();
    return pt::wrap_results(rv, uq_vxc);
}

} // namespace scf::xc
