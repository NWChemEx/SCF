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
#include <scf/property_types.hpp>
#include <simde/simde.hpp>

namespace scf::xc {
namespace {

const auto desc = R"(
Exchange-Correlation Energy via Quadrature
==========================================

Computes the exchange-correlation energy of a restricted wavefunction by
quadrature. This module builds the integration grid for the molecule the AO
basis set is centered on, and then passes the grid, along with its inputs, to
the "Quadrature" submodule, which evaluates the quadrature on that grid. The
choice of functional, and the checks on the inputs, are made by the
"Quadrature" submodule.
)";

} // namespace

using wf_type = simde::type::rscf_wf;
using pt      = simde::eval_braket<wf_type, simde::type::XC_e_type, wf_type>;
using grid_pt = simde::MolecularGrid;
using quad_pt = scf::property_types::xc_energy_quadrature;

MODULE_CTOR(XCEnergyViaQuadrature) {
    satisfies_property_type<pt>();
    description(desc);
    add_submodule<grid_pt>("Integration grid");
    add_submodule<quad_pt>("Quadrature");
}

MODULE_RUN(XCEnergyViaQuadrature) {
    const auto& [braket] = pt::unwrap_inputs(inputs);

    const auto& P   = braket.op().get_rhs_particle();
    const auto& aos = P.basis_set().ao_basis_set();

    auto& grid_mod   = submods.at("Integration grid");
    const auto& grid = grid_mod.run_as<grid_pt>(libxc::aos2molecule(aos));

    auto& quad_mod = submods.at("Quadrature");
    const auto exc = quad_mod.run_as<quad_pt>(braket, grid);

    auto rv = results();
    return pt::wrap_results(rv, exc);
}

} // namespace scf::xc
