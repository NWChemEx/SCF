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

#include "../libxc/libxc.hpp"
#include "native.hpp"
#include <simde/simde.hpp>

namespace scf::xc::native {
namespace {

const auto desc = R"(
SVWN3 Exchange-Correlation Energy
=================================

Computes the SVWN3 (Slater exchange plus VWN3 correlation) exchange-correlation
energy of a restricted wavefunction by numerical integration. The functional is
implemented natively, so this module works with any of the floating-point types
supported by the grid, AOs-on-a-grid, and density-on-a-grid submodules
(including the UQ types). For a spin-unpolarized density SVWN3 is identical to
LibXC's XC_LDA_X plus XC_LDA_C_VWN_3. Unlike LibXC, no density threshold is
applied.
)";

} // namespace

using XC_e_t = simde::type::XC_e_type;

template<typename WFType>
using pt = simde::eval_braket<WFType, XC_e_t, WFType>;

using grid_pt     = simde::MolecularGrid;
using rho2grid_pt = simde::EDensityCollocationMatrix;

MODULE_CTOR(SVWN3Energy) {
    using wf_type = simde::type::rscf_wf;
    satisfies_property_type<pt<wf_type>>();
    description(desc);
    add_submodule<grid_pt>("Integration grid");
    add_submodule<rho2grid_pt>("Density on a grid");
}

MODULE_RUN(SVWN3Energy) {
    using wf_type        = simde::type::rscf_wf;
    const auto& [braket] = pt<wf_type>::unwrap_inputs(inputs);

    const auto& bra_wf = braket.bra();
    const auto& xc_op  = braket.op();
    const auto& ket_wf = braket.ket();

    if(bra_wf != ket_wf)
        throw std::runtime_error("Expected the same basis set");

    if(xc_op.get_functional_name() !=
       chemist::qm_operator::xc_functional::SVWN3)
        throw std::runtime_error("SVWN3 Energy only supports SVWN3");

    const auto& P   = xc_op.get_rhs_particle();
    const auto& aos = P.basis_set().ao_basis_set();

    // Get grid
    auto& grid_mod   = submods.at("Integration grid");
    const auto& grid = grid_mod.run_as<grid_pt>(libxc::aos2molecule(aos));

    // Get density on grid
    auto& rho_mod   = submods.at("Density on a grid");
    const auto& rho = rho_mod.run_as<rho2grid_pt>(grid, P);

    auto e_xc    = svwn3_energy_density(rho);
    auto weights = tensorify_weights(grid, rho);

    simde::type::tensor exc;
    exc("") = weights("i") * e_xc("i");

    auto rv = results();
    return pt<wf_type>::wrap_results(rv, exc);
}

} // namespace scf::xc::native
