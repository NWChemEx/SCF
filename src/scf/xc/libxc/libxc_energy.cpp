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
#include <simde/simde.hpp>

namespace scf::xc::libxc {
namespace {

const auto desc = R"()";

} // namespace

using XC_e_t = simde::type::XC_e_type;

template<typename WFType>
using pt = simde::eval_braket<WFType, XC_e_t, WFType>;

using grid_pt     = simde::MolecularGrid;
using rho2grid_pt = simde::EDensityCollocationMatrix;

MODULE_CTOR(LibXCEnergy) {
    using wf_type = simde::type::rscf_wf;
    satisfies_property_type<pt<wf_type>>();
    description(desc);
    add_submodule<grid_pt>("Integration grid");
    add_submodule<rho2grid_pt>("Density on a grid");
}

MODULE_RUN(LibXCEnergy) {
    using wf_type        = simde::type::rscf_wf;
    const auto& [braket] = pt<wf_type>::unwrap_inputs(inputs);

    const auto& bra_wf = braket.bra();
    const auto& xc_op  = braket.op();
    const auto& ket_wf = braket.ket();

    if(bra_wf != ket_wf)
        throw std::runtime_error("Expected the same basis set");

    const auto func_name = xc_op.functional_name();
    const auto& P        = xc_op.rhs_particle();
    const auto& aos      = P.basis_set().ao_basis_set();

    // Molecule from AOs
    const auto mol = libxc::aos2molecule(aos);

    // Get grid
    auto& grid_mod   = submods.at("Integration grid");
    const auto& grid = grid_mod.run_as<simde::MolecularGrid>(mol);

    // Get density on grid
    auto rho_mod    = submods.at("Density on a grid");
    const auto& rho = rho_mod.run_as<rho2grid_pt>(grid, P);

    auto x_xc = libxc::libxc_lda_energy_density(func_name, rho);

    // Put weights into a tensor
    auto weight = libxc::tensorify_weights(grid, get_runtime());

    simde::type::tensor exc_xc;
    exc_xc("") = weight("i") * x_xc("i");

    simde::type::tensor n_electrons;
    n_electrons("") = weight("i") * rho("i");

    auto rv = results();
    return pt<wf_type>::wrap_results(rv, exc_xc);
}

} // namespace scf::xc::libxc
