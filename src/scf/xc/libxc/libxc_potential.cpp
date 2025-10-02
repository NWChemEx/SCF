/*
 * Copyright 2022 NWChemEx-Project
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
const auto desc = R"(
Exchange-Correlation Potential via LibXC
========================================
)";
}

using pt            = simde::aos_xc_e_aos;
using grid_pt       = simde::MolecularGrid;
using ao_on_grid_pt = simde::AOCollocationMatrix;
using rho2grid_pt   = simde::EDensityCollocationMatrix;

MODULE_CTOR(LibXCPotential) {
    satisfies_property_type<pt>();
    description(desc);
    add_submodule<grid_pt>("Integration grid");
    add_submodule<ao_on_grid_pt>("AOs on a grid");
    add_submodule<rho2grid_pt>("Density on a grid");
}

MODULE_RUN(LibXCPotential) {
    const auto& [braket] = pt::unwrap_inputs(inputs);

    const auto& bra_aos = braket.bra();
    const auto& xc_op   = braket.op();
    const auto& ket_aos = braket.ket();

    if(bra_aos != ket_aos)
        throw std::runtime_error("Expected the same basis set!");

    const auto func = xc_op.functional_name();
    const auto& P   = xc_op.rhs_particle();
    const auto& aos = bra_aos.ao_basis_set();

    // Get grid
    auto& grid_mod   = submods.at("Integration grid");
    const auto& grid = grid_mod.run_as<grid_pt>(libxc::aos2molecule(aos));

    // Get AOs on the grid
    auto& aos_on_grid_mod   = submods.at("AOs on a grid");
    const auto& aos_on_grid = aos_on_grid_mod.run_as<ao_on_grid_pt>(grid, aos);

    // Get density on grid
    auto rho_mod    = submods.at("Density on a grid");
    const auto& rho = rho_mod.run_as<rho2grid_pt>(grid, P);

    auto dx_xc = libxc::libxc_lda_energy_density_derivative(func, rho);

    // Put weights into a tensor
    auto weights = tensorify_weights(grid, get_runtime());

    simde::type::tensor vxc, x;
    x("i") = dx_xc("i") * weights("i");

    auto temp  = libxc::weight_a_matrix(x, aos_on_grid);
    vxc("m,n") = temp("m,i") * aos_on_grid("n,i");

    auto rv = results();
    return pt::wrap_results(rv, vxc);
}

} // namespace scf::xc::libxc
