// /*
//  * Copyright 2022 NWChemEx-Project
//  *
//  * Licensed under the Apache License, Version 2.0 (the "License");
//  * you may not use this file except in compliance with the License.
//  * You may obtain a copy of the License at
//  *
//  * http://www.apache.org/licenses/LICENSE-2.0
//  *
//  * Unless required by applicable law or agreed to in writing, software
//  * distributed under the License is distributed on an "AS IS" BASIS,
//  * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  * See the License for the specific language governing permissions and
//  * limitations under the License.
//  */

#include "gauxc.hpp"

#include <gauxc/molecular_weights.hpp>
#include <gauxc/molgrid/defaults.hpp>
#include <gauxc/runtime_environment.hpp>

namespace scf::xc::gauxc {

using xc_grid_type = GauXC::AtomicGridSizeDefault;

/**
 *  All default grids are product quadratures consisting of NR Mura-Knowles
 (MK)
 *  radial quadratures and NA Lebedev (L) spherical quadratures.
 *
 *  Fine(Grid)
 *    - 75 (MK) x 302 (L) grid for each atom
 *    - Target Accuracy 1e-6
 *  UltraFine(Grid)
 *    - 99 (MK) x 590 (L) grid for each atom
 *    - Target Accuracy 1e-8
 *  SuperFine(Grid)
 *    - 175 (MK) x 974 (L) grid for atoms with Z <= 2
 *    - 250 (MK) x 974 (L) grid for atoms with Z >= 3
 *    - Target Accuracy 1e-10
 */
static utilities::CaseInsensitiveMap<xc_grid_type> atomic_grid_types = {
  {"Fine", xc_grid_type::FineGrid},
  {"UltraFine", xc_grid_type::UltraFineGrid},
  {"SuperFine", xc_grid_type::SuperFineGrid},
  {"FineGrid", xc_grid_type::FineGrid},
  {"UltraFineGrid", xc_grid_type::UltraFineGrid},
  {"SuperFineGrid", xc_grid_type::SuperFineGrid}};

MODULE_CTOR(QuadratureBatches) {
    satisfies_property_type<XCQuadratureBatches>();

    add_submodule<gauxc_basis_conversion_t>("GauXC Basis Converter")
      .set_description("Converts NWX Basis -> GauXC Basis");
    add_submodule<basis_to_gauxc_molecule_conversion_t>(
      "GauXC Molecule Converter")
      .set_description("Converts NWX Basis -> GauXC Molecule");

    add_input<std::string>("Grid Type").set_default("UltraFine");

    add_input<GauXC::PruningScheme>("Pruning Scheme")
      .set_default(GauXC::PruningScheme::Unpruned);

    add_input<GauXC::RadialQuad>("Radial Quadrature Type")
      .set_default(GauXC::RadialQuad::MuraKnowles);

    /// TODO: Replace with information from the runtime
    add_input<bool>("On GPU").set_default(false);
}

MODULE_RUN(QuadratureBatches) {
    const auto& [obs] = XCQuadratureBatches::unwrap_inputs(inputs);

    // Unpack module-specific inputs
    auto grid_spec =
      atomic_grid_types.at(inputs.at("Grid Type").value<std::string>());

    auto pruning_scheme =
      inputs.at("Pruning Scheme").value<GauXC::PruningScheme>();

    auto radial_quad =
      inputs.at("Radial Quadrature Type").value<GauXC::RadialQuad>();

    auto on_gpu = inputs.at("On GPU").value<bool>();

    auto ex_space = (on_gpu) ?
                      GauXC::ExecutionSpace::Device :
                      GauXC::ExecutionSpace::Host; // TODO handle CUDA/HIP
    auto comm     = get_runtime().mpi_comm();
    std::shared_ptr<GauXC::RuntimeEnvironment> gauxc_rt =
#ifdef GAUXC_ENABLE_DEVICE
      (on_gpu) ? std::make_shared<GauXC::DeviceRuntimeEnvironment>(comm, 0.8) :
#endif
                 std::make_shared<GauXC::RuntimeEnvironment>(comm);

    // Convert NWX Molecule -> GauXC Molecule
    auto& mol_submod = submods.at("GauXC Molecule Converter");
    const auto& gauxc_mol =
      mol_submod.run_as<basis_to_gauxc_molecule_conversion_t>(obs);

    // Convert NWX Basis -> GauXC Basis
    auto& basis_submod = submods.at("GauXC Basis Converter");
    const auto& gauxc_basis =
      basis_submod.run_as<gauxc_basis_conversion_t>(obs);

    // Set up integration grid
    const auto bsz = GauXC::BatchSize(512);
    auto mol_grid  = GauXC::MolGridFactory::create_default_molgrid(
      gauxc_mol, pruning_scheme, bsz, radial_quad, grid_spec);

    // Generate LoadBalancer
    GauXC::LoadBalancerFactory lb_factory(ex_space, "Default");
    auto lb =
      lb_factory.get_instance(*gauxc_rt, gauxc_mol, mol_grid, gauxc_basis);

    // Apply Molecular Partition Weights
    GauXC::MolecularWeightsFactory mw_factory(
      ex_space, "Default", GauXC::MolecularWeightsSettings{});
    mw_factory.get_instance().modify_weights(lb);

    auto rv = results();
    return XCQuadratureBatches::wrap_results(rv, lb);
}

} // namespace scf::xc::gauxc
