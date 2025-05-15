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

#pragma once
#include <gauxc/basisset.hpp>
#include <gauxc/grid_factory.hpp>
#include <gauxc/load_balancer.hpp>
#include <gauxc/molecule.hpp>
#include <simde/simde.hpp>

namespace scf::xc::gauxc {

using xc_grid_type = GauXC::AtomicGridSizeDefault;

using gauxc_basis_conversion_t =
  simde::Convert<GauXC::BasisSet<double>, simde::type::ao_basis_set>;

using basis_to_gauxc_molecule_conversion_t =
  simde::Convert<GauXC::Molecule, simde::type::ao_basis_set>;

DECLARE_PROPERTY_TYPE(XCQuadratureBatches);

PROPERTY_TYPE_INPUTS(XCQuadratureBatches) {
    using basis_type          = simde::type::ao_basis_set;
    using pruning_scheme_type = GauXC::PruningScheme;
    using radial_quad_type    = GauXC::RadialQuad;

    auto rv = pluginplay::declare_input()
                .add_field<const basis_type&>("AO Basis")
                .add_field<const pruning_scheme_type&>("Pruning Scheme")
                .add_field<const radial_quad_type&>("Radial Quadrature")
                .add_field<const xc_grid_type&>("XC Grid Specification");
    return rv;
}

PROPERTY_TYPE_RESULTS(XCQuadratureBatches) {
    using xc_quad_batch_type = GauXC::LoadBalancer;
    return pluginplay::declare_result().add_field<xc_quad_batch_type>(
      "Quadrature Batches");
}

} // namespace scf::xc::gauxc
