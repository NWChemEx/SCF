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

#pragma once
#include "gauxc_property_types.hpp"
#include <simde/simde.hpp>

namespace scf::xc::gauxc {

DECLARE_MODULE(GauXCDriver);
DECLARE_MODULE(XCPotential);
DECLARE_MODULE(XCEnergy);

DECLARE_MODULE(snLinK);

// Module to convert Chemist Basis -> GauXC Basis
DECLARE_MODULE(BasisConversion);

// Module to convert Chemist AO Basis -> GauXC Molecule
DECLARE_MODULE(MoleculeConversion);

// Modules to generate XC Quadrature Batches
DECLARE_MODULE(QuadratureBatches);

void set_defaults(pluginplay::ModuleManager& mm);
void load_modules(pluginplay::ModuleManager& mm);

} // namespace scf::xc::gauxc