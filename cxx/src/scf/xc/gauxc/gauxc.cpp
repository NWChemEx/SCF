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

#include "gauxc.hpp"

namespace scf::xc::gauxc {

void load_modules(pluginplay::ModuleMap& modules) {
    // Conversion Utilities
    modules.add_module<BasisConversion>("GauXC Basis Converter");
    modules.add_module<MoleculeConversion>("GauXC Molecule Converter");

    // Data geneation utilities
    modules.add_module<QuadratureBatches>("GauXC Quadrature Batches");

    // XC Integration
    modules.add_module<GauXCDriver>("GauXC Driver");
    modules.add_module<XCPotential>("GauXC XC Potential");
    modules.add_module<XCEnergy>("GauXC XC Energy");

    // sn-LinK Integration
    modules.add_module<snLinK>("snLinK");
}

void set_defaults(pluginplay::ModuleMap& modules) {
    // Data geneation utilities
    modules.change_submod("GauXC Quadrature Batches", "GauXC Basis Converter",
                          "GauXC Basis Converter");
    modules.change_submod("GauXC Quadrature Batches",
                          "GauXC Molecule Converter",
                          "GauXC Molecule Converter");

    // XC Integration
    modules.change_submod("GauXC Driver", "Quadrature Batches",
                          "GauXC Quadrature Batches");

    modules.change_submod("GauXC XC Potential", "XC Driver", "GauXC Driver");
    modules.change_submod("GauXC XC Energy", "XC Driver", "GauXC Driver");

    // sn-LinK Integration
    modules.change_submod("snLinK", "Quadrature Batches",
                          "GauXC Quadrature Batches");
}

} // namespace scf::xc::gauxc
