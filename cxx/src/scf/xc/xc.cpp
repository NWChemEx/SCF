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

#ifdef BUILD_GAUXC
#include "gauxc/gauxc.hpp"
#endif
#include "libxc/libxc.hpp"
#include "native/native.hpp"
#include "xc.hpp"

namespace scf::xc {
void load_modules(pluginplay::ModuleMap& modules) {
#ifdef BUILD_GAUXC
    gauxc::load_modules(modules);
#endif
    libxc::load_modules(modules);
    native::load_modules(modules);
    modules.add_module<AOsOnGrid>("AOs on a grid");
    modules.add_module<Gau2Grid>("Gau2Grid");
    modules.add_module<GridFromFile>("Grid From File");
    modules.add_module<GridFromIntegratorXX>("Grid From IntegratorXX");
    modules.add_module<Density2Grid>("Density2Grid");
    modules.add_module<UQXCEnergy>("UQ XC Energy");
    modules.add_module<UQXCPotential>("UQ XC Potential");
}

void set_defaults(pluginplay::ModuleMap& modules) {
#ifdef BUILD_GAUXC
    gauxc::set_defaults(modules);
#endif
    libxc::set_defaults(modules);
    native::set_defaults(modules);
    modules.change_submod("Density2Grid", "AOs on a grid", "AOs on a grid");

    // The UQ XC modules default to the native SVWN3 quadratures, since they
    // support the UQ types. "Incompleteness Error Model" has no default yet.
    const auto grid_key = "Grid From IntegratorXX";
    modules.change_submod("UQ XC Energy", "Integration grid", grid_key);
    modules.change_submod("UQ XC Potential", "Integration grid", grid_key);
    modules.change_submod("UQ XC Energy", "XC Energy",
                          "SVWN3 Energy Quadrature");
    modules.change_submod("UQ XC Potential", "XC Potential",
                          "SVWN3 Potential Quadrature");
}

} // namespace scf::xc
