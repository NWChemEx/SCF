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

#include "gauxc/gauxc.hpp"
#include "libxc/libxc.hpp"
#include "xc.hpp"

namespace scf::xc {
void load_modules(pluginplay::ModuleManager& mm) {
    gauxc::load_modules(mm);
    libxc::load_modules(mm);
    mm.add_module<AOsOnGrid>("AOs on a grid");
    mm.add_module<Gau2Grid>("Gau2Grid");
    mm.add_module<GridFromFile>("Grid From File");
    mm.add_module<Density2Grid>("Density2Grid");
}

void set_defaults(pluginplay::ModuleManager& mm) {
    gauxc::set_defaults(mm);
    libxc::set_defaults(mm);
    mm.change_submod("Density2Grid", "AOs on a grid", "AOs on a grid");
}

} // namespace scf::xc
