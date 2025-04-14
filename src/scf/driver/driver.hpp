/*
 * Copyright 2024 NWChemEx-Project
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
#include <simde/simde.hpp>

namespace scf::driver {

DECLARE_MODULE(SCFDriver);
DECLARE_MODULE(SCFLoop);

inline void load_modules(pluginplay::ModuleManager& mm) {
    mm.add_module<SCFDriver>("SCF Driver");
    mm.add_module<SCFLoop>("Loop");
}

inline void set_defaults(pluginplay::ModuleManager& mm) {
    mm.change_submod("Loop", "Electronic energy", "Electronic energy");
    mm.change_submod("Loop", "Density matrix", "Density matrix builder");
    mm.change_submod("Loop", "Guess update", "Diagonalization Fock update");
    mm.change_submod("Loop", "One-electron Fock operator",
                     "Restricted One-Electron Fock op");
    mm.change_submod("Loop", "Fock operator", "Restricted Fock Op");
    mm.change_submod("Loop", "Charge-charge", "Coulomb's Law");
    mm.change_submod("Loop", "Fock matrix builder", "Fock matrix builder");

    mm.change_submod("SCF Driver", "Guess", "Core guess");
    mm.change_submod("SCF Driver", "Optimizer", "Loop");
}

} // namespace scf::driver
