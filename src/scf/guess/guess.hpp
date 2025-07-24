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

namespace scf::guess {

DECLARE_MODULE(Core);
DECLARE_MODULE(SAD);

inline void load_modules(pluginplay::ModuleManager& mm) {
    mm.add_module<Core>("Core guess");
    mm.add_module<SAD>("SAD guess");
}

inline void set_defaults(pluginplay::ModuleManager& mm) {
    mm.change_submod("Core guess", "Build Fock operator",
                     "Restricted One-Electron Fock Op");
    mm.change_submod("Core guess", "Guess updater",
                     "Diagonalization Fock update");

    mm.change_submod("SAD guess", "Build Fock operator",
                     "Restricted One-Electron Fock Op");
    mm.change_submod("SAD guess", "Guess updater",
                     "Diagonalization Fock update");
}

} // namespace scf::guess