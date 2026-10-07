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
#include <pluginplay/pluginplay.hpp>

namespace scf {
/** @brief Loads the modules contained in the SCF module collection into the
 *         provided ModuleManager instance.
 */
void load_modules(pluginplay::ModuleManager& mm);

/** @brief Connects SCF's modules to the modules other plugins provide.
 *
 *  Several SCF modules need a submodule that SCF itself does not provide (the
 *  Hamiltonian, the AO integrals, the overlap matrix, and the SAD density).
 *  This function sets those submodules to the default implementations from
 *  the NUX, Integrals, and ChemCache plugins.
 *
 *  Each submodule is only set if the module providing it has been loaded into
 *  @p mm, so this function must be called after those plugins' load_modules.
 *  Submodules whose provider is missing are left unset, and SCF does not
 *  depend on any of those plugins.
 *
 *  @param[in,out] mm The ModuleManager to set the submodules in. Must already
 *                    contain SCF's modules (see load_modules).
 */
void set_defaults(pluginplay::ModuleManager& mm);

} // namespace scf
