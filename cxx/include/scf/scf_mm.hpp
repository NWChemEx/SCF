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
/** @brief The SCF plugin.
 *
 *  Load it into a ModuleManager with `mm.load_plugin<scf::SCF>()`
 *  (`mm.load_plugin(scf.SCF())` from Python). The plugin's name is "scf".
 *
 *  Several SCF modules need a submodule that SCF itself does not provide (the
 *  Hamiltonian, the AO integrals, the overlap matrix, and the SAD density).
 *  The plugin records defaults for those submodules from the NUX, Integrals,
 *  and ChemCache plugins, which ModuleManager::wire_defaults sets once the
 *  plugins are loaded. Submodules whose provider is not loaded are left unset,
 *  and SCF does not depend on any of those plugins.
 */
DECLARE_PLUGIN(SCF);

} // namespace scf
