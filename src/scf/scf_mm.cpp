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

#include "fock_operator/fock_operator.hpp"
#include "scf_modules.hpp"
#include <scf/scf_mm.hpp>

namespace scf {

void load_modules(pluginplay::ModuleManager& mm) {
    fock_operator::load_modules(mm);
    mm.add_module<CoulombsLaw>("Coulomb's Law");
#ifdef BUILD_TAMM_SCF
    mm.add_module<TAMMEnergy>("SCF Energy via TAMM");
#endif
}

} // namespace scf
