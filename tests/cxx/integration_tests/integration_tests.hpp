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
#include "../test_scf.hpp"
#include <integrals/integrals.hpp>
#include <scf/scf.hpp>
#include <simde/simde.hpp>

namespace test_scf {

/// Factors out setting submodules for SCF plugin from other plugins
inline auto load_modules() {
    pluginplay::ModuleManager mm;
    scf::load_modules(mm);
    integrals::load_modules(mm);

    mm.change_submod("Four center J builder", "Four-center ERI", "ERI4");
    mm.change_submod("Four center K builder", "Four-center ERI", "ERI4");

    const auto ao_driver_key = "AO integral driver";
    mm.change_submod(ao_driver_key, "Kinetic", "Kinetic");
    mm.change_submod(ao_driver_key, "Electron-Nuclear attraction", "Nuclear");

    mm.change_submod("Diagonalization Fock update", "Overlap matrix builder",
                     "Overlap");

    return mm;
}

} // namespace test_scf