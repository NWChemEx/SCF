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
#include <chemcache/chemcache.hpp>
#include <integrals/integrals.hpp>
#include <nux/nux.hpp>
#include <scf/scf.hpp>
#include <simde/simde.hpp>

namespace test_scf {

/// Factors out setting submodules for SCF plugin from other plugins
template<typename FloatType>
pluginplay::ModuleManager load_modules() {
    auto rv = std::make_shared<parallelzone::runtime::RuntimeView>();
    pluginplay::ModuleManager mm(rv);
    scf::load_modules(mm);
    integrals::load_modules(mm);
    nux::load_modules(mm);
    chemcache::load_modules(mm);

    mm.change_submod("SCF Driver", "Hamiltonian",
                     "Born-Oppenheimer approximation");

    const auto ao_driver_key = "SCF integral driver";
    mm.change_submod(ao_driver_key, "Fundamental matrices",
                     "AO integral driver");

    mm.change_submod("Diagonalization Fock update", "Overlap matrix builder",
                     "Overlap");

    mm.change_submod("Loop", "Overlap matrix builder", "Overlap");

    mm.change_submod("SAD guess", "SAD Density", "sto-3g SAD density");

    if constexpr(!std::is_same_v<FloatType, double>) {
        mm.change_input("Evaluate 2-Index BraKet", "With UQ?", true);
        mm.change_input("Evaluate 4-Index BraKet", "With UQ?", true);
        mm.change_input("Overlap", "With UQ?", true);
        mm.change_input("ERI4", "With UQ?", true);
        mm.change_input("Kinetic", "With UQ?", true);
        mm.change_input("Nuclear", "With UQ?", true);
        mm.change_input("sto-3g atomic density matrix", "With UQ?", true);
    }

    return mm;
}

} // namespace test_scf