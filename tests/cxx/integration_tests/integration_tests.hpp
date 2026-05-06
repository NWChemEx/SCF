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
    integrals::set_defaults(mm);
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

    if constexpr(std::is_same_v<FloatType, tensorwrapper::types::udouble>) {
        std::string key = "uncertain";
        mm.change_input("Evaluate 2-Index BraKet", "UQ Type", key);
        mm.change_input("Evaluate 4-Index BraKet", "UQ Type", key);
        mm.change_input("Overlap", "UQ Type", key);
        mm.change_input("ERI4", "UQ Type", key);
        mm.change_input("Kinetic", "UQ Type", key);
        mm.change_input("Nuclear", "UQ Type", key);
        mm.change_input("sto-3g atomic density matrix", "With UQ?", true);
    } else if constexpr(tensorwrapper::types::is_interval_v<FloatType>) {
        std::string key = "interval";
        mm.at("Generalized eigensolve via Eigen").turn_off_memoization();
        mm.at("Density matrix builder").turn_off_memoization();
        mm.at("Electronic energy").turn_off_memoization();
        mm.at("Fock matrix builder").turn_off_memoization();
        mm.at("Restricted One-Electron Fock Op").turn_off_memoization();
        mm.at("Restricted Fock Op").turn_off_memoization();
        mm.at("Four center J builder").turn_off_memoization();
        mm.at("Four center K builder").turn_off_memoization();
        mm.change_input("Evaluate 2-Index BraKet", "UQ Type", key);
        mm.change_input("Evaluate 4-Index BraKet", "UQ Type", key);
        mm.change_input("Overlap", "UQ Type", key);
        mm.change_input("ERI4", "UQ Type", key);
        mm.change_input("Kinetic", "UQ Type", key);
        mm.change_input("Nuclear", "UQ Type", key);
        mm.change_submod("Loop", "Diagonalizer",
                         "Generalized eigensolve via Ball Arithmetic");
        mm.change_submod("Diagonalization Fock Update", "Diagonalizer",
                         "Generalized eigensolve via Ball Arithmetic");
    }

    return mm;
}

} // namespace test_scf
