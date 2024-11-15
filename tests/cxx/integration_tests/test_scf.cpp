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

#include <catch2/catch_approx.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <chemcache/chemcache.hpp>
#include <iostream>
#include <scf/scf.hpp>

TEST_CASE("SCF") {

    // Populate modules
    pluginplay::ModuleManager mm;
    chemcache::load_modules(mm);
    scf::load_modules(mm);

    // Create ChemicalSystem
    std::string mol_name = "water";
    auto mol = mm.at("NWX Molecules").run_as<simde::MoleculeFromString>(mol_name);
    simde::type::chemical_system cs(mol);

    // Create BasisSet
    std::string basis_name = "sto-3g"; // This is the only supported basis in ChemCache
    auto aos = mm.at(basis_name).run_as<simde::MolecularBasisSet>(mol);

    // Run module
    mm.change_input("SCF Energy", "molecule_name", mol_name);
    auto E = mm.at("SCF Energy").run_as<simde::AOEnergy>(aos, cs);
    std::cout << "SCF Energy = " << E << " Hartree" << std::endl;
    
    REQUIRE(E == Catch::Approx(-74.3670617803483).margin(1.0e-6));
}

TEST_CASE("DFT") {

    // Populate modules
    pluginplay::ModuleManager mm;
    chemcache::load_modules(mm);
    scf::load_modules(mm);

    // Create ChemicalSystem
    std::string mol_name = "water";
    auto mol = mm.at("NWX Molecules").run_as<simde::MoleculeFromString>(mol_name);
    simde::type::chemical_system cs(mol);

    // Create BasisSet
    std::string basis_name = "sto-3g"; // This is the only supported basis in ChemCache
    auto aos = mm.at(basis_name).run_as<simde::MolecularBasisSet>(mol);

    // Run module
    std::vector<std::string> xc_type = {"pbe0"};
    mm.change_input("SCF Energy", "xc_type", xc_type);
    mm.change_input("SCF Energy", "molecule_name", mol_name);
    auto E = mm.at("SCF Energy").run_as<simde::AOEnergy>(aos, cs);
    std::cout << "SCF Energy = " << E << " Hartree" << std::endl;
        
    REQUIRE(E == Catch::Approx(-74.81168986385825).margin(1.0e-6));
}