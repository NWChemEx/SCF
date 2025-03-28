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

using tensorwrapper::operations::approximately_equal;

TEST_CASE("TAMM SCF") {
    // Populate modules
    pluginplay::ModuleManager mm;
    chemcache::load_modules(mm);
    scf::load_modules(mm);

    bool test = mm.count("SCF Energy via TAMM") != 0;
    if(test) {
        // Create ChemicalSystem
        std::string mol_name = "water";
        auto mol =
          mm.at("NWX Molecules").run_as<simde::MoleculeFromString>(mol_name);
        simde::type::chemical_system cs(mol);

        // Create BasisSet
        std::string basis_name =
          "sto-3g"; // This is the only supported basis in ChemCache
        auto aos = mm.at(basis_name).run_as<simde::MolecularBasisSet>(mol);

        // Run module
        mm.change_input("SCF Energy via TAMM", "molecule_name", mol_name);
        auto E = mm.at("SCF Energy via TAMM").run_as<simde::AOEnergy>(aos, cs);
        std::cout << "SCF Energy = " << E << " Hartree" << std::endl;

        simde::type::tensor corr{-74.3670617803483};
        REQUIRE(approximately_equal(corr, E, 1E-6));
    }
}

TEST_CASE("TAMM DFT") {
    // Populate modules
    pluginplay::ModuleManager mm;
    chemcache::load_modules(mm);
    scf::load_modules(mm);

    bool test = mm.count("SCF Energy via TAMM") != 0;
    if(test) {
        // Create ChemicalSystem
        std::string mol_name = "water";
        auto mol =
          mm.at("NWX Molecules").run_as<simde::MoleculeFromString>(mol_name);
        simde::type::chemical_system cs(mol);

        // Create BasisSet
        std::string basis_name =
          "sto-3g"; // This is the only supported basis in ChemCache
        auto aos = mm.at(basis_name).run_as<simde::MolecularBasisSet>(mol);

        // Run module
        std::vector<std::string> xc_type = {"pbe0"};
        mm.change_input("SCF Energy via TAMM", "xc_type", xc_type);
        mm.change_input("SCF Energy via TAMM", "molecule_name", mol_name);
        auto E = mm.at("SCF Energy via TAMM").run_as<simde::AOEnergy>(aos, cs);
        std::cout << "SCF Energy = " << E << " Hartree" << std::endl;

        simde::type::tensor corr{-74.81168986385825};
        REQUIRE(approximately_equal(corr, E, 1E-6));
    }
}
