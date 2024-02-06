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

#include <iostream>

#include <fastscf/fastscf.hpp>

int main(int argc, char** argv) {

    // Populate modules
    pluginplay::ModuleManager mm;
    fastscf::load_modules(mm);

    // TODO Create ChemicalSystem
    simde::type::chemical_system cs;

    // TODO Create BasisSet
    simde::type::ao_basis_set aos;

    // Run module
    auto E = mm.at("FastSCF Energy").run_as<simde::AOEnergy>(aos, cs);
    std::cout << "SCF Energy = " << E << " Eh" << std::endl;
    

    return 0;
}
