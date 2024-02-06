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

    // Create ChemicalSystem
    simde::type::atom h1("H", 1ul, 0.0, 0.0, 0.0, 0.0);
    simde::type::atom h2("H", 1ul, 0.0, 0.0, 0.0, 1.0);
    simde::type::molecule mol({h1, h2});
    simde::type::chemical_system cs(mol);

    // Create BasisSet
    std::vector<double> h_sto_alpha = {3.42525091, 0.62391373, 0.16885540};
    std::vector<double> h_sto_coeff = {0.15432897, 0.53532814, 0.44463454};
    simde::type::contracted_gaussian h1_cg(h_sto_coeff.begin(), h_sto_coeff.end(), h_sto_alpha.begin(), h_sto_alpha.end(), 0, 0, 0);
    simde::type::contracted_gaussian h2_cg(h_sto_coeff.begin(), h_sto_coeff.end(), h_sto_alpha.begin(), h_sto_alpha.end(), 0, 0, 1);
    simde::type::atomic_basis_set h1_bs("STO-3G", 1, 0.0, 0.0, 0.0);
    simde::type::atomic_basis_set h2_bs("STO-3G", 1, 0.0, 0.0, 1.0);
    h1_bs.add_shell(chemist::ShellType::pure, 0, h1_cg);
    h2_bs.add_shell(chemist::ShellType::pure, 0, h2_cg);

    simde::type::ao_basis_set aos;
    aos.add_center(h1_bs);
    aos.add_center(h2_bs);

    // Run module
    auto E = mm.at("FastSCF Energy").run_as<simde::AOEnergy>(aos, cs);
    std::cout << "SCF Energy = " << E << " Eh" << std::endl;
    

    return 0;
}
