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

#include "driver/driver.hpp"
#include "eigen_solver/eigen_solver.hpp"
#include "fock_operator/fock_operator.hpp"
#include "guess/guess.hpp"
#include "matrix_builder/matrix_builder.hpp"
#include "scf_modules.hpp"
#include "update/update.hpp"
#include "xc/xc.hpp"
#include <scf/scf_mm.hpp>

namespace scf {

PLUGIN_CTOR(SCF) {
    set_name("scf");

    // The defaults refer to modules from more than one subdirectory, so all
    // modules are added before any defaults are set
    driver::load_modules(modules());
    eigen_solver::load_modules(modules());
    fock_operator::load_modules(modules());
    guess::load_modules(modules());
    matrix_builder::load_modules(modules());
    update::load_modules(modules());
    xc::load_modules(modules());

    add_module<CoulombsLaw>("Coulomb's Law");

#ifdef BUILD_TAMM_SCF
    add_module<TAMMEnergy>("SCF Energy via TAMM");
#endif

    driver::set_defaults(modules());
    eigen_solver::set_defaults(modules());
    guess::set_defaults(modules());
    matrix_builder::set_defaults(modules());
    update::set_defaults(modules());
    xc::set_defaults(modules());

    // Submodules SCF does not provide come from these plugins, if loaded
    set_preferred_providers("nux", "integrals", "chemcache");

    default_submod_from("nux", "SCF Driver", "Hamiltonian",
                        "Born-Oppenheimer approximation");

    default_submod_from("integrals", "SCF integral driver",
                        "Fundamental matrices", "AO integral driver");
    default_submod_from("integrals", "Diagonalization Fock update",
                        "Overlap matrix builder", "Overlap");
    default_submod_from("integrals", "Loop", "Overlap matrix builder",
                        "Overlap");

    default_submod_from("chemcache", "SAD guess", "SAD Density",
                        "sto-3g SAD density");
}

} // namespace scf
