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

void load_modules(pluginplay::ModuleManager& mm) {
    driver::load_modules(mm);
    eigen_solver::load_modules(mm);
    fock_operator::load_modules(mm);
    guess::load_modules(mm);
    matrix_builder::load_modules(mm);
    update::load_modules(mm);
    xc::load_modules(mm);

    mm.add_module<CoulombsLaw>("Coulomb's Law");

    // mm.add_module<SCFDriver>("Driver");
#ifdef BUILD_TAMM_SCF
    mm.add_module<TAMMEnergy>("SCF Energy via TAMM");
#endif

    driver::set_defaults(mm);
    guess::set_defaults(mm);
    matrix_builder::set_defaults(mm);
    update::set_defaults(mm);
    xc::set_defaults(mm);
}

} // namespace scf
