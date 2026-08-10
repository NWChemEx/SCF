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
#include <simde/simde.hpp>

namespace scf::eigen_solver {

DECLARE_MODULE(GeneralizedEigenSolver);
DECLARE_MODULE(EigenSolveDriver);
DECLARE_MODULE(EigenGeneralized);
DECLARE_MODULE(EigenNormal);
DECLARE_MODULE(JacobiNormal);

inline void set_defaults(pluginplay::ModuleManager& mm) {
    mm.change_submod("Eigen Solve", "none", "Eigen Solve via Eigen");
    mm.change_submod("Eigen Solve", "uncertain", "Eigen Solve via Jacobi");
    mm.change_submod("Eigen Solve", "interval", "Eigen Solve via Jacobi");
    mm.change_submod("Generalized eigensolve", "Eigen Solve", "Eigen Solve");
}

inline void load_modules(pluginplay::ModuleManager& mm) {
    mm.add_module<EigenSolveDriver>("Eigen Solve");
    mm.add_module<EigenNormal>("Eigen Solve via Eigen");
    mm.add_module<JacobiNormal>("Eigen Solve via Jacobi");
    mm.add_module<EigenGeneralized>("Generalized eigensolve via Eigen");
    mm.add_module<GeneralizedEigenSolver>("Generalized eigensolve");
    set_defaults(mm);
}

} // namespace scf::eigen_solver
