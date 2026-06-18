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

#include "h2_dimer_pencil.hpp"
#include "test_eigen_solver.hpp"

using types = std::tuple<float, double>;
using namespace test_eigen_solver;

TEMPLATE_LIST_TEST_CASE("GeneralizedEigenSolver H2 dimer", "", types) {
    using pt = simde::GeneralizedEigenSolve;
    pluginplay::ModuleManager mm;
    scf::load_modules(mm);

    auto rtol = std::is_same_v<TestType, float> ? 5e-4 : 1e-5;
    auto A    = h2_dimer_fock_as<TestType>();
    auto B    = h2_dimer_overlap_as<TestType>();

    auto& mod              = mm.at("Generalized eigensolve");
    auto [values, vectors] = mod.run_as<pt>(A, B);
    auto eval_corr         = h2_dimer_evals<TestType>();
    require_eigenvalues_approx(values, eval_corr, rtol);
    require_eigenpair_residual(A, values, vectors, rtol);
}
