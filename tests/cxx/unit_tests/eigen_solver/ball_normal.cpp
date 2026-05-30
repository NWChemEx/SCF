/*
 * Copyright 2026 NWChemEx-Project
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

#include "test_eigen_solver.hpp"
#ifdef ENABLE_SIGMA
using types = std::tuple<tensorwrapper::types::idouble>;

TEMPLATE_LIST_TEST_CASE("BallNormal", "", types) {
    using uq_type = TestType;
    pluginplay::ModuleManager mm;
    scf::load_modules(mm);

    auto& mod = mm.at("Eigen Solve via Ball arithmetic");

    SECTION("classic 2 by 2 with noise") {
        auto system = test_eigen_solver::classic_2x2();
        auto A      = test_eigen_solver::noisy_matrix<uq_type>(system, 0.001);
        auto [values, vectors] = mod.run_as<simde::EigenSolve>(A);
        test_eigen_solver::require_uq_eigenvalues_contain<uq_type>(
          values, test_eigen_solver::eigenvalues_vector(system));
    }
}
#endif
