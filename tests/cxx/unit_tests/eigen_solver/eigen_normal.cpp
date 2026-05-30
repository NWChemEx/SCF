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

#include "test_eigen_solver.hpp"

using types = std::tuple<float, double>;
using namespace test_eigen_solver;
using namespace tensorwrapper::generate;

TEMPLATE_LIST_TEST_CASE("EigenNormal", "", types) {
    pluginplay::ModuleManager mm;
    scf::load_modules(mm);

    auto& mod = mm.at("Eigen Solve via Eigen");
    auto rtol = std::is_same_v<TestType, float> ? 5e-4 : 1e-5;
    using pt  = simde::EigenSolve;
    SECTION("classic 2 by 2") {
        auto system            = test_eigen_solver::classic_2x2<TestType>();
        auto [values, vectors] = mod.run_as<pt>(system.matrix);

        require_eigenvalues_approx(values, system.eigenvalues, rtol);
        require_eigenpair_residual(system.matrix, values, vectors, rtol);
    }

    SECTION("generated n=4 condition number 1e3") {
        test_eigen_solver::SymmetricMatrixSpec spec;
        spec.n                 = 4;
        spec.condition_number  = 1e3;
        spec.spacing           = test_eigen_solver::EigenvalueSpacing::Linear;
        spec.seed              = 11;
        auto system            = generate_eigen_system<TestType>(spec);
        auto [values, vectors] = mod.run_as<pt>(system.matrix);
        require_eigenvalues_approx(values, system.eigenvalues, rtol);
        require_eigenpair_residual(system.matrix, values, vectors, rtol);
    }

    SECTION("generated clustered n=6") {
        test_eigen_solver::SymmetricMatrixSpec spec;
        spec.n                = 6;
        spec.condition_number = 100.0;
        spec.spacing          = test_eigen_solver::EigenvalueSpacing::Clustered;
        spec.n_clusters       = 3;
        spec.cluster_width    = 1e-8;
        spec.seed             = 23;
        auto system           = generate_eigen_system<TestType>(spec);
        auto [values, vectors] = mod.run_as<pt>(system.matrix);
        require_eigenvalues_approx(values, system.eigenvalues, rtol);
        require_eigenpair_residual(system.matrix, values, vectors, rtol);
    }
}
