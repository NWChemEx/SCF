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

using types = std::tuple<float, double, tensorwrapper::types::idouble>;
using namespace test_eigen_solver;
using namespace tensorwrapper::generate;

TEMPLATE_LIST_TEST_CASE("QRNormal", "", types) {
    pluginplay::ModuleManager mm;
    scf::load_modules(mm);

    auto& mod = mm.at("Eigen Solve via QR");
    auto rtol = std::is_same_v<TestType, float> ? 5e-4 : 1e-5;
    using pt  = simde::EigenSolve;

    SECTION("classic 2 by 2") {
        auto system            = classic_2x2<TestType>();
        auto [values, vectors] = mod.run_as<pt>(system.matrix);

        require_eigenvalues_approx(values, system.eigenvalues, rtol);
        require_eigenpair_residual(system.matrix, values, vectors, rtol);
    }

    SECTION("generated n=4 condition number 1e3") {
        SymmetricMatrixSpec spec;
        spec.n                 = 4;
        spec.condition_number  = 1e3;
        spec.spacing           = EigenvalueSpacing::Linear;
        spec.seed              = 11;
        auto system            = generate_eigen_system<TestType>(spec);
        auto [values, vectors] = mod.run_as<pt>(system.matrix);
        require_eigenvalues_approx(values, system.eigenvalues, 10 * rtol);
        require_eigenpair_residual(system.matrix, values, vectors, rtol);
    }

    SECTION("generated clustered n=6") {
        SymmetricMatrixSpec spec;
        spec.n                 = 6;
        spec.condition_number  = 100.0;
        spec.spacing           = EigenvalueSpacing::Clustered;
        spec.n_clusters        = 3;
        spec.cluster_width     = 1e-3;
        spec.seed              = 23;
        auto system            = generate_eigen_system<TestType>(spec);
        auto [values, vectors] = mod.run_as<pt>(system.matrix);
        std::cout << "values: " << values << std::endl;
        std::cout << "eigenvalues: " << system.eigenvalues << std::endl;
        require_eigenvalues_approx(values, system.eigenvalues, rtol);
        require_eigenpair_residual(system.matrix, values, vectors, rtol);
    }
}

#ifdef ENABLE_SIGMA
using types2 = std::tuple<tensorwrapper::types::idouble>;
TEMPLATE_LIST_TEST_CASE("QRNormal with noise", "", types2) {
    pluginplay::ModuleManager mm;
    scf::load_modules(mm);

    auto& mod        = mm.at("Eigen Solve via QR");
    auto noise_level = 1e-6;
    using pt         = simde::EigenSolve;
    SECTION("classic 2 by 2") {
        auto system       = classic_2x2<TestType>();
        auto noisy_matrix = add_noise<TestType>(system.matrix, noise_level);
        auto [values, vectors] = mod.run_as<pt>(noisy_matrix);
        require_uq_eigenvalues_contain<TestType>(values, system.eigenvalues,
                                                 noise_level);
        require_uq_eigenpair_residual<TestType>(noisy_matrix, values, vectors,
                                                noise_level);
    }

    SECTION("generated n=4 condition number 1e3") {
        SymmetricMatrixSpec spec;
        spec.n                = 4;
        spec.condition_number = 1e3;
        spec.spacing          = EigenvalueSpacing::Linear;
        spec.seed             = 11;
        auto system           = generate_eigen_system<TestType>(spec);
        auto noisy_matrix     = add_noise<TestType>(system.matrix, noise_level);
        auto [values, vectors] = mod.run_as<pt>(noisy_matrix);
        require_uq_eigenvalues_contain<TestType>(values, system.eigenvalues,
                                                 noise_level);
        require_uq_eigenpair_residual<TestType>(noisy_matrix, values, vectors,
                                                noise_level);
    }

    SECTION("generated clustered n=6") {
        SymmetricMatrixSpec spec;
        spec.n                = 6;
        spec.condition_number = 100.0;
        spec.spacing          = EigenvalueSpacing::Clustered;
        spec.n_clusters       = 3;
        spec.cluster_width    = 1e-8;
        spec.seed             = 23;
        auto system           = generate_eigen_system<TestType>(spec);
        auto noisy_matrix     = add_noise<TestType>(system.matrix, noise_level);
        auto [values, vectors] = mod.run_as<pt>(noisy_matrix);
        require_uq_eigenvalues_contain<TestType>(values, system.eigenvalues,
                                                 noise_level);
        require_uq_eigenpair_residual<TestType>(noisy_matrix, values, vectors,
                                                noise_level);
    }
}
#endif
