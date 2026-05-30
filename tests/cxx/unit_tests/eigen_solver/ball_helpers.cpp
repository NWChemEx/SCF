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

#include "scf/eigen_solver/ball_helpers.hpp"
#include "test_eigen_solver.hpp"
#include <scf/scf.hpp>
#include <simde/simde.hpp>

#ifdef ENABLE_SIGMA
#include <sigma/sigma.hpp>

using namespace scf::eigen_solver;

using types = std::tuple<tensorwrapper::types::idouble>;
TEMPLATE_LIST_TEST_CASE("compute_residual", "", types) {
    using uq_type = TestType;
    SECTION("2 by 2 test case") {
        auto system = test_eigen_solver::classic_2x2();
        auto A      = test_eigen_solver::noisy_matrix<uq_type>(system, 1e-7);
        auto C      = test_eigen_solver::eigenvectors_as<uq_type>(system);
        auto L      = test_eigen_solver::eigenvalues_as<uq_type>(system);

        auto H = compute_residual<uq_type>(A, C, L);

        auto A_exact  = test_eigen_solver::noisy_matrix<uq_type>(system, 0.0);
        auto H_ref    = compute_residual<uq_type>(A_exact, C, L);
        auto H_norm   = ball_matrix_norm<uq_type>(H);
        auto ref_norm = ball_matrix_norm<uq_type>(H_ref);
        REQUIRE(H_norm == Catch::Approx(ref_norm).margin(1e-4));
    }
}

TEMPLATE_LIST_TEST_CASE("ball_matrix_norm", "", types) {
    using uq_type = TestType;
    SECTION("2 by 2 test case") {
        auto system = test_eigen_solver::classic_2x2();
        auto A      = test_eigen_solver::noisy_matrix<uq_type>(system, 0.001);
        auto norm   = ball_matrix_norm<uq_type>(A);
        REQUIRE(norm == Catch::Approx(4.24452635755745966).margin(5e-3));
    }
}

TEST_CASE("compute_k", "") {
    using float_type = double;
    SECTION("2 by 2 test case") {
        auto system = test_eigen_solver::classic_2x2();
        auto L      = system.eigenvalues;
        auto k      = compute_k<float_type>(L);
        REQUIRE(k == Catch::Approx(0.23606793581038488).margin(1e-10));
    }
}

#endif
