// /*
//  * Copyright 2026 NWChemEx-Project
//  *
//  * Licensed under the Apache License, Version 2.0 (the "License");
//  * you may not use this file except in compliance with the License.
//  * You may obtain a copy of the License at
//  *
//  * http://www.apache.org/licenses/LICENSE-2.0
//  *
//  * Unless required by applicable law or agreed to in writing, software
//  * distributed under the License is distributed on an "AS IS" BASIS,
//  * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  * See the License for the specific language governing permissions and
//  * limitations under the License.
//  */

// #include "scf/eigen_solver/ball_helpers.hpp"
// #include "test_eigen_solver.hpp"
// #include <scf/scf.hpp>
// #include <simde/simde.hpp>

// #ifdef ENABLE_SIGMA
// #include <sigma/sigma.hpp>

// using namespace scf::eigen_solver;

// using types = std::tuple<tensorwrapper::types::idouble>;
// TEMPLATE_LIST_TEST_CASE("compute_residual", "", types) {
//     using uq_type = TestType;
//     SECTION("2 by 2 test case") {
//         auto system = test_eigen_solver::classic_2x2();
//         auto A      = test_eigen_solver::noisy_matrix<uq_type>(system, 1e-7);
//         auto C      = test_eigen_solver::eigenvectors_as<uq_type>(system);
//         auto L      = test_eigen_solver::eigenvalues_as<uq_type>(system);

//         auto H = compute_residual<uq_type>(A, C, L);

//         auto A_exact  = test_eigen_solver::noisy_matrix<uq_type>(system,
//         0.0); auto H_ref    = compute_residual<uq_type>(A_exact, C, L); auto
//         H_norm   = ball_matrix_norm<uq_type>(H); auto ref_norm =
//         ball_matrix_norm<uq_type>(H_ref); REQUIRE(H_norm ==
//         Catch::Approx(ref_norm).margin(1e-4));
//     }
// }

// TEMPLATE_LIST_TEST_CASE("ball_matrix_norm", "", types) {
//     using uq_type = TestType;
//     SECTION("2 by 2 test case") {
//         auto system = test_eigen_solver::classic_2x2();
//         auto A      = test_eigen_solver::noisy_matrix<uq_type>(system,
//         0.001); auto norm   = ball_matrix_norm<uq_type>(A); REQUIRE(norm ==
//         Catch::Approx(7.07376242552457146).margin(1e-6));
//     }
// }

// TEST_CASE("compute_kappa", "") {
//     using float_type = double;
//     SECTION("2 by 2 test case") {
//         auto system = test_eigen_solver::classic_2x2();
//         auto L      = system.eigenvalues;
//         auto kappa  = compute_kappa<float_type>(L);
//         REQUIRE(kappa == Catch::Approx(4.23606793581038488).margin(1e-10));
//     }
// }

// TEST_CASE("certification_condition_met", "") {
//     using float_type = double;
//     SECTION("small residual passes for classic 2 by 2") {
//         auto system      = test_eigen_solver::classic_2x2();
//         auto kappa       = compute_kappa<float_type>(system.eigenvalues);
//         auto lambda_norm =
//         diagonal_vector_norm<float_type>(system.eigenvalues);
//         REQUIRE(certification_condition_met(1e-12, 2, kappa, lambda_norm));
//     }
// }

// TEST_CASE("l_ball contains its center", "") {
//     using uq_type = tensorwrapper::types::idouble;
//     SECTION("negative eigenvalue") {
//         tensorwrapper::shape::Smooth shape{1};
//         auto buf = tensorwrapper::buffer::make_contiguous<uq_type>(shape);
//         buf.set_elem({0}, uq_type(-0.2360680401325226));
//         simde::type::tensor L(shape, std::move(buf));

//         auto L_ball = l_ball<uq_type>(L, 1e-6);
//         auto out    =
//         tensorwrapper::buffer::make_contiguous(L_ball.buffer()); auto ball =
//         wtf::fp::float_cast<uq_type>(out.get_elem({0}));
//         REQUIRE(ball.contains(-0.2360680401325226));
//     }
// }

// #endif
