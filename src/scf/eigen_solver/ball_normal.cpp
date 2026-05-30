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

#include "ball_helpers.hpp"
#include "eigen_solver.hpp"

#include <simde/simde.hpp>
#include <tensorwrapper/tensorwrapper.hpp>
#ifdef ENABLE_SIGMA
#include <sigma/sigma.hpp>

using pt = simde::EigenSolve;

namespace scf::eigen_solver {
namespace {

const auto desc = R"(
 Eigen Solve via Ball Arithmetic
 -------------------------------------------
 https://www.texmacs.org/joris/ball/ball.html section 6.3

 H = T^T M T - Λ
 κ = max[max(1/|Λ_jj - Λ_ii|), max(1/|Λ_ii|)]
 certification when ||H|| <= 1 / (80 n κ^2 ||Λ||)
 η = 6 n sqrt(κ) ||H||
 T_ball = T (1 + η Ω_n)
 L_ball = B(1, η) Λ
 )";

constexpr int kMaxNewtonIter = 50;
constexpr double kNewtonTol  = 1e-14;

} // namespace

MODULE_CTOR(BallNormal) {
    description(desc);
    satisfies_property_type<pt>();

    add_submodule<pt>("Eigen Solve");
}

MODULE_RUN(BallNormal) {
    auto&& [A] = pt::unwrap_inputs(inputs);

    using uq_type    = tensorwrapper::types::idouble;
    using float_type = uq_type::value_t;
    using namespace tensorwrapper::buffer;
    using namespace tensorwrapper::shape;

    auto A_buf    = make_contiguous(A.buffer());
    auto shape    = A_buf.shape().make_smooth();
    auto n        = shape.extent(0);
    auto A_median = median_matrix<uq_type>(A_buf);

    auto eigen_solver_mod = submods.at("Eigen Solve");
    auto [values, vectors] =
      wrap_subdiagonalization<uq_type>(A_buf, eigen_solver_mod);

    auto clusters = finest_cluster(n);
    auto T        = vectors;

    for(int it = 0; it < kMaxNewtonIter; ++it) {
        auto [T_new, lambda] =
          fundamental_iteration_step(T, A_median, clusters);

        tensorwrapper::Tensor TA, M;
        TA("i,k") = T("j,i") * A_median("j,k");
        M("i,k")  = TA("i,j") * T("j,k");

        auto off_norm = off_block_max_norm_double(M, clusters);
        T             = T_new;
        values        = lambda;
        if(off_norm < kNewtonTol) break;
    }

    auto [uq_values, uq_vectors] = convert_to_uq<uq_type>(values, T);

    auto H           = compute_residual<uq_type>(A, uq_vectors, uq_values);
    auto h_norm      = ball_matrix_norm<uq_type>(H);
    auto kappa       = compute_kappa<float_type>(values);
    auto lambda_norm = diagonal_vector_norm<float_type>(values);

    if(!certification_condition_met(h_norm, n, kappa, lambda_norm)) {
        throw std::runtime_error(
          "BallNormal: certification condition (32) not satisfied");
    }

    auto eta_val = compute_eta(h_norm, n, kappa);

    auto C_ball = t_ball<uq_type>(uq_vectors, eta_val);
    auto L_ball = l_ball<uq_type>(uq_values, eta_val);

    auto rv = results();
    return pt::wrap_results(rv, L_ball, C_ball);
}
} // namespace scf::eigen_solver
#else
namespace scf::eigen_solver {
using pt = simde::EigenSolve;
MODULE_CTOR(BallNormal) {
    description("Sigma was not enabled.");
    satisfies_property_type<pt>();

    add_submodule<pt>("Eigen Solve");
}

MODULE_RUN(BallNormal) { throw std::runtime_error("Sigma was not enabled."); }
} // namespace scf::eigen_solver
#endif
