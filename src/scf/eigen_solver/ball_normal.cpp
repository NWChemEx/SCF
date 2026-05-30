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
 https://www.texmacs.org/joris/ball/ball.html
 TODO: Write me!!!
 H = TMT - L
 n = number of rows (or columns)
 k = max[max(1/(|L_jj -Lii|)), max(1/L_ii)]
 Omega_n = n by n matrix 0+/-1
 eta = 6 sqrt(n) k ||H||
 T_ball = T(1 + eta Omega_n)
 L_ball = B(1, eta) * L
 )";
} // namespace

MODULE_CTOR(BallNormal) {
    description(desc);
    satisfies_property_type<pt>();

    add_submodule<pt>("Eigen Solve");
}

MODULE_RUN(BallNormal) {
    auto&& [A] = pt::unwrap_inputs(inputs);

    using uq_type = tensorwrapper::types::idouble;
    using namespace tensorwrapper::buffer;
    using namespace tensorwrapper::shape;

    auto A_buf = make_contiguous(A.buffer());
    auto shape = A_buf.shape().make_smooth();
    auto n     = shape.extent(0);

    auto eigen_solver_mod = submods.at("Eigen Solve");
    auto [values, vectors] =
      wrap_subdiagonalization<uq_type>(A_buf, eigen_solver_mod);

    auto k = compute_k<uq_type::value_t>(values);

    auto [uq_values, uq_vectors] = convert_to_uq<uq_type>(values, vectors);

    auto H = compute_residual<uq_type>(A, uq_vectors, uq_values);

    using tensorwrapper::utilities::diagonal_matrix;
    auto H_norm      = ball_matrix_norm<uq_type>(H);
    auto eval_matrix = diagonal_matrix(uq_values);
    auto eval_norm   = ball_matrix_norm<uq_type>(eval_matrix);
    auto eta_val     = 6.0 * std::sqrt(n) * k * H_norm;

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
