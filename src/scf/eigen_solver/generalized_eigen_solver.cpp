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

#include "eigen_solver.hpp"
#include <wtf/fp/float_view.hpp>

#include <simde/simde.hpp>
#include <tensorwrapper/tensorwrapper.hpp>

using pt        = simde::GeneralizedEigenSolve;
using pt_normal = simde::EigenSolve;
namespace scf::eigen_solver {
namespace {
const auto desc = R"(
Generalized Eigen Solve
-----------------------

TODO: Write me!!!
)";
}

MODULE_CTOR(GeneralizedEigenSolver) {
    description(desc);
    satisfies_property_type<pt>();

    add_submodule<pt_normal>("Eigen Solve");
}

MODULE_RUN(GeneralizedEigenSolver) {
    auto&& [A, B] = pt::unwrap_inputs(inputs);

    auto eigen_solver_mod = submods.at("Eigen Solve");

    // Step 1: Diagonalize B to get B_values and B_vectors
    auto [B_values, B_vectors] = eigen_solver_mod.run_as<pt_normal>(B);

    // Step 2: Compute (B_values)**-1/2
    using tensorwrapper::operations::power;
    using tensorwrapper::utilities::diagonal_matrix;
    auto B_values_inv_sqrt = power(B_values, -0.5);
    auto B_values_matrix   = diagonal_matrix(B_values_inv_sqrt);

    // Step 3: Compute C = B_vectors * (B_values)**-1/2
    simde::type::tensor C;
    C("i,k") = B_vectors("i,j") * B_values_matrix("j,k");

    // Step 4: A' = C^T * A * C
    simde::type::tensor CA, A_prime;
    CA("i,k")      = C("j,i") * A("j,k");
    A_prime("i,k") = CA("i,j") * C("j,k");

    // Step 5: Diagonalize A' to get A_values and A'_vectors
    auto [A_values, A_vectors] = eigen_solver_mod.run_as<pt_normal>(A_prime);

    // Step 6: A_vectors = C * A'_vectors
    simde::type::tensor evectors;
    evectors("i,k") = C("i,j") * A_vectors("j,k");

    auto rv = results();
    return pt::wrap_results(rv, A_values, evectors);
}
} // namespace scf::eigen_solver
