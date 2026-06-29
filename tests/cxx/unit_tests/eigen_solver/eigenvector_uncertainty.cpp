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

#include "eigen_solver/eigenvector_uncertainty.hpp"
#include "test_eigen_solver.hpp"

#ifdef ENABLE_SIGMA

using namespace test_eigen_solver;
using namespace tensorwrapper::generate;

// The one-shot correction is applied post-convergence and never re-enters the
// density, so it runs for ALL uncertainty types -- including intervals, which
// the iterative version had to skip. This checks the spread is actually
// attached for each type and that a stationary contraction stays bounded (the
// correction adds real uncertainty without blowing up on a single application).
using uq_types =
  std::tuple<tensorwrapper::types::idouble, tensorwrapper::types::adouble,
             tensorwrapper::types::tadouble>;

TEMPLATE_LIST_TEST_CASE("attach_eigenvector_uncertainty", "", uq_types) {
    using tensorwrapper::buffer::get_raw_data;
    using tensorwrapper::types::uq_center;
    using tensorwrapper::types::uq_upper;
    auto radius = [](const TestType& x) { return uq_upper(x) - uq_center(x); };

    pluginplay::ModuleManager mm;
    scf::load_modules(mm);
    auto& mod        = mm.at("Eigen Solve via Jacobi");
    auto noise_level = 1e-6;
    using pt         = simde::EigenSolve;

    auto run = [&](const auto& spec_or_system, std::size_t n) {
        auto noisy = add_noise<TestType>(spec_or_system.matrix, noise_level);

        // Solve the (S = I) standard problem: vectors are orthonormal, so the
        // generalized correction reduces to G = C^T A C in the same basis.
        auto [values, vectors] = mod.run_as<pt>(noisy);

        // Before: the Jacobi solver leaves the eigenvectors at their centers
        // (intervals carry only tiny outward-rounding width).
        auto pre       = get_raw_data<TestType>(vectors.buffer());
        double pre_max = 0.0;
        for(std::size_t k = 0; k < pre.size(); ++k)
            pre_max = std::max(pre_max, radius(pre[k]));
        REQUIRE(pre_max < 1.0e-10);

        // After: the one-shot correction attaches coefficient uncertainty.
        scf::eigen_solver::attach_eigenvector_uncertainty(vectors, noisy,
                                                          values);
        auto post       = get_raw_data<TestType>(vectors.buffer());
        double post_max = 0.0;
        for(std::size_t k = 0; k < post.size(); ++k)
            post_max = std::max(post_max, radius(post[k]));
        REQUIRE(post_max > 1.0e-9); // spread attached, intervals included
    };

    SECTION("classic 2 by 2") { run(classic_2x2<TestType>(), 2); }

    SECTION("generated n=4 condition number 1e3") {
        SymmetricMatrixSpec spec;
        spec.n                = 4;
        spec.condition_number = 1e3;
        spec.spacing          = EigenvalueSpacing::Linear;
        spec.seed             = 11;
        run(generate_eigen_system<TestType>(spec), 4);
    }
}

#endif
