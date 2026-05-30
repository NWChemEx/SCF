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

using types = std::tuple<float, double>;
using namespace test_eigen_solver;

TEMPLATE_LIST_TEST_CASE("EigenSolveDriver", "", types) {
    pluginplay::ModuleManager mm;
    scf::load_modules(mm);
    auto& mod = mm.at("Eigen Solve");
    using pt  = simde::EigenSolve;

    auto rtol              = std::is_same_v<TestType, float> ? 5e-4 : 1e-5;
    auto system            = classic_2x2<TestType>();
    auto [values, vectors] = mod.run_as<pt>(system.matrix);
    require_eigenvalues_approx(values, system.eigenvalues, rtol);
    require_eigenpair_residual(system.matrix, values, vectors, rtol);
}
