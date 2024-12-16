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

#include "../../test_scf.hpp"
#include <scf/scf.hpp>
#include <simde/simde.hpp>

TEST_CASE("EigenGeneralized") {
    pluginplay::ModuleManager mm;
    scf::load_modules(mm);

    using pt = simde::GeneralizedEigenSolve;

    auto& mod = mm.at("Generalized eigensolve via Eigen");

    simde::type::tensor A({{1.0, 2.0}, {2.0, 3.0}});
    simde::type::tensor B({{1.0, 0.0}, {0.0, 1.0}});

    auto&& [values, vector] = mod.run_as<pt>(A, B);

    using value_alloc_t      = tensorwrapper::allocator::Eigen<double, 1>;
    const auto& eigen_values = value_alloc_t::rebind(values.buffer());

    using Catch::Matchers::WithinAbs;
    REQUIRE_THAT(eigen_values.value()(0), WithinAbs(-0.236068, 1E-6));
    REQUIRE_THAT(eigen_values.value()(1), WithinAbs(4.236068, 1E-6));
}