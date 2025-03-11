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

TEMPLATE_LIST_TEST_CASE("EigenGeneralized", "", test_scf::float_types) {
    using float_type = TestType;
    pluginplay::ModuleManager mm;
    scf::load_modules(mm);

    using pt = simde::GeneralizedEigenSolve;

    auto& mod = mm.at("Generalized eigensolve via Eigen");

    tensorwrapper::allocator::Eigen<float_type> alloc(mm.get_runtime());
    tensorwrapper::shape::Smooth shape{2, 2};
    tensorwrapper::layout::Physical l(shape);
    auto A_buffer      = alloc.allocate(l);
    A_buffer->at(0, 0) = 1.0;
    A_buffer->at(0, 1) = 2.0;
    A_buffer->at(1, 0) = 2.0;
    A_buffer->at(1, 1) = 3.0;

    auto B_buffer      = alloc.allocate(l);
    B_buffer->at(0, 0) = 1.0;
    B_buffer->at(0, 1) = 0.0;
    B_buffer->at(1, 0) = 0.0;
    B_buffer->at(1, 1) = 1.0;

    simde::type::tensor A(shape, std::move(A_buffer));
    simde::type::tensor B(shape, std::move(B_buffer));

    auto&& [values, vector] = mod.run_as<pt>(A, B);

    auto corr_buffer = alloc.construct({-0.236068, 4.236068});
    tensorwrapper::shape::Smooth corr_shape{2};
    simde::type::tensor corr(corr_shape, std::move(corr_buffer));

    REQUIRE(tensorwrapper::operations::approximately_equal(corr, values, 1E-6));
}