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

TEST_CASE("Commutator Function") {
    SECTION("Commutator Test") {
        simde::type::tensor A{{1.0, 2.0}, {3.0, 4.0}};
        simde::type::tensor B{{2.0, 3.0}, {4.0, 5.0}};
        simde::type::tensor S{{3.0, 4.0}, {5.0, 6.0}};

        simde::type::tensor test_grad{{-14, -42}, {42, 14}};
        auto grad = scf::driver::commutator(A, B, S);

        REQUIRE(grad == test_grad);
    }

    SECTION("Input not Matrix Rank") {
        simde::type::tensor A{{1.0, 2.0}, {3.0, 4.0}};
        simde::type::tensor B{4.0, 5.0};
        simde::type::tensor S{{3.0, 4.0}, {5.0, 6.0}};
        REQUIRE_THROWS(scf::driver::commutator(A, B, S));
    }
}
