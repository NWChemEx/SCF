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
    SECTION("Commutator") {
        simde::type::tensor Fock_Matrix{{1.0, 2.0}, {3.0, 4.0}};
        simde::type::tensor Density_Matrix{{2.0, 3.0}, {4.0, 5.0}};
        simde::type::tensor Overlap_Matrix{{3.0, 4.0}, {5.0, 6.0}};

        simde::type::tensor test_grad{{-14, -42}, {42, 14}};
        auto grad =
          scf::driver::commutator(Fock_Matrix, Density_Matrix, Overlap_Matrix);

        REQUIRE(grad == test_grad);
    }

    SECTION("Empty Tensor") {
        simde::type::tensor Fock_Matrix{{1.0, 2.0}, {3.0, 4.0}};
        simde::type::tensor Density_Matrix{{1.0, 2.0, 3.0}, {4.0, 5.0}};
        simde::type::tensor Overlap_Matrix{{3.0, 4.0}, {5.0, 6.0}};

        auto grad =
          scf::driver::commutator(Fock_Matrix, Density_Matrix, Overlap_Matrix);

        REQUIRE_THROWS_AS("Not smooth", std::runtime_error);
    }
}
