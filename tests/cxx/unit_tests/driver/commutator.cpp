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
 simde::type::tensor Fock_Matrix{{1.0, 2.0}, {3.0, 4.0}};
 simde::type::tensor Density_Matrix{{2.0, 3.0}, {4.0, 5.0}};
 simde::type::tensor Overlap_Matrix{{3.0, 4.0}, {5.0, 6.0}};

 simde::type::tensor FP, PF, FPS, SPF, test_grad;
 FP("m, l") = Fock_Matrix("m,n") * Density_Matrix("n,l");
 FPS("m, l") = FP("m,n") * Overlap_Matrix("n,l");
 PF("m,l") = Density_Matrix("m,n") * Fock_Matrix("n,l");
 SPF("m,l") = Overlap_Matrix("m,n") * PF("n,l");
 test_grad("m,n") = FPS("m,n") - SPF("m,n");

 auto grad = commutator(Fock_Matrix, Density_Matrix, Overlap_Matrix);

 REQUIRE(grad == test_grad);
}
