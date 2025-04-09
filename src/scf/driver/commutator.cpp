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

#include "driver.hpp"

simde::type::tensor commutator(simde::type::tensor Fock_Matrix, simde::type::tensor Density_Matrix, simde::type::tensor Overlap_Matrix) {
    
    simde::type::tensor FPS;
    FPS("m,l") = Fock_Matrix("m,n") * Density_Matrix("n,l");
    FPS("m,l") = FPS("m,n") * Overlap_Matrix("n,l");

    simde::type::tensor SPF;
    SPF("m,l") = Density_Matrix("m,n") * Fock_Matrix("n,l");
    SPF("m,l") = Overlap_Matrix("m,n") * SPF("n,l");

    simde::type::tensor grad;
    grad("m,n")   = FPS("m,n") - SPF("m,n");

    return grad;
}
