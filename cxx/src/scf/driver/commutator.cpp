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

#include <scf/driver/commutator.hpp>

namespace scf::driver {

simde::type::tensor commutator(simde::type::tensor A, simde::type::tensor B,
                               simde::type::tensor S) {
    if(A.rank() != 2 || B.rank() != 2 || S.rank() != 2) {
        std::stringstream ss;
        ss << "Matrix 1 Rank: " << A.rank() << "\nMatrix 2 Rank: " << B.rank()
           << "\nMatrix 3 Rank: " << S.rank() << std::endl;
        throw std::runtime_error("Input matrix rank not 2!\n" + ss.str());
    }
    simde::type::tensor AB, BA, ABS, SBA;
    AB("m,l")  = A("m,n") * B("n,l");
    ABS("m,l") = AB("m,n") * S("n,l");

    BA("m,l")  = B("m,n") * A("n,l");
    SBA("m,l") = S("m,n") * BA("n,l");

    simde::type::tensor grad;
    grad("m,n") = ABS("m,n") - SBA("m,n");

    return grad;
}
} // namespace scf::driver
