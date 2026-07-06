/*
 * Copyright 2025 NWChemEx-Project
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

#pragma once
#include <simde/simde.hpp>

namespace scf::driver {

/** @brief Evaluates the commutator of matrix A, matrix B, with the Overlap S.
 *
 *  The commutator function evaluates the commutator between matrix A, matrix B,
 *  and the Overlap matrix S. This requires that all of the matrices are square.
 *
 *
 *  @param[in] A: matrix A, a tensor of rank 2.
 *  @param[in] B: matrix B, a tensor of rank 2.
 *  @param[in] S: Overlap matrix S, a tensor of rank 2.
 *
 *  @return simde::type::tensor
 */

simde::type::tensor commutator(simde::type::tensor A, simde::type::tensor B,
                               simde::type::tensor S);
} // namespace scf::driver
