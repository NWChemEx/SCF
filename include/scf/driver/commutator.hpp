#pragma once
#include <simde/simde.hpp>

namespace scf::driver {

/** @brief Evaluates the commutator of matrix A, matrix B, with the Overlap S.
 *
 *  The commutator function evaluates the commutator between matrix A, matrix B,
 *  and the Overlap matrix S. This requires that all of the matrices are square.
 *
 *  @tparam simde::type::tensor
 *
 *  @param[in] A: matrix A, a tensor of rank 2.
 *  @param[in] B: matrix B, a tensor of rank 2.
 *  @param[in] S: matrix S, a tensor of rank 2. The overlap matrix of A and B.
 *
 *  @return simde::type::tensor
 */

simde::type::tensor commutator(simde::type::tensor A, simde::type::tensor B,
                               simde::type::tensor S);
} // namespace scf::driver
