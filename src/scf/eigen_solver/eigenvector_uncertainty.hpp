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

#pragma once
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <simde/simde.hpp>
#include <span>
#include <tensorwrapper/tensorwrapper.hpp>
#include <utility>
#include <vector>

namespace scf::eigen_solver {
namespace detail {

/** @brief Applies the first-order eigenvector correction to @p C in place.
 *
 *  Operates on the raw element spans of the eigenvectors @p C (columns,
 *  row-major n by n), the MO-basis Fock @p G = C^T F C (row-major n by n), and
 *  the eigenvalues @p eps (length n). For each column i it adds
 *
 *      dc_i = sum_{j : |e_i - e_j| > 2*sigma} G(j,i) / (e_i - e_j) * c_j ,
 *
 *  where G(j,i) is the off-diagonal coupling v_j^T (delta F) v_i evaluated in
 *  native UQ arithmetic (so the spread shares F's error symbols) and the gap
 *  divisor uses centers only. sigma is the max abs row-sum of the off-diagonal
 *  coupling radii; couplings across gaps <= 2*sigma are dropped
 * (near-degenerate block: gauge freedom + avoids an ill-conditioned ~0
 * divisor).
 *
 *  Runs for every UQ type, intervals included. This is sound here because the
 *  correction is applied ONCE to a converged solution and only reported -- it
 * is not propagated back through the SCF iterations, so the interval dependency
 *  problem cannot compound it.
 */
struct EigenvectorUncertaintyKernel {
    std::size_t m_n;
    const tensorwrapper::buffer::Contiguous& m_C;
    const tensorwrapper::buffer::Contiguous& m_G;
    const tensorwrapper::buffer::Contiguous& m_eps;

    // Dispatches on the (writable) output buffer; the inputs C, G, eps share
    // its runtime type and are pulled with the resolved type. This avoids a
    // cartesian-product instantiation over every buffer's type.
    template<typename FloatType>
    void operator()(std::span<FloatType> out) {
        using clean_t = std::decay_t<FloatType>;
        // The dispatch instantiates const-qualified spans too; the writing body
        // is only valid (and only selected at runtime) for the mutable buffer.
        if constexpr(!std::is_const_v<FloatType> &&
                     tensorwrapper::types::is_uq_type_v<clean_t>) {
            using tensorwrapper::buffer::get_raw_data;
            using tensorwrapper::types::uq_center;
            using tensorwrapper::types::uq_upper;
            using value_t = decltype(uq_center(std::declval<clean_t>()));
            const auto n  = m_n;
            auto radius   = [](const clean_t& x) {
                return uq_upper(x) - uq_center(x);
            };

            auto C   = get_raw_data<clean_t>(m_C);
            auto G   = get_raw_data<clean_t>(m_G);
            auto eps = get_raw_data<clean_t>(m_eps);

            // Degeneracy threshold: max abs row-sum of off-diagonal coupling
            // radii.
            value_t sigma(0);
            for(std::size_t i = 0; i < n; ++i) {
                value_t row(0);
                for(std::size_t j = 0; j < n; ++j) {
                    if(i != j) { row += radius(G[j * n + i]); }
                }
                sigma = std::max(sigma, row);
            }

            // Eigenvalue centers.
            std::vector<value_t> eps_c(n);
            for(std::size_t i = 0; i < n; ++i) { eps_c[i] = uq_center(eps[i]); }

            // out = C + correction, accumulated from the ORIGINAL columns so
            // one column's correction never feeds into another.
            for(std::size_t k = 0; k < n * n; ++k) { out[k] = C[k]; }
            for(std::size_t i = 0; i < n; ++i) {
                for(std::size_t j = 0; j < n; ++j) {
                    if(i == j) { continue; }
                    const value_t gap = eps_c[i] - eps_c[j];
                    if(std::abs(gap) <= value_t(2) * sigma) { continue; }
                    const clean_t coeff = G[j * n + i] / clean_t(gap);
                    for(std::size_t r = 0; r < n; ++r) {
                        out[r * n + i] += coeff * C[r * n + j];
                    }
                }
            }
        }
    }
};

} // namespace detail

/** @brief Attaches first-order MO-coefficient uncertainty to converged
 * orbitals.
 *
 *  Given the converged Fock matrix @p F (carrying the input uncertainty), the
 *  S-orthonormal eigenvectors @p C (columns), and the eigenvalues @p eps, this
 *  replaces each column with its first-order perturbed value
 *  c_i + sum_{j} G(j,i)/(e_i - e_j) c_j, where G = C^T F C. The off-diagonal
 *  G(j,i) = c_j^T (delta F) c_i is the projected SCF residual and is the source
 *  of the coefficient uncertainty.
 *
 *  Roothaan-Hall is the generalized problem F c = e S c, whose rigorous
 *  first-order coupling is c_j^T (delta F - e_i delta S) c_i and which also
 *  carries an S-renormalization term -1/2 (c_i^T delta S c_i) c_i. Here S is
 *  exact (delta S = 0): the uncertainty lives at the Fock level, not in the AO
 *  overlap integrals. Both S-dependent terms therefore vanish and the result
 *  reduces to the ordinary projected-delta-F correction below; S still enters
 *  through the S-orthonormality of C that the C^T F C projection bakes in. If
 *  UQ is ever injected upstream (geometry, basis exponents) so delta S != 0,
 *  the metric terms must be restored.
 *
 *  This is a ONE-SHOT, post-convergence step. Propagating the eigenvector
 * spread through every SCF diagonalization compounds it (the energy's
 * parametric uncertainty is the Hellmann-Feynman explicit term and needs no
 * eigenvector spread); applying the correction once to the converged solution
 * avoids that feedback. Because the result does not re-enter the density, it is
 * applied for ALL uncertainty types, intervals included. A no-op for plain
 * floats.
 *
 *  @param[in,out] C The converged eigenvectors (columns); corrected in place.
 *  @param[in] F The converged, uncertainty-bearing Fock matrix.
 *  @param[in] eps The converged eigenvalues.
 */
inline void attach_eigenvector_uncertainty(simde::type::tensor& C,
                                           const simde::type::tensor& F,
                                           const simde::type::tensor& eps) {
    using tensorwrapper::buffer::make_contiguous;
    using tensorwrapper::buffer::visit_contiguous_buffer;

    // MO-basis Fock G(i,k) = c_i^T F c_k. The diagonal is the eigenvalues; the
    // off-diagonals are the couplings that drive the correction.
    simde::type::tensor CF, G;
    CF("i,k") = C("j,i") * F("j,k");
    G("i,k")  = CF("i,j") * C("j,k");

    // Read-only views of the actual data (single-arg form returns a reference;
    // the 2-arg form would allocate a fresh, default-initialized buffer).
    const auto& c_in = make_contiguous(C.buffer());
    const auto& g_in = make_contiguous(G.buffer());
    const auto& e_in = make_contiguous(eps.buffer());

    const auto n = e_in.shape().extent(0);
    tensorwrapper::shape::Smooth mat_shape{n, n};

    // Writable output buffer of the right shape/type; the kernel fills it.
    auto out_buf = make_contiguous(C.buffer(), mat_shape);

    detail::EigenvectorUncertaintyKernel kernel{n, c_in, g_in, e_in};
    visit_contiguous_buffer(kernel, out_buf);

    C = simde::type::tensor(mat_shape, std::move(out_buf));
}

} // namespace scf::eigen_solver
