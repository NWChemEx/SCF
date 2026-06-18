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
#include <Eigen/Jacobi>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

template<typename T>
using dynamic_matrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

namespace scf::eigen_solver::detail {

template<typename T, bool = tensorwrapper::types::is_uq_type_v<T>>
struct value_type_impl {
    using type = T;
};

template<typename T>
struct value_type_impl<T, true> {
    using type = typename T::value_t;
};

template<typename T>
using value_type = typename value_type_impl<T>::type;

// Computes the Frobenius norm of the off-diagonal elements of n by n matrix S.
template<typename T>
T off_diagonal_frobenius(const dynamic_matrix<T>& S, std::size_t n) {
    using tensorwrapper::types::uq_center;
    value_type<T> frob(0);
    for(std::size_t i = 0; i < n; ++i) {
        for(std::size_t j = 0; j < n; ++j) {
            if(i != j) {
                auto x = uq_center(S(i, j));
                frob += x * x;
            }
        }
    }
    return T(std::sqrt(frob));
}

/** @brief Computes the cosine and sine of a Jacobi rotation for the p, q plane.
 *
 *  @param[in] s_pp_uq The (p, p) element of the matrix to be diagonalized.
 *  @param[in] s_pq_uq The (p, q) element of the matrix to be diagonalized.
 *  @param[in] s_qq_uq The (q, q) element of the matrix to be diagonalized.
 *  @param[out] c_out Output parameter for the cosine of the rotation.
 *  @param[out] s_out Output parameter for the sine of the rotation.
 *
 *  @return True if the rotation was successfully computed, false if the off-
 *          diagonal element is too small.
 */
template<typename T>
bool make_jacobi(const T& s_pp_uq, const T& s_pq_uq, const T& s_qq_uq, T& c_out,
                 T& s_out) {
    using value_t = value_type<T>;
    using tensorwrapper::types::uq_center;

    const value_t s_pp = uq_center(s_pp_uq);
    const value_t s_pq = uq_center(s_pq_uq);
    const value_t s_qq = uq_center(s_qq_uq);

    const value_t deno = value_t(2) * std::abs(s_pq);
    // If the off-diagonal element is already almost zero, the angle will be
    // zero, so we can skip the rotation by just returning cos(0) = 1 and
    // sin(0)=0.
    if(deno < std::numeric_limits<value_t>::min()) {
        c_out = T(1);
        s_out = T(0);
        return false;
    }

    // tau = cot(2*theta)
    const value_t tau = (s_pp - s_qq) / deno;

    // w = |csc(2*theta)| = sqrt(tau^2 + 1)
    const value_t w = std::sqrt(tau * tau + value_t(1));

    // t = tan(theta)
    value_t t;
    if(tau > value_t(0)) {
        t = value_t(1) / (tau + w);
    } else {
        t = value_t(1) / (tau - w);
    }
    const value_t sign_t = t > value_t(0) ? value_t(1) : value_t(-1);
    const value_t n      = value_t(1) / std::sqrt(t * t + value_t(1));
    const value_t y_sign = s_pq >= value_t(0) ? value_t(1) : value_t(-1);
    const value_t s      = -sign_t * y_sign * std::abs(t) * n;
    const value_t c      = n;
    c_out                = T(c);
    s_out                = T(s);
    return true;
}

/** @brief Applies the adjoint of the Jacobi rotation to the left of a matrix.
 *
 *  @param[in,out] mat The matrix to which the rotation is applied.
 *  @param[in] p The first index of the rotation plane.
 *  @param[in] q The second index of the rotation plane.
 *  @param[in] c The cosine of the rotation.
 *  @param[in] s The sine of the rotation.
 */
template<typename T>
void apply_on_the_left_adjoint(dynamic_matrix<T>& mat, Eigen::Index p,
                               Eigen::Index q, const T& c, const T& s) {
    const Eigen::Index n_cols = mat.cols();
    for(Eigen::Index j = 0; j < n_cols; ++j) {
        const T xp    = mat(p, j);
        const T xq    = mat(q, j);
        const T new_p = c * xp - s * xq;
        const T new_q = s * xp + c * xq;
        mat(p, j)     = new_p;
        mat(q, j)     = new_q;
    }
}

/** @brief Applies the Jacobi rotation to the right of a matrix.
 *
 *  @param[in,out] mat The matrix to which the rotation is applied.
 *  @param[in] p The first index of the rotation plane.
 *  @param[in] q The second index of the rotation plane.
 *  @param[in] c The cosine of the rotation.
 *  @param[in] s The sine of the rotation.
 */
template<typename T>
void apply_on_the_right(dynamic_matrix<T>& mat, Eigen::Index p, Eigen::Index q,
                        const T& c, const T& s) {
    const Eigen::Index n_rows = mat.rows();
    for(Eigen::Index i = 0; i < n_rows; ++i) {
        const T xp    = mat(i, p);
        const T xq    = mat(i, q);
        const T new_p = c * xp - s * xq;
        const T new_q = s * xp + c * xq;
        mat(i, p)     = new_p;
        mat(i, q)     = new_q;
    }
}

// Diagaonalizes a matrix via V^T A V, extracts diagonal as a vector.
template<typename T>
std::vector<T> rayleigh_diagonal(const dynamic_matrix<T>& V,
                                 const dynamic_matrix<T>& A_orig,
                                 std::size_t n) {
    std::vector<T> diag(n);
    const dynamic_matrix<T> rayleigh = V.transpose() * A_orig * V;
    for(std::size_t j = 0; j < n; ++j) { diag[j] = rayleigh(j, j); }
    return diag;
}

/** @brief Performs one sweep of Jacobi rotations.
 *
 *  A sweep consists of applying a Jacobi rotation to every pair of indices
 *  (p, q).
 *
 *  @param[in,out] S The matrix to be diagonalized. Gets modified in-place.
 *  @param[in,out] V The matrix that accumulates the rotations. Gets modified
 *                   in-place.
 *  @param[in] n The size of the matrix S (assumed to be n by n).
 */
template<typename T>
void jacobi_sweep(dynamic_matrix<T>& S, dynamic_matrix<T>& V, std::size_t n) {
    for(std::size_t p = 0; p < n; ++p) {
        for(std::size_t q = p + 1; q < n; ++q) {
            T c;
            T s;
            if(make_jacobi(S(p, p), S(p, q), S(q, q), c, s)) {
                apply_on_the_left_adjoint(S, p, q, c, s);
                apply_on_the_right(S, p, q, c, s);
                apply_on_the_right(V, p, q, c, s);
                S(p, q) = T(0);
                S(q, p) = T(0);
            }
        }
    }
}

template<typename T>
inline std::pair<std::vector<T>, std::vector<T>> symmetric_jacobi_eigen(
  std::span<const T> A, std::size_t n, double tol, std::size_t max_sweeps) {
    using matrix_type = dynamic_matrix<T>;

    // Copy of A, used to get eigenvalues via V^T A V at the end. We can't just
    // use S because the Jacobi rotations are applied in-place to S, so S
    // doesn't equal V^T A V until convergence.
    matrix_type A_orig(n, n);

    // Iterative approximation to the eigenvalues of A. Converges to a diagonal
    // matrix.
    matrix_type S(n, n);
    for(std::size_t i = 0; i < n; ++i) {
        for(std::size_t j = 0; j < n; ++j) {
            A_orig(i, j) = A[i * n + j];
            S(i, j)      = A[i * n + j];
        }
    }

    // Accumulates the rotations. Converges to the eigenvectors of A.
    matrix_type V = matrix_type::Identity(n, n);
    using tensorwrapper::types::strictly_less;
    using tensorwrapper::types::uq_center;
    for(std::size_t sweep = 0; sweep < max_sweeps; ++sweep) {
        const auto frob_sqrt = off_diagonal_frobenius(S, n);
        if(strictly_less(frob_sqrt, tol)) { break; }
        if(sweep == max_sweeps - 1) {
            throw std::runtime_error("Jacobi algorithm did not converge");
        }

        jacobi_sweep(S, V, n);
        // Symmetrize S to prevent numerical issues from destroying symmetry.
        S = (S + S.transpose()) / T(2.0);
    }

    const auto diag = rayleigh_diagonal(V, A_orig, n);

    // Order eigenvalue in increasing order
    std::vector<std::size_t> order(n);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](std::size_t i, std::size_t j) {
        return uq_center(diag[i]) < uq_center(diag[j]);
    });

    // Fill in return values/vectors
    std::vector<T> values(n);
    std::vector<T> vectors(n * n);
    for(std::size_t k = 0; k < n; ++k) {
        const auto idx = order[k];
        values[k]      = diag[idx];
        for(std::size_t i = 0; i < n; ++i) { vectors[i * n + k] = V(i, idx); }
    }
    return {std::move(values), std::move(vectors)};
}

} // namespace scf::eigen_solver::detail
