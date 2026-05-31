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
#include <Eigen/Eigen>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <span>
#include <tensorwrapper/types/floating_point.hpp>
#include <utility>
#include <vector>
#ifdef ENABLE_SIGMA
#include <sigma/operations/exponents.hpp>
#endif

namespace sigma {

template<typename T1, typename T2>
bool operator<(const Interval<T1>& lhs, const Interval<T2>& rhs) {
    return lhs.median() < rhs.median();
}

template<typename T1, typename T2>
bool operator>(const Interval<T1>& lhs, const Interval<T2>& rhs) {
    return rhs < lhs;
}

template<typename T1, typename T2>
bool operator<=(const Interval<T1>& lhs, const Interval<T2>& rhs) {
    return (lhs == rhs) || (lhs < rhs);
}

template<typename T1, typename T2>
bool operator>=(const Interval<T1>& lhs, const Interval<T2>& rhs) {
    return rhs <= lhs;
}

} // namespace sigma

namespace scf::eigen_solver::detail {

namespace fp {

template<typename T>
T zero() {
    return T(0);
}

template<typename T>
T one() {
    return T(1);
}

template<typename T>
T two() {
    return T(2);
}

template<typename T>
T abs2(const T& x) {
    return x * x;
}

template<typename T>
T sqrt(const T& x) {
    if constexpr(tensorwrapper::types::is_interval_v<T>) {
        return sigma::sqrt(x);
    } else if constexpr(tensorwrapper::types::is_uncertain_v<T>) {
        return sigma::sqrt(x);
    } else {
        return std::sqrt(x);
    }
}

template<typename T>
T min_magnitude() {
    if constexpr(tensorwrapper::types::is_interval_v<T>) {
        using value_t = typename T::value_t;
        return T(std::numeric_limits<value_t>::min());
    } else if constexpr(tensorwrapper::types::is_uncertain_v<T>) {
        using value_t = typename T::value_t;
        return T(std::numeric_limits<value_t>::min());
    } else {
        return std::numeric_limits<T>::min();
    }
}

template<typename T>
bool is_nearly_zero(const T& x) {
    if constexpr(tensorwrapper::types::is_interval_v<T>) {
        using value_t = typename T::value_t;
        return x.lower() >= value_t(0) && x.upper() <= value_t(0);
    } else {
        return x == zero<T>();
    }
}

template<typename T>
bool is_below_or_equal(const T& lhs, const T& rhs) {
    if constexpr(tensorwrapper::types::is_interval_v<T>) {
        return lhs.upper() <= rhs.lower();
    } else {
        return lhs <= rhs;
    }
}

template<typename T>
bool is_above_or_equal(const T& lhs, const T& rhs) {
    if constexpr(tensorwrapper::types::is_interval_v<T>) {
        return lhs.lower() >= rhs.upper();
    } else {
        return lhs >= rhs;
    }
}

/// Stable Householder pivot u0 = x0 +/- ||x|| avoiding x0 - beta cancellation.
template<typename T>
T stable_householder_pivot(const T& x0, const T& nrm) {
    if constexpr(tensorwrapper::types::is_interval_v<T>) {
        if(is_below_or_equal(x0, zero<T>())) { return x0 - nrm; }
        if(is_above_or_equal(x0, zero<T>())) { return x0 + nrm; }
        const T minus_branch = x0 - nrm;
        const T plus_branch  = x0 + nrm;
        return T(std::min(minus_branch.lower(), plus_branch.lower()),
                 std::max(minus_branch.upper(), plus_branch.upper()));
    } else {
        if(x0 <= zero<T>()) { return x0 - nrm; }
        return x0 + nrm;
    }
}

template<typename T>
T reflector_beta(const T& x0, const T& nrm) {
    if constexpr(tensorwrapper::types::is_interval_v<T>) {
        if(is_below_or_equal(x0, zero<T>())) { return nrm; }
        if(is_above_or_equal(x0, zero<T>())) { return -nrm; }
        return T(-nrm.upper(), nrm.upper());
    } else if constexpr(tensorwrapper::types::is_uncertain_v<T>) {
        return (x0 <= zero<T>()) ? nrm : -nrm;
    } else {
        return (x0 <= zero<T>()) ? nrm : -nrm;
    }
}

template<typename T>
T as_nonnegative(const T& x) {
    if constexpr(tensorwrapper::types::is_interval_v<T>) {
        using value_t = typename T::value_t;
        return T(std::max(x.lower(), value_t(0)),
                 std::max(x.upper(), value_t(0)));
    } else {
        return x;
    }
}

template<typename T>
bool is_at_most_tolerance(const T& value, double tol) {
    if constexpr(tensorwrapper::types::is_interval_v<T>) {
        return value.upper() <= typename T::value_t(tol);
    } else {
        return value <= T(tol);
    }
}

/// Builds a reflector for interval types using median values to avoid
/// dependency-induced pivots that straddle zero, then stores the result as
/// point intervals for use with interval matrix arithmetic in the apply step.
template<typename Derived>
void make_householder_interval(Eigen::MatrixBase<Derived>& x,
                               typename Derived::Scalar& tau,
                               typename Derived::Scalar& beta) {
    using T              = typename Derived::Scalar;
    using value_t        = typename T::value_t;
    const Eigen::Index m = x.size();

    const value_t c0 = x(0).median();
    value_t tail_sq  = value_t(0);
    for(Eigen::Index i = 1; i < m; ++i) {
        const value_t xi = x(i).median();
        tail_sq += xi * xi;
    }

    const value_t tol = std::numeric_limits<value_t>::min();
    if(tail_sq <= tol && c0 * c0 <= tol) {
        tau  = T(0);
        beta = x(0);
        x.tail(m - 1).setZero();
        return;
    }

    const value_t nrm = std::sqrt(c0 * c0 + tail_sq);
    value_t u0        = (c0 <= value_t(0)) ? c0 - nrm : c0 + nrm;
    const value_t u0_tol =
      std::numeric_limits<value_t>::epsilon() * (value_t(1) + nrm);
    if(std::abs(u0) < u0_tol) { u0 = (c0 <= value_t(0)) ? -nrm : nrm; }
    beta = T((c0 <= value_t(0)) ? nrm : -nrm);

    value_t denom = value_t(1);
    for(Eigen::Index i = 1; i < m; ++i) {
        const value_t essential_i = x(i).median() / u0;
        x(i)                      = T(essential_i);
        denom += essential_i * essential_i;
    }
    tau = T(value_t(2) / denom);
}

} // namespace fp

/// Computes a Householder reflector in Eigen's compact storage format.
///
/// On input, @p x holds the vector to reflect. On output, @p x(0) is @p beta,
/// @p x.tail(m-1) holds the essential part of the reflector, and @p tau is the
/// scaling factor for H = I - tau v v^T with v = [1, essential].
template<typename Derived>
void make_householder(Eigen::MatrixBase<Derived>& x,
                      typename Derived::Scalar& tau,
                      typename Derived::Scalar& beta) {
    using T              = typename Derived::Scalar;
    const Eigen::Index m = x.size();
    if(m == 1) {
        tau  = fp::zero<T>();
        beta = x(0);
        return;
    }

    if constexpr(tensorwrapper::types::is_interval_v<T>) {
        fp::make_householder_interval(x, tau, beta);
        return;
    }

    const T c0 = x(0);
    T tail_sq  = fp::zero<T>();
    for(Eigen::Index i = 1; i < m; ++i) { tail_sq += fp::abs2(x(i)); }

    const T tol = fp::min_magnitude<T>();
    if(fp::is_below_or_equal(tail_sq, tol) &&
       fp::is_below_or_equal(fp::abs2(c0), tol)) {
        tau  = fp::zero<T>();
        beta = c0;
        x.tail(m - 1).setZero();
        return;
    }

    const T nrm = fp::sqrt(fp::as_nonnegative(fp::abs2(c0) + tail_sq));
    const T u0  = fp::stable_householder_pivot(c0, nrm);
    beta        = fp::reflector_beta(c0, nrm);

    x.tail(m - 1) = x.tail(m - 1) / u0;

    T denom = fp::one<T>();
    for(Eigen::Index i = 1; i < m; ++i) { denom += fp::abs2(x(i)); }
    tau = fp::two<T>() / denom;
}

template<typename Derived, typename EssentialPart>
void apply_householder_on_the_right(Eigen::MatrixBase<Derived>& mat,
                                    const EssentialPart& essential,
                                    const typename Derived::Scalar& tau,
                                    typename Derived::Scalar* workspace) {
    using T = typename Derived::Scalar;
    if(mat.cols() == 1) {
        mat *= fp::one<T>() - tau;
        return;
    }
    if(fp::is_nearly_zero(tau)) { return; }

    if constexpr(tensorwrapper::types::is_interval_v<T>) {
        using value_t             = typename T::value_t;
        const value_t tau_m       = tau.median();
        const Eigen::Index n_rows = mat.rows();
        const Eigen::Index n_cols = mat.cols() - 1;

        std::vector<value_t> tmp(n_rows);
        for(Eigen::Index i = 0; i < n_rows; ++i) {
            value_t sum = mat(i, 0).median();
            for(Eigen::Index j = 0; j < n_cols; ++j) {
                sum += mat(i, j + 1).median() * essential(j).median();
            }
            tmp[static_cast<std::size_t>(i)] = sum;
        }

        for(Eigen::Index i = 0; i < n_rows; ++i) {
            const value_t new_v =
              mat(i, 0).median() - tau_m * tmp[static_cast<std::size_t>(i)];
            mat(i, 0) = T(new_v);
        }
        for(Eigen::Index j = 0; j < n_cols; ++j) {
            const value_t e_j = essential(j).median();
            for(Eigen::Index i = 0; i < n_rows; ++i) {
                const value_t new_v =
                  mat(i, j + 1).median() -
                  tau_m * tmp[static_cast<std::size_t>(i)] * e_j;
                mat(i, j + 1) = T(new_v);
            }
        }
        return;
    }

    using ColType = Eigen::Matrix<T, Eigen::Dynamic, 1>;
    Eigen::Map<ColType> tmp(workspace, mat.rows());
    auto right    = mat.rightCols(mat.cols() - 1);
    tmp.noalias() = right * essential;
    tmp += mat.col(0);
    mat.col(0) -= tau * tmp;
    right.noalias() -= tau * tmp * essential.adjoint();
}

template<typename Derived, typename EssentialPart>
void apply_householder_on_the_left(Eigen::MatrixBase<Derived>& mat,
                                   const EssentialPart& essential,
                                   const typename Derived::Scalar& tau,
                                   typename Derived::Scalar* workspace) {
    using T = typename Derived::Scalar;
    if(mat.rows() == 1) {
        mat *= fp::one<T>() - tau;
        return;
    }
    if(fp::is_nearly_zero(tau)) { return; }

    if constexpr(tensorwrapper::types::is_interval_v<T>) {
        using value_t             = typename T::value_t;
        const value_t tau_m       = tau.median();
        const Eigen::Index n_rows = mat.rows() - 1;
        const Eigen::Index n_cols = mat.cols();

        std::vector<value_t> tmp(static_cast<std::size_t>(n_cols));
        for(Eigen::Index j = 0; j < n_cols; ++j) {
            value_t sum = mat(0, j).median();
            for(Eigen::Index i = 0; i < n_rows; ++i) {
                sum += essential(i).median() * mat(i + 1, j).median();
            }
            tmp[static_cast<std::size_t>(j)] = sum;
        }

        for(Eigen::Index j = 0; j < n_cols; ++j) {
            const value_t new_v =
              mat(0, j).median() - tau_m * tmp[static_cast<std::size_t>(j)];
            mat(0, j) = T(new_v);
        }
        for(Eigen::Index i = 0; i < n_rows; ++i) {
            const value_t e_i = essential(i).median();
            for(Eigen::Index j = 0; j < n_cols; ++j) {
                const value_t new_v =
                  mat(i + 1, j).median() -
                  tau_m * e_i * tmp[static_cast<std::size_t>(j)];
                mat(i + 1, j) = T(new_v);
            }
        }
        return;
    }

    using RowType = Eigen::Matrix<T, 1, Eigen::Dynamic>;
    Eigen::Map<RowType> tmp(workspace, mat.cols());
    auto bottom   = mat.bottomRows(mat.rows() - 1);
    tmp.noalias() = essential.adjoint() * bottom;
    tmp += mat.row(0);
    mat.row(0) -= tau * tmp;
    bottom.noalias() -= tau * essential * tmp;
}

template<typename MatrixType, typename HCoeffsType>
void householder_qr_inplace(MatrixType& mat, HCoeffsType& h_coeffs,
                            typename MatrixType::Scalar* temp_data) {
    using T          = typename MatrixType::Scalar;
    using Index      = typename MatrixType::Index;
    const Index rows = mat.rows();
    const Index cols = mat.cols();
    const Index size = (std::min)(rows, cols);

    h_coeffs.resize(size);

    for(Index k = 0; k < size; ++k) {
        const Index remaining_rows = rows - k;
        const Index remaining_cols = cols - k - 1;

        auto col_tail = mat.col(k).tail(remaining_rows);
        T beta;
        make_householder(col_tail, h_coeffs.coeffRef(k), beta);
        mat.coeffRef(k, k) = beta;

        if(remaining_cols > 0) {
            auto block = mat.bottomRightCorner(remaining_rows, remaining_cols);
            apply_householder_on_the_left(
              block, mat.col(k).tail(remaining_rows - 1), h_coeffs.coeffRef(k),
              temp_data + k + 1);
        }
    }
}

template<typename T>
using dynamic_matrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

template<typename T>
dynamic_matrix<T> extract_upper_triangular(const dynamic_matrix<T>& factored) {
    const auto n = factored.rows();
    dynamic_matrix<T> R(n, n);
    R.setZero();
    for(std::size_t i = 0; i < static_cast<std::size_t>(n); ++i) {
        for(std::size_t j = i; j < static_cast<std::size_t>(n); ++j) {
            R(i, j) = factored(i, j);
        }
    }
    return R;
}

template<typename T>
dynamic_matrix<T> householder_q_from_compact(
  const dynamic_matrix<T>& factored,
  const Eigen::Matrix<T, Eigen::Dynamic, 1>& h_coeffs) {
    using Index         = Eigen::Index;
    const Index n       = factored.rows();
    dynamic_matrix<T> Q = dynamic_matrix<T>::Identity(n, n);
    apply_householder_q_on_the_right(Q, factored, h_coeffs);
    return Q;
}

template<typename MatrixType, typename HCoeffsType>
void apply_householder_q_on_the_right(MatrixType& mat,
                                      const MatrixType& factored,
                                      const HCoeffsType& h_coeffs) {
    using Index   = typename MatrixType::Index;
    const Index n = mat.cols();
    std::vector<typename MatrixType::Scalar> workspace(
      static_cast<std::size_t>(n));

    for(Index k = 0; k < n; ++k) {
        const Index remaining = n - k;
        auto trailing         = mat.rightCols(remaining);
        const auto essential  = factored.col(k).tail(remaining - 1);
        apply_householder_on_the_right(trailing, essential, h_coeffs.coeff(k),
                                       workspace.data());
    }
}

template<typename T>
inline std::pair<std::vector<T>, std::vector<T>> symmetric_qr_eigen(
  std::span<const T> A, std::size_t n, double tol, std::size_t max_iter) {
    using matrix_type = dynamic_matrix<T>;

    matrix_type A_orig(n, n);
    matrix_type A_mat(n, n);
    for(std::size_t i = 0; i < n; ++i) {
        for(std::size_t j = 0; j < n; ++j) {
            A_orig(i, j) = A[i * n + j];
            A_mat(i, j)  = A[i * n + j];
        }
    }

    matrix_type V = matrix_type::Identity(n, n);
    Eigen::Matrix<T, Eigen::Dynamic, 1> h_coeffs;
    std::vector<T> workspace(2 * n);

    for(std::size_t iter = 0; iter < max_iter; ++iter) {
        T frob = fp::zero<T>();
        for(std::size_t i = 0; i < n; ++i) {
            for(std::size_t j = 0; j < n; ++j) {
                if(i != j) { frob += fp::abs2(A_mat(i, j)); }
            }
        }
        const auto frob_sqrt =
          tensorwrapper::types::pow(fp::as_nonnegative(frob), 0.5);
        if(fp::is_at_most_tolerance(frob_sqrt, tol)) { break; }

        matrix_type factored = A_mat;
        householder_qr_inplace(factored, h_coeffs, workspace.data());

        const matrix_type R = extract_upper_triangular(factored);
        A_mat               = R;
        apply_householder_q_on_the_right(A_mat, factored, h_coeffs);
        A_mat = (A_mat + A_mat.transpose()) / fp::two<T>();
        apply_householder_q_on_the_right(V, factored, h_coeffs);
    }

    const matrix_type rayleigh = V.transpose() * A_orig * V;

    std::vector<std::size_t> order(n);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](std::size_t i, std::size_t j) {
        return rayleigh(i, i) < rayleigh(j, j);
    });

    std::vector<T> values(n);
    std::vector<T> vectors(n * n);
    for(std::size_t k = 0; k < n; ++k) {
        const auto idx = order[k];
        values[k]      = rayleigh(idx, idx);
        for(std::size_t i = 0; i < n; ++i) { vectors[i * n + k] = V(i, idx); }
    }
    return {std::move(values), std::move(vectors)};
}

} // namespace scf::eigen_solver::detail
