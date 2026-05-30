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
#include <pluginplay/submodule_request.hpp>
#include <simde/simde.hpp>
#include <tensorwrapper/buffer/contiguous.hpp>
#include <tensorwrapper/tensor/tensor.hpp>
#include <tensorwrapper/utilities/diagonal_matrix.hpp>
#include <wtf/fp/float_view.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace scf::eigen_solver {

/// Similarity transform T^T A T (no subtraction of eigenvalues).
template<typename UQType>
auto compute_similarity_transform(const tensorwrapper::Tensor& A,
                                  const tensorwrapper::Tensor& T) {
    tensorwrapper::Tensor TA, TAT;
    TA("i,k")  = T("j,i") * A("j,k");
    TAT("i,k") = TA("i,j") * T("j,k");
    return TAT;
}

// Similarity transform plus eigenvalue subtraction.
template<typename UQType>
auto compute_residual(const tensorwrapper::Tensor& A,
                      const tensorwrapper::Tensor& C,
                      const tensorwrapper::Tensor& L) {
    auto TAC = compute_similarity_transform<UQType>(A, C);

    using tensorwrapper::utilities::diagonal_matrix;
    auto L_diag = diagonal_matrix(L);

    tensorwrapper::Tensor H;
    H("i,j") = TAC("i,j") - L_diag("i,j");
    return H;
}

/// Upper bound on the Euclidean operator norm ||A||_2 from ball.html §2.3.
///
/// For interval/ball matrices we majorize each entry by its magnitude interval
/// and use ||A||_2 <= sqrt(n) max_i sum_j |A_ij|, which is a valid upper bound
/// on the default operator norm used in §6.3.
template<typename UQType>
auto ball_matrix_norm(const tensorwrapper::Tensor& A) {
    auto A_buffer = tensorwrapper::buffer::make_contiguous(A.buffer());
    auto shape    = A_buffer.shape();
    auto n_rows   = shape.extent(0);
    assert(n_rows == shape.extent(1));

    using float_type = typename UQType::value_t;
    using wtf::fp::float_cast;
    float_type max_row_sum = 0.0;
    for(std::size_t i = 0; i < n_rows; ++i) {
        float_type row_sum = 0.0;
        for(std::size_t j = 0; j < n_rows; ++j) {
            auto Aij      = A_buffer.get_elem({i, j});
            auto Aij_uq   = float_cast<UQType>(Aij);
            auto Aij_low  = std::abs(Aij_uq.lower());
            auto Aij_high = std::abs(Aij_uq.upper());
            row_sum += std::max(Aij_low, Aij_high);
        }
        max_row_sum = std::max(max_row_sum, row_sum);
    }
    return std::sqrt(static_cast<float_type>(n_rows)) * max_row_sum;
}

/// Operator norm of a diagonal matrix given by a vector of eigenvalues.
template<typename FloatType>
auto diagonal_vector_norm(const tensorwrapper::Tensor& values) {
    auto values_buf =
      tensorwrapper::buffer::get_raw_data<FloatType>(values.buffer());
    FloatType norm = 0.0;
    for(const auto& li : values_buf) { norm = std::max(norm, std::abs(li)); }
    return norm;
}

/// Separation number kappa(Λ) from ball.html §6.3.
template<typename FloatType>
auto compute_kappa(const tensorwrapper::Tensor& values) {
    auto values_buf =
      tensorwrapper::buffer::get_raw_data<FloatType>(values.buffer());
    auto n_rows         = values_buf.size();
    FloatType max_value = 0.0;
    for(std::size_t i = 0; i < n_rows; ++i) {
        auto Lii = values_buf[i];
        if(Lii == 0.0) { return std::numeric_limits<FloatType>::infinity(); }
        max_value = std::max(max_value, 1.0 / std::abs(Lii));
        for(std::size_t j = i + 1; j < n_rows; ++j) {
            auto Ljj  = values_buf[j];
            auto diff = std::abs(Ljj - Lii);
            if(diff == 0.0) {
                return std::numeric_limits<FloatType>::infinity();
            }
            max_value = std::max(max_value, 1.0 / diff);
        }
    }
    return max_value;
}

/// ball.html §6.3 condition (32): ||H|| <= 1 / (80 n kappa^2 ||Λ||).
template<typename FloatType>
bool certification_condition_met(const FloatType& h_norm, std::size_t n,
                                 const FloatType& kappa,
                                 const FloatType& lambda_norm) {
    if(!std::isfinite(h_norm) || !std::isfinite(kappa) ||
       !std::isfinite(lambda_norm)) {
        return false;
    }
    if(lambda_norm == 0.0) return false;
    const auto threshold =
      1.0 / (80.0 * static_cast<FloatType>(n) * kappa * kappa * lambda_norm);
    return h_norm <= threshold;
}

/// ball.html §6.3 certified inflation radius: eta <= 6 n sqrt(kappa) ||H||.
template<typename FloatType>
auto compute_eta(const FloatType& h_norm, std::size_t n,
                 const FloatType& kappa) {
    return 6.0 * static_cast<FloatType>(n) * std::sqrt(kappa) * h_norm;
}

/// Finest clustering: each index is its own cluster.
inline std::vector<std::size_t> finest_cluster(std::size_t n) {
    std::vector<std::size_t> clusters(n);
    std::iota(clusters.begin(), clusters.end(), 0);
    return clusters;
}

inline bool same_cluster(std::size_t i, std::size_t j,
                         const std::vector<std::size_t>& clusters) {
    return clusters[i] == clusters[j];
}

/// Max-norm of off-diagonal part of M (double precision).
inline double off_block_max_norm_double(
  const tensorwrapper::Tensor& M, const std::vector<std::size_t>& clusters) {
    using tensorwrapper::buffer::make_contiguous;
    using wtf::fp::float_cast;
    auto M_buf         = make_contiguous(M.buffer());
    auto n             = M_buf.shape().extent(0);
    double max_row_sum = 0.0;
    for(std::size_t i = 0; i < n; ++i) {
        double row_sum = 0.0;
        for(std::size_t j = 0; j < n; ++j) {
            if(same_cluster(i, j, clusters)) continue;
            row_sum += std::abs(float_cast<double>(M_buf.get_elem({i, j})));
        }
        max_row_sum = std::max(max_row_sum, row_sum);
    }
    return max_row_sum;
}

/// One fundamental Newton step (ball.html §6.3) in double precision.
inline std::pair<tensorwrapper::Tensor, tensorwrapper::Tensor>
fundamental_iteration_step(const tensorwrapper::Tensor& T,
                           const tensorwrapper::Tensor& A,
                           const std::vector<std::size_t>& clusters) {
    using tensorwrapper::buffer::get_raw_data;
    using tensorwrapper::buffer::make_contiguous;
    using tensorwrapper::shape::Smooth;

    auto T_buf = make_contiguous(T.buffer());
    auto A_buf = make_contiguous(A.buffer());
    auto n     = T_buf.shape().extent(0);
    assert(n == A_buf.shape().extent(0));

    auto T_data = get_raw_data<double>(T_buf);
    auto A_data = get_raw_data<double>(A_buf);

    std::vector<double> M(n * n, 0.0);
    for(std::size_t i = 0; i < n; ++i) {
        for(std::size_t k = 0; k < n; ++k) {
            double sum = 0.0;
            for(std::size_t j = 0; j < n; ++j) {
                for(std::size_t l = 0; l < n; ++l) {
                    sum +=
                      T_data[j * n + i] * A_data[j * n + l] * T_data[l * n + k];
                }
            }
            M[i * n + k] = sum;
        }
    }

    std::vector<double> lambda(n);
    for(std::size_t i = 0; i < n; ++i) lambda[i] = M[i * n + i];

    std::vector<double> E(n * n, 0.0);
    for(std::size_t i = 0; i < n; ++i) {
        for(std::size_t j = 0; j < n; ++j) {
            if(same_cluster(i, j, clusters)) continue;
            auto denom = lambda[j] - lambda[i];
            if(std::abs(denom) < 1e-300) {
                throw std::runtime_error(
                  "fundamental_iteration_step: zero eigenvalue separation");
            }
            E[i * n + j] = M[i * n + j] / denom;
        }
    }

    std::vector<double> T_new(n * n, 0.0);
    for(std::size_t i = 0; i < n; ++i) {
        for(std::size_t k = 0; k < n; ++k) {
            double sum = 0.0;
            for(std::size_t j = 0; j < n; ++j) {
                double ipe = (j == k ? 1.0 : 0.0) + E[j * n + k];
                sum += T_data[i * n + j] * ipe;
            }
            T_new[i * n + k] = sum;
        }
    }

    Smooth vector_shape{n};
    Smooth matrix_shape{n, n};
    auto T_new_buffer  = make_contiguous<double>(matrix_shape);
    auto lambda_buffer = make_contiguous<double>(vector_shape);
    for(std::size_t i = 0; i < n; ++i) {
        lambda_buffer.set_elem({i}, lambda[i]);
        for(std::size_t j = 0; j < n; ++j) {
            T_new_buffer.set_elem({i, j}, T_new[i * n + j]);
        }
    }

    tensorwrapper::Tensor T_new_tensor(matrix_shape, std::move(T_new_buffer));
    tensorwrapper::Tensor lambda_tensor(vector_shape, std::move(lambda_buffer));
    return {T_new_tensor, lambda_tensor};
}

/// Build a double matrix from the medians of a UQ buffer.
template<typename UQType>
auto median_matrix(const tensorwrapper::buffer::Contiguous& A_buf) {
    using namespace tensorwrapper::buffer;
    using namespace tensorwrapper::shape;

    auto A_data  = get_raw_data<UQType>(A_buf);
    auto A_shape = A_buf.shape().make_smooth();
    assert(A_shape.rank() == 2);
    Smooth shape({A_shape.extent(0), A_shape.extent(1)});

    std::vector<double> A_median(shape.size());
    for(std::size_t i = 0; i < shape.size(); ++i) {
        A_median[i] = A_data[i].median();
    }

    Contiguous A_median_buffer(std::move(A_median), shape);
    return tensorwrapper::Tensor(shape, std::move(A_median_buffer));
}

template<typename UQType>
auto wrap_subdiagonalization(const tensorwrapper::buffer::Contiguous& A_buf,
                             pluginplay::SubmoduleRequest& eigen_solver_mod) {
    auto A_tensor = median_matrix<UQType>(A_buf);
    using pt      = simde::EigenSolve;
    return eigen_solver_mod.run_as<pt>(A_tensor);
}

template<typename UQType>
auto convert_to_uq(const tensorwrapper::Tensor& values,
                   const tensorwrapper::Tensor& vectors) {
    using float_type = typename UQType::value_t;
    using tensorwrapper::buffer::make_contiguous;
    auto values_buf  = make_contiguous(values.buffer());
    auto vectors_buf = make_contiguous(vectors.buffer());
    auto n           = values_buf.shape().extent(0);
    tensorwrapper::shape::Smooth vector_shape{n};
    tensorwrapper::shape::Smooth matrix_shape{n, n};
    auto uq_values_buffer  = make_contiguous<UQType>(vector_shape);
    auto uq_vectors_buffer = make_contiguous<UQType>(matrix_shape);
    using wtf::fp::float_cast;
    for(std::size_t i = 0; i < n; ++i) {
        auto vi = float_cast<float_type>(values_buf.get_elem({i}));
        uq_values_buffer.set_elem({i}, UQType(vi));
        for(std::size_t j = 0; j < n; ++j) {
            auto vij = float_cast<float_type>(vectors_buf.get_elem({i, j}));
            uq_vectors_buffer.set_elem({i, j}, UQType(vij));
        }
    }
    tensorwrapper::Tensor uq_values(vector_shape, std::move(uq_values_buffer));
    tensorwrapper::Tensor uq_vectors(matrix_shape,
                                     std::move(uq_vectors_buffer));
    return std::make_pair(uq_values, uq_vectors);
}

template<typename UQType>
auto t_ball(const tensorwrapper::Tensor& C,
            const typename UQType::value_t& eta) {
    using namespace tensorwrapper::buffer;
    UQType pm1(-1.0, 1.0);
    auto shape = make_contiguous(C.buffer()).shape().make_smooth();
    std::vector<UQType> eta_omega(shape.size(), eta * pm1);
    Contiguous eta_omega_buffer(std::move(eta_omega), shape);
    tensorwrapper::Tensor eta_omega_tensor(shape, std::move(eta_omega_buffer));
    tensorwrapper::Tensor T_ball, COmega;
    COmega("i,k") = C("i,j") * eta_omega_tensor("j,k");
    T_ball("i,j") = C("i,j") + COmega("i,j");
    return T_ball;
}

template<typename UQType>
auto l_ball(const tensorwrapper::Tensor& L,
            const typename UQType::value_t& eta) {
    using namespace tensorwrapper::buffer;
    UQType one_pm_eta(1.0 - eta, 1.0 + eta);
    auto shape = make_contiguous(L.buffer()).shape().make_smooth();
    std::vector<UQType> one_pm_eta_data(shape.size(), one_pm_eta);
    Contiguous one_pm_eta_buffer(std::move(one_pm_eta_data), shape);
    tensorwrapper::Tensor one_pm_eta_tensor(shape,
                                            std::move(one_pm_eta_buffer));
    tensorwrapper::Tensor L_ball;
    L_ball("i") = L("i") * one_pm_eta_tensor("i");
    return L_ball;
}

} // namespace scf::eigen_solver
