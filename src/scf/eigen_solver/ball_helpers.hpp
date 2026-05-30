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

template<typename UQType>
auto ball_matrix_norm(const tensorwrapper::Tensor& A) {
    auto A_buffer = tensorwrapper::buffer::make_contiguous(A.buffer());
    auto shape    = A_buffer.shape();
    auto n_rows   = shape.extent(0);
    assert(n_rows == shape.extent(1)); // Assume square matrix for now

    using float_type = typename UQType::value_t;
    using wtf::fp::float_cast;
    float_type sum = 0.0;
    for(std::size_t i = 0; i < n_rows; ++i) {
        for(std::size_t j = 0; j < n_rows; ++j) {
            auto Aij      = A_buffer.get_elem({i, j});
            auto Aij_uq   = float_cast<UQType>(Aij);
            auto Aij_low  = std::abs(Aij_uq.lower());
            auto Aij_high = std::abs(Aij_uq.upper());
            auto mag      = std::max(Aij_low, Aij_high);
            sum += mag * mag;
        }
    }
    return std::sqrt(sum);
}

template<typename FloatType>
auto compute_k(const tensorwrapper::Tensor& values) {
    // max[max(1/(L_jj -Lii)), max(1/L_ii)]
    auto values_buf =
      tensorwrapper::buffer::get_raw_data<FloatType>(values.buffer());
    auto n_rows         = values_buf.size();
    FloatType max_diff  = 0.0;
    FloatType max_value = 0.0;
    for(std::size_t i = 0; i < n_rows; ++i) {
        auto Lii  = values_buf[i];
        max_value = std::max(max_value, 1.0 / Lii);
        for(std::size_t j = i + 1; j < n_rows; ++j) {
            auto Ljj  = values_buf[j];
            auto diff = std::abs(Ljj - Lii);
            max_value = std::max(max_value, 1.0 / diff);
        }
    }
    return std::max(max_diff, max_value);
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
