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
#include "../../test_scf.hpp"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <pluginplay/pluginplay.hpp>
#include <random>
#include <scf/scf.hpp>
#include <simde/simde.hpp>
#include <tensorwrapper/buffer/contiguous.hpp>
#include <tensorwrapper/generate/generate.hpp>
#include <tensorwrapper/tensorwrapper.hpp>
#include <tuple>
#include <vector>
#include <wtf/fp/float_view.hpp>

namespace test_eigen_solver {

using EigenSystem         = tensorwrapper::generate::EigenSystem;
using EigenvalueSpacing   = tensorwrapper::generate::EigenvalueSpacing;
using SymmetricMatrixSpec = tensorwrapper::generate::SymmetricMatrixSpec;
constexpr std::size_t kMaxMatrixDim = tensorwrapper::generate::kMaxMatrixDim;

// EigenSystem consistent with A = [[1,2],[2,1]]
template<typename FloatType>
inline EigenSystem classic_2x2() {
    using tensorwrapper::utilities::make_tensor;
    std::vector evals{static_cast<FloatType>(-0.2360680401325226),
                      static_cast<FloatType>(4.2360687255859375)};
    std::vector evecs{static_cast<FloatType>(-0.8506508469581604),
                      static_cast<FloatType>(-0.5257311463356018),
                      static_cast<FloatType>(0.5257311463356018),
                      static_cast<FloatType>(-0.8506508469581604)};
    std::vector mat{static_cast<FloatType>(1.0), static_cast<FloatType>(2.0),
                    static_cast<FloatType>(2.0), static_cast<FloatType>(3.0)};
    EigenSystem rv;
    rv.n            = 2;
    rv.eigenvalues  = make_tensor({2}, evals);
    rv.matrix       = make_tensor({2, 2}, mat);
    rv.eigenvectors = make_tensor({2, 2}, evecs);
    return rv;
}

// Checks that the eigenvalues are in the same order and approximately equal.
inline void require_eigenvalues_approx(const simde::type::tensor& values,
                                       const simde::type::tensor& expected,
                                       double rtol = 1e-6) {
    using tensorwrapper::operations::approximately_equal;
    REQUIRE(approximately_equal(values, expected, rtol));
}

// Ensures residual H = C^T A C - diag(L) is small.
inline void require_eigenpair_residual(const simde::type::tensor& A,
                                       const simde::type::tensor& L,
                                       const simde::type::tensor& C,
                                       double rtol = 1e-6) {
    using tensorwrapper::operations::approximately_equal;
    using tensorwrapper::utilities::diagonal_matrix;
    auto L_diag = diagonal_matrix(L);

    tensorwrapper::Tensor VA, VAV;
    VA("i,k")  = C("j,i") * A("j,k");
    VAV("i,k") = VA("i,j") * C("j,k");
    REQUIRE(approximately_equal(VAV, L_diag, rtol));
}

// #ifdef ENABLE_SIGMA

// template<typename UQType>
// inline void require_uq_eigenvalues_contain(
//   const simde::type::tensor& values, const std::vector<double>& expected) {
//     using tensorwrapper::buffer::make_contiguous;
//     using wtf::fp::float_cast;

//     auto value_buffer = make_contiguous(values.buffer());
//     REQUIRE(value_buffer.shape().extent(0) == expected.size());
//     using value_type = typename UQType::value_t;

//     std::vector<UQType> balls;
//     balls.reserve(expected.size());
//     for(std::size_t i = 0; i < expected.size(); ++i) {
//         const auto elem = value_buffer.get_elem({i});
//         try {
//             balls.push_back(float_cast<UQType>(elem));
//         } catch(const std::runtime_error&) {
//             const auto median = float_cast<value_type>(elem);
//             balls.emplace_back(median);
//         }
//     }

//     auto expected_sorted = expected;
//     std::sort(expected_sorted.begin(), expected_sorted.end());
//     std::sort(balls.begin(), balls.end(), [](const UQType& a, const UQType&
//     b) {
//         return a.median() < b.median();
//     });

//     for(std::size_t i = 0; i < expected_sorted.size(); ++i) {
//         REQUIRE(balls[i].contains(expected_sorted[i]));
//     }
// }

// template<typename UQType>
// inline void require_uq_eigenvalues_contain_balls(
//   const simde::type::tensor& values, const std::vector<UQType>& expected) {
//     using tensorwrapper::buffer::make_contiguous;
//     using wtf::fp::float_cast;

//     auto value_buffer = make_contiguous(values.buffer());
//     REQUIRE(value_buffer.shape().extent(0) == expected.size());
//     for(std::size_t i = 0; i < expected.size(); ++i) {
//         auto value_uq = float_cast<UQType>(value_buffer.get_elem({i}));
//         REQUIRE(value_uq.contains(expected[i].median()));
//     }
// }

// template<typename UQType>
// inline void require_clustered_eigenvalues(const simde::type::tensor& values,
//                                           double center, double rtol = 1e-6)
//                                           {
//     using tensorwrapper::buffer::get_raw_data;
//     using tensorwrapper::buffer::make_contiguous;

//     using wtf::fp::float_cast;
//     auto value_buffer = make_contiguous(values.buffer());
//     REQUIRE(value_buffer.shape().extent(0) >= 2);
//     auto read_uq = [](const auto& elem) -> UQType {
//         using wtf::fp::float_cast;
//         using value_type = typename UQType::value_t;
//         try {
//             return float_cast<UQType>(elem);
//         } catch(const std::runtime_error&) {
//             const auto m = float_cast<value_type>(elem);
//             return UQType(m, m);
//         }
//     };
//     auto v0 = read_uq(value_buffer.get_elem({0}));
//     auto v1 = read_uq(value_buffer.get_elem({1}));
//     // Each certified eigenvalue ball should bracket the (shared) cluster
//     // center, and the two medians should themselves be clustered together to
//     // within the requested tolerance.
//     const auto m0  = v0.median();
//     const auto m1  = v1.median();
//     const auto tol = rtol * (1.0 + std::abs(center));
//     REQUIRE(std::abs(m0 - center) <= tol);
//     REQUIRE(std::abs(m1 - center) <= tol);
//     REQUIRE(std::abs(m0 - m1) <= tol);
// }

// template<typename UQType>
// inline void require_uq_eigenvalues_approx(const simde::type::tensor& values,
//                                           const std::vector<double>&
//                                           expected, double rtol = 1e-6) {
//     using value_type = typename UQType::value_t;
//     using tensorwrapper::buffer::make_contiguous;
//     using tensorwrapper::utilities::make_tensor;
//     using wtf::fp::float_cast;

//     auto value_buffer = make_contiguous(values.buffer());
//     REQUIRE(value_buffer.shape().extent(0) == expected.size());
//     std::vector<double> computed(value_buffer.shape().extent(0));
//     for(std::size_t i = 0; i < computed.size(); ++i) {
//         const auto elem = value_buffer.get_elem({i});
//         try {
//             computed[i] =
//               static_cast<double>(float_cast<UQType>(elem).median());
//         } catch(const std::runtime_error&) {
//             computed[i] = static_cast<double>(float_cast<value_type>(elem));
//         }
//     }
//     auto corr = computed;
//     auto ref  = expected;
//     std::sort(corr.begin(), corr.end());
//     std::sort(ref.begin(), ref.end());
//     auto corr_tensor = make_tensor({corr.size()}, corr);
//     auto ref_tensor  = make_tensor({ref.size()}, ref);
//     REQUIRE(tensorwrapper::operations::approximately_equal(corr_tensor,
//                                                            ref_tensor,
//                                                            rtol));
// }

// template<typename UQType>
// inline simde::type::tensor uq_tensor_as_double(const simde::type::tensor& M)
// {
//     using value_type = typename UQType::value_t;
//     using tensorwrapper::buffer::make_contiguous;
//     using tensorwrapper::shape::Smooth;
//     using wtf::fp::float_cast;

//     auto in          = make_contiguous(M.buffer());
//     auto shape       = in.shape().make_smooth();
//     const auto n0    = shape.extent(0);
//     const auto n1    = shape.rank() == 1 ? 1 : shape.extent(1);
//     Smooth out_shape = shape.rank() == 1 ? Smooth{n0} : Smooth{n0, n1};
//     auto out         = make_contiguous<double>(out_shape);
//     if(shape.rank() == 1) {
//         for(std::size_t i = 0; i < n0; ++i) {
//             const auto elem = in.get_elem({i});
//             try {
//                 out.set_elem(
//                   {i},
//                   static_cast<double>(float_cast<UQType>(elem).median()));
//             } catch(const std::runtime_error&) {
//                 out.set_elem({i},
//                              static_cast<double>(float_cast<value_type>(elem)));
//             }
//         }
//     } else {
//         for(std::size_t i = 0; i < n0; ++i) {
//             for(std::size_t j = 0; j < n1; ++j) {
//                 const auto elem = in.get_elem({i, j});
//                 try {
//                     out.set_elem({i, j}, static_cast<double>(
//                                            float_cast<UQType>(elem).median()));
//                 } catch(const std::runtime_error&) {
//                     out.set_elem({i, j}, static_cast<double>(
//                                            float_cast<value_type>(elem)));
//                 }
//             }
//         }
//     }
//     return simde::type::tensor(out_shape, std::move(out));
// }

// #endif

} // namespace test_eigen_solver
