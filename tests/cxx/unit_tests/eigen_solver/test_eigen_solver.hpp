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

#ifdef ENABLE_SIGMA
template<typename UQType>
inline void require_uq_eigenvalues_contain(const simde::type::tensor& values,
                                           const simde::type::tensor& expected,
                                           double noise_level) {
    using tensorwrapper::buffer::make_contiguous;
    using wtf::fp::float_cast;

    auto values_data   = get_raw_data<UQType>(values.buffer());
    auto expected_data = get_raw_data<UQType>(expected.buffer());
    REQUIRE(values_data.size() == expected_data.size());

    for(std::size_t i = 0; i < values_data.size(); ++i) {
        REQUIRE(values_data[i].contains(expected_data[i].median()));
        // Ensure width did not grow too much
        REQUIRE(values_data[i].width() <= 50 * noise_level);
    }
}

template<typename UQType>
inline void require_uq_eigenpair_residual(const simde::type::tensor& A,
                                          const simde::type::tensor& L,
                                          const simde::type::tensor& C,
                                          double noise_level) {
    using tensorwrapper::operations::approximately_equal;
    using tensorwrapper::utilities::diagonal_matrix;
    auto L_diag = diagonal_matrix(L);

    tensorwrapper::Tensor VA, VAV, H;
    VA("i,k")  = C("j,i") * A("j,k");
    VAV("i,k") = VA("i,j") * C("j,k");
    H("i,j")   = VAV("i,j") - L_diag("i,j");

    auto H_data = get_raw_data<UQType>(H.buffer());
    for(std::size_t i = 0; i < H_data.size(); ++i) {
        REQUIRE(H_data[i].contains(0.0));
        REQUIRE(H_data[i].width() <= 50 * noise_level);
    }
}

#endif

} // namespace test_eigen_solver
