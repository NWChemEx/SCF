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

inline double buffer_elem_as_double(
  const tensorwrapper::buffer::Contiguous::const_reference& elem) {
    using wtf::fp::float_cast;
    try {
        return float_cast<double>(elem);
    } catch(const std::runtime_error&) {
        return static_cast<double>(float_cast<float>(elem));
    }
}

inline std::vector<double> eigenvalues_vector(const EigenSystem& system) {
    using tensorwrapper::buffer::get_raw_data;
    using tensorwrapper::buffer::make_contiguous;
    auto data =
      get_raw_data<const double>(make_contiguous(system.eigenvalues.buffer()));
    return std::vector<double>(data.begin(), data.end());
}

inline EigenSystem identity_system(std::size_t n) {
    using tensorwrapper::utilities::make_tensor;

    REQUIRE_NOTHROW(tensorwrapper::generate::require_valid_n(n));
    EigenSystem rv;
    rv.n            = n;
    rv.eigenvalues  = make_tensor({n}, std::vector<double>(n, 1.0));
    rv.matrix       = tensorwrapper::generate::identity_matrix(n);
    rv.eigenvectors = rv.matrix;
    return rv;
}

inline EigenSystem classic_2x2() {
    using tensorwrapper::utilities::make_tensor;

    EigenSystem rv;
    rv.n           = 2;
    rv.eigenvalues = make_tensor(
      {2}, std::vector<double>{-0.2360680401325226, 4.2360687255859375});
    rv.matrix = make_tensor({2, 2}, std::vector<double>{1.0, 2.0, 2.0, 3.0});
    rv.eigenvectors = make_tensor(
      {2, 2}, std::vector<double>{-0.8506508469581604, -0.5257311463356018,
                                  0.5257311463356018, -0.8506508469581604});
    return rv;
}

inline EigenSystem identity_2x2_degenerate() {
    SymmetricMatrixSpec spec;
    spec.n                = 2;
    spec.condition_number = 1.0;
    spec.min_eigenvalue   = 1.0;
    spec.spacing          = EigenvalueSpacing::Degenerate;
    spec.n_clusters       = 1;
    spec.seed             = 7;
    return tensorwrapper::generate::generate_eigen_system(spec);
}

template<typename FloatType>
inline simde::type::tensor matrix_as(const EigenSystem& system) {
    using tensorwrapper::buffer::get_raw_data;
    using tensorwrapper::buffer::make_contiguous;
    using tensorwrapper::utilities::make_tensor;

    const auto n = system.n;
    auto data =
      get_raw_data<const double>(make_contiguous(system.matrix.buffer()));
    REQUIRE(data.size() == n * n);
    std::vector<FloatType> elements(n * n);
    for(std::size_t i = 0; i < n * n; ++i) {
        elements[i] = static_cast<FloatType>(data[i]);
    }
    return make_tensor({n, n}, elements);
}

#ifdef ENABLE_SIGMA
template<typename UQType>
inline simde::type::tensor matrix_as_uq(const EigenSystem& system) {
    using value_type = typename UQType::value_t;
    using tensorwrapper::buffer::get_raw_data;
    using tensorwrapper::buffer::make_contiguous;
    using tensorwrapper::shape::Smooth;

    const auto n = system.n;
    auto exact =
      get_raw_data<const double>(make_contiguous(system.matrix.buffer()));
    REQUIRE(exact.size() == n * n);
    Smooth shape{n, n};
    auto buffer = make_contiguous<UQType>(shape);
    for(std::size_t i = 0; i < n; ++i) {
        for(std::size_t j = 0; j < n; ++j) {
            const auto x = exact[i * n + j];
            buffer.set_elem({i, j}, UQType(static_cast<value_type>(x)));
        }
    }
    return simde::type::tensor(shape, std::move(buffer));
}

template<typename UQType>
inline simde::type::tensor eigenvectors_as(const EigenSystem& system) {
    using value_type = typename UQType::value_t;
    using tensorwrapper::buffer::get_raw_data;
    using tensorwrapper::buffer::make_contiguous;
    using tensorwrapper::shape::Smooth;

    const auto n = system.n;
    auto exact =
      get_raw_data<const double>(make_contiguous(system.eigenvectors.buffer()));
    REQUIRE(exact.size() == n * n);
    Smooth shape{n, n};
    auto buffer = make_contiguous<UQType>(shape);
    for(std::size_t i = 0; i < n; ++i) {
        for(std::size_t j = 0; j < n; ++j) {
            const auto x = exact[i * n + j];
            buffer.set_elem({i, j}, UQType(static_cast<value_type>(x)));
        }
    }
    return simde::type::tensor(shape, std::move(buffer));
}

template<typename UQType>
inline simde::type::tensor eigenvalues_as(const EigenSystem& system) {
    using value_type = typename UQType::value_t;
    using tensorwrapper::buffer::get_raw_data;
    using tensorwrapper::buffer::make_contiguous;
    using tensorwrapper::shape::Smooth;

    const auto n = system.n;
    auto exact =
      get_raw_data<const double>(make_contiguous(system.eigenvalues.buffer()));
    REQUIRE(exact.size() == n);
    Smooth shape{n};
    auto buffer = make_contiguous<UQType>(shape);
    for(std::size_t i = 0; i < n; ++i) {
        buffer.set_elem({i}, UQType(static_cast<value_type>(exact[i])));
    }
    return simde::type::tensor(shape, std::move(buffer));
}
#endif

#ifdef ENABLE_SIGMA

template<typename UQType>
inline simde::type::tensor noisy_matrix(const EigenSystem& system, double t,
                                        std::uint64_t seed = 99) {
    using tensorwrapper::buffer::get_raw_data;
    using tensorwrapper::buffer::make_contiguous;
    using tensorwrapper::shape::Smooth;
    using value_type = typename UQType::value_t;

    const auto n = system.n;
    auto matrix_data =
      get_raw_data<const double>(make_contiguous(system.matrix.buffer()));
    REQUIRE(matrix_data.size() == n * n);

    auto noisy = tensorwrapper::generate::add_noise(system.matrix, t, seed);

    auto noisy_buf = make_contiguous(noisy.buffer());
    Smooth shape{n, n};
    auto buffer = make_contiguous<UQType>(shape);
    for(std::size_t i = 0; i < n; ++i) {
        for(std::size_t j = 0; j < n; ++j) {
            const auto v = static_cast<value_type>(
              buffer_elem_as_double(noisy_buf.get_elem({i, j})));
            buffer.set_elem({i, j}, UQType(v, v));
        }
    }
    return simde::type::tensor(shape, std::move(buffer));
}

#endif

/** @brief Builds a "fuzzy" identity matrix for generalized eigen tests.
 *
 *  Returns an @p n by @p n metric matrix whose every entry carries a symmetric
 *  absolute interval of half-width @p halfwidth: the diagonal entries are
 *  `[1 - hw, 1 + hw]` and the off-diagonals are `[-hw, hw]`. The small
 *  off-diagonal width lifts the otherwise exact eigenvalue degeneracy of the
 *  identity, which keeps the ball-arithmetic eigensolver (and the subsequent
 *  `B^{-1/2}` step of the generalized solver) numerically well posed.
 */
template<typename UQType>
inline simde::type::tensor noisy_identity(std::size_t n, double halfwidth) {
    using tensorwrapper::buffer::make_contiguous;
    using tensorwrapper::shape::Smooth;
    using value_type = typename UQType::value_t;

    Smooth shape{n, n};
    auto buffer   = make_contiguous<UQType>(shape);
    const auto hw = static_cast<value_type>(halfwidth);
    for(std::size_t i = 0; i < n; ++i) {
        for(std::size_t j = 0; j < n; ++j) {
            const auto center = (i == j) ? value_type{1} : value_type{0};
            buffer.set_elem({i, j}, UQType(center - hw, center + hw));
        }
    }
    return simde::type::tensor(shape, std::move(buffer));
}

inline std::vector<double> tensor_eigenvalues(
  const simde::type::tensor& values) {
    using tensorwrapper::buffer::make_contiguous;
    auto buf = make_contiguous(values.buffer());
    std::vector<double> rv(buf.shape().extent(0));
    for(std::size_t i = 0; i < rv.size(); ++i) {
        rv[i] = buffer_elem_as_double(buf.get_elem({i}));
    }
    return rv;
}

template<typename FloatType = double>
inline void require_eigenvalues_approx(const simde::type::tensor& values,
                                       const std::vector<double>& expected,
                                       double rtol = 1e-6) {
    using tensorwrapper::utilities::make_tensor;
    auto computed = tensor_eigenvalues(values);
    auto corr     = computed;
    auto ref      = expected;
    std::sort(corr.begin(), corr.end());
    std::sort(ref.begin(), ref.end());

    auto corr_tensor = make_tensor({corr.size()}, corr);
    auto ref_tensor  = make_tensor({ref.size()}, ref);
    REQUIRE(tensorwrapper::operations::approximately_equal(corr_tensor,
                                                           ref_tensor, rtol));
}

template<typename FloatType = double>
inline void require_eigenpair_residual(const simde::type::tensor& A,
                                       const simde::type::tensor& values,
                                       const simde::type::tensor& vectors,
                                       double rtol = 1e-6) {
    using tensorwrapper::buffer::make_contiguous;
    const auto n        = make_contiguous(A.buffer()).shape().extent(0);
    const auto vals     = tensor_eigenvalues(values);
    const auto A_buf    = make_contiguous(A.buffer());
    const auto V_buf    = make_contiguous(vectors.buffer());
    double max_residual = 0.0;

    for(std::size_t k = 0; k < n; ++k) {
        const auto lambda = vals[k];
        for(std::size_t i = 0; i < n; ++i) {
            double Av_i = 0.0;
            for(std::size_t j = 0; j < n; ++j) {
                Av_i += buffer_elem_as_double(A_buf.get_elem({i, j})) *
                        buffer_elem_as_double(V_buf.get_elem({j, k}));
            }
            const auto residual = std::abs(
              Av_i - lambda * buffer_elem_as_double(V_buf.get_elem({i, k})));
            max_residual = std::max(max_residual, residual);
        }
    }
    REQUIRE(max_residual <= rtol * (1.0 + std::abs(vals.back())));
}

#ifdef ENABLE_SIGMA

template<typename UQType>
inline void require_uq_eigenvalues_contain(
  const simde::type::tensor& values, const std::vector<double>& expected) {
    using tensorwrapper::buffer::make_contiguous;
    using wtf::fp::float_cast;

    auto value_buffer = make_contiguous(values.buffer());
    REQUIRE(value_buffer.shape().extent(0) == expected.size());
    using value_type = typename UQType::value_t;
    for(std::size_t i = 0; i < expected.size(); ++i) {
        const auto elem = value_buffer.get_elem({i});
        try {
            auto value_uq = float_cast<UQType>(elem);
            REQUIRE(value_uq.contains(expected[i]));
        } catch(const std::runtime_error&) {
            const auto median = float_cast<value_type>(elem);
            REQUIRE(median == Catch::Approx(expected[i]).margin(1e-2));
        }
    }
}

template<typename UQType>
inline void require_uq_eigenvalues_contain_balls(
  const simde::type::tensor& values, const std::vector<UQType>& expected) {
    using tensorwrapper::buffer::make_contiguous;
    using wtf::fp::float_cast;

    auto value_buffer = make_contiguous(values.buffer());
    REQUIRE(value_buffer.shape().extent(0) == expected.size());
    for(std::size_t i = 0; i < expected.size(); ++i) {
        auto value_uq = float_cast<UQType>(value_buffer.get_elem({i}));
        REQUIRE(value_uq.contains(expected[i].median()));
    }
}

template<typename UQType>
inline void require_clustered_eigenvalues(const simde::type::tensor& values,
                                          double center, double rtol = 1e-6) {
    using tensorwrapper::buffer::get_raw_data;
    using tensorwrapper::buffer::make_contiguous;

    using wtf::fp::float_cast;
    auto value_buffer = make_contiguous(values.buffer());
    REQUIRE(value_buffer.shape().extent(0) >= 2);
    auto read_uq = [](const auto& elem) -> UQType {
        using wtf::fp::float_cast;
        using value_type = typename UQType::value_t;
        try {
            return float_cast<UQType>(elem);
        } catch(const std::runtime_error&) {
            const auto m = float_cast<value_type>(elem);
            return UQType(m, m);
        }
    };
    auto v0 = read_uq(value_buffer.get_elem({0}));
    auto v1 = read_uq(value_buffer.get_elem({1}));
    // Each certified eigenvalue ball should bracket the (shared) cluster
    // center, and the two medians should themselves be clustered together to
    // within the requested tolerance.
    const auto m0  = v0.median();
    const auto m1  = v1.median();
    const auto tol = rtol * (1.0 + std::abs(center));
    REQUIRE(std::abs(m0 - center) <= tol);
    REQUIRE(std::abs(m1 - center) <= tol);
    REQUIRE(std::abs(m0 - m1) <= tol);
}

template<typename UQType>
inline void require_uq_eigenvalues_approx(const simde::type::tensor& values,
                                          const std::vector<double>& expected,
                                          double rtol = 1e-6) {
    using value_type = typename UQType::value_t;
    using tensorwrapper::buffer::make_contiguous;
    using tensorwrapper::utilities::make_tensor;
    using wtf::fp::float_cast;

    auto value_buffer = make_contiguous(values.buffer());
    REQUIRE(value_buffer.shape().extent(0) == expected.size());
    std::vector<double> computed(value_buffer.shape().extent(0));
    for(std::size_t i = 0; i < computed.size(); ++i) {
        const auto elem = value_buffer.get_elem({i});
        try {
            computed[i] =
              static_cast<double>(float_cast<UQType>(elem).median());
        } catch(const std::runtime_error&) {
            computed[i] = static_cast<double>(float_cast<value_type>(elem));
        }
    }
    auto corr = computed;
    auto ref  = expected;
    std::sort(corr.begin(), corr.end());
    std::sort(ref.begin(), ref.end());
    auto corr_tensor = make_tensor({corr.size()}, corr);
    auto ref_tensor  = make_tensor({ref.size()}, ref);
    REQUIRE(tensorwrapper::operations::approximately_equal(corr_tensor,
                                                           ref_tensor, rtol));
}

#endif

} // namespace test_eigen_solver
