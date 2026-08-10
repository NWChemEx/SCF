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
#include <cmath>
#include <cstddef>
#include <simde/simde.hpp>
#include <span>
#include <stdexcept>
#include <tensorwrapper/tensorwrapper.hpp>

namespace scf::eigen_solver {
namespace detail {

/** @brief Adds a fresh, independent uncertainty of radius @p delta to every
 *         element of a buffer, in place.
 *
 *  For each element x the kernel forms x + e, where e is a brand-new error
 *  term of radius @p delta centered at zero (built with construct_uq_type). A
 *  new error symbol is minted per element, so the added uncertainties are
 *  mutually independent.
 *
 *  The body is only instantiated/selected for a writable buffer whose element
 *  type is an actual UQ type. For plain float/double it is a no-op, so the
 *  inflation only happens "when UQ is enabled". (construct_uq_type on a non-UQ
 *  type returns center + radius, which would SHIFT rather than widen -- another
 *  reason to gate on is_uq_type_v.)
 */
struct InflateUncertaintyKernel {
    double m_delta;
    explicit InflateUncertaintyKernel(double delta) : m_delta(delta) {}

    template<typename FloatType>
    void operator()(std::span<FloatType> data) {
        using clean_t = std::decay_t<FloatType>;
        if constexpr(!std::is_const_v<FloatType> &&
                     tensorwrapper::types::is_uq_type_v<clean_t>) {
            using tensorwrapper::types::construct_uq_type;
            using tensorwrapper::types::uq_center;
            using value_t = decltype(uq_center(std::declval<clean_t>()));
            const value_t delta(m_delta);
            for(auto& x : data) {
                x = x + construct_uq_type<clean_t>(value_t(0), delta);
            }
        }
    }
};

/** @brief Adds a fresh, independent uncertainty whose per-element radius is
 *         |center(source)|, in place.
 *
 *  data[k] += e_k, where e_k is a new, independent error term of radius
 *  |center(source[k])|. Using the magnitude of the (real) center keeps the
 *  radius real and non-negative regardless of the sign of the source.
 *
 *  Applied to the MO coefficients with source = the converged density change
 *  dp, the density P = sum_i C_mi C_ni picks up a first-order (linear) change
 *  of order |C| * radius per element. With radius = |center(dp_mn)| (~dP) and
 *  |C| ~ O(1), each density element's uncertainty lands at ~dP. (A sqrt(dP)
 *  radius would instead make the linear term ~sqrt(dP), overshooting by orders
 *  of magnitude for small dP.) A no-op for non-UQ types.
 */
struct InflateFromMagnitudeKernel {
    std::size_t m_n;
    const tensorwrapper::buffer::Contiguous& m_src;

    InflateFromMagnitudeKernel(std::size_t n,
                               const tensorwrapper::buffer::Contiguous& src) :
      m_n(n), m_src(src) {}

    template<typename FloatType>
    void operator()(std::span<FloatType> data) {
        using clean_t = std::decay_t<FloatType>;
        if constexpr(!std::is_const_v<FloatType> &&
                     tensorwrapper::types::is_uq_type_v<clean_t>) {
            using tensorwrapper::buffer::get_raw_data;
            using tensorwrapper::types::construct_uq_type;
            using tensorwrapper::types::uq_center;
            using value_t = decltype(uq_center(std::declval<clean_t>()));
            auto src      = get_raw_data<clean_t>(m_src);
            for(std::size_t k = 0; k < m_n; ++k) {
                const value_t radius = std::abs(value_t(uq_center(src[k])));
                data[k] =
                  data[k] + construct_uq_type<clean_t>(value_t(0), radius);
            }
        }
    }
};

/** @brief Returns the absolute value of the center of element 0.
 *
 *  For a convergence residual (e.g. de = e - e_old) this is the quantity the
 *  convergence check bounds by the tolerance -- the incomplete-convergence
 *  error -- NOT the propagated uncertainty (which is the radius and is of the
 *  same order as the energy's own uncertainty). Dispatched on the element type;
 *  for non-UQ types the center is just the value.
 */
struct CenterMagnitudeKernel {
    template<typename FloatType>
    double operator()(std::span<FloatType> data) const {
        using tensorwrapper::types::uq_center;
        return std::abs(double(uq_center(data[0])));
    }
};

} // namespace detail

/** @brief Inflates the uncertainty of every element of @p t by @p delta.
 *
 *  Adds an independent error term of radius @p delta to each element of @p t,
 *  in place. Existing uncertainty on each element is preserved. A no-op for
 *  non-UQ element types and for @p delta <= 0.
 */
inline void inflate_uncertainty(simde::type::tensor& t, double delta) {
    if(delta <= 0.0) return;
    using tensorwrapper::buffer::make_contiguous;
    using tensorwrapper::buffer::visit_contiguous_buffer;

    auto& buf = make_contiguous(t.buffer());
    detail::InflateUncertaintyKernel kernel(delta);
    visit_contiguous_buffer(kernel, buf);
}

/** @brief Inflates each element of @p t by an independent error term of radius
 *         |center(source)|.
 *
 *  Element-wise: t[k] gains an independent uncertainty of radius
 *  |center(source[k])|. @p t and @p source must have the same number of
 *  elements. A no-op for non-UQ element types.
 *
 *  @param[in,out] t The tensor whose element uncertainties are inflated.
 *  @param[in] source The tensor supplying the per-element radii (the absolute
 *                    values of its element centers).
 *
 *  @throw std::runtime_error if @p t and @p source have different sizes.
 */
inline void inflate_uncertainty_from(simde::type::tensor& t,
                                     const simde::type::tensor& source) {
    using tensorwrapper::buffer::make_contiguous;
    using tensorwrapper::buffer::visit_contiguous_buffer;

    auto& buf        = make_contiguous(t.buffer());
    const auto& sbuf = make_contiguous(source.buffer());
    if(buf.size() != sbuf.size())
        throw std::runtime_error(
          "inflate_uncertainty_from: tensor and source must have equal sizes");

    detail::InflateFromMagnitudeKernel kernel(buf.size(), sbuf);
    visit_contiguous_buffer(kernel, buf);
}

/** @brief Returns |center| of a scalar tensor's element.
 *
 *  For a converged residual this is the incomplete-convergence error (bounded
 * by the convergence tolerance), as opposed to the propagated uncertainty
 * radius.
 */
inline double uq_center_magnitude(const simde::type::tensor& t) {
    using tensorwrapper::buffer::make_contiguous;
    using tensorwrapper::buffer::visit_contiguous_buffer;

    const auto& buf = make_contiguous(t.buffer());
    detail::CenterMagnitudeKernel kernel;
    return visit_contiguous_buffer(kernel, buf);
}

} // namespace scf::eigen_solver
