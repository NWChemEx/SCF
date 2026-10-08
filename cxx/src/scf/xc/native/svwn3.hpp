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
#include <numbers>
#include <stdexcept>

/** @file svwn3.hpp
 *
 *  Templated, closed-shell implementations of Slater exchange and VWN3
 *  correlation. These are written so that they work for any floating-point
 *  type T which supports arithmetic with doubles and has log, atan, sqrt, and
 *  pow (with a double exponent) findable by ADL (e.g., double and Sigma's
 *  Uncertain/Interval types).
 *
 *  The formulas follow LibXC's XC_LDA_X and XC_LDA_C_VWN_3 (see
 *  maple/lda_exc/lda_x.mpl, maple/vwn.mpl, and maple/lda_exc/lda_c_vwn_3.mpl
 *  in the LibXC source). For a spin-unpolarized density the spin polarization
 *  is zero, which removes every term of VWN3 except the paramagnetic
 *  interpolation formula. Unlike LibXC, no density threshold is applied.
 *
 *  As the density goes to zero, the energy density and the potential both go
 *  to zero; however, the expressions below evaluate to NaN at exactly zero
 *  (r_s diverges). A density which is exactly zero (with no uncertainty) is
 *  therefore special cased to return the limit. The potential's derivative
 *  diverges at zero, so a density whose value could be zero or negative
 *  (e.g., an uncertain value whose mean is zero, or an interval with a
 *  non-positive lower bound) has no meaningful result and is rejected.
 */

namespace scf::xc::native::svwn3 {
namespace detail_ {

/// Parameters of the paramagnetic VWN interpolation formula
inline constexpr double A  = 0.0310907;
inline constexpr double b  = 3.72744;
inline constexpr double c  = 12.9352;
inline constexpr double x0 = -0.10498;

/// Q = sqrt(4c - b^2)
inline double Q() { return std::sqrt(4.0 * c - b * b); }

/// b * x0 / X(x0) where X(x) = x^2 + b*x + c
inline constexpr double f2() { return b * x0 / (x0 * x0 + b * x0 + c); }

/// Coefficient of the arctangent term
inline double atan_coef() {
    const auto f3 = 2.0 * (2.0 * x0 + b) / Q();
    return 2.0 * b / Q() - f2() * f3;
}

/// (3/pi)^(1/3)
inline double cx() { return std::cbrt(3.0 / std::numbers::pi); }

/// (3/(4 pi))^(1/3), i.e., r_s * rho^(1/3)
inline double crs() { return std::cbrt(3.0 / (4.0 * std::numbers::pi)); }

/** @brief Is @p rho exactly zero?
 *
 *  @param[in] rho The density to check.
 *
 *  @return True if @p rho is exactly zero (with no uncertainty) and false if
 *          @p rho is strictly positive.
 *
 *  @throw std::domain_error if @p rho could be zero or negative, but is not
 *                           exactly zero. Strong throw guarantee.
 */
template<typename T>
bool is_zero_density(const T& rho) {
    bool is_zero     = false;
    bool is_positive = false;
    if constexpr(requires {
                     rho.lower();
                     rho.upper();
                 }) { // Intervals
        is_zero     = rho.lower() == 0.0 && rho.upper() == 0.0;
        is_positive = rho.lower() > 0.0;
    } else if constexpr(requires {
                            rho.mean();
                            rho.sd();
                        }) { // Uncertain
        is_zero     = rho.mean() == 0.0 && rho.sd() == 0.0;
        is_positive = rho.mean() > 0.0;
    } else {
        is_zero     = rho == 0.0;
        is_positive = rho > 0.0;
    }
    if(!is_zero && !is_positive)
        throw std::domain_error("SVWN3: density must be exactly zero or "
                                "strictly positive");
    return is_zero;
}

/// rho^(1/3)
template<typename T>
T cube_root(const T& rho) {
    using std::pow;
    return pow(rho, 1.0 / 3.0);
}

/// x = sqrt(r_s) from rho^(1/3)
template<typename T>
T x_from_rho13(const T& rho13) {
    using std::sqrt;
    return sqrt(crs() / rho13);
}

/// VWN correlation energy per particle as a function of x = sqrt(r_s)
template<typename T>
T vwn_eps(const T& x) {
    using std::atan;
    using std::log;
    const T X      = x * x + b * x + c;
    const T x_m_x0 = x - x0;
    const T term1  = log(x * x / X);
    const T term2  = atan_coef() * atan(Q() / (2.0 * x + b));
    const T term3  = f2() * log(x_m_x0 * x_m_x0 / X);
    const T eps    = A * (term1 + term2 - term3);
    return eps;
}

/// Derivative of vwn_eps with respect to x
template<typename T>
T vwn_deps_dx(const T& x) {
    const T X       = x * x + b * x + c;
    const T two_x_b = 2.0 * x + b;
    const T dterm1  = 2.0 / x - two_x_b / X;
    const T dterm2 =
      atan_coef() * (-2.0 * Q()) / (two_x_b * two_x_b + Q() * Q());
    const T dterm3 = f2() * (2.0 / (x - x0) - two_x_b / X);
    return A * (dterm1 + dterm2 - dterm3);
}

} // namespace detail_

/// Slater exchange energy per particle, -(3/4)(3/pi)^(1/3) rho^(1/3)
template<typename T>
T slater_eps(const T& rho) {
    return -0.75 * detail_::cx() * detail_::cube_root(rho);
}

/// Slater exchange potential, d(rho * eps_x)/d rho = -(3/pi)^(1/3) rho^(1/3)
template<typename T>
T slater_v(const T& rho) {
    return -detail_::cx() * detail_::cube_root(rho);
}

/// VWN3 correlation energy per particle
template<typename T>
T vwn3_eps(const T& rho) {
    const auto x = detail_::x_from_rho13(detail_::cube_root(rho));
    return detail_::vwn_eps(x);
}

/** @brief VWN3 correlation potential, d(rho * eps_c)/d rho
 *
 *  Since r_s is proportional to rho^(-1/3), rho d/d rho = -(r_s/3) d/d r_s and
 *  r_s d/d r_s = (x/2) d/dx. Thus v_c = eps_c - (x/6) d eps_c / dx.
 */
template<typename T>
T vwn3_v(const T& rho) {
    const auto x = detail_::x_from_rho13(detail_::cube_root(rho));
    return detail_::vwn_eps(x) - x * detail_::vwn_deps_dx(x) / 6.0;
}

/** @brief SVWN3 energy density, rho * (eps_x + eps_c)
 *
 *  @param[in] rho The density at the point of interest.
 *
 *  @return The energy density at @p rho. If @p rho is exactly zero, this is
 *          the zero-density limit, i.e., zero.
 *
 *  @throw std::domain_error if @p rho could be zero or negative, but is not
 *                           exactly zero. Strong throw guarantee.
 */
template<typename T>
T energy_density(const T& rho) {
    if(detail_::is_zero_density(rho)) return T(0.0);
    return rho * (slater_eps(rho) + vwn3_eps(rho));
}

/** @brief SVWN3 potential, v_x + v_c
 *
 *  @param[in] rho The density at the point of interest.
 *
 *  @return The potential at @p rho. If @p rho is exactly zero, this is the
 *          zero-density limit, i.e., zero.
 *
 *  @throw std::domain_error if @p rho could be zero or negative, but is not
 *                           exactly zero. Strong throw guarantee.
 */
template<typename T>
T potential(const T& rho) {
    if(detail_::is_zero_density(rho)) return T(0.0);
    return slater_v(rho) + vwn3_v(rho);
}

} // namespace scf::xc::native::svwn3
