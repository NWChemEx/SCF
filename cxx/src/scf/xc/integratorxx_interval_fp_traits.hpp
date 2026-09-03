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
#include <integratorxx/util/fp_traits.hpp>
#include <limits>
#include <sigma/interval/interval.hpp>

namespace IntegratorXX {

/// Specializes fp_traits<T> for sigma::Interval<T>. IntegratorXX::fp_traits
/// is a full-struct customization point (see fp_traits.hpp), so this
/// reimplements every method the code actually exercises rather than
/// overriding just one -- everything below except clamp() and infinity() is
/// a faithful copy of the primary template's own (ADL-based) default
/// behavior, kept here only because specializing the struct at all replaces
/// the whole thing.
template<typename T>
struct fp_traits<sigma::Interval<T>, void> {
    using value_type = sigma::Interval<T>;

    static value_type log(const value_type& x) {
        using std::log;
        return log(x);
    }

    static value_type exp(const value_type& x) {
        using std::exp;
        return exp(x);
    }

    static value_type sqrt(const value_type& x) {
        using std::sqrt;
        return sqrt(x);
    }

    static value_type cos(const value_type& x) {
        using std::cos;
        return cos(x);
    }

    static value_type sin(const value_type& x) {
        using std::sin;
        return sin(x);
    }

    static value_type abs(const value_type& x) {
        using std::abs;
        return abs(x);
    }

    template<typename U>
    static value_type pow(const value_type& x, const U& p) {
        using std::pow;
        return pow(x, p);
    }

    /// Intersects with [lo, hi] instead of comparing: sigma::Interval's
    /// relational operators are a pragmatic, median-based total order meant
    /// for sorting (see grid_from_integratorxx.cpp), not a containment/
    /// enclosure test, so a std::clamp-style compare-and-select would give a
    /// semantically wrong answer here. set_intersection is the
    /// mathematically correct narrowing operation: it can only discard
    /// values that are provably not in [lo, hi], never a value that might be
    /// the true one, so intersecting a sound enclosure with a provably-exact
    /// range keeps it sound while removing the excess width naive interval
    /// arithmetic can introduce (see the callers in partition_weights.hpp).
    static value_type clamp(const value_type& x, const value_type& lo,
                            const value_type& hi) {
        return x.set_intersection(value_type(lo.lower(), hi.upper()));
    }

    static value_type from_integer(ixx_int v) { return value_type(T(v)); }

    static value_type from_real(ixx_real v) {
#ifdef ENABLE_STRING_REALS
        double d{};
        std::from_chars(v.data(), v.data() + v.size(), d);
        return value_type(T(d));
#else
        return value_type(T(v));
#endif
    }

    static value_type divide_integer(ixx_int num, ixx_int den) {
        return from_integer(num) / from_integer(den);
    }

    /// Unlike the primary template's std::numeric_limits<T>::infinity()
    /// default, sigma has no std::numeric_limits<Interval<T>> specialization
    /// -- that default would silently return a default-constructed (empty)
    /// Interval here, not an infinite one.
    static value_type infinity() {
        return value_type(std::numeric_limits<T>::infinity());
    }
};

} // namespace IntegratorXX
