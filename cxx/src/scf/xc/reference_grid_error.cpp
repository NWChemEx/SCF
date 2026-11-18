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

#include "libxc/libxc.hpp"
#include "xc.hpp"
#include <cmath>
#include <scf/property_types/xc_uncertain.hpp>
#include <simde/simde.hpp>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace scf::xc {
namespace {

const auto desc = R"(
Reference Grid Error
====================

Estimates the error of an exchange-correlation (XC) quadrature from the
incompleteness of its grid, and attaches it to the result of the quadrature.

The quadrature is evaluated a second time, by the "Reference XC Module"
submodule, on the much finer grid built by the "Reference grid" submodule. Let
d be the reference value minus the given value. The sign of d is taken as the
direction of the error, so the given value becomes an uncertain value whose
nominal value is the given value, whose deviation towards the reference value
is |d|, and whose deviation away from it is "Opposite Side Fraction" times
|d|. This is done element-wise, e.g., for each element of the XC potential.

The kind of uncertain value comes from the "UQ Initializer" submodule. Only the
triplex kind stores both deviations exactly; the other kinds store the
tightest value they can which contains them.

This is an estimate, not a bound: refinement need not converge monotonically,
so the reference value can lie on the wrong side of the exact integral, and a
one-sided deviation then excludes it. Increasing "Opposite Side Fraction"
guards against this at the cost of a wider interval.

The "Reference XC Module" must evaluate the same functional as the module
whose result is being given. Each call costs one quadrature on the reference
grid, including the density and AOs on that grid.

The given value must hold doubles.
)";

/// Returns the AO basis set of an XC energy or XC potential braket
template<typename BraKetType>
const auto& get_aos(const BraKetType& braket) {
    using bra_type = std::decay_t<decltype(braket.bra())>;
    if constexpr(std::is_same_v<bra_type, simde::type::aos>) {
        return braket.bra().ao_basis_set();
    } else {
        return braket.op().get_rhs_particle().basis_set().ao_basis_set();
    }
}

/// Returns the elements of @p t, which must hold doubles, in row-major order
std::vector<double> to_doubles(const simde::type::tensor& t,
                               const std::string& what) {
    using tensorwrapper::buffer::make_contiguous;
    using tensorwrapper::buffer::visit_contiguous_buffer;
    auto visitor = [&](auto span) -> std::vector<double> {
        using value_type = typename decltype(span)::value_type;
        if constexpr(std::is_same_v<std::remove_cv_t<value_type>, double>) {
            return std::vector<double>(span.begin(), span.end());
        } else {
            throw std::runtime_error("Reference Grid Error: the " + what +
                                     " must hold doubles");
        }
    };
    return visit_contiguous_buffer(visitor, make_contiguous(t.buffer()));
}

/// Returns the shape of @p t
tensorwrapper::shape::Smooth get_shape(const simde::type::tensor& t) {
    using tensorwrapper::buffer::make_contiguous;
    const auto shape = make_contiguous(t.buffer()).shape();
    std::vector<std::size_t> extents;
    for(std::size_t i = 0; i < shape.rank(); ++i)
        extents.push_back(shape.extent(i));
    return tensorwrapper::shape::Smooth(extents.begin(), extents.end());
}

} // namespace

template<typename QuadPT>
using error_pt = scf::property_types::XCUncertain<QuadPT>;

using grid_pt = simde::MolecularGrid;
using uq_pt   = simde::UQInitializer;

template<typename QuadPT>
TEMPLATED_MODULE_CTOR(ReferenceGridError, QuadPT) {
    satisfies_property_type<error_pt<QuadPT>>();
    description(desc);

    add_submodule<grid_pt>("Reference grid")
      .set_description("Builds the finer grid the reference value is "
                       "evaluated on");
    add_submodule<QuadPT>("Reference XC Module")
      .set_description("Evaluates the reference value. Must use the same "
                       "functional as the given value");
    add_submodule<uq_pt>("UQ Initializer")
      .set_description("Determines the kind of uncertain value returned");

    add_input<double>("Opposite Side Fraction")
      .set_default(0.0)
      .set_description("The deviation away from the reference value, as a "
                       "fraction of |reference - value|. Must be "
                       "non-negative");
}

template<typename QuadPT>
TEMPLATED_MODULE_RUN(ReferenceGridError, QuadPT) {
    const auto& [braket, grid, value] = error_pt<QuadPT>::unwrap_inputs(inputs);
    const auto f = inputs.at("Opposite Side Fraction").value<double>();
    if(f < 0.0)
        throw std::invalid_argument(
          "Reference Grid Error: Opposite Side Fraction must be non-negative");

    const auto& aos      = get_aos(braket);
    auto& grid_mod       = submods.at("Reference grid");
    const auto& ref_grid = grid_mod.run_as<grid_pt>(libxc::aos2molecule(aos));

    auto& ref_mod  = submods.at("Reference XC Module");
    const auto ref = ref_mod.run_as<QuadPT>(braket, ref_grid);

    const auto v = to_doubles(value, "given value");
    const auto r = to_doubles(ref, "reference value");
    if(v.size() != r.size())
        throw std::runtime_error(
          "Reference Grid Error: the given and reference values differ in "
          "size");

    const auto& factory = submods.at("UQ Initializer").run_as<uq_pt>();

    auto shape = get_shape(value);
    auto rv_t  = simde::visit_uq_kind(factory.kind(), [&](auto type) {
        using float_type = typename decltype(type)::type;
        std::vector<float_type> elements;
        elements.reserve(v.size());
        for(std::size_t i = 0; i < v.size(); ++i) {
            const auto d     = r[i] - v[i];
            const auto plus  = d >= 0.0 ? std::fabs(d) : f * std::fabs(d);
            const auto minus = d >= 0.0 ? f * std::fabs(d) : std::fabs(d);
            auto element     = factory(v[i], minus, plus);
            elements.push_back(wtf::fp::float_cast<float_type>(element));
        }
        tensorwrapper::buffer::Contiguous buffer(std::move(elements), shape);
        return simde::type::tensor(shape, std::move(buffer));
    });

    auto rv = results();
    return error_pt<QuadPT>::wrap_results(rv, std::move(rv_t));
}

template struct ReferenceGridError<scf::property_types::xc_energy_quadrature>;
template struct ReferenceGridError<
  scf::property_types::xc_potential_quadrature>;

} // namespace scf::xc
