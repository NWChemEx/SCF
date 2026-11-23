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
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <integratorxx/quadratures/s2/lebedev_laikov.hpp>
#include <optional>
#include <scf/property_types/xc_uncertain.hpp>
#include <simde/simde.hpp>
#include <stdexcept>
#include <string>
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
submodule, on a grid which is finer than the given grid. Let d be the
largest distance between the center of the given value and the range of the
reference value, i.e., the largest of |upper(ref) - center(value)| and
|lower(ref) - center(value)|. For a reference which holds doubles, this is
|ref - value|. The given value becomes an uncertain value whose nominal value
is the given value and whose deviation is d on both sides. This is done element-wise, e.g., for each element of the XC
potential. The deviation is symmetric because the side the reference is on
is not a reliable estimate of the direction of the error: as a grid is refined, the error of
the quadrature can change sign.

The finer grid is chosen from the radial and angular sizes the given grid
records (see chemist::Grid), so the given grid must record both. The radial
size is multiplied by 1.5 (rounding up). The angular size, which must be the
size of a Lebedev-Laikov grid, is replaced by the size of the smallest
Lebedev-Laikov grid with at least twice as many points, or by the largest
Lebedev-Laikov grid (5810 points) if there is no such grid. The "Reference
grid" submodule builds the finer grid. It is run with its "Radial Size" and
"Angular Size" inputs set to the finer sizes, so it must have these inputs,
e.g., "Grid From IntegratorXX". Its other inputs are used as they are, and
should match those used to make the given grid.

If the given value holds doubles, the kind of uncertain value comes from the
"UQ Initializer" submodule. If the given value is already uncertain, e.g.,
because the density is, it keeps its type and its uncertainty, and each
element gets an extra, independent uncertainty of d.

This is an estimate, not a bound: the difference between two grids can be
smaller than the error of the coarser one.

The "Reference XC Module" must evaluate the same functional as the module
whose result is being given. Each call costs one quadrature on the finer grid,
including the density and AOs on that grid.

The given value must hold doubles or an uncertain type.
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

/// The centers, lower bounds, and upper bounds of the elements of a tensor
struct Bounds {
    std::vector<double> center;
    std::vector<double> lower;
    std::vector<double> upper;
};

/// Returns the bounds of the elements of @p t, in row-major order
Bounds get_bounds(const simde::type::tensor& t) {
    using tensorwrapper::buffer::make_contiguous;
    using tensorwrapper::buffer::visit_contiguous_buffer;
    using tensorwrapper::types::uq_center;
    using tensorwrapper::types::uq_lower;
    using tensorwrapper::types::uq_upper;
    auto visitor = [](auto span) {
        Bounds rv;
        for(const auto& x : span) {
            rv.center.push_back(double(uq_center(x)));
            rv.lower.push_back(double(uq_lower(x)));
            rv.upper.push_back(double(uq_upper(x)));
        }
        return rv;
    };
    return visit_contiguous_buffer(visitor, make_contiguous(t.buffer()));
}

/// Does @p t hold doubles?
bool holds_doubles(const simde::type::tensor& t) {
    using tensorwrapper::buffer::make_contiguous;
    using tensorwrapper::buffer::visit_contiguous_buffer;
    auto visitor = [](auto span) {
        using value_type = typename decltype(span)::value_type;
        return std::is_same_v<std::remove_cv_t<value_type>, double>;
    };
    return visit_contiguous_buffer(visitor, make_contiguous(t.buffer()));
}

/** @brief Adds an independent uncertainty of radius @p d[i] to element i of
 *         @p t, which must hold an uncertain type, in place.
 *
 *  @throw std::runtime_error if @p t does not hold an uncertain type.
 */
void add_errors(simde::type::tensor& t, const std::vector<double>& d) {
    using tensorwrapper::buffer::make_contiguous;
    using tensorwrapper::buffer::visit_contiguous_buffer;
    auto visitor = [&](auto span) {
        using element_t = typename decltype(span)::element_type;
        using clean_t   = std::remove_cv_t<element_t>;
        if constexpr(!std::is_const_v<element_t> &&
                     tensorwrapper::types::is_uq_type_v<clean_t>) {
            using tensorwrapper::types::construct_uq_type;
            using tensorwrapper::types::uq_center;
            using center_t = decltype(uq_center(std::declval<clean_t>()));
            for(std::size_t i = 0; i < span.size(); ++i)
                span[i] = span[i] + construct_uq_type<clean_t>(center_t(0),
                                                               center_t(d[i]));
        } else {
            throw std::runtime_error("Reference Grid Error: the given value "
                                     "must hold doubles or an uncertain type");
        }
    };
    visit_contiguous_buffer(visitor, make_contiguous(t.buffer()));
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

/// The radial size of the finer grid is this times the given radial size
constexpr double radial_factor = 1.5;

/// The angular size of the finer grid is at least this times the given size
constexpr std::size_t angular_factor = 2;

/// Returns the radial size of the grid finer than one with @p radial_size
std::size_t finer_radial_size(std::size_t radial_size) {
    if(radial_size == 0)
        throw std::invalid_argument(
          "Reference Grid Error: the radial size of the given grid must be "
          "positive");
    const auto r = std::ceil(radial_factor * static_cast<double>(radial_size));
    return static_cast<std::size_t>(r);
}

/** @brief Returns the angular size of the grid finer than one with
 *         @p angular_size.
 *
 *  This is the size of the smallest Lebedev-Laikov grid with at least
 *  angular_factor times @p angular_size points, or the size of the largest
 *  Lebedev-Laikov grid if there is no such grid.
 *
 *  @throw std::invalid_argument if @p angular_size is not the size of a
 *                               Lebedev-Laikov grid.
 */
std::size_t finer_angular_size(std::size_t angular_size) {
    using traits_type =
      IntegratorXX::quadrature_traits<IntegratorXX::LebedevLaikov<double>>;
    auto order = traits_type::algebraic_order_by_npts(
      static_cast<std::int64_t>(angular_size));
    if(order < 0)
        throw std::invalid_argument(
          "Reference Grid Error: the angular size of the given grid (" +
          std::to_string(angular_size) +
          ") is not the size of a Lebedev-Laikov grid");

    const auto target = angular_factor * angular_size;
    auto n_points     = angular_size;
    while(n_points < target) {
        const auto next_order = traits_type::next_algebraic_order(order + 1);
        if(next_order <= order) break; // Already the largest grid
        order = next_order;
        n_points =
          static_cast<std::size_t>(traits_type::npts_by_algebraic_order(order));
    }
    return n_points;
}

/// Returns the @p what size of the given grid, throwing if it is not set
std::size_t get_size(const std::optional<std::size_t>& size,
                     const std::string& what) {
    if(!size.has_value())
        throw std::invalid_argument("Reference Grid Error: the given grid "
                                    "does not record its " +
                                    what + " size");
    return *size;
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
      .set_description("Builds the finer grid. Must have \"Radial Size\" and "
                       "\"Angular Size\" inputs");
    add_submodule<QuadPT>("Reference XC Module")
      .set_description("Evaluates the reference value. Must use the same "
                       "functional as the given value");
    add_submodule<uq_pt>("UQ Initializer")
      .set_description("Determines the kind of uncertain value returned");
}

template<typename QuadPT>
TEMPLATED_MODULE_RUN(ReferenceGridError, QuadPT) {
    const auto& [braket, grid, value] = error_pt<QuadPT>::unwrap_inputs(inputs);

    const auto radial            = get_size(grid.get_radial_size(), "radial");
    const auto angular           = get_size(grid.get_angular_size(), "angular");
    const std::size_t ref_radial = finer_radial_size(radial);
    const std::size_t ref_angular = finer_angular_size(angular);

    // Runs the "Reference grid" submodule with the finer sizes
    const auto& aos = get_aos(braket);
    auto& grid_mod  = submods.at("Reference grid").value();
    auto ref_inputs = grid_mod.inputs();
    ref_inputs = grid_pt::wrap_inputs(ref_inputs, libxc::aos2molecule(aos));
    for(const auto* key : {"Radial Size", "Angular Size"})
        if(!ref_inputs.count(key))
            throw std::invalid_argument(
              "Reference Grid Error: the \"Reference grid\" submodule must "
              "have a \"" +
              std::string(key) + "\" input");
    ref_inputs.at("Radial Size").change(ref_radial);
    ref_inputs.at("Angular Size").change(ref_angular);
    const auto ref_grid =
      std::get<0>(grid_pt::unwrap_results(grid_mod.run(ref_inputs)));

    auto& ref_mod  = submods.at("Reference XC Module");
    const auto ref = ref_mod.run_as<QuadPT>(braket, ref_grid);

    const auto v = get_bounds(value).center;
    const auto r = get_bounds(ref);
    if(v.size() != r.center.size())
        throw std::runtime_error(
          "Reference Grid Error: the given and reference values differ in "
          "size");

    std::vector<double> d(v.size());
    for(std::size_t i = 0; i < v.size(); ++i)
        d[i] =
          std::max(std::fabs(r.upper[i] - v[i]), std::fabs(r.lower[i] - v[i]));

    // An uncertain value keeps its type and gets the error added to it
    if(!holds_doubles(value)) {
        auto rv_t = value;
        add_errors(rv_t, d);
        auto rv = results();
        return error_pt<QuadPT>::wrap_results(rv, std::move(rv_t));
    }

    const auto& factory = submods.at("UQ Initializer").run_as<uq_pt>();

    auto shape = get_shape(value);
    auto rv_t  = simde::visit_uq_kind(factory.kind(), [&](auto type) {
        using float_type = typename decltype(type)::type;
        std::vector<float_type> elements;
        elements.reserve(v.size());
        for(std::size_t i = 0; i < v.size(); ++i) {
            auto element = factory(v[i], d[i], d[i]);
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
