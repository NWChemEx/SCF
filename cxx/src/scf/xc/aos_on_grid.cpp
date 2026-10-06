/*
 * Copyright 2025 NWChemEx-Project
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

#include "xc.hpp"
#include <chemist/experimental/basis_set/ao_shell_view.hpp>
#include <chemist/experimental/basis_set/cartesian_ao.hpp>
#include <chemist/experimental/basis_set/molecular_basis_set.hpp>
#include <chemist/experimental/basis_set/spherical_ao.hpp>
#include <simde/integration_grids/collocation_matrix.hpp>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace scf::xc {
namespace {
const auto desc = R"(

AO CollocationMatrix
--------------------

Evaluates each AO of the basis set at each point of the grid. Element (m, i)
of the result is the value of AO m at grid point i.

The AOs are normalized following the convention of the integral libraries
(e.g., libint and gau2grid): each primitive and the contraction as a whole are
normalized, but the per-component Cartesian factor (e.g., sqrt(3) for d_xy) is
not applied. Spherical AOs are fully normalized.
)";

using legacy_basis_type = simde::type::ao_basis_set;
using basis_set_type    = chemist::experimental::MolecularBasisSet;
using fp_types          = chemist::types::floating_point_types;
using buffer_type       = wtf::buffer::FloatBuffer;

/** @brief Converts the legacy AO basis set to the experimental one.
 *
 *  The legacy basis set only holds doubles, whereas the experimental one can
 *  only be evaluated at points holding the same type as its parameters. The
 *  parameters are therefore converted to @p T, the type of the grid.
 */
template<typename T>
basis_set_type to_experimental(const legacy_basis_type& ao_basis) {
    using chemist::experimental::AtomicBasisSet;
    using chemist::experimental::Point;
    using chemist::experimental::ShellPurity;

    std::vector<AtomicBasisSet> atoms;
    for(const auto& atomic_basis : ao_basis) {
        const auto& center = atomic_basis.center();
        auto purity        = ShellPurity::cartesian;
        if(atomic_basis.size() > 0 &&
           atomic_basis.at(0).pure() == chemist::ShellType::pure)
            purity = ShellPurity::pure;

        // Name and atomic number are optional in the legacy class
        Point r(static_cast<T>(center.x()), static_cast<T>(center.y()),
                static_cast<T>(center.z()));
        AtomicBasisSet atom(atomic_basis.basis_set_name().value_or(""),
                            atomic_basis.atomic_number().value_or(0),
                            std::move(r), purity);
        for(const auto& shell_i : atomic_basis) {
            const bool is_pure = shell_i.pure() == chemist::ShellType::pure;
            if(is_pure != (purity == ShellPurity::pure))
                throw std::runtime_error(
                  "AOsOnGrid: shells on a center must all be pure or all be "
                  "Cartesian.");

            const auto& cg = shell_i.contracted_gaussian();
            std::vector<T> coefs, exps;
            for(const auto& prim : cg) {
                coefs.push_back(static_cast<T>(prim.coefficient()));
                exps.push_back(static_cast<T>(prim.exponent()));
            }
            atom.add_shell(shell_i.l(), coefs.begin(), coefs.end(),
                           exps.begin(), exps.end());
        }
        atoms.push_back(std::move(atom));
    }
    return basis_set_type(atoms.begin(), atoms.end());
}

/// Converts @p ao_basis to the experimental basis in the type of @p points
template<typename PointSetType>
basis_set_type to_experimental(const legacy_basis_type& ao_basis,
                               const PointSetType& points) {
    // An empty grid has no type to match, so any type will do
    if(points.size() == 0) return to_experimental<double>(ao_basis);
    auto visitor = [&]<typename T>(std::span<T>) {
        return to_experimental<std::remove_const_t<T>>(ao_basis);
    };
    return wtf::buffer::visit_contiguous_buffer_view<fp_types>(
      visitor, points.get_x_buffer());
}

/** @brief Evaluates the AOs of @p shell at @p points.
 *
 *  Cartesian AOs leave out N^AO_ijk, which is the convention of the integral
 *  libraries. Spherical AOs are fully normalized (and agree with the integral
 *  libraries by construction).
 *
 *  @return A buffer holding the values in row-major order, i.e., all of the
 *          points for the first AO, then all of the points for the second AO,
 *          etc.
 */
template<typename ShellType, typename PointSetType>
buffer_type evaluate_shell(const ShellType& shell, const PointSetType& points) {
    std::vector<buffer_type> aos;
    if(shell.is_cartesian()) {
        const auto& cart = chemist::experimental::as_cartesian_shell(shell);
        for(std::size_t i = 0; i < cart.size(); ++i)
            aos.push_back(wtf::buffer::make_float_buffer<fp_types>(
              cart.at(i)->cg_normalized_evaluate(points)));
    } else {
        const auto& sph = chemist::experimental::as_spherical_shell(shell);
        for(std::size_t i = 0; i < sph.size(); ++i)
            aos.push_back(wtf::buffer::make_float_buffer<fp_types>(
              sph.at(i)->normalized_evaluate(points)));
    }
    return wtf::buffer::concatenate<fp_types>(aos);
}

} // namespace

using pt = simde::AOCollocationMatrix;

MODULE_CTOR(AOsOnGrid) {
    satisfies_property_type<pt>();
    description(desc);
}

MODULE_RUN(AOsOnGrid) {
    const auto& [grid, ao_basis] = pt::unwrap_inputs(inputs);
    const auto points            = grid.get_points();
    const auto basis             = to_experimental(ao_basis, points);

    std::vector<buffer_type> shells;
    for(std::size_t s = 0; s < basis.n_shells(); ++s)
        shells.push_back(evaluate_shell(*basis.shell(s), points));

    tensorwrapper::shape::Smooth matrix_shape{basis.n_aos(), grid.size()};
    tensorwrapper::buffer::Contiguous buffer(
      wtf::buffer::concatenate<fp_types>(shells), matrix_shape);
    simde::type::tensor collocation_matrix(matrix_shape, std::move(buffer));
    auto rv = results();
    return pt::wrap_results(rv, std::move(collocation_matrix));
};

} // namespace scf::xc
