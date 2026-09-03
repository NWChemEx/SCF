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

#include "xc.hpp"
#include <algorithm>
#include <concepts>
#include <integratorxx/molecular_grid/defaults.hpp>
#include <integratorxx/molecular_grid/molecular_grid.hpp>
#include <simde/integration_grids/molecular_grid.hpp>
#include <tensorwrapper/types/floating_point.hpp>
#include <wtf/enums/enums.hpp>
#ifdef ENABLE_SIGMA
#include "integratorxx_interval_fp_traits.hpp"
#endif

namespace scf::xc {
namespace {

const auto desc = R"(

GridFromIntegratorXX
---------------------

Generates a molecular integration grid directly from the molecule using
IntegratorXX: a per-element radial x angular product quadrature (optionally
pruned), recentered on each atom and combined with fuzzy-cell partition
weights. Unlike GridFromFile, this module actually consumes the "Molecule"
property-type input.

The "Float Type" input selects the concrete floating-point type used to
generate the grid (points and weights alike) before it is type-erased into
the returned chemist::Grid, by name as registered with WTF (see
wtf::enums::from_string) -- e.g. "float", "double", "udouble", "idouble".
"udouble" and "idouble" additionally require Sigma support (ENABLE_SIGMA) to
be compiled in: without it, tensorwrapper's sigma-backed types are never
registered with WTF, so wtf::enums::from_string throws for those names
exactly as it would for any other unrecognized Float Type. "adouble",
"tadouble", and "tmdouble" are NOT currently supported -- sigma::Affine
(and, transitively, ThresholdedAffine/TaylorModel, which wrap
Affine/Interval) doesn't implement cos, which two of IntegratorXX's four
radial quadratures call -- and always throw regardless of ENABLE_SIGMA, so
as not to require those (currently-uncompilable) template instantiations at
all.

"Partition Scheme" (Becke/SSF/LKO fuzzy-cell weighting) is applied only for
Float Types that satisfy std::totally_ordered ("double", "float", "udouble",
and "idouble"). For "idouble" this order is sigma::Interval's median-based
comparison, a pragmatic ordering for sorting rather than a rigorous
enclosure relation -- see sigma::Interval's relational operators for
details.
)";

IntegratorXX::PruningScheme pruning_scheme_from_string(
  const std::string& pruning_scheme) {
    if(pruning_scheme == "Unpruned")
        return IntegratorXX::PruningScheme::Unpruned;
    if(pruning_scheme == "Robust") return IntegratorXX::PruningScheme::Robust;
    if(pruning_scheme == "Treutler")
        return IntegratorXX::PruningScheme::Treutler;
    throw std::runtime_error("GridFromIntegratorXX: Invalid Pruning Scheme " +
                             pruning_scheme);
}

IntegratorXX::PartitionScheme partition_scheme_from_string(
  const std::string& partition_scheme) {
    if(partition_scheme == "Becke") return IntegratorXX::PartitionScheme::Becke;
    if(partition_scheme == "SSF") return IntegratorXX::PartitionScheme::SSF;
    if(partition_scheme == "LKO") return IntegratorXX::PartitionScheme::LKO;
    throw std::runtime_error("GridFromIntegratorXX: Invalid Partition Scheme " +
                             partition_scheme);
}

/// Generates the whole molecular grid for concrete type @p T and flattens
/// it into a type-erased chemist::Grid.
template<typename T, typename MoleculeType>
chemist::Grid build_grid(const MoleculeType& molecule,
                         IntegratorXX::RadialQuad radial_quad,
                         std::size_t radial_size, std::size_t angular_size,
                         IntegratorXX::PruningScheme pruning,
                         IntegratorXX::PartitionScheme partition,
                         std::size_t max_batch_sz) {
    std::vector<IntegratorXX::AtomicId> atomic_ids;
    std::vector<IntegratorXX::cartesian_pt_t<T>> positions;
    atomic_ids.reserve(molecule.size());
    positions.reserve(molecule.size());
    for(const auto& atom : molecule) {
        atomic_ids.push_back(static_cast<IntegratorXX::AtomicId>(atom.Z()));
        positions.push_back({T(atom.x()), T(atom.y()), T(atom.z())});
    }

    std::vector<IntegratorXX::AtomicId> unique_ids(atomic_ids);
    std::sort(unique_ids.begin(), unique_ids.end());
    unique_ids.erase(std::unique(unique_ids.begin(), unique_ids.end()),
                     unique_ids.end());

    using defaults_type = IntegratorXX::MolecularGridDefaults<T>;
    auto spec_map       = defaults_type::create_default_grid_spec_map(
      unique_ids, pruning, radial_quad, radial_size, angular_size);
    auto element_grids = defaults_type::generate_gridmap(spec_map);
    auto atoms = IntegratorXX::make_atom_instances<T>(atomic_ids, positions,
                                                      element_grids);

    IntegratorXX::MolecularGrid<T> mg(std::move(atoms), max_batch_sz);
    // IntegratorXX::reference_lko_partition_weights<T> std::sorts on
    // inter-atom distances of type T, which requires a strict total order;
    // instantiating MolecularGrid<T>::apply_partition_weights at all pulls
    // in that (and the Becke/SSF) branch together (a switch's branches are
    // all compiled, not just the one taken at runtime), so it only compiles
    // for a T that is totally ordered. sigma::Interval supplies a
    // pragmatic, median-based set of relational operators for exactly this
    // purpose (true interval/enclosure comparison isn't a total order, so
    // it isn't what's used here); sigma::Affine/ThresholdedAffine/
    // TaylorModel do not yet expose one, so partition weighting is skipped
    // for those, leaving raw (unpartitioned) atomic-quadrature weights
    // rather than failing to compile entirely.
    if constexpr(std::totally_ordered<T>) {
        mg.apply_partition_weights(partition);
    }

    const auto& points  = mg.points();
    const auto& weights = mg.weights();

    std::vector<chemist::GridPoint> grid_points;
    grid_points.reserve(mg.npts());
    for(std::size_t i = 0; i < mg.npts(); ++i) {
        grid_points.emplace_back(weights[i], points[i][0], points[i][1],
                                 points[i][2]);
    }
    return chemist::Grid(grid_points.begin(), grid_points.end());
}

} // namespace

using pt = simde::MolecularGrid;

MODULE_CTOR(GridFromIntegratorXX) {
    satisfies_property_type<pt>();

    add_input<std::string>("Float Type").set_default("double");
    add_input<std::string>("Radial Quadrature Type").set_default("MuraKnowles");
    add_input<std::size_t>("Radial Size").set_default(std::size_t(150));
    add_input<std::size_t>("Angular Size").set_default(std::size_t(194));
    add_input<std::string>("Pruning Scheme").set_default("Unpruned");
    add_input<std::string>("Partition Scheme").set_default("Becke");
    add_input<std::size_t>("Max Batch Size").set_default(std::size_t(512));

    description(desc);
}

MODULE_RUN(GridFromIntegratorXX) {
    const auto& [molecule] = pt::unwrap_inputs(inputs);

    const auto& float_type = inputs.at("Float Type").value<std::string>();

    // Deliberately rejected here, before any WTF/build_grid machinery runs
    // (not even behind #ifdef ENABLE_SIGMA): sigma::Interval now provides
    // the cos/pow and the (median-based, sort-only) total order that
    // IntegratorXX's radial quadratures and box-partitioning/LKO weighting
    // need (see sigma::Interval's relational operators), but Affine (and,
    // transitively, ThresholdedAffine/TaylorModel, which wrap
    // Affine/Interval) doesn't implement cos, which IntegratorXX's Becke and
    // TreutlerAhlrichs radial quadratures call by way of their underlying
    // Gauss-Chebyshev primitive quadrature -- and, because
    // IntegratorXX::make_radial_traits<T> switches on RadialQuad at runtime,
    // ALL FOUR radial quadrature types get instantiated for T regardless of
    // which one is actually selected, so this can't be dodged by picking a
    // different "Radial Quadrature Type" input. Attempting the instantiation
    // anyway would make *any* ENABLE_SIGMA=ON build of this project fail to
    // compile, not just a call with this Float Type.
    if(float_type == "adouble" || float_type == "tadouble" ||
       float_type == "tmdouble") {
        throw std::runtime_error(
          "GridFromIntegratorXX: Float Type '" + float_type +
          "' is not yet supported (blocked on IntegratorXX gaps for this "
          "type -- see grid_from_integratorxx.cpp).");
    }

    auto kind        = wtf::enums::from_string(float_type);
    auto radial_quad = IntegratorXX::radial_from_string(
      inputs.at("Radial Quadrature Type").value<std::string>());
    auto radial_size  = inputs.at("Radial Size").value<std::size_t>();
    auto angular_size = inputs.at("Angular Size").value<std::size_t>();
    auto pruning      = pruning_scheme_from_string(
      inputs.at("Pruning Scheme").value<std::string>());
    auto partition = partition_scheme_from_string(
      inputs.at("Partition Scheme").value<std::string>());
    auto max_batch_sz = inputs.at("Max Batch Size").value<std::size_t>();

    // Only float/double/udouble/idouble reach here (adouble and friends are
    // rejected above); this tuple is what wtf::enums::detail_::
    // dispatch_by_kind searches to find the concrete type identified by
    // kind. tensorwrapper::types::udouble/idouble are plain double when
    // ENABLE_SIGMA is off, but that's harmless here: without ENABLE_SIGMA
    // those names are never registered with WTF, so wtf::enums::from_string
    // above already throws before kind could resolve to either entry --
    // kind can only actually resolve to them when Sigma support -- and thus
    // the real sigma-backed types -- is available.
    using supported_types =
      std::tuple<float, double, tensorwrapper::types::udouble,
                 tensorwrapper::types::idouble>;

    auto grid = wtf::enums::detail_::dispatch_by_kind<supported_types>(
      kind, [&]<typename T>() {
          return build_grid<T>(molecule, radial_quad, radial_size, angular_size,
                               pruning, partition, max_batch_sz);
      });

    auto rv = results();
    return pt::wrap_results(rv, std::move(grid));
}

} // namespace scf::xc
