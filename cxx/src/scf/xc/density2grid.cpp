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

#include "libxc/libxc.hpp"
#include "xc.hpp"
#include <simde/integration_grids/collocation_matrix.hpp>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace scf::xc {
namespace {
const auto desc = R"(

DensityCollocationMatrix
-----------------

The AOs on the grid are converted to the floating-point type of the density,
e.g., an uncertain type, since tensors of different types can not be
contracted.
)";

} // namespace

simde::type::tensor to_type_of(const simde::type::tensor& t,
                               const simde::type::tensor& like) {
    using tensorwrapper::buffer::make_contiguous;
    using tensorwrapper::buffer::visit_contiguous_buffer;
    const auto& buffer = make_contiguous(t.buffer());

    auto like_visitor = [&](auto like_span) -> simde::type::tensor {
        using float_type =
          std::remove_cv_t<typename decltype(like_span)::value_type>;
        auto t_visitor = [&](auto t_span) -> simde::type::tensor {
            using t_type =
              std::remove_cv_t<typename decltype(t_span)::value_type>;
            if constexpr(std::is_same_v<t_type, float_type>) {
                return t;
            } else if constexpr(std::is_same_v<t_type, double>) {
                std::vector<float_type> elements;
                elements.reserve(t_span.size());
                for(const auto x : t_span) elements.push_back(float_type(x));
                const auto& view = buffer.shape();
                std::vector<std::size_t> extents;
                for(std::size_t i = 0; i < view.rank(); ++i)
                    extents.push_back(view.extent(i));
                tensorwrapper::shape::Smooth shape(extents.begin(),
                                                   extents.end());
                tensorwrapper::buffer::Contiguous rv(std::move(elements),
                                                     shape);
                return simde::type::tensor(shape, std::move(rv));
            } else {
                throw std::runtime_error(
                  "to_type_of: the tensor must hold doubles or the type to "
                  "convert to");
            }
        };
        return visit_contiguous_buffer(t_visitor, buffer);
    };
    return visit_contiguous_buffer(like_visitor,
                                   make_contiguous(like.buffer()));
}

using pt         = simde::EDensityCollocationMatrix;
using ao2grid_pt = simde::AOCollocationMatrix;

MODULE_CTOR(Density2Grid) {
    satisfies_property_type<pt>();
    description(desc);

    add_submodule<ao2grid_pt>("AOs on a grid");
}

MODULE_RUN(Density2Grid) {
    const auto& [grid, density] = pt::unwrap_inputs(inputs);

    const auto& rho = density.value();
    const auto& aos = density.basis_set().ao_basis_set();

    auto& ao2grid_mod = submods.at("AOs on a grid");
    const auto aos_on_grid =
      to_type_of(ao2grid_mod.run_as<ao2grid_pt>(grid, aos), rho);

    simde::type::tensor X, rho2;
    rho2("m,n") = rho("m,n") * 2.0; // Assumes restricted orbitals
    X("m,i")    = rho2("m,n") * aos_on_grid("n,i");

    auto rho_on_grid = libxc::batched_dot(aos_on_grid, X);

    auto rv = results();
    return pt::wrap_results(rv, std::move(rho_on_grid));
}

} // namespace scf::xc
