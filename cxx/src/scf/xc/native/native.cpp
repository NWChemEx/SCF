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

#include "native.hpp"
#include "svwn3.hpp"
#include <span>
#include <type_traits>
#include <vector>

namespace scf::xc::native {
namespace {

/** @brief Is @p T a floating-point type the native functionals support?
 *
 *  The visitors are instantiated for every type in
 *  tensorwrapper::types::floating_point_types, but the functionals need log,
 *  atan, and fractional powers, as well as arithmetic with double. Only
 *  double and the double-based Uncertain and Interval types provide all of
 *  those (the affine types lack atan and fractional powers).
 */
template<typename T>
constexpr bool is_supported_v =
  std::is_same_v<T, float> || std::is_same_v<T, double> ||
  std::is_same_v<T, tensorwrapper::types::udouble> ||
  std::is_same_v<T, tensorwrapper::types::idouble>;

/// Wraps a vector of values in a rank 1 tensor
template<typename T>
simde::type::tensor make_vector(std::vector<T> values) {
    tensorwrapper::shape::Smooth shape{values.size()};
    tensorwrapper::buffer::Contiguous buffer(std::move(values), shape);
    return simde::type::tensor(shape, std::move(buffer));
}

/// Applies @p fxn to each element of the visited buffer
template<typename FunctionType>
struct ElementwiseKernel {
    template<typename FloatType>
    simde::type::tensor operator()(const std::span<FloatType>& rho) {
        using clean_type = std::remove_cv_t<FloatType>;
        if constexpr(is_supported_v<clean_type>) {
            std::vector<clean_type> rv(rho.size());
            for(std::size_t i = 0; i < rho.size(); ++i) rv[i] = fxn(rho[i]);
            return make_vector(std::move(rv));
        } else {
            throw std::runtime_error(
              "Native XC: floating-point type not supported");
        }
    }

    FunctionType fxn;
};

/// Extracts the grid weights as the type of the visited buffer
struct WeightsKernel {
    template<typename FloatType>
    simde::type::tensor operator()(const std::span<FloatType>&) {
        using clean_type = std::remove_cv_t<FloatType>;
        std::vector<clean_type> weights(m_pgrid->size());
        for(std::size_t i = 0; i < weights.size(); ++i)
            weights[i] = m_pgrid->at(i).get_weight().value<clean_type>();
        return make_vector(std::move(weights));
    }

    const chemist::Grid* m_pgrid;
};

template<typename KernelType>
simde::type::tensor visit(KernelType&& k, const simde::type::tensor& t) {
    using tensorwrapper::buffer::make_contiguous;
    using tensorwrapper::buffer::visit_contiguous_buffer;
    const auto& buffer = make_contiguous(t.buffer());
    return visit_contiguous_buffer(std::forward<KernelType>(k), buffer);
}

} // namespace

void load_modules(pluginplay::ModuleManager& mm) {
    mm.add_module<SVWN3Energy>("SVWN3 Energy");
    mm.add_module<SVWN3Potential>("SVWN3 Potential");
}

void set_defaults(pluginplay::ModuleManager& mm) {
    mm.change_submod("SVWN3 Energy", "Integration grid",
                     "Grid From IntegratorXX");
    mm.change_submod("SVWN3 Potential", "Integration grid",
                     "Grid From IntegratorXX");
    mm.change_submod("SVWN3 Energy", "Density on a grid", "Density2Grid");
    mm.change_submod("SVWN3 Potential", "Density on a grid", "Density2Grid");
    mm.change_submod("SVWN3 Potential", "AOs on a grid", "AOs on a grid");
}

simde::type::tensor svwn3_energy_density(
  const simde::type::tensor& rho_on_grid) {
    auto fxn = [](const auto& rho) { return svwn3::energy_density(rho); };
    return visit(ElementwiseKernel<decltype(fxn)>{fxn}, rho_on_grid);
}

simde::type::tensor svwn3_potential(const simde::type::tensor& rho_on_grid) {
    auto fxn = [](const auto& rho) { return svwn3::potential(rho); };
    return visit(ElementwiseKernel<decltype(fxn)>{fxn}, rho_on_grid);
}

simde::type::tensor tensorify_weights(const chemist::Grid& grid,
                                      const simde::type::tensor& like) {
    return visit(WeightsKernel{&grid}, like);
}

} // namespace scf::xc::native
