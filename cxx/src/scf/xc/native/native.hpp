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
#include <simde/simde.hpp>

/** @file native.hpp
 *
 *  XC functionals implemented natively (i.e., without LibXC or GauXC). Since
 *  they are templated C++ they work with the uncertainty-quantification (UQ)
 *  floating-point types, not just double.
 */

namespace scf::xc::native {

DECLARE_MODULE(SVWN3Energy);
DECLARE_MODULE(SVWN3Potential);

void set_defaults(pluginplay::ModuleManager& mm);
void load_modules(pluginplay::ModuleManager& mm);

/** @brief Computes the SVWN3 energy density, rho * (eps_x + eps_c).
 *
 *  @param[in] rho_on_grid The (total, spin-unpolarized) density evaluated on
 *                         the grid points.
 *
 *  @return A tensor with the same shape and floating-point type as
 *          @p rho_on_grid.
 *
 *  @throw std::runtime_error if @p rho_on_grid's floating-point type is not
 *                            supported. Strong throw guarantee.
 */
simde::type::tensor svwn3_energy_density(
  const simde::type::tensor& rho_on_grid);

/** @brief Computes the SVWN3 potential, v_x + v_c.
 *
 *  @param[in] rho_on_grid The (total, spin-unpolarized) density evaluated on
 *                         the grid points.
 *
 *  @return A tensor with the same shape and floating-point type as
 *          @p rho_on_grid.
 *
 *  @throw std::runtime_error if @p rho_on_grid's floating-point type is not
 *                            supported. Strong throw guarantee.
 */
simde::type::tensor svwn3_potential(const simde::type::tensor& rho_on_grid);

/** @brief Extracts the weights from @p grid and puts them in a tensor.
 *
 *  Unlike libxc::tensorify_weights, the weights are not assumed to be
 *  doubles. Instead they are extracted as the floating-point type of @p like.
 *
 *  @param[in] grid The grid whose weights are being extracted.
 *  @param[in] like A tensor holding the floating-point type to use.
 *
 *  @throw std::runtime_error if the weights of @p grid can not be retrieved as
 *                            @p like's floating-point type. Strong throw
 *                            guarantee.
 */
simde::type::tensor tensorify_weights(const chemist::Grid& grid,
                                      const simde::type::tensor& like);

} // namespace scf::xc::native
