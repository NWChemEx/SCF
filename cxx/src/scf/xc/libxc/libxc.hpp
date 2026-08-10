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

#pragma once
#include <simde/simde.hpp>
#include <utility>

namespace scf::xc::libxc {

DECLARE_MODULE(LibXCEnergy);
DECLARE_MODULE(LibXCPotential);

void set_defaults(pluginplay::ModuleManager& mm);
void load_modules(pluginplay::ModuleManager& mm);

/** @brief Wraps the process of going from NWChemEx's XC specification to
 *         LibXC's.
 */
std::pair<int, int> to_libxc_codes(chemist::qm_operator::xc_functional func);

/// Wraps converting an AOBasisSet to a simde::type::molecule
simde::type::molecule aos2molecule(const simde::type::ao_basis_set& aos);

/** @brief Computes the energy density.
 *
 *  N.b. LibXC computes the energy per particle. This function will weight that
 *  by @p rho_on_grid to get the energy density.
 */
simde::type::tensor libxc_lda_energy_density(
  chemist::qm_operator::xc_functional func,
  const simde::type::tensor& rho_on_grid);

simde::type::tensor libxc_lda_energy_density_derivative(
  chemist::qm_operator::xc_functional func,
  const simde::type::tensor& rho_on_grid);

/** @brief Extracts the weights from a Grid object and puts them in a tensor.
 *
 *  This function is a stop gap. When we have type-erased floats and the
 *  ability to create a tensor from said type-erased floats this function can
 *  go away.
 */
simde::type::tensor tensorify_weights(const chemist::Grid& grid,
                                      parallelzone::runtime::RuntimeView rv);

/// Wraps the process of doing A(m,i) = w(i) * B(m,i);
simde::type::tensor weight_a_matrix(const simde::type::tensor& w,
                                    const simde::type::tensor& b);

/// Wraps the process of doing rho(i) = X(n,i) * aos(n,i);
simde::type::tensor batched_dot(const simde::type::tensor& aos_on_grid,
                                const simde::type::tensor& X,
                                bool sum_row = true);

simde::type::tensor normalize_row(const simde::type::tensor& w,
                                  const simde::type::tensor& b);

} // namespace scf::xc::libxc
