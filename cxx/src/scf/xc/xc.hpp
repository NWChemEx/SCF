/*
 * Copyright 2022 NWChemEx-Project
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
#include <pluginplay/pluginplay.hpp>
#include <scf/property_types/quadrature.hpp>

namespace scf::xc {

DECLARE_MODULE(AOsOnGrid);
DECLARE_MODULE(Density2Grid);
DECLARE_MODULE(Gau2Grid);
DECLARE_MODULE(GridFromFile);
DECLARE_MODULE(GridFromIntegratorXX);
DECLARE_MODULE(UQXCEnergy);
DECLARE_MODULE(UQXCPotential);
DECLARE_MODULE(XCEnergyViaQuadrature);
DECLARE_MODULE(XCPotentialViaQuadrature);

/** @brief Estimates the error of a quadrature from a finer grid.
 *
 *  @tparam QuadPT The property type of the quadrature, e.g.,
 *                 scf::property_types::xc_potential_quadrature.
 */
template<typename QuadPT>
DECLARE_MODULE(ReferenceGridError);

extern template struct ReferenceGridError<
  scf::property_types::xc_energy_quadrature>;
extern template struct ReferenceGridError<
  scf::property_types::xc_potential_quadrature>;

void set_defaults(pluginplay::ModuleMap& modules);
void load_modules(pluginplay::ModuleMap& modules);

} // namespace scf::xc
