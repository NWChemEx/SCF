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
#include <simde/simde.hpp>
#include <string>
#include <vector>

namespace scf::xc {

/// The property type for the XC energy of a restricted SCF wavefunction
using xc_energy_pt =
  simde::eval_braket<simde::type::rscf_wf, simde::type::XC_e_type,
                     simde::type::rscf_wf>;

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

/** @brief Forwards an XC calculation to the module for its functional.
 *
 *  The functional is read from the XC operator of the input braket. Each
 *  functional in xc_dispatch_table() is a submodule, keyed by
 *  XCDispatchEntry::key.
 *
 *  @tparam PT The property type of the XC calculation, i.e.,
 *             simde::aos_xc_e_aos or xc_energy_pt.
 */
template<typename PT>
DECLARE_MODULE(XCDispatcher);

extern template struct XCDispatcher<simde::aos_xc_e_aos>;
extern template struct XCDispatcher<xc_energy_pt>;

/// An XC functional and the modules which, by default, evaluate it
struct XCDispatchEntry {
    /// The functional
    chemist::qm_operator::xc_functional functional;

    /// The key of the functional's submodule in the XC dispatchers
    std::string key;

    /// The key of the module computing the functional's XC potential
    std::string potential;

    /// The key of the module computing the functional's XC energy
    std::string energy;
};

/** @brief The functionals the XC dispatchers support in this build.
 *
 *  SVWN3 always uses the native modules. The remaining functionals use GauXC,
 *  if SCF was built with it, and otherwise LibXC, if SCF was built with it
 *  and LibXC supports the functional.
 *
 *  @return One entry per supported functional.
 *
 *  @throw std::bad_alloc if there is a problem allocating the table. Strong
 *                        throw guarantee.
 */
std::vector<XCDispatchEntry> xc_dispatch_table();

/** @brief Returns @p t converted to the floating-point type of @p like.
 *
 *  Tensors of different floating-point types can not be contracted, so the
 *  quadratures use this to convert, e.g., the AOs on a grid, which hold
 *  doubles, to the type of an uncertain density. @p t is returned as is if it
 *  already has that type. Otherwise it must hold doubles, each of which
 *  becomes an element of the new type with no uncertainty.
 *
 *  @throw std::runtime_error if @p t has neither the type of @p like nor
 *                            double. Strong throw guarantee.
 */
simde::type::tensor to_type_of(const simde::type::tensor& t,
                               const simde::type::tensor& like);

void set_defaults(pluginplay::ModuleMap& modules);
void load_modules(pluginplay::ModuleMap& modules);

} // namespace scf::xc
