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
#include <pluginplay/pluginplay.hpp>
#include <scf/property_types/quadrature.hpp>

namespace scf::property_types {

/** @brief Attaches uncertainty to the result of another property type.
 *
 *  Modules satisfying this property type are error models for the XC energy
 *  and potential, *e.g.*, for the error from the incompleteness of the grid.
 *  They are given the inputs of @p BasePT and the results of calling a module
 *  satisfying @p BasePT with those inputs. The results are those of
 *  @p BasePT, but with uncertainty attached: each value is a UQ type whose
 *  nominal value is the corresponding value of the given results. The module
 *  satisfying this property type gets the UQ type from its own
 *  ``UQ Initializer`` submodule.
 *
 *  This property type is a stop gap. It should eventually replace
 *  SimDE's ``Uncertainty`` property type.
 *
 *  @tparam BasePT The property type of the module whose result is having
 *                 uncertainty attached, *e.g.*, xc_potential_quadrature.
 */
template<typename BasePT>
DECLARE_TEMPLATED_PROPERTY_TYPE(XCUncertain, BasePT);

template<typename BasePT>
TEMPLATED_PROPERTY_TYPE_INPUTS(XCUncertain, BasePT) {
    // Appends each result of BasePT, with the same key and type, as an input
    auto rv = BasePT::inputs() + BasePT::results();
    for(const auto& [key, value] : BasePT::results())
        rv[key].set_description("The value to attach uncertainty to");
    return rv;
}

template<typename BasePT>
TEMPLATED_PROPERTY_TYPE_RESULTS(XCUncertain, BasePT) {
    return BasePT::results();
}

/// Property type for attaching uncertainty to the XC energy from quadrature
using xc_energy_uncertainty = XCUncertain<xc_energy_quadrature>;

/// Property type for attaching uncertainty to the XC potential from quadrature
using xc_potential_uncertainty = XCUncertain<xc_potential_quadrature>;

} // namespace scf::property_types
