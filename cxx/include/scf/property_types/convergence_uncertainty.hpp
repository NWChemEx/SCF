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
#include <simde/types.hpp>

namespace scf::property_types {

/** @brief Attaches the uncertainty from incomplete convergence to an energy.
 *
 *  Called by an iterative method once it has converged. The result is the
 *  energy with uncertainty attached, *e.g.*, from the change in the energy
 *  in the final iteration.
 */
DECLARE_PROPERTY_TYPE(ConvergenceUncertainty);

PROPERTY_TYPE_INPUTS(ConvergenceUncertainty) {
    using tensor_type = const simde::type::tensor&;
    auto rv           = pluginplay::declare_input()
                .add_field<tensor_type>("Energy")
                .template add_field<tensor_type>("Energy change");
    rv["Energy"].set_description("The converged energy");
    rv["Energy change"].set_description(
      "The change in the energy in the final iteration");
    return rv;
}

PROPERTY_TYPE_RESULTS(ConvergenceUncertainty) {
    auto rv =
      pluginplay::declare_result().add_field<simde::type::tensor>("Energy");
    rv["Energy"].set_description("The energy with uncertainty attached");
    return rv;
}

} // namespace scf::property_types
