/*
 * Copyright 2024 NWChemEx-Project
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

#include "fastscf_modules.hpp"
#include <simde/simde.hpp>

namespace fastscf {

using energy_pt = simde::AOEnergy;

MODULE_CTOR(FastSCFEnergy) {
    satisfies_property_type<energy_pt>();
}

MODULE_RUN(FastSCFEnergy) {
    const auto& [aos, cs] = energy_pt::unwrap_inputs(inputs);

    double E0 = 0; /// This is a total energy
    auto rv = results();
    return energy_pt::wrap_results(rv, E0);
}

} // namespace fastscf
