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

#include "uncertainty.hpp"

namespace scf::uncertainty {
namespace {

const auto desc = R"(
No Convergence Uncertainty
--------------------------

Returns the energy unchanged, i.e., attaches no uncertainty from incomplete
convergence to it.
)";

} // namespace

using pt = scf::property_types::ConvergenceUncertainty;

MODULE_CTOR(NoConvergenceUncertainty) {
    description(desc);
    satisfies_property_type<pt>();
}

MODULE_RUN(NoConvergenceUncertainty) {
    const auto& [E, dE] = pt::unwrap_inputs(inputs);
    auto rv             = results();
    return pt::wrap_results(rv, E);
}

} // namespace scf::uncertainty
