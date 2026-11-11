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
No MO Uncertainty
-----------------

Returns the orbital coefficients unchanged, i.e., attaches no uncertainty to
them.
)";

} // namespace

using pt = scf::property_types::MOUncertainty;

MODULE_CTOR(NoMOUncertainty) {
    description(desc);
    satisfies_property_type<pt>();
}

MODULE_RUN(NoMOUncertainty) {
    const auto& [F, S, eps, C, dP] = pt::unwrap_inputs(inputs);
    auto rv                        = results();
    return pt::wrap_results(rv, C);
}

} // namespace scf::uncertainty
