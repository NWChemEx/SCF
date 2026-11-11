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

#include "../eigen_solver/inflate_uncertainty.hpp"
#include "uncertainty.hpp"

namespace scf::uncertainty {
namespace {

const auto desc = R"(
Energy Change Uncertainty
-------------------------

Attaches the uncertainty from incomplete convergence to the converged energy.
The energy gets an extra, independent uncertainty whose radius is the
magnitude of the center of the change in the energy in the final iteration.
The radius of that change is not used, since it is roughly the energy's own
uncertainty and would count it twice. Plain floating-point energies are
returned unchanged.
)";

} // namespace

using pt = scf::property_types::ConvergenceUncertainty;

MODULE_CTOR(EnergyChangeUncertainty) {
    description(desc);
    satisfies_property_type<pt>();
}

MODULE_RUN(EnergyChangeUncertainty) {
    const auto& [E, dE] = pt::unwrap_inputs(inputs);

    auto E_uq = E;
    eigen_solver::inflate_uncertainty(E_uq,
                                      eigen_solver::uq_center_magnitude(dE));

    auto rv = results();
    return pt::wrap_results(rv, E_uq);
}

} // namespace scf::uncertainty
