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

#include "../eigen_solver/eigenvector_uncertainty.hpp"
#include "../eigen_solver/inflate_uncertainty.hpp"
#include "uncertainty.hpp"

namespace scf::uncertainty {
namespace {

const auto desc = R"(
First-Order MO Uncertainty
--------------------------

Attaches uncertainty to the converged orbital coefficients from two sources:

- The uncertainty in the Fock matrix. Each orbital gets the first-order
  correction from its coupling, through the Fock matrix, to the orbitals it is
  not (nearly) degenerate with. The coupling is evaluated in the arithmetic of
  the uncertain type, so the coefficients share the Fock matrix's error terms.
- Incomplete convergence. Each coefficient is inflated by the magnitude of the
  matching element of the change in the density matrix in the final
  iteration.

This is a one-shot correction. The coefficients it returns must not be fed
back into the SCF iterations, where the uncertainty would compound. The
overlap matrix is assumed to be exact, so the metric terms of the generalized
eigenproblem are omitted. Plain floating-point coefficients are returned
unchanged.
)";

} // namespace

using pt = scf::property_types::MOUncertainty;

MODULE_CTOR(FirstOrderMOUncertainty) {
    description(desc);
    satisfies_property_type<pt>();
}

MODULE_RUN(FirstOrderMOUncertainty) {
    const auto& [F, S, eps, C, dP] = pt::unwrap_inputs(inputs);

    auto C_uq = C;
    eigen_solver::attach_eigenvector_uncertainty(C_uq, F, eps);
    eigen_solver::inflate_uncertainty_from(C_uq, dP);

    auto rv = results();
    return pt::wrap_results(rv, C_uq);
}

} // namespace scf::uncertainty
