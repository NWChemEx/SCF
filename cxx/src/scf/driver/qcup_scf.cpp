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

#include "driver.hpp"

namespace scf::driver {
namespace {

const auto desc = R"(
QCUP-SCF
--------

An uncertainty-quantified Hartree-Fock calculation. Numerical approximations,
most notably the screening of the electron repulsion integrals (ERIs), are
tracked as uncertainties: the calculation uses the uncertain type chosen with
"UQ Type", and the result is an energy with a bound or estimate of its error.

This driver runs the "SCF" submodule, which is "SCF Driver" by default, and
only changes how the modules below it are wired:

- "UQ Type" configures the Integrals plugin's "UQ Initializer", which every
  module that creates uncertain values from doubles uses.
- "ERIs" is used for the four-center ERIs of the Coulomb and exchange
  matrices. By default it attaches the error of "ERI Error Model" to the
  screened ERIs.
- "MO Uncertainty Model" and "Convergence Uncertainty Model" attach
  uncertainty to the orbitals and the energy after convergence.
- The SCF uses "Diagonalizer", which must support the uncertain type.
- Memoization is turned off for every module except "ERI4".

Settings of modules which are not in the call graph, e.g., "ERI Error
Estimate" when "ERI Error Model" is not "Primitive Error Model", have no
effect and are reported when the driver is locked.
)";

} // namespace

using pt           = simde::AOEnergy;
using eri_pt       = simde::ERI4;
using eri_error_pt = simde::Uncertainty<eri_pt>;
using mo_uq_pt     = scf::property_types::MOUncertainty;
using conv_uq_pt   = scf::property_types::ConvergenceUncertainty;

DRIVER_CTOR(QCUPSCF) {
    using integrals_key = pluginplay::QualifiedKey;

    description(desc);
    satisfies_property_type<pt>();

    add_submodule<pt>("SCF")
      .set_description("The SCF calculation which is wired for UQ")
      .set_default("SCF Driver");

    add_submodule<simde::GeneralizedEigenSolve>("Diagonalizer")
      .set_description("Diagonalizes the Fock matrix. Must support the "
                       "uncertain type")
      .set_default("Generalized eigensolve")
      .set_required_target("Loop", "Diagonalizer")
      .set_target("Diagonalization Fock update", "Diagonalizer");

    // -- Where the uncertain type enters ------------------------------------
    const integrals_key uq_initializer{"integrals", "UQ Initializer"};
    add_input<std::string>("UQ Type")
      .set_description("The uncertain type: uncertain, interval, affine, "
                       "or thresholded affine")
      .set_default(std::string("interval"))
      .set_required_target(uq_initializer, "UQ Type");

    // -- The ERIs and their error model -------------------------------------
    const integrals_key uq_eris{"integrals", "UQ Atom Symm Blocked Driver"};
    add_submodule<eri_pt>("ERIs")
      .set_description("The four-center ERIs, with uncertainty attached")
      .set_required_target({"integrals", "Four center J builder"},
                           "Four-center ERI")
      .set_required_target({"integrals", "Four center K builder"},
                           "Four-center ERI");

    add_submodule<eri_error_pt>("ERI Error Model")
      .set_description("Estimates the error of the screened ERIs")
      .set_required_target(uq_eris, "ERI Error");

    add_submodule<eri_pt>("Benchmark ERIs")
      .set_description("The reference ERIs of the analytic error model")
      .set_target({"integrals", "Analytic Error"}, "ERI4s");

    add_input<std::string>("ERI Error Estimate")
      .set_description("The error estimate of the primitive error model")
      .set_default(std::string("Fine"))
      .set_target({"integrals", "Primitive Error Model"}, "Error estimate");

    add_input<std::string>("Mean Type")
      .set_description("How the error is averaged over each block of ERIs")
      .set_default(std::string("none"))
      .set_target(uq_eris, "Mean Type");

    add_input<double>("ERI Threshold")
      .set_description("The screening threshold of the ERIs")
      .set_default(1.0E-6)
      .set_target({"integrals", "ERI4"}, "Threshold");

    // -- Uncertainty attached after convergence -----------------------------
    add_submodule<mo_uq_pt>("MO Uncertainty Model")
      .set_description("Attaches uncertainty to the converged orbitals")
      .set_default("First-order MO uncertainty")
      .set_required_target("Loop", "MO uncertainty");

    add_submodule<conv_uq_pt>("Convergence Uncertainty Model")
      .set_description("Attaches the uncertainty from incomplete convergence "
                       "to the converged energy")
      .set_default("Energy change uncertainty")
      .set_required_target("Loop", "Convergence uncertainty");

    wire_memoization(false, {integrals_key{"integrals", "ERI4"}});
}

DRIVER_RUN(QCUPSCF) {
    const auto& [aos, sys] = pt::unwrap_inputs(inputs);
    const auto e           = submods.at("SCF").run_as<pt>(aos, sys);
    auto rv                = results();
    return pt::wrap_results(rv, e);
}

} // namespace scf::driver
