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

/** @file property_types.hpp
 *
 *  The property types defined by SCF. They are candidates for SimDE once
 *  their APIs have settled.
 */
#pragma once
#include <pluginplay/pluginplay.hpp>
#include <simde/types.hpp>

/** @namespace scf::property_types
 *
 *  @brief The namespace for property types defined by SCF
 */
namespace scf::property_types {

/** @brief Attaches uncertainty to converged molecular orbitals.
 *
 *  Called by an SCF once it has converged. The orbital energies and
 *  coefficients are the eigenpairs of the converged Fock matrix in the
 *  metric of the overlap matrix. The result is the coefficients with
 *  uncertainty attached, *e.g.*, from the uncertainty in the Fock matrix and
 *  from incomplete convergence.
 */
DECLARE_PROPERTY_TYPE(MOUncertainty);

PROPERTY_TYPE_INPUTS(MOUncertainty) {
    using tensor_type = const simde::type::tensor&;
    auto rv           = pluginplay::declare_input()
                .add_field<tensor_type>("Fock matrix")
                .template add_field<tensor_type>("Overlap matrix")
                .template add_field<tensor_type>("Orbital energies")
                .template add_field<tensor_type>("Orbital coefficients")
                .template add_field<tensor_type>("Density change");
    rv["Fock matrix"].set_description("The converged Fock matrix");
    rv["Overlap matrix"].set_description("The AO overlap matrix");
    rv["Orbital energies"].set_description(
      "The eigenvalues of the converged Fock matrix");
    rv["Orbital coefficients"].set_description(
      "The eigenvectors of the converged Fock matrix, as columns");
    rv["Density change"].set_description(
      "The change in the density matrix in the final iteration");
    return rv;
}

PROPERTY_TYPE_RESULTS(MOUncertainty) {
    auto rv = pluginplay::declare_result().add_field<simde::type::tensor>(
      "Orbital coefficients");
    rv["Orbital coefficients"].set_description(
      "The orbital coefficients with uncertainty attached");
    return rv;
}

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
