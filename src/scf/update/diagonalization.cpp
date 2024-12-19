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

#include "update.hpp"

namespace scf::update {

using rscf_wf         = simde::type::rscf_wf;
using fock_matrix_pt  = simde::aos_f_e_aos;
using pt              = simde::UpdateGuess<rscf_wf>;
using diagonalizer_pt = simde::GeneralizedEigenSolve;
using s_pt            = simde::aos_s_e_aos;

const auto desc = R"(
)";

MODULE_CTOR(Diagonalization) {
    description(desc);
    satisfies_property_type<pt>();
    add_submodule<fock_matrix_pt>("Fock matrix builder");
    add_submodule<diagonalizer_pt>("Diagonalizer");
    add_submodule<s_pt>("Overlap matrix builder");
}

MODULE_RUN(Diagonalization) {
    const auto&& [f, old_guess] = pt::unwrap_inputs(inputs);
    const auto& aos             = old_guess.orbitals().from_space();

    // Comput F
    chemist::braket::BraKet f_mn(aos, f, aos);
    auto& fock_matrix_mod = submods.at("Fock matrix builder");
    const auto& f_matrix  = fock_matrix_mod.run_as<fock_matrix_pt>(f_mn);

    // Compute S
    chemist::braket::BraKet s_mn(aos, simde::type::s_e_type{}, aos);
    auto& s_mod          = submods.at("Overlap matrix builder");
    const auto& s_matrix = s_mod.run_as<s_pt>(s_mn);

    // Diagonalize
    auto& diagonalizer_mod = submods.at("Diagonalizer");
    const auto&& [evalues, evectors] =
      diagonalizer_mod.run_as<diagonalizer_pt>(f_matrix, s_matrix);

    // Create new guess
    simde::type::cmos cmos(evalues, aos, evectors);
    rscf_wf new_guess(old_guess.orbital_indices(), cmos);

    auto rv = results();
    return pt::wrap_results(rv, new_guess);
}

} // namespace scf::update