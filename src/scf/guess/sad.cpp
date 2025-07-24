/*
 * Copyright 2025 NWChemEx-Project
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

#include "guess.hpp"

namespace scf::guess {
namespace {
const auto desc = R"(
SAD Guess
---------

TODO: Write me!!!
)";
}

using rscf_wf        = simde::type::rscf_wf;
using density_t      = simde::type::decomposable_e_density;
using pt             = simde::InitialGuess<rscf_wf>;
using fock_op_pt     = simde::FockOperator<density_t>;
using update_pt      = simde::UpdateGuess<rscf_wf>;
using initial_rho_pt = simde::InitialDensity;

using simde::type::tensor;

// TODO: move to chemist?
struct NElectronCounter : public chemist::qm_operator::OperatorVisitor {
    NElectronCounter() : chemist::qm_operator::OperatorVisitor(false) {}

    void run(const simde::type::T_e_type& T_e) { set_n(T_e.particle().size()); }

    void run(const simde::type::V_en_type& V_en) {
        set_n(V_en.lhs_particle().size());
    }

    void run(const simde::type::V_ee_type& V_ee) {
        set_n(V_ee.lhs_particle().size());
        set_n(V_ee.rhs_particle().size());
    }

    void set_n(unsigned int n) {
        if(n_electrons == 0)
            n_electrons = n;
        else if(n_electrons != n) {
            throw std::runtime_error("Deduced a different number of electrons");
        }
    }

    unsigned int n_electrons = 0;
};

MODULE_CTOR(SAD) {
    description(desc);
    satisfies_property_type<pt>();
    add_submodule<fock_op_pt>("Build Fock operator");
    add_submodule<update_pt>("Guess updater");
    add_submodule<initial_rho_pt>("SAD Density");
}

MODULE_RUN(SAD) {
    const auto&& [H, aos] = pt::unwrap_inputs(inputs);

    // Step 1: Build Fock Operator with zero density
    auto& initial_rho_mod = submods.at("SAD Density");
    const auto& rho       = initial_rho_mod.run_as<initial_rho_pt>(H);
    auto& fock_op_mod     = submods.at("Build Fock operator");
    const auto& f         = fock_op_mod.run_as<fock_op_pt>(H, rho);

    // Step 2: Get number of electrons and occupations
    simde::type::cmos cmos(tensor{}, aos, tensor{});
    NElectronCounter visitor;
    H.visit(visitor);
    auto n_electrons = visitor.n_electrons;
    if(n_electrons % 2 != 0)
        throw std::runtime_error("Assumed even number of electrons");

    typename rscf_wf::orbital_index_set_type occs;
    using value_type = typename rscf_wf::orbital_index_set_type::value_type;
    for(value_type i = 0; i < n_electrons / 2; ++i) occs.insert(i);

    rscf_wf zero_guess(occs, cmos);
    auto& update_mod = submods.at("Guess updater");
    const auto& Psi0 = update_mod.run_as<update_pt>(f, zero_guess);

    auto rv = results();
    return pt::wrap_results(rv, Psi0);
}

} // namespace scf::guess