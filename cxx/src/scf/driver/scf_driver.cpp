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

#include "../xc/xc.hpp"
#include "driver.hpp"
#include <algorithm>
#include <string>

namespace scf::driver {
namespace {

const auto desc = R"(
SCF Driver
----------

Computes the restricted SCF energy of a chemical system: "Hamiltonian" builds
the Hamiltonian, "Guess" the initial wavefunction, and "Optimizer" converges
it.

"DFT Functional" chooses between Hartree-Fock and Kohn-Sham DFT. With "None",
the default, the modules are used as they are wired, i.e., Hartree-Fock by
default. With the name of an XC functional, e.g., "SVWN3", this driver wires
the Kohn-Sham operators, with that functional, into "Loop" and "Core guess".
Which functionals are available depends on how SCF was built; the input's
check lists them.
)";

/// The value of "DFT Functional" which leaves the wiring as it is
const std::string no_functional = "None";

/// The values "DFT Functional" may have, separated by commas
std::string dft_functional_names() {
    std::string rv = no_functional;
    for(const auto& entry : xc::xc_dispatch_table()) rv += ", " + entry.key;
    return rv;
}

/// Is @p name a value "DFT Functional" may have?
bool is_dft_functional(const std::string& name) {
    if(name == no_functional) return true;
    const auto table = xc::xc_dispatch_table();
    return std::any_of(table.begin(), table.end(),
                       [&](const auto& entry) { return entry.key == name; });
}

} // namespace

using simde::type::hamiltonian;

using pt = simde::AOEnergy;

using ham_pt = simde::Convert<hamiltonian, simde::type::chemical_system>;

template<typename WfType>
using guess_pt = simde::InitialGuess<WfType>;

template<typename WfType>
using braket_t = chemist::braket::BraKet<WfType, hamiltonian, WfType>;

template<typename WfType>
using egy_pt = simde::EvaluateBraKet<braket_t<WfType>>;

template<typename WfType>
using opt_pt = simde::Optimize<egy_pt<WfType>, WfType>;

DRIVER_CTOR(SCFDriver) {
    using wf_type = simde::type::rscf_wf;

    description(desc);
    satisfies_property_type<pt>();
    add_submodule<ham_pt>("Hamiltonian");
    add_submodule<guess_pt<wf_type>>("Guess");
    add_submodule<opt_pt<wf_type>>("Optimizer");

    const std::string dft_key = "DFT Functional";
    pluginplay::ModuleInput::validity_check<std::string> check{
      is_dft_functional};
    add_input<std::string>(dft_key)
      .set_description("The XC functional, or \"None\" for Hartree-Fock")
      .set_default(no_functional)
      .add_check(check, "Is one of: " + dft_functional_names());

    // One branch per functional. There is no branch for "None", so that
    // Hartree-Fock leaves the wiring as it is.
    const auto RKS_op = "Restricted Kohn-Sham Op";
    const auto rks_op = "Restricted One-Electron Kohn-Sham Op";
    pluginplay::ConditionalWiring* chain = nullptr;
    for(const auto& entry : xc::xc_dispatch_table()) {
        chain = chain ? &chain->otherwise(dft_key, entry.key) :
                        &when(dft_key, entry.key);
        chain->set_input(RKS_op, "XC Potential", entry.functional)
          .set_input(rks_op, "XC Potential", entry.functional)
          .set_submod("Loop", "One-electron Fock operator", rks_op)
          .set_submod("Loop", "Fock operator", RKS_op)
          .set_submod("Core guess", "Build Fock operator", rks_op)
          .required();
    }
}

DRIVER_RUN(SCFDriver) {
    using wf_type = simde::type::rscf_wf;

    const auto& [ao_params, sys] = pt::unwrap_inputs(inputs);
    simde::type::aos aos(ao_params);

    auto& ham_mod = submods.at("Hamiltonian");
    const auto& H = ham_mod.run_as<ham_pt>(sys);

    auto& guess_mod  = submods.at("Guess");
    const auto& Psi0 = guess_mod.run_as<guess_pt<wf_type>>(H, aos);

    auto& opt_mod = submods.at("Optimizer");
    braket_t<wf_type> H_00(Psi0, H, Psi0);
    const auto&& [e, Psi] = opt_mod.run_as<opt_pt<wf_type>>(H_00, Psi0);

    auto rv = results();
    return pt::wrap_results(rv, e);
}

} // namespace scf::driver
