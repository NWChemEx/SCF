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

#include "xc.hpp"
#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

namespace scf::xc {
namespace {

const auto desc = R"(
Exchange-Correlation Dispatcher
===============================

Forwards an exchange-correlation (XC) calculation to the module for the XC
functional of the operator. Each supported functional has a submodule, named
after the functional (e.g., "SVWN3" or "PBE"), and the inputs are passed to it
unchanged. Which functionals have a submodule depends on the XC backends SCF
was built with.
)";

} // namespace

std::vector<XCDispatchEntry> xc_dispatch_table() {
    using chemist::qm_operator::xc_functional;

    // The native SVWN3 modules are always built and support the UQ types
    std::vector<XCDispatchEntry> table{
      {xc_functional::SVWN3, "SVWN3", "SVWN3 Potential", "SVWN3 Energy"}};
#ifdef BUILD_GAUXC
    const std::string v_key = "GauXC XC Potential";
    const std::string e_key = "GauXC XC Energy";
    table.push_back({xc_functional::SVWN5, "SVWN5", v_key, e_key});
    table.push_back({xc_functional::BLYP, "BLYP", v_key, e_key});
    table.push_back({xc_functional::B3LYP, "B3LYP", v_key, e_key});
    table.push_back({xc_functional::PBE, "PBE", v_key, e_key});
    table.push_back({xc_functional::revPBE, "revPBE", v_key, e_key});
    table.push_back({xc_functional::PBE0, "PBE0", v_key, e_key});
#elif defined(BUILD_LIBXC)
    table.push_back(
      {xc_functional::SVWN5, "SVWN5", "LibXC Potential", "LibXC Energy"});
#endif
    return table;
}

template<typename PT>
TEMPLATED_MODULE_CTOR(XCDispatcher, PT) {
    satisfies_property_type<PT>();
    description(desc);

    for(const auto& entry : xc_dispatch_table())
        add_submodule<PT>(entry.key).set_description("Evaluates the " +
                                                     entry.key + " functional");
}

template<typename PT>
TEMPLATED_MODULE_RUN(XCDispatcher, PT) {
    const auto& [braket] = PT::unwrap_inputs(inputs);
    const auto func      = braket.op().get_functional_name();

    const auto table = xc_dispatch_table();
    auto is_func     = [func](const auto& entry) {
        return entry.functional == func;
    };
    auto it = std::find_if(table.begin(), table.end(), is_func);
    if(it == table.end()) {
        std::string available;
        for(const auto& entry : table)
            available += (available.empty() ? "" : ", ") + entry.key;
        throw std::runtime_error("XC Dispatcher: the requested functional is "
                                 "not supported by this build of SCF. "
                                 "Supported functionals: " +
                                 available);
    }

    const auto& value = submods.at(it->key).template run_as<PT>(braket);

    auto rv = results();
    return PT::wrap_results(rv, value);
}

template struct XCDispatcher<simde::aos_xc_e_aos>;
template struct XCDispatcher<xc_energy_pt>;

} // namespace scf::xc
