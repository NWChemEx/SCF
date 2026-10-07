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

#include "../test_scf.hpp"
#include <array>
#include <scf/scf.hpp>
#include <simde/simde.hpp>
#include <string>
#include <vector>

namespace {

// Stands in for a module from another plugin. It only needs to satisfy the
// property type; set_defaults never runs it.
template<typename PropertyType>
struct StubProvider : pluginplay::ModuleBase {
    StubProvider() : pluginplay::ModuleBase(this) {
        satisfies_property_type<PropertyType>();
    }

private:
    result_map_type run_(input_map_type, submodule_map_type) const override {
        return results();
    }
};

using ham_pt =
  simde::Convert<simde::type::hamiltonian, simde::type::chemical_system>;

// (module, submodule, provider) for every submodule set_defaults sets
const std::vector<std::array<std::string, 3>> wiring{
  {"SCF Driver", "Hamiltonian", "Born-Oppenheimer approximation"},
  {"SCF integral driver", "Fundamental matrices", "AO integral driver"},
  {"Diagonalization Fock update", "Overlap matrix builder", "Overlap"},
  {"Loop", "Overlap matrix builder", "Overlap"},
  {"SAD guess", "SAD Density", "sto-3g SAD density"}};

} // namespace

TEST_CASE("set_defaults") {
    pluginplay::ModuleManager mm;
    scf::load_modules(mm);

    SECTION("No providers loaded") {
        scf::set_defaults(mm);
        for(const auto& [module, submod, provider] : wiring) {
            REQUIRE_FALSE(mm.at(module).submods().at(submod).has_module());
        }
    }

    SECTION("All providers loaded") {
        mm.add_module<StubProvider<ham_pt>>("Born-Oppenheimer approximation");
        mm.add_module<StubProvider<simde::aos_op_base_aos>>(
          "AO integral driver");
        mm.add_module<StubProvider<simde::aos_s_e_aos>>("Overlap");
        mm.add_module<StubProvider<simde::InitialDensity>>(
          "sto-3g SAD density");

        scf::set_defaults(mm);
        for(const auto& [module, submod, provider] : wiring) {
            const auto& request = mm.at(module).submods().at(submod);
            REQUIRE(request.has_module());
            REQUIRE(request.uuid() == mm.at(provider).uuid());
        }
    }
}
