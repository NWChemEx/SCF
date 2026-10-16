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
// property type; wiring never runs it.
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

// Stand-ins for the plugins SCF gets default submodules from
class FakeNUX : public pluginplay::Plugin {
public:
    FakeNUX() {
        set_name("nux");
        add_module<StubProvider<ham_pt>>("Born-Oppenheimer approximation");
    }
};

class FakeIntegrals : public pluginplay::Plugin {
public:
    FakeIntegrals() {
        set_name("integrals");
        add_module<StubProvider<simde::aos_op_base_aos>>("AO integral driver");
        add_module<StubProvider<simde::aos_s_e_aos>>("Overlap");
    }
};

class FakeChemCache : public pluginplay::Plugin {
public:
    FakeChemCache() {
        set_name("chemcache");
        add_module<StubProvider<simde::InitialDensity>>("sto-3g SAD density");
    }
};

// (module, submodule, provider plugin, provider) for every default SCF records
const std::vector<std::array<std::string, 4>> wiring{
  {"SCF Driver", "Hamiltonian", "nux", "Born-Oppenheimer approximation"},
  {"SCF integral driver", "Fundamental matrices", "integrals",
   "AO integral driver"},
  {"Diagonalization Fock update", "Overlap matrix builder", "integrals",
   "Overlap"},
  {"Loop", "Overlap matrix builder", "integrals", "Overlap"},
  {"SAD guess", "SAD Density", "chemcache", "sto-3g SAD density"}};

} // namespace

TEST_CASE("SCF plugin") {
    scf::SCF plugin;
    const auto& modules = plugin.get_modules();

    SECTION("Identity") {
        REQUIRE(plugin.name() == "scf");
        using provider_list = scf::SCF::provider_list_type;
        provider_list corr{"nux", "integrals", "chemcache"};
        REQUIRE(plugin.preferred_providers() == corr);
    }

    SECTION("Defaults set across subdirectories") {
        REQUIRE(modules.has_submod("SCF Driver", "Optimizer"));
        REQUIRE(modules.has_submod("Generalized eigensolve", "Eigen Solve"));
        REQUIRE(modules.has_submod("Density2Grid", "AOs on a grid"));
        REQUIRE(modules.has_submod("Loop", "Charge-charge"));
    }

    pluginplay::ModuleManager mm;
    mm.load_plugin<scf::SCF>();

    SECTION("No providers loaded") {
        mm.wire_defaults();
        for(const auto& [module, submod, plugin_name, provider] : wiring) {
            REQUIRE_FALSE(mm.at(module).submods().at(submod).has_module());
        }
    }

    SECTION("All providers loaded") {
        mm.load_plugin<FakeNUX>();
        mm.load_plugin<FakeIntegrals>();
        mm.load_plugin<FakeChemCache>();
        mm.wire_defaults();
        for(const auto& [module, submod, plugin_name, provider] : wiring) {
            const auto& request = mm.at(module).submods().at(submod);
            REQUIRE(request.has_module());
            REQUIRE(request.uuid() == mm.at({plugin_name, provider}).uuid());
        }
    }
}
