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
#include <scf/property_types.hpp>

using namespace scf::property_types;

namespace {

using string_list = std::initializer_list<std::string>;

// Checks that the property type T has exactly the given inputs and results
template<typename T>
void test_property_type(string_list input_fields, string_list result_fields) {
    auto inputs = T::inputs();
    REQUIRE(inputs.size() == input_fields.size());
    for(const auto& field : input_fields) REQUIRE(inputs.count(field) == 1);

    auto results = T::results();
    REQUIRE(results.size() == result_fields.size());
    for(const auto& field : result_fields) REQUIRE(results.count(field) == 1);
}

} // namespace

TEST_CASE("MOUncertainty") {
    test_property_type<MOUncertainty>(
      {"Fock matrix", "Overlap matrix", "Orbital energies",
       "Orbital coefficients", "Density change"},
      {"Orbital coefficients"});
}

TEST_CASE("ConvergenceUncertainty") {
    test_property_type<ConvergenceUncertainty>({"Energy", "Energy change"},
                                               {"Energy"});
}
