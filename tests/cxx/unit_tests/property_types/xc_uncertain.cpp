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

#include "test_property_types.hpp"
#include <scf/property_types/xc_uncertain.hpp>

using namespace scf::property_types;

TEST_CASE("XCUncertain") {
    SECTION("xc_energy_uncertainty") {
        test_scf::test_property_type<xc_energy_uncertainty>(
          {"BraKet", "Grid", "tensor representation"},
          {"tensor representation"});
    }

    SECTION("xc_potential_uncertainty") {
        test_scf::test_property_type<xc_potential_uncertainty>(
          {"BraKet", "Grid", "tensor representation"},
          {"tensor representation"});
    }
}
