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

#include "../integration_tests.hpp"

using Catch::Matchers::WithinAbs;

using pt = simde::AOEnergy;

TEST_CASE("SCFDriver") {
    using float_type = double;
    auto mm          = test_scf::load_modules<float_type>();
    auto h2          = test_scf::make_h2<simde::type::chemical_system>();
    auto aos         = test_scf::h2_aos().ao_basis_set();

    const auto e = mm.run_as<pt>("SCF Driver", aos, h2);
    REQUIRE_THAT(e, WithinAbs(-1.1167592336, 1E-6));
}