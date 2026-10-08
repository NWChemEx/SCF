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

#include "../../../test_scf.hpp"
#include <cmath>
#include <scf/xc/native/svwn3.hpp>

using namespace scf::xc::native;

TEST_CASE("svwn3 zero density") {
    SECTION("double") {
        REQUIRE(svwn3::energy_density(0.0) == 0.0);
        REQUIRE(svwn3::potential(0.0) == 0.0);

        // The limit is approached continuously
        REQUIRE(std::isfinite(svwn3::energy_density(1.0e-30)));
        REQUIRE(std::abs(svwn3::energy_density(1.0e-30)) < 1.0e-30);
        REQUIRE(std::abs(svwn3::potential(1.0e-30)) < 1.0e-9);

        REQUIRE_THROWS_AS(svwn3::energy_density(-1.0e-10), std::domain_error);
        REQUIRE_THROWS_AS(svwn3::potential(-1.0e-10), std::domain_error);
    }

    SECTION("float") {
        REQUIRE(svwn3::energy_density(0.0f) == 0.0f);
        REQUIRE(svwn3::potential(0.0f) == 0.0f);
        REQUIRE_THROWS_AS(svwn3::potential(-1.0f), std::domain_error);
    }

#ifdef ENABLE_SIGMA
    using udouble = tensorwrapper::types::udouble;
    using idouble = tensorwrapper::types::idouble;

    SECTION("udouble") {
        REQUIRE(svwn3::energy_density(udouble(0.0)).mean() == 0.0);
        REQUIRE(svwn3::potential(udouble(0.0)).mean() == 0.0);

        udouble zero_mean(0.0, 1.0e-6);
        REQUIRE_THROWS_AS(svwn3::energy_density(zero_mean), std::domain_error);
        REQUIRE_THROWS_AS(svwn3::potential(zero_mean), std::domain_error);
    }

    SECTION("idouble") {
        REQUIRE(svwn3::energy_density(idouble(0.0)) == idouble(0.0));
        REQUIRE(svwn3::potential(idouble(0.0)) == idouble(0.0));

        idouble spans_zero(0.0, 1.0e-6);
        REQUIRE_THROWS_AS(svwn3::energy_density(spans_zero), std::domain_error);
        REQUIRE_THROWS_AS(svwn3::potential(spans_zero), std::domain_error);
    }
#endif
}
