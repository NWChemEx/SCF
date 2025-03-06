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

#include "integration_tests.hpp"

using pt = simde::charge_charge_interaction;
using Catch::Matchers::WithinAbs;

TEMPLATE_LIST_TEST_CASE("CoulombsLaw", "", test_scf::float_types) {
    using float_type = double;
    auto mm          = test_scf::load_modules<float_type>();
    auto& mod        = mm.at("Coulomb's Law");

    tensorwrapper::allocator::Eigen<float_type> alloc(mm.get_runtime());
    tensorwrapper::shape::Smooth shape_corr{};
    auto pcorr = alloc.allocate(tensorwrapper::layout::Physical(shape_corr));
    using tensorwrapper::operations::approximately_equal;

    auto h2_nuclei = test_scf::make_h2<simde::type::nuclei>();
    // TODO: Conversions are missing in Chemist. Use those when they're in place
    simde::type::charges qs;
    for(const auto& nucleus : h2_nuclei)
        qs.push_back(nucleus.as_point_charge());

    auto e_nuclear = mod.run_as<pt>(qs, qs);

    pcorr->at() = 0.71510297482837526;
    tensorwrapper::Tensor corr(shape_corr, std::move(pcorr));
    REQUIRE(approximately_equal(corr, e_nuclear, 1E-6));
}