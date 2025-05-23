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

using pt = simde::aos_f_e_aos;

using simde::type::t_e_type;
using simde::type::v_en_type;

TEMPLATE_LIST_TEST_CASE("Fock Matrix Builder", "", test_scf::float_types) {
    using float_type = TestType;
    auto mm          = test_scf::load_modules<float_type>();

    auto& mod = mm.at("Fock Matrix Builder");
    auto aos  = test_scf::h2_aos();

    using tensorwrapper::operations::approximately_equal;
    tensorwrapper::allocator::Eigen<float_type> alloc(mm.get_runtime());
    tensorwrapper::shape::Smooth shape_corr{2, 2};
    tensorwrapper::layout::Physical l(shape_corr);
    auto pcorr = alloc.allocate(l);

    SECTION("No J or K") {
        auto h2 = test_scf::make_h2<simde::type::nuclei>();
        simde::type::electron e;

        simde::type::fock f_e;
        f_e.emplace_back(1.0, std::make_unique<t_e_type>(e));
        f_e.emplace_back(1.0, std::make_unique<v_en_type>(e, h2));
        chemist::braket::BraKet f_mn(aos, f_e, aos);
        const auto& F = mod.template run_as<pt>(f_mn);

        pcorr->set_elem({0, 0}, -1.120958);
        pcorr->set_elem({0, 1}, -0.959374);
        pcorr->set_elem({1, 0}, -0.959374);
        pcorr->set_elem({1, 1}, -1.120958);

        tensorwrapper::Tensor corr(shape_corr, std::move(pcorr));

        REQUIRE(approximately_equal(F, corr, 1E-6));
    }

    SECTION("With J and K") {
        auto f_e = test_scf::h2_fock<simde::type::electron, float_type>();
        chemist::braket::BraKet f_mn(aos, f_e, aos);
        const auto& F = mod.template run_as<pt>(f_mn);

        pcorr->set_elem({0, 0}, -0.319459);
        pcorr->set_elem({0, 1}, -0.571781);
        pcorr->set_elem({1, 0}, -0.571781);
        pcorr->set_elem({1, 1}, -0.319459);

        tensorwrapper::Tensor corr(shape_corr, std::move(pcorr));

        REQUIRE(approximately_equal(F, corr, 1E-6));
    }
}
