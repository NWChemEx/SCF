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

#include "../../test_scf.hpp"

using pt = simde::aos_rho_e_aos<simde::type::cmos>;

TEMPLATE_LIST_TEST_CASE("Density Matrix Builder", "", test_scf::float_types) {
    using float_type = TestType;
    pluginplay::ModuleManager mm;
    scf::load_modules(mm);
    auto& mod = mm.at("Density matrix builder");
    auto aos  = test_scf::h2_aos();
    auto cmos = test_scf::h2_cmos<float_type>();
    std::vector<int> occs{1, 0};
    simde::type::rho_e<simde::type::cmos> rho_hat(cmos, occs);

    chemist::braket::BraKet p_mn(aos, rho_hat, aos);
    const auto& P = mod.run_as<pt>(p_mn);
    tensorwrapper::allocator::Eigen<float_type> alloc(mm.get_runtime());
    tensorwrapper::shape::Smooth corr_shape{2, 2};
    tensorwrapper::layout::Physical l(corr_shape);
    auto corr_buffer = alloc.construct(l, 0.31980835);
    tensorwrapper::Tensor corr(corr_shape, std::move(corr_buffer));

    using tensorwrapper::operations::approximately_equal;
    REQUIRE(approximately_equal(P, corr, 1E-6));
}
