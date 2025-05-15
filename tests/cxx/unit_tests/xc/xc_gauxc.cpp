/*
 * Copyright 2022 NWChemEx-Project
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
#include <iostream>
#include <pluginplay/pluginplay.hpp>

using namespace scf;

TEST_CASE("XC") {
    // pluginplay::ModuleManager mm;
    // scf::load_modules(mm);

    // auto mod = mm.at("XC");

    // auto aos = test_scf::h2_aos();

    // std::vector ao_bases{aos, aos};
    // auto ot_shape  = testing::make_one_tile_shape(ao_bases);
    // auto shape_mod = pluginplay::make_facade<simde::IntegralShape>(ot_shape);
    // mod.change_submod("Tensor Shape", shape_mod);

    // // Build density and j operator
    // auto occ = get_space(property::occupied, name, bs);

    // simde::type::tensor rho;
    // rho("mu,nu") = occ.C()("mu, i") * occ.C()("nu, i");
    // simde::type::el_density P(rho, aos);
    // simde::type::el e;
    // using xc_functional_type =
    //   chemist::operators::ExchangeCorrelationFunctional;
    // auto func = xc_functional_type::PBE0;
    // simde::type::el_scf_xc vxc(func, e, P);

    // SECTION("As Energy Module") {
    //     auto EXC = mod.run_as<exc_pt>(aos, vxc, aos);
    //     REQUIRE(EXC == Approx(-7.03384).margin(1.0e-5));
    // }

    // SECTION("As Potential Module") {
    //     auto corr_vxc = get_ao_data(
    //       name, {bs, bs}, property::pbe0_exchange_correlation_potential);
    //     auto VXC = mod.run_as<vxc_pt>(aos, vxc, aos);
    //     REQUIRE(tensorwrapper::tensor::allclose(VXC, corr_vxc));
    // }
}

TEST_CASE("snLinK") {
    // pluginplay::ModuleManager mm;
    // scf::load_modules(mm);

    // auto mod = mm.at("snLinK");

    // const auto name = molecule::h2o;
    // const auto bs   = basis_set::sto3g;
    // auto mol        = get_molecule(name);
    // auto aos        = get_bases(name, bs);

    // // Build density and j operator
    // auto occ = get_space(property::occupied, name, bs);

    // simde::type::tensor rho;
    // rho("mu,nu") = occ.C()("mu, i") * occ.C()("nu, i");
    // simde::type::el_density P(rho, aos);
    // simde::type::el e;
    // simde::type::el_scf_k k(e, P);

    // auto corr_k = get_ao_data(name, {bs, bs}, property::exchange);
    // auto K      = mod.run_as<k_pt>(aos, k, aos);
    // REQUIRE(tensorwrapper::tensor::allclose(K, corr_k));
}
