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

#include "../../../test_scf.hpp"

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