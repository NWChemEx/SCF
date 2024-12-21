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

template<typename WFType>
using pt =
  simde::eval_braket<WFType, simde::type::electronic_hamiltonian, WFType>;

TEST_CASE("ElectronicEnergy") {
    auto mm  = test_scf::load_modules();
    auto mod = mm.at("Electronic energy");

    using wf_type   = simde::type::rscf_wf;
    using index_set = typename wf_type::orbital_index_set_type;

    wf_type psi(index_set{0}, test_scf::h2_cmos());
    simde::type::many_electrons es{2};

    simde::type::T_e_type T_e(es);

    auto h2_nuclei = test_scf::make_h2<simde::type::nuclei>();
    simde::type::V_en_type V_en(es, h2_nuclei);

    auto rho = test_scf::h2_density();
    simde::type::J_e_type J_e(es, rho);
    simde::type::K_e_type K_e(es, rho);
    simde::type::electronic_hamiltonian H_e(T_e * 2.0 + V_en * 2.0 + J_e * 2.0 -
                                            K_e);
    chemist::braket::BraKet braket(psi, H_e, psi);

    const auto& E_elec = mod.run_as<pt<wf_type>>(braket);
    REQUIRE_THAT(E_elec, WithinAbs(-1.90066758625308307, 1E-6));
}