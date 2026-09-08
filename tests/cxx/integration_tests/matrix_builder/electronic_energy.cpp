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

TEMPLATE_LIST_TEST_CASE("ElectronicEnergy", "", test_scf::float_types) {
    using float_type = TestType;
    auto mm          = test_scf::load_modules<float_type>();
    auto mod         = mm.at("Electronic energy");

    using wf_type   = simde::type::rscf_wf;
    using index_set = typename wf_type::orbital_index_set_type;

    using tensorwrapper::buffer::make_contiguous;
    tensorwrapper::shape::Smooth shape_corr{};
    auto pcorr = make_contiguous<float_type>(shape_corr);
    using tensorwrapper::operations::approximately_equal;

    wf_type psi(index_set{0}, test_scf::h2_cmos<float_type>());
    simde::type::many_electrons es{2};

    simde::type::T_e_type T_e(es);

    auto h2_nuclei = test_scf::make_h2<simde::type::nuclei>();
    simde::type::V_en_type V_en(es, h2_nuclei);

    auto rho = test_scf::h2_density<float_type>();
    simde::type::J_e_type J_e(es, rho);
    SECTION("RHF") {
        simde::type::K_e_type K_e(es, rho);
        simde::type::electronic_hamiltonian H_e(T_e * 2.0 + V_en * 2.0 +
                                                J_e * 2.0 - K_e);
        chemist::braket::BraKet braket(psi, H_e, psi);

        const auto& E_elec = mod.template run_as<pt<wf_type>>(braket);

        pcorr.set_elem({}, float_type{-1.90066758625308307});
        tensorwrapper::Tensor corr(shape_corr, std::move(pcorr));
        REQUIRE(approximately_equal(corr, E_elec, 1E-6));
    }
    SECTION("RKS (GauXC)") {
#ifdef BUILD_GAUXC
        if constexpr(std::is_same_v<float_type, double>) {
            auto func = chemist::qm_operator::xc_functional::PBE;
            simde::type::XC_e_type XC_e(func, es, rho);
            simde::type::electronic_hamiltonian H_e(T_e * 2.0 + V_en * 2.0 +
                                                    J_e * 2.0 + XC_e);
            chemist::braket::BraKet braket(psi, H_e, psi);
            const auto& E_elec = mod.template run_as<pt<wf_type>>(braket);
            tensorwrapper::Tensor corr(-1.90692);
            REQUIRE(approximately_equal(corr, E_elec, 1E-5));
        }
#endif
    }
    SECTION("RKS (LibXC)") {
#ifdef BUILD_LIBXC
        if constexpr(std::is_same_v<float_type, double>) {
            // Wired explicitly so this stays a LibXC test even when
            mm.change_submod("Electronic energy", "XC Energy", "LibXC Energy");
            // mod above is a copy taken before the rewiring, so re-fetch.
            auto& xc_mod = mm.at("Electronic energy");
            auto func    = chemist::qm_operator::xc_functional::SVWN3;
            simde::type::XC_e_type XC_e(func, es, rho);
            simde::type::electronic_hamiltonian H_e(T_e * 2.0 + V_en * 2.0 +
                                                    J_e * 2.0 + XC_e);
            chemist::braket::BraKet braket(psi, H_e, psi);
            const auto& E_elec = xc_mod.template run_as<pt<wf_type>>(braket);

            // Psi4 1.11, STO-3G, at this same (fixed) density:
            // 2*Tr(PT) + 2*Tr(PV) + 2*Tr(PJ) + E_xc, with E_xc from
            // {x: LDA_X, c: LDA_C_VWN_3} on a 150x194 unpruned Becke grid.
            // The same script reproduces this test's RHF reference above and
            // the LibXCEnergy/LibXCPotential references in
            // tests/cxx/integration_tests/xc/.
            tensorwrapper::Tensor corr(-1.87456183);
            REQUIRE(approximately_equal(corr, E_elec, 1E-5));
        }
#endif
    }
}
