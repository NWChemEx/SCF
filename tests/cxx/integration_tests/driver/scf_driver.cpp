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

using pt = simde::AOEnergy;
using tensorwrapper::operations::approximately_equal;

using types = std::tuple<double, tensorwrapper::types::udouble,
                         tensorwrapper::types::interval_type<double>,
                         tensorwrapper::types::thresholded_affine_type<double>>;

TEMPLATE_LIST_TEST_CASE("SCFDriver", "", types) {
    using float_type = TestType;
    auto mm          = test_scf::load_modules<float_type>();
    using tensorwrapper::buffer::make_contiguous;
    tensorwrapper::shape::Smooth shape_corr{};
    auto pcorr = make_contiguous<float_type>(shape_corr);

    SECTION("He") {
        auto he  = test_scf::make_he<simde::type::chemical_system>();
        auto aos = test_scf::he_aos().ao_basis_set();

        // Correct energy is validated against Psi4
        pcorr.set_elem({}, float_type{-2.8077839566141960});
        simde::type::tensor corr(shape_corr, std::move(pcorr));
        const auto e = mm.template run_as<pt>("SCF Driver", aos, he);
        REQUIRE(approximately_equal(corr, e, 1E-6));
    }

    SECTION("H2") {
        auto h2  = test_scf::make_h2<simde::type::chemical_system>();
        auto aos = test_scf::h2_aos().ao_basis_set();

        SECTION("SCF") {
            // Correct energy is validated against Psi4
            pcorr.set_elem({}, float_type{-1.1167592336});
            simde::type::tensor corr(shape_corr, std::move(pcorr));
            const auto e = mm.template run_as<pt>("SCF Driver", aos, h2);
            REQUIRE(approximately_equal(corr, e, 1E-6));
        }
        SECTION("DFT (GauXC)") {
#ifdef BUILD_GAUXC
            if constexpr(!tensorwrapper::types::is_uq_type_v<float_type>) {
                auto func         = chemist::qm_operator::xc_functional::PBE;
                const auto RKS_op = "Restricted Kohn-Sham Op";
                const auto rks_op = "Restricted One-Electron Kohn-Sham Op";
                mm.change_input(RKS_op, "XC Potential", func);
                mm.change_input(rks_op, "XC Potential", func);
                mm.change_submod("Loop", "One-electron Fock operator", rks_op);
                mm.change_submod("Loop", "Fock operator", RKS_op);
                mm.change_submod("Core guess", "Build Fock Operator", rks_op);
                const auto e = mm.template run_as<pt>("SCF Driver", aos, h2);
                pcorr.set_elem({}, float_type{-1.15207});
                simde::type::tensor corr(shape_corr, std::move(pcorr));
                REQUIRE(approximately_equal(corr, e, 1E-5));
            }
#endif
        }
        SECTION("DFT (LibXC)") {
#ifdef BUILD_LIBXC

            // LibXC's kernels take a raw double* out of the buffer and
            // NormalizeKernel throws for anything else, so the uncertain
            // types are skipped exactly as they are for GauXC above.
            if constexpr(!tensorwrapper::types::is_uq_type_v<float_type>) {
                auto func         = chemist::qm_operator::xc_functional::SVWN3;
                const auto RKS_op = "Restricted Kohn-Sham Op";
                const auto rks_op = "Restricted One-Electron Kohn-Sham Op";
                mm.change_submod("SCF integral driver", "XC Potential",
                                 "LibXC Potential");
                mm.change_submod("Electronic energy", "XC Energy",
                                 "LibXC Energy");
                mm.change_input(RKS_op, "XC Potential", func);
                mm.change_input(rks_op, "XC Potential", func);
                mm.change_submod("Loop", "One-electron Fock operator", rks_op);
                mm.change_submod("Loop", "Fock operator", RKS_op);
                mm.change_submod("Core guess", "Build Fock Operator", rks_op);
                const auto e = mm.template run_as<pt>("SCF Driver", aos, h2);

                // Validated against Psi4 1.11: RKS/STO-3G, scf_type pk,
                // dft_functional = {x: LDA_X, c: LDA_C_VWN_3} (spelled out
                // because Psi4's "SVWN" alias is Slater+VWN3RPA, a different
                // parameterization than the XC_LDA_C_VWN_3 that
                // to_libxc_codes selects). Grid-converged: 75x302, 99x590,
                // 150x194, and 200x974 unpruned Becke grids all agree to
                // 1e-10.
                pcorr.set_elem({}, float_type{-1.121206107});
                simde::type::tensor corr(shape_corr, std::move(pcorr));
                REQUIRE(approximately_equal(corr, e, 1E-5));
            }
#endif
        }
    }

    SECTION("H2 Dimer") {
        // Correct energy is validated against Psi4
        simde::type::nucleus h0("H", 1ul, 1836.15, 0.0, 0.0, 0.0);
        simde::type::nucleus h1("H", 1ul, 1836.15, 0.0, 0.0, 1.39839);
        simde::type::nucleus h2("H", 1ul, 1836.15, 0.0, 0.0, 4.39839);
        simde::type::nucleus h3("H", 1ul, 1836.15, 0.0, 0.0, 5.79678);
        simde::type::nuclei h2_dimer_nuclei{h0, h1, h2, h3};
        auto ao_bs = test_scf::h_basis(h2_dimer_nuclei);
        simde::type::molecule h2_dimer_mol(0, 1, h2_dimer_nuclei);
        simde::type::chemical_system h2_dimer_sys(h2_dimer_mol);
        const auto e =
          mm.template run_as<pt>("SCF Driver", ao_bs, h2_dimer_sys);
        pcorr.set_elem({}, float_type{-2.2260535919670001});
        simde::type::tensor corr(shape_corr, std::move(pcorr));
        REQUIRE(approximately_equal(corr, e, 1E-6));
    }

    SECTION("Water") {
        using atom_t      = simde::type::atom;
        using molecule_t  = simde::type::molecule;
        const auto a2b    = 1.8897259886; // Angstroms to bohrs
        const auto H_mass = 1822.877;     // Hydrogen mass in atomic units
        const auto O_mass = 29166.037;    // Oxygen mass in atomic units
        atom_t H0("H", 1ul, H_mass, -1.958940 * a2b, -0.032063 * a2b,
                  0.725554 * a2b);
        atom_t H1("H", 1ul, H_mass, -0.607485 * a2b, 0.010955 * a2b,
                  0.056172 * a2b);
        atom_t O0("O", 8ul, O_mass, -1.538963 * a2b, 0.004548 * a2b,
                  -0.117331 * a2b);
        molecule_t water{H0, H1, O0};
        auto aos =
          mm.template run_as<simde::MolecularBasisSet>("STO-3G", water);

        simde::type::chemical_system water_cs(water);
        auto e = mm.template run_as<pt>("SCF Driver", aos, water_cs);
        pcorr.set_elem({}, float_type{-74.9602586404361944});
        simde::type::tensor corr(shape_corr, std::move(pcorr));
        REQUIRE(approximately_equal(corr, e, 1E-6));
    }
}
