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

/* These tests use QCUP
 */
TEST_CASE("QCUP-SCFDriver") {
    using float_type = tensorwrapper::types::udouble;
    auto mm          = test_scf::load_modules<float_type>();
    auto key         = "UQ Atom Symm Blocked Driver";
    mm.change_submod("Four center J builder", "Four-center ERI", key);
    mm.change_submod("Four center K builder", "Four-center ERI", key);
    mm.change_input(key, "Mean Type", "max");
    mm.change_input("ERI4", "With UQ?", false);
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
        std::cout << e << " " << corr << std::endl;
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
