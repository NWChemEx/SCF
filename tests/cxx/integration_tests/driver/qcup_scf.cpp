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

#include "../integration_tests.hpp"
#include <algorithm>
#include <cmath>
#include <string>
#include <tuple>
#include <vector>

#ifdef ENABLE_SIGMA

using pt       = simde::AOEnergy;
using tensor_t = simde::type::tensor;
using Catch::Matchers::WithinAbs;

namespace {

/// The center and radius of the scalar @p e, whose type is FloatType
template<typename FloatType>
std::pair<double, double> center_and_radius(const tensor_t& e) {
    using tensorwrapper::buffer::get_raw_data;
    using tensorwrapper::buffer::make_contiguous;
    using tensorwrapper::types::uq_center;
    using tensorwrapper::types::uq_upper;
    const auto x = get_raw_data<FloatType>(make_contiguous(e.buffer()))[0];
    return {double(uq_center(x)), double(uq_upper(x) - uq_center(x))};
}

/// The centers and radii of the elements of the tensor @p t of FloatType
template<typename FloatType>
auto centers_and_radii(const tensor_t& t) {
    using tensorwrapper::buffer::get_raw_data;
    using tensorwrapper::buffer::make_contiguous;
    using tensorwrapper::types::uq_center;
    using tensorwrapper::types::uq_upper;
    std::vector<double> centers, radii;
    for(const auto& x : get_raw_data<FloatType>(make_contiguous(t.buffer()))) {
        centers.push_back(double(uq_center(x)));
        radii.push_back(double(uq_upper(x) - uq_center(x)));
    }
    return std::pair{centers, radii};
}

/// Water, so that both screening and incomplete convergence matter
simde::type::molecule make_water() {
    using atom_t      = simde::type::atom;
    const auto a2b    = 1.8897259886; // Angstroms to bohrs
    const auto H_mass = 1822.877;     // Hydrogen mass in atomic units
    const auto O_mass = 29166.037;    // Oxygen mass in atomic units
    atom_t H0("H", 1ul, H_mass, -1.958940 * a2b, -0.032063 * a2b,
              0.725554 * a2b);
    atom_t H1("H", 1ul, H_mass, -0.607485 * a2b, 0.010955 * a2b,
              0.056172 * a2b);
    atom_t O0("O", 8ul, O_mass, -1.538963 * a2b, 0.004548 * a2b,
              -0.117331 * a2b);
    return simde::type::molecule{H0, H1, O0};
}

/// True if any of the inactive settings of @p mod mentions @p what
bool reports(const pluginplay::Module& mod, const std::string& what) {
    const auto inactive = mod.inactive_settings();
    return std::any_of(inactive.begin(), inactive.end(), [&](const auto& msg) {
        return msg.find(what) != std::string::npos;
    });
}

} // namespace

TEST_CASE("QCUP SCF") {
    using idouble  = tensorwrapper::types::idouble;
    using tadouble = tensorwrapper::types::tadouble;

    auto mm = test_scf::load_modules<double>();

    const auto water = make_water();
    auto aos         = mm.run_as<simde::MolecularBasisSet>("STO-3G", water);
    simde::type::chemical_system h2o(water);

    // Correct energy, without screening, is validated against Psi4
    const double corr = -74.9602586404361944;

    // The energy of an ordinary SCF which screens the ERIs as QCUP SCF does
    // by default
    const double screened = -74.9602557626054;

    // The uncertainty must account for screening: [c - r, c + r] contains corr
    auto check = [&](double c, double radius) {
        REQUIRE_THAT(c, WithinAbs(screened, 1.0E-9));
        REQUIRE(std::fabs(c - corr) <= radius);
    };

    SECTION("Defaults: intervals and the analytic error model") {
        const auto e           = mm.run_as<pt>("QCUP SCF", aos, h2o);
        const auto [c, radius] = center_and_radius<idouble>(e);
        check(c, radius);

        // The primitive error model is not in the call graph, so its
        // setting has no effect
        REQUIRE(reports(mm.at("QCUP SCF"), "Primitive Error Model"));
    }

    SECTION("Ordinary SCF is unaffected") {
        mm.run_as<pt>("QCUP SCF", aos, h2o);

        // The ModuleManager's modules are not changed...
        const auto& uq_type = mm.at("UQ Initializer").inputs().at("UQ Type");
        REQUIRE(uq_type.value<std::string>() == "none");
        const auto& eri =
          mm.at("Four center J builder").submods().at("Four-center ERI");
        REQUIRE(&eri.value() == &mm.at("ERI4"));

        // ... so the SCF Driver still runs with doubles, without screening
        const auto e = mm.run_as<pt>("SCF Driver", aos, h2o);
        using tensorwrapper::buffer::get_raw_data;
        using tensorwrapper::buffer::make_contiguous;
        const auto x = get_raw_data<double>(make_contiguous(e.buffer()))[0];
        REQUIRE_THAT(x, WithinAbs(corr, 1.0E-6));
    }

    SECTION("The screened reference energy") {
        mm.change_input("ERI4", "Threshold", 1.0E-6);
        const auto e = mm.run_as<pt>("SCF Driver", aos, h2o);
        using tensorwrapper::buffer::get_raw_data;
        using tensorwrapper::buffer::make_contiguous;
        const auto x = get_raw_data<double>(make_contiguous(e.buffer()))[0];
        REQUIRE_THAT(x, WithinAbs(screened, 1.0E-9));
    }

    SECTION("Primitive error model") {
        mm.change_submod("QCUP SCF", "ERI Error Model",
                         "Primitive Error Model");
        const auto e           = mm.run_as<pt>("QCUP SCF", aos, h2o);
        const auto [c, radius] = center_and_radius<idouble>(e);
        check(c, radius);

        // The benchmark ERIs are only used by the analytic error model
        REQUIRE(reports(mm.at("QCUP SCF"), "Analytic Error"));
    }

    SECTION("A copy with another uncertain type") {
        mm.copy_module("QCUP SCF", "QCUP SCF taffine");
        mm.change_input("QCUP SCF taffine", "UQ Type",
                        std::string("thresholded affine"));

        const auto e_ta        = mm.run_as<pt>("QCUP SCF taffine", aos, h2o);
        const auto [c, radius] = center_and_radius<tadouble>(e_ta);
        check(c, radius);

        // The original still uses intervals
        const auto e = mm.run_as<pt>("QCUP SCF", aos, h2o);
        REQUIRE_NOTHROW(center_and_radius<idouble>(e));
    }

    SECTION("Uncertainty attached after convergence") {
        // Converge loosely, so that incomplete convergence matters
        for(const auto* tol :
            {"energy tolerance", "density tolerance", "gradient tolerance"})
            mm.change_input("Loop", tol, 1.0E-3);

        mm.copy_module("QCUP SCF", "QCUP SCF no post");
        mm.change_submod("QCUP SCF no post", "MO Uncertainty Model",
                         "No MO uncertainty");
        mm.change_submod("QCUP SCF no post", "Convergence Uncertainty Model",
                         "No convergence uncertainty");

        const auto e            = mm.run_as<pt>("QCUP SCF", aos, h2o);
        const auto e_no         = mm.run_as<pt>("QCUP SCF no post", aos, h2o);
        const auto [c, r]       = center_and_radius<idouble>(e);
        const auto [c_no, r_no] = center_and_radius<idouble>(e_no);
        REQUIRE(c == c_no);
        REQUIRE(r > r_no);
    }
}

// "SCF Driver" only returns the energy, which does not depend on the MO
// coefficients "Loop" returns. So "Loop" is run directly, wired by hand as
// "QCUP SCF" wires it, to check that the MO uncertainty reaches the orbitals.
TEST_CASE("Loop with MO uncertainty") {
    using idouble = tensorwrapper::types::idouble;
    using wf_type = simde::type::rscf_wf;
    using ham_pt =
      simde::Convert<simde::type::hamiltonian, simde::type::chemical_system>;
    using guess_pt = simde::InitialGuess<wf_type>;
    using braket_t =
      chemist::braket::BraKet<wf_type, simde::type::hamiltonian, wf_type>;
    using egy_pt  = simde::EvaluateBraKet<braket_t>;
    using loop_pt = simde::Optimize<egy_pt, wf_type>;

    auto mm = test_scf::load_modules<double>();

    // Intervals, with the error of the screened ERIs from the analytic model
    mm.change_input("UQ Initializer", "UQ Type", std::string("interval"));
    mm.change_submod("Loop", "Diagonalizer", "Generalized eigensolve");
    mm.change_submod("Diagonalization Fock update", "Diagonalizer",
                     "Generalized eigensolve");
    for(const auto* key : {"Four center J builder", "Four center K builder"})
        mm.change_submod(key, "Four-center ERI", "UQ Atom Symm Blocked Driver");
    mm.change_submod("UQ Atom Symm Blocked Driver", "ERI Error",
                     "Analytic Error");
    mm.change_submod("Analytic Error", "ERI4s", "Analytic Benchmark ERI4");
    mm.change_input("ERI4", "Threshold", 1.0E-6);

    mm.copy_module("Loop", "Loop with MO uncertainty");
    mm.change_submod("Loop with MO uncertainty", "MO uncertainty",
                     "First-order MO uncertainty");

    const auto water = make_water();
    simde::type::chemical_system h2o(water);
    simde::type::aos aos(mm.run_as<simde::MolecularBasisSet>("STO-3G", water));
    const auto H    = mm.run_as<ham_pt>("Born-Oppenheimer approximation", h2o);
    const auto psi0 = mm.run_as<guess_pt>("Core guess", H, aos);
    braket_t H_00(psi0, H, psi0);

    const auto psi =
      std::get<1>(mm.run_as<loop_pt>("Loop with MO uncertainty", H_00, psi0));
    const auto psi_no = std::get<1>(mm.run_as<loop_pt>("Loop", H_00, psi0));

    const auto [c, r] = centers_and_radii<idouble>(psi.orbitals().transform());
    const auto [c_no, r_no] =
      centers_and_radii<idouble>(psi_no.orbitals().transform());

    // The orbitals diagonalize the center of the Fock matrix, so the couplings,
    // and thus the correction, are centered at zero. The centers are unchanged
    // and every radius grows.
    REQUIRE(c.size() == c_no.size());
    for(std::size_t i = 0; i < c.size(); ++i) {
        REQUIRE_THAT(c[i], WithinAbs(c_no[i], 1.0E-10));
        REQUIRE(r[i] > r_no[i]);
    }
}

#endif
