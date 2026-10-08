/*
 * Copyright 2025 NWChemEx-Project
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
#include <integrals/integrals.hpp>
#include <string>
#include <tuple>
#include <type_traits>
using grid_pt = simde::MolecularGrid;
using rho_pt  = simde::EDensityCollocationMatrix;
using wf_type = simde::type::rscf_wf;
using pt      = simde::eval_braket<wf_type, simde::type::XC_e_type, wf_type>;

TEST_CASE("LibXCEnergy") {
    using float_type = double;
    using tensorwrapper::operations::approximately_equal;

    pluginplay::ModuleManager mm;
    scf::load_modules(mm);
    mm.change_submod("LibXC Energy", "Integration grid",
                     "Grid From IntegratorXX");
    auto& mod = mm.at("LibXC Energy");

    // Matches the (radial x angular) point count of the reference
    // he_grid.txt/h2_grid.txt files (150*194 = 29100, 2*150*194 = 58200)
    // generated with IntegratorXX's own defaults (MuraKnowles radial,
    // Lebedev-Laikov angular, unpruned, Becke partition), so this
    // cross-validates the generated grid against the known-good file-based
    // one -- not a bit-identical grid, hence the looser tolerance below
    // (vs. the GridFromFile-based unit tests' 1e-5).
    mm.change_input("Grid From IntegratorXX", "Radial Size", std::size_t(150));
    mm.change_input("Grid From IntegratorXX", "Angular Size", std::size_t(194));

    SECTION("He") {
        auto rho    = test_scf::he_density<float_type>();
        auto aos    = test_scf::he_aos();
        auto psi    = test_scf::he_wave_function<float_type>();
        auto func   = chemist::qm_operator::xc_functional::SVWN3;
        auto xc_hat = test_scf::he_xc<float_type>(func);
        chemist::braket::BraKet braket(psi, xc_hat, psi);

#ifdef BUILD_LIBXC
        auto exc = mod.run_as<pt>(braket);
        simde::type::tensor corr(-1.01982);
        REQUIRE(approximately_equal(exc, corr, 1e-3));
#else
        REQUIRE_THROWS_AS(mod.run_as<pt>(braket), std::runtime_error);
#endif
    }

    SECTION("H2") {
        auto rho    = test_scf::h2_density<float_type>();
        auto aos    = test_scf::h2_aos();
        auto psi    = test_scf::h2_wave_function<float_type>();
        auto func   = chemist::qm_operator::xc_functional::SVWN3;
        auto xc_hat = test_scf::h2_xc<float_type>(func);
        chemist::braket::BraKet braket(psi, xc_hat, psi);

#ifdef BUILD_LIBXC
        auto exc = mod.run_as<pt>(braket);
        simde::type::tensor corr(-0.734458);
        REQUIRE(approximately_equal(exc, corr, 1e-3));
#else
        REQUIRE_THROWS_AS(mod.run_as<pt>(braket), std::runtime_error);
#endif
    }
}

namespace {

#ifdef ENABLE_SIGMA
using svwn3_types = std::tuple<double, tensorwrapper::types::idouble>;
#else
using svwn3_types = std::tuple<double>;
#endif

/// The name the "Float Type" input of the grid module uses for @p T
template<typename T>
std::string float_type_name() {
    if constexpr(std::is_same_v<T, double>) {
        return "double";
    } else {
        return "idouble";
    }
}

/// Converts the scalar @p t, which holds a double, to a scalar holding a T
template<typename T>
simde::type::tensor to_float_type(const simde::type::tensor& t) {
    using tensorwrapper::buffer::make_contiguous;
    const auto& buffer = make_contiguous(t.buffer());
    auto value         = wtf::fp::float_cast<double>(buffer.get_elem({}));
    tensorwrapper::shape::Smooth shape{};
    auto pbuffer = make_contiguous<T>(shape);
    pbuffer.set_elem({}, T{value});
    return simde::type::tensor(shape, std::move(pbuffer));
}

/// Makes a ModuleManager whose IntegratorXX grid matches the LibXC tests
pluginplay::ModuleManager make_svwn3_mm(const std::string& float_type) {
    pluginplay::ModuleManager mm;
    scf::load_modules(mm);
    const auto grid_key = "Grid From IntegratorXX";
    mm.change_input(grid_key, "Radial Size", std::size_t(150));
    mm.change_input(grid_key, "Angular Size", std::size_t(194));
    mm.change_input(grid_key, "Float Type", float_type);
    return mm;
}

} // namespace

// Tests the native (LibXC-free) SVWN3 against the same reference values as
// the LibXCEnergy test above and, when LibXC is available, against LibXC's
// result itself. The latter comparison uses a much tighter tolerance since
// both modules use the same grid and the same functional.
TEMPLATE_LIST_TEST_CASE("SVWN3Energy", "", svwn3_types) {
    using float_type = TestType;
    using tensorwrapper::operations::approximately_equal;

    auto mm   = make_svwn3_mm(float_type_name<float_type>());
    auto& mod = mm.at("SVWN3 Energy");
    auto func = chemist::qm_operator::xc_functional::SVWN3;

    SECTION("He") {
        auto psi    = test_scf::he_wave_function<float_type>();
        auto xc_hat = test_scf::he_xc<float_type>(func);
        chemist::braket::BraKet braket(psi, xc_hat, psi);
        auto exc = mod.template run_as<pt>(braket);

        simde::type::tensor corr(-1.01982);
        REQUIRE(
          approximately_equal(exc, to_float_type<float_type>(corr), 1e-3));

#ifdef BUILD_LIBXC
        auto libxc_mm = make_svwn3_mm("double");
        auto psi_d    = test_scf::he_wave_function<double>();
        auto xc_hat_d = test_scf::he_xc<double>(func);
        chemist::braket::BraKet braket_d(psi_d, xc_hat_d, psi_d);
        auto libxc_exc =
          libxc_mm.at("LibXC Energy").template run_as<pt>(braket_d);
        auto libxc_corr = to_float_type<float_type>(libxc_exc);
        REQUIRE(approximately_equal(exc, libxc_corr, 1e-8));
#endif
    }

    SECTION("H2") {
        auto psi    = test_scf::h2_wave_function<float_type>();
        auto xc_hat = test_scf::h2_xc<float_type>(func);
        chemist::braket::BraKet braket(psi, xc_hat, psi);
        auto exc = mod.template run_as<pt>(braket);

        simde::type::tensor corr(-0.734458);
        REQUIRE(
          approximately_equal(exc, to_float_type<float_type>(corr), 1e-3));

#ifdef BUILD_LIBXC
        auto libxc_mm = make_svwn3_mm("double");
        auto psi_d    = test_scf::h2_wave_function<double>();
        auto xc_hat_d = test_scf::h2_xc<double>(func);
        chemist::braket::BraKet braket_d(psi_d, xc_hat_d, psi_d);
        auto libxc_exc =
          libxc_mm.at("LibXC Energy").template run_as<pt>(braket_d);
        auto libxc_corr = to_float_type<float_type>(libxc_exc);
        REQUIRE(approximately_equal(exc, libxc_corr, 1e-8));
#endif
    }
}
