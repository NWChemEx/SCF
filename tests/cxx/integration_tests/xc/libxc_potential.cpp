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
using ao_pt   = simde::AOCollocationMatrix;
using rho_pt  = simde::EDensityCollocationMatrix;
using pt      = simde::aos_xc_e_aos;

TEST_CASE("LibXCPotential") {
    using float_type = double;
    using tensorwrapper::operations::approximately_equal;

    pluginplay::ModuleManager mm;
    scf::load_modules(mm);
    mm.change_submod("LibXC Potential", "Integration grid",
                     "Grid From IntegratorXX");
    auto& mod = mm.at("LibXC Potential");

    // See libxc_energy.cpp's LibXCEnergy test for why these sizes and the
    // loosened tolerance below.
    mm.change_input("Grid From IntegratorXX", "Radial Size", std::size_t(150));
    mm.change_input("Grid From IntegratorXX", "Angular Size", std::size_t(194));

    SECTION("He") {
        auto rho  = test_scf::he_density<float_type>();
        auto aos  = test_scf::he_aos();
        auto func = chemist::qm_operator::xc_functional::SVWN3;
        simde::type::xc_e_type xc_op(func, simde::type::electron{}, rho);
        chemist::braket::BraKet braket(aos, xc_op, aos);

#ifdef BUILD_LIBXC
        auto vxc = mod.run_as<pt>(braket);
        typename simde::type::tensor::matrix_il_type il{{-0.668319}};
        simde::type::tensor corr(il);
        REQUIRE(approximately_equal(vxc, corr, 1e-3));
#else
        REQUIRE_THROWS_AS(mod.run_as<pt>(braket), std::runtime_error);
#endif
    }

    SECTION("H2") {
        auto rho  = test_scf::h2_density<float_type>();
        auto aos  = test_scf::h2_aos();
        auto func = chemist::qm_operator::xc_functional::SVWN3;
        simde::type::xc_e_type xc_op(func, simde::type::electron{}, rho);
        chemist::braket::BraKet braket(aos, xc_op, aos);

#ifdef BUILD_LIBXC
        auto vxc = mod.run_as<pt>(braket);
        simde::type::tensor corr{{-0.453301, -0.296985},
                                 {-0.296985, -0.453301}};
        REQUIRE(approximately_equal(vxc, corr, 1e-3));
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

/// Converts the matrix @p t, which holds doubles, to a matrix holding Ts
template<typename T>
simde::type::tensor to_float_type(const simde::type::tensor& t) {
    using tensorwrapper::buffer::make_contiguous;
    const auto& buffer = make_contiguous(t.buffer());
    const auto n_rows  = buffer.shape().extent(0);
    const auto n_cols  = buffer.shape().extent(1);
    tensorwrapper::shape::Smooth shape{n_rows, n_cols};
    auto pbuffer = make_contiguous<T>(shape);
    for(std::size_t i = 0; i < n_rows; ++i) {
        for(std::size_t j = 0; j < n_cols; ++j) {
            auto value = wtf::fp::float_cast<double>(buffer.get_elem({i, j}));
            pbuffer.set_elem({i, j}, T{value});
        }
    }
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
// the LibXCPotential test above and, when LibXC is available, against LibXC's
// result itself. The latter comparison uses a much tighter tolerance since
// both modules use the same grid and the same functional.
TEMPLATE_LIST_TEST_CASE("SVWN3Potential", "", svwn3_types) {
    using float_type = TestType;
    using tensorwrapper::operations::approximately_equal;

    auto mm   = make_svwn3_mm(float_type_name<float_type>());
    auto& mod = mm.at("SVWN3 Potential");
    auto func = chemist::qm_operator::xc_functional::SVWN3;

    SECTION("He") {
        auto rho = test_scf::he_density<float_type>();
        auto aos = test_scf::he_aos();
        simde::type::xc_e_type xc_op(func, simde::type::electron{}, rho);
        chemist::braket::BraKet braket(aos, xc_op, aos);
        auto vxc = mod.template run_as<pt>(braket);

        typename simde::type::tensor::matrix_il_type il{{-0.668319}};
        simde::type::tensor corr(il);
        REQUIRE(
          approximately_equal(vxc, to_float_type<float_type>(corr), 1e-3));

#ifdef BUILD_LIBXC
        auto libxc_mm = make_svwn3_mm("double");
        auto rho_d    = test_scf::he_density<double>();
        simde::type::xc_e_type xc_op_d(func, simde::type::electron{}, rho_d);
        chemist::braket::BraKet braket_d(aos, xc_op_d, aos);
        auto libxc_vxc =
          libxc_mm.at("LibXC Potential").template run_as<pt>(braket_d);
        auto libxc_corr = to_float_type<float_type>(libxc_vxc);
        REQUIRE(approximately_equal(vxc, libxc_corr, 1e-8));
#endif
    }

    SECTION("H2") {
        auto rho = test_scf::h2_density<float_type>();
        auto aos = test_scf::h2_aos();
        simde::type::xc_e_type xc_op(func, simde::type::electron{}, rho);
        chemist::braket::BraKet braket(aos, xc_op, aos);
        auto vxc = mod.template run_as<pt>(braket);

        simde::type::tensor corr{{-0.453301, -0.296985},
                                 {-0.296985, -0.453301}};
        REQUIRE(
          approximately_equal(vxc, to_float_type<float_type>(corr), 1e-3));

#ifdef BUILD_LIBXC
        auto libxc_mm = make_svwn3_mm("double");
        auto rho_d    = test_scf::h2_density<double>();
        simde::type::xc_e_type xc_op_d(func, simde::type::electron{}, rho_d);
        chemist::braket::BraKet braket_d(aos, xc_op_d, aos);
        auto libxc_vxc =
          libxc_mm.at("LibXC Potential").template run_as<pt>(braket_d);
        auto libxc_corr = to_float_type<float_type>(libxc_vxc);
        REQUIRE(approximately_equal(vxc, libxc_corr, 1e-8));
#endif
    }
}
