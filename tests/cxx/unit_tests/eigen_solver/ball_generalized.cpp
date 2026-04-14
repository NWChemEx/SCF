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
#include <scf/scf.hpp>
#include <simde/simde.hpp>

namespace {

// 7x7 Fock / overlap-style fixtures and expected generalized eigenvalues
// (L_ball) transcribed from console output: each printed "m+/-r" token is
// stored as idouble(lower, upper) with lower = m - r, upper = m + r.

inline simde::type::tensor uncertain_A() {
    using uq_type = tensorwrapper::types::idouble;
    tensorwrapper::shape::Smooth shape{7, 7};
    auto buf = tensorwrapper::buffer::make_contiguous<uq_type>(shape);
    buf.set_elem({0, 0},
                 uq_type{-5.9761941965020215e-01, -5.9761890608307133e-01});
    buf.set_elem({0, 1},
                 uq_type{-4.0307029666324207e-01, -4.0306978309613190e-01});
    buf.set_elem({0, 2},
                 uq_type{-1.2374285072777236e+00, -1.2374279937106170e+00});
    buf.set_elem({0, 3},
                 uq_type{-1.0494685522547358e+00, -1.0494680386876154e+00});
    buf.set_elem({0, 4},
                 uq_type{2.1281939270258100e-02, 2.1282452837359501e-02});
    buf.set_elem({0, 5},
                 uq_type{-5.0779585063444665e-01, -5.0779533706733271e-01});
    buf.set_elem({0, 6},
                 uq_type{1.9650251601842519e-01, 1.9650302958553259e-01});
    buf.set_elem({1, 0},
                 uq_type{-4.0307029666324190e-01, -4.0306978309613173e-01});
    buf.set_elem({1, 1},
                 uq_type{-5.9047593004366428e-01, -5.9047541647653523e-01});
    buf.set_elem({1, 2},
                 uq_type{-1.2227376779992984e+00, -1.2227371644321914e+00});
    buf.set_elem({1, 3},
                 uq_type{-1.0408334164996009e+00, -1.0408329029324810e+00});
    buf.set_elem({1, 4},
                 uq_type{-1.7104644978368001e-03, -1.7099509307360000e-03});
    buf.set_elem({1, 5},
                 uq_type{-1.4859304502705128e-01, -1.4859253145994750e-01});
    buf.set_elem({1, 6},
                 uq_type{-5.2200750954699371e-01, -5.2200699597987887e-01});
    buf.set_elem({2, 0},
                 uq_type{-1.2374285072777242e+00, -1.2374279937106172e+00});
    buf.set_elem({2, 1},
                 uq_type{-1.2227376779992984e+00, -1.2227371644321918e+00});
    buf.set_elem({2, 2},
                 uq_type{-2.0229432721235455e+01, -2.0229432207668175e+01});
    buf.set_elem({2, 3},
                 uq_type{-5.1628583778067823e+00, -5.1628578642396556e+00});
    buf.set_elem({2, 4},
                 uq_type{7.9385631994439998e-04, 7.9436988704520005e-04});
    buf.set_elem({2, 5},
                 uq_type{-2.6603689531534101e-02, -2.6603175964432901e-02});
    buf.set_elem({2, 6},
                 uq_type{-1.3130169293639199e-02, -1.3129655726538200e-02});
    buf.set_elem({3, 0},
                 uq_type{-1.0494685522547353e+00, -1.0494680386876150e+00});
    buf.set_elem({3, 1},
                 uq_type{-1.0408334164996000e+00, -1.0408329029324805e+00});
    buf.set_elem({3, 2},
                 uq_type{-5.1628583778067814e+00, -5.1628578642396548e+00});
    buf.set_elem({3, 3},
                 uq_type{-2.4488877696649971e+00, -2.4488872560978492e+00});
    buf.set_elem({3, 4},
                 uq_type{3.5224479263691001e-03, 3.5229614934699002e-03});
    buf.set_elem({3, 5},
                 uq_type{-1.1753004101803359e-01, -1.1752952745093080e-01});
    buf.set_elem({3, 6},
                 uq_type{-5.6952104344840299e-02, -5.6951590777738301e-02});
    buf.set_elem({4, 0},
                 uq_type{2.1281939270258100e-02, 2.1282452837359501e-02});
    buf.set_elem({4, 1},
                 uq_type{-1.7104644978368001e-03, -1.7099509307360000e-03});
    buf.set_elem({4, 2},
                 uq_type{7.9385631994439998e-04, 7.9436988704520005e-04});
    buf.set_elem({4, 3},
                 uq_type{3.5224479263691001e-03, 3.5229614934699002e-03});
    buf.set_elem({4, 4},
                 uq_type{-3.9094388879395919e-01, -3.9094337522681283e-01});
    buf.set_elem({4, 5},
                 uq_type{-1.5239204150682001e-03, -1.5234068479674000e-03});
    buf.set_elem({4, 6},
                 uq_type{1.1589217834762000e-03, 1.1594353505770000e-03});
    buf.set_elem({5, 0},
                 uq_type{-5.0779585063444688e-01, -5.0779533706733360e-01});
    buf.set_elem({5, 1},
                 uq_type{-1.4859304502705112e-01, -1.4859253145994750e-01});
    buf.set_elem({5, 2},
                 uq_type{-2.6603689531534198e-02, -2.6603175964433200e-02});
    buf.set_elem({5, 3},
                 uq_type{-1.1753004101803349e-01, -1.1752952745093090e-01});
    buf.set_elem({5, 4},
                 uq_type{-1.5239204150682001e-03, -1.5234068479674000e-03});
    buf.set_elem({5, 5},
                 uq_type{-3.5326421481592180e-01, -3.5326370124877843e-01});
    buf.set_elem({5, 6},
                 uq_type{-1.0335912128416299e-02, -1.0335398561315099e-02});
    buf.set_elem({6, 0},
                 uq_type{1.9650251601842419e-01, 1.9650302958553159e-01});
    buf.set_elem({6, 1},
                 uq_type{-5.2200750954699304e-01, -5.2200699597987910e-01});
    buf.set_elem({6, 2},
                 uq_type{-1.3130169293640000e-02, -1.3129655726539001e-02});
    buf.set_elem({6, 3},
                 uq_type{-5.6952104344839098e-02, -5.6951590777737100e-02});
    buf.set_elem({6, 4},
                 uq_type{1.1589217834762000e-03, 1.1594353505770000e-03});
    buf.set_elem({6, 5},
                 uq_type{-1.0335912128416299e-02, -1.0335398561315099e-02});
    buf.set_elem({6, 6},
                 uq_type{-3.3401382547702690e-01, -3.3401331190988387e-01});
    return simde::type::tensor(shape, std::move(buf));
}

inline simde::type::tensor uncertain_B() {
    using uq_type = tensorwrapper::types::idouble;
    tensorwrapper::shape::Smooth shape{7, 7};
    auto buf = tensorwrapper::buffer::make_contiguous<uq_type>(shape);
    buf.set_elem({0, 0},
                 uq_type{9.9999999999999978e-01, 1.0000000000000002e+00});
    buf.set_elem({0, 1},
                 uq_type{2.5358257133619461e-01, 2.5358257133619483e-01});
    buf.set_elem({0, 2},
                 uq_type{5.5900804544533403e-02, 5.5900804544533597e-02});
    buf.set_elem({0, 3},
                 uq_type{4.8454762713240557e-01, 4.8454762713240579e-01});
    buf.set_elem({0, 4},
                 uq_type{-1.5482109990821101e-02, -1.5482109990820900e-02});
    buf.set_elem({0, 5},
                 uq_type{3.5644036709221749e-01, 3.5644036709221771e-01});
    buf.set_elem({0, 6},
                 uq_type{-1.7760045089221962e-01, -1.7760045089221940e-01});
    buf.set_elem({1, 0},
                 uq_type{2.5358257133619461e-01, 2.5358257133619483e-01});
    buf.set_elem({1, 1},
                 uq_type{9.9999999999999978e-01, 1.0000000000000002e+00});
    buf.set_elem({1, 2},
                 uq_type{5.5230922941791903e-02, 5.5230922941792097e-02});
    buf.set_elem({1, 3},
                 uq_type{4.8117972220632926e-01, 4.8117972220632949e-01});
    buf.set_elem({1, 4},
                 uq_type{2.6830766568749997e-03, 2.6830766568752000e-03});
    buf.set_elem({1, 5},
                 uq_type{7.2658318900858498e-02, 7.2658318900858693e-02});
    buf.set_elem({1, 6},
                 uq_type{3.9007755239467901e-01, 3.9007755239467923e-01});
    buf.set_elem({2, 0},
                 uq_type{5.5900804544533403e-02, 5.5900804544533597e-02});
    buf.set_elem({2, 1},
                 uq_type{5.5230922941791903e-02, 5.5230922941792097e-02});
    buf.set_elem({2, 2},
                 uq_type{1.0000000000000002e+00, 1.0000000000000007e+00});
    buf.set_elem({2, 3},
                 uq_type{2.3670392057272610e-01, 2.3670392057272632e-01});
    buf.set_elem({2, 4},
                 uq_type{-9.9999999999999998e-17, 9.9999999999999998e-17});
    buf.set_elem({2, 5},
                 uq_type{-9.9999999999999998e-17, 9.9999999999999998e-17});
    buf.set_elem({2, 6},
                 uq_type{-2.0000000000000000e-16, 0.0000000000000000e+00});
    buf.set_elem({3, 0},
                 uq_type{4.8454762713240557e-01, 4.8454762713240579e-01});
    buf.set_elem({3, 1},
                 uq_type{4.8117972220632926e-01, 4.8117972220632949e-01});
    buf.set_elem({3, 2},
                 uq_type{2.3670392057272618e-01, 2.3670392057272641e-01});
    buf.set_elem({3, 3},
                 uq_type{9.9999999999999967e-01, 9.9999999999999989e-01});
    buf.set_elem({3, 4},
                 uq_type{-9.9999999999999998e-17, 9.9999999999999998e-17});
    buf.set_elem({3, 5},
                 uq_type{-9.9999999999999998e-17, 9.9999999999999998e-17});
    buf.set_elem({3, 6},
                 uq_type{-7.9999999999999998e-16, -6.0000000000000009e-16});
    buf.set_elem({4, 0},
                 uq_type{-1.5482109990821101e-02, -1.5482109990820900e-02});
    buf.set_elem({4, 1},
                 uq_type{2.6830766568749997e-03, 2.6830766568752000e-03});
    buf.set_elem({4, 2},
                 uq_type{-9.9999999999999998e-17, 9.9999999999999998e-17});
    buf.set_elem({4, 3},
                 uq_type{-9.9999999999999998e-17, 9.9999999999999998e-17});
    buf.set_elem({4, 4},
                 uq_type{1.0000000000000000e+00, 1.0000000000000004e+00});
    buf.set_elem({4, 5},
                 uq_type{-9.9999999999999998e-17, 9.9999999999999998e-17});
    buf.set_elem({4, 6},
                 uq_type{-9.9999999999999998e-17, 9.9999999999999998e-17});
    buf.set_elem({5, 0},
                 uq_type{3.5644036709221738e-01, 3.5644036709221760e-01});
    buf.set_elem({5, 1},
                 uq_type{7.2658318900858498e-02, 7.2658318900858693e-02});
    buf.set_elem({5, 2},
                 uq_type{-9.9999999999999998e-17, 9.9999999999999998e-17});
    buf.set_elem({5, 3},
                 uq_type{-9.9999999999999998e-17, 9.9999999999999998e-17});
    buf.set_elem({5, 4},
                 uq_type{-9.9999999999999998e-17, 9.9999999999999998e-17});
    buf.set_elem({5, 5},
                 uq_type{9.9999999999999978e-01, 1.0000000000000002e+00});
    buf.set_elem({5, 6},
                 uq_type{-9.9999999999999998e-17, 9.9999999999999998e-17});
    buf.set_elem({6, 0},
                 uq_type{-1.7760045089221932e-01, -1.7760045089221910e-01});
    buf.set_elem({6, 1},
                 uq_type{3.9007755239467889e-01, 3.9007755239467912e-01});
    buf.set_elem({6, 2},
                 uq_type{-2.9999999999999999e-16, -9.9999999999999998e-17});
    buf.set_elem({6, 3},
                 uq_type{-9.0000000000000003e-16, -6.9999999999999994e-16});
    buf.set_elem({6, 4},
                 uq_type{-9.9999999999999998e-17, 9.9999999999999998e-17});
    buf.set_elem({6, 5},
                 uq_type{-9.9999999999999998e-17, 9.9999999999999998e-17});
    buf.set_elem({6, 6},
                 uq_type{1.0000000000000000e+00, 1.0000000000000004e+00});
    return simde::type::tensor(shape, std::move(buf));
}

// Reference L_ball from pasted console output (m+/-0 -> idouble(m,m)).
inline simde::type::tensor l_ball_corr() {
    using uq_type = tensorwrapper::types::idouble;
    tensorwrapper::shape::Smooth shape{7};
    auto buf = tensorwrapper::buffer::make_contiguous<uq_type>(shape);
    buf.set_elem({0}, uq_type{-20.2383215358238289, -20.2383215358238289});
    buf.set_elem({1}, uq_type{-1.2730342062134834, -1.2730342062134834});
    buf.set_elem({2}, uq_type{-0.6267577709783282, -0.6267577709783282});
    buf.set_elem({3}, uq_type{-0.4510422475343768, -0.4510422475343768});
    buf.set_elem({4}, uq_type{-0.3910154399987634, -0.3910154399987634});
    buf.set_elem({5}, uq_type{0.6153219563024107, 0.6153219563024107});
    buf.set_elem({6}, uq_type{0.7613504962941617, 0.7613504962941617});
    return simde::type::tensor(shape, std::move(buf));
}

} // namespace

using types = std::tuple<tensorwrapper::types::idouble>;
TEMPLATE_LIST_TEST_CASE("BallGeneralized", "", types) {
    using uq_type    = TestType;
    using float_type = typename uq_type::value_t;
    pluginplay::ModuleManager mm;
    scf::load_modules(mm);

    using pt = simde::GeneralizedEigenSolve;

    auto& mod = mm.at("Generalized eigensolve via Ball Arithmetic");

    SECTION("2 by 2 test case") {
        tensorwrapper::shape::Smooth shape{2, 2};
        using tensorwrapper::buffer::make_contiguous;
        auto A_buffer    = make_contiguous<uq_type>(shape);
        float_type noise = 0.001;

        A_buffer.set_elem({0, 0}, uq_type{1.0 - noise, 1.0 + noise});
        A_buffer.set_elem({0, 1}, uq_type{2.0 - noise, 2.0 + noise});
        A_buffer.set_elem({1, 0}, uq_type{2.0 - noise, 2.0 + noise});
        A_buffer.set_elem({1, 1}, uq_type{3.0 - noise, 3.0 + noise});

        auto B_buffer = make_contiguous<uq_type>(shape);
        B_buffer.set_elem({0, 0}, uq_type{1.0 - noise, 1 + noise});
        B_buffer.set_elem({0, 1}, uq_type{0.0 - noise, 0.0 + noise});
        B_buffer.set_elem({1, 0}, uq_type{0.0 - noise, 0.0 + noise});
        B_buffer.set_elem({1, 1}, uq_type{1.0 - noise, 1.0 + noise});

        simde::type::tensor A(shape, std::move(A_buffer));
        simde::type::tensor B(shape, std::move(B_buffer));

        auto&& [values, vector] = mod.run_as<pt>(A, B);

        std::vector<uq_type> expected_values{-0.236068, 4.236068};
        tensorwrapper::shape::Smooth corr_shape{2};
        tensorwrapper::buffer::Contiguous corr_buffer(expected_values,
                                                      corr_shape);
        auto value_buffer = make_contiguous(values.buffer());
        for(std::size_t i = 0; i < corr_buffer.shape().size(); ++i) {
            auto corr_value = corr_buffer.get_elem({i});
            auto value      = value_buffer.get_elem({i});
            auto corr_uq    = wtf::fp::float_cast<uq_type>(corr_value);
            auto value_uq   = wtf::fp::float_cast<uq_type>(value);
            REQUIRE(value_uq.contains(corr_uq.median()));
        }
    }

    SECTION("7 by 7 test case") {
        auto A                  = uncertain_A();
        auto B                  = uncertain_B();
        auto&& [values, vector] = mod.run_as<pt>(A, B);
        auto L_ball             = l_ball_corr();

        auto value_buffer = make_contiguous(values.buffer());
        auto corr_buffer  = make_contiguous(L_ball.buffer());
        for(std::size_t i = 0; i < corr_buffer.shape().size(); ++i) {
            auto corr_value = corr_buffer.get_elem({i});
            auto value      = value_buffer.get_elem({i});
            auto corr_uq    = wtf::fp::float_cast<uq_type>(corr_value);
            auto value_uq   = wtf::fp::float_cast<uq_type>(value);
            REQUIRE(value_uq.contains(corr_uq.median()));
        }
    }
}
