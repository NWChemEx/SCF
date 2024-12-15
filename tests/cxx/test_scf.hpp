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
#pragma once
#include <catch2/catch_approx.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <simde/simde.hpp>

namespace test_scf {

/// Makes a H nucleus at the point @p x, @p y, @p z
inline auto h_nucleus(double x, double y, double z) {
    return simde::type::nucleus("H", 1ul, 1836.15, x, y, z);
}

template<typename ResultType>
ResultType make_h2() {
    using simde::type::chemical_system;
    using simde::type::molecule;
    using simde::type::nuclei;
    if constexpr(std::is_same_v<ResultType, nuclei>) {
        auto h0 = h_nucleus(0.0, 0.0, 0.0);
        auto h1 = h_nucleus(0.0, 0.0, 1.3984);
        return nuclei{h0, h1};
    } else if constexpr(std::is_same_v<ResultType, molecule>) {
        return molecule(0, 1, make_h2<nuclei>());
    } else if constexpr(std::is_same_v<ResultType, chemical_system>) {
        return chemical_system(make_h2<molecule>());
    } else {
        // We know this assert fails if we're in the else statement
        // Getting here means you provided a bad type.
        static_assert(std::is_same_v<ResultType, nuclei>);
    }
}

// Applies the STO-3G basis to a Nuclei object filled with hydrogens
template<typename NucleiType>
inline auto h_basis(NucleiType& hydrogens) {
    using ao_basis_type            = chemist::basis_set::AOBasisSetD;
    using atomic_basis_type        = typename ao_basis_type::value_type;
    using shell_type               = typename atomic_basis_type::value_type;
    using contracted_gaussian_type = typename shell_type::cg_type;
    using center_type = typename atomic_basis_type::shell_traits::center_type;

    std::vector<double> h_coefs{0.1543289673, 0.5353281423, 0.4446345422};
    std::vector<double> h_exps{3.425250914, 0.6239137298, 0.1688554040};
    auto cartesian = shell_type::pure_type::cartesian;
    shell_type::angular_momentum_type l0{0};

    ao_basis_type rv;
    for(std::size_t i = 0; i < hydrogens.size(); ++i) {
        const auto& hi = hydrogens[i];
        center_type coords(hi.x(), hi.y(), hi.z());
        contracted_gaussian_type h_cg(h_coefs.begin(), h_coefs.end(),
                                      h_exps.begin(), h_exps.end(), coords);
        atomic_basis_type h_basis("STO-3G", 1, coords);
        h_basis.add_shell(cartesian, l0, h_cg);
        rv.add_center(h_basis);
    }

    return rv;
}

inline auto h2_hamiltonian() {
    simde::type::many_electrons es(2);
    auto h2 = make_h2<simde::type::nuclei>();
    simde::type::T_e_type T_e(es);
    simde::type::V_en_type V_en(es, h2);
    simde::type::V_ee_type V_ee(es, es);
    simde::type::V_nn_type V_nn(h2, h2);
    return simde::type::hamiltonian(T_e + V_en + V_ee + V_nn);
}

inline auto h2_aos() {
    auto h2 = make_h2<simde::type::nuclei>();
    return simde::type::aos(h_basis(h2));
}

inline auto h2_mos() {
    using mos_type    = simde::type::mos;
    using tensor_type = typename mos_type::transform_type;
    tensor_type c({{-0.565516, -1.07019}, {-0.565516, 1.07019}});
    return mos_type(h2_aos(), std::move(c));
}

inline auto h2_density() {
    using density_type = simde::type::decomposable_e_density;
    typename density_type::value_type rho(
      {{0.92501791, -0.54009707}, {-0.28540122, 1.7505162}});
    return density_type(rho, h2_mos());
}

/// The Fock matrix consistent with h2_hamiltonian and h2_density
inline auto h2_fock() {
    simde::type::many_electrons es(2);
    auto h2  = make_h2<simde::type::nuclei>();
    auto rho = h2_density();
    simde::type::fock F;
    F.emplace_back(1.0, std::make_unique<simde::type::T_e_type>(es));
    F.emplace_back(1.0, std::make_unique<simde::type::V_en_type>(es, h2));
    F.emplace_back(2.0, std::make_unique<simde::type::J_e_type>(es, rho));
    F.emplace_back(-1.0, std::make_unique<simde::type::K_e_type>(es, rho));
    return F;
}

} // namespace test_scf