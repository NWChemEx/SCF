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
#include <scf/scf.hpp>
#include <simde/simde.hpp>

namespace test_scf {

using float_types = std::tuple<double, tensorwrapper::types::udouble>;

/// Makes a H nucleus at the point @p x, @p y, @p z
inline auto h_nucleus(double x, double y, double z) {
    return simde::type::nucleus("H", 1ul, 1836.15, x, y, z);
}

inline auto he_nucleus(double x, double y, double z) {
    return simde::type::nucleus("He", 2ul, 7344.61, x, y, z);
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

template<typename ResultType>
ResultType make_he() {
    using simde::type::chemical_system;
    using simde::type::molecule;
    using simde::type::nuclei;
    if constexpr(std::is_same_v<ResultType, nuclei>) {
        auto he0 = he_nucleus(0.0, 0.0, 0.0);
        return nuclei{he0};
    } else if constexpr(std::is_same_v<ResultType, molecule>) {
        return molecule(0, 1, make_he<nuclei>());
    } else if constexpr(std::is_same_v<ResultType, chemical_system>) {
        return chemical_system(make_he<molecule>());
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

template<typename NucleiType>
inline auto he_basis(NucleiType& heliums) {
    using ao_basis_type            = chemist::basis_set::AOBasisSetD;
    using atomic_basis_type        = typename ao_basis_type::value_type;
    using shell_type               = typename atomic_basis_type::value_type;
    using contracted_gaussian_type = typename shell_type::cg_type;
    using center_type = typename atomic_basis_type::shell_traits::center_type;

    std::vector<double> he_coefs{1.5432896730e-01, 5.3532814230e-01,
                                 4.4463454220e-01};
    std::vector<double> he_exps{6.3624213940e+00, 1.1589229990e+00,
                                3.1364979150e-01};
    auto cartesian = shell_type::pure_type::cartesian;
    shell_type::angular_momentum_type l0{0};

    ao_basis_type rv;
    for(std::size_t i = 0; i < heliums.size(); ++i) {
        const auto& hi = heliums[i];
        center_type coords(hi.x(), hi.y(), hi.z());
        contracted_gaussian_type he_cg(he_coefs.begin(), he_coefs.end(),
                                       he_exps.begin(), he_exps.end(), coords);
        atomic_basis_type he_basis("STO-3G", 2, coords);
        he_basis.add_shell(cartesian, l0, he_cg);
        rv.add_center(he_basis);
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

inline auto he_hamiltonian() {
    simde::type::many_electrons es(2);
    auto he = make_he<simde::type::nuclei>();
    simde::type::T_e_type T_e(es);
    simde::type::V_en_type V_en(es, he);
    simde::type::V_ee_type V_ee(es, es);
    return simde::type::hamiltonian(T_e + V_en + V_ee);
}

inline auto h2_aos() {
    auto h2 = make_h2<simde::type::nuclei>();
    return simde::type::aos(h_basis(h2));
}

inline auto he_aos() {
    auto he = make_he<simde::type::nuclei>();
    return simde::type::aos(he_basis(he));
}

template<typename FloatType>
inline auto h2_mos() {
    using mos_type       = simde::type::mos;
    using tensor_type    = typename mos_type::transform_type;
    using allocator_type = tensorwrapper::allocator::Eigen<FloatType>;
    allocator_type alloc(parallelzone::runtime::RuntimeView{});
    tensorwrapper::shape::Smooth shape{2, 2};
    tensorwrapper::layout::Physical l(shape);
    auto c_buffer = alloc.allocate(l);
    c_buffer->set_elem({0, 0}, -0.565516);
    c_buffer->set_elem({0, 1}, -1.07019);
    c_buffer->set_elem({1, 0}, -0.565516);
    c_buffer->set_elem({1, 1}, 1.07019);
    tensor_type t(shape, std::move(c_buffer));
    return mos_type(h2_aos(), std::move(t));
}

template<typename FloatType>
inline auto he_mos() {
    using mos_type       = simde::type::mos;
    using tensor_type    = typename mos_type::transform_type;
    using allocator_type = tensorwrapper::allocator::Eigen<FloatType>;
    allocator_type alloc(parallelzone::runtime::RuntimeView{});
    tensorwrapper::shape::Smooth shape{1, 1};
    tensorwrapper::layout::Physical l(shape);
    auto c_buffer = alloc.allocate(l);
    c_buffer->set_elem({0, 0}, 1.0000);
    tensor_type t(shape, std::move(c_buffer));
    return mos_type(he_aos(), std::move(t));
}

template<typename FloatType>
inline auto h2_cmos() {
    using cmos_type      = simde::type::cmos;
    using tensor_type    = typename cmos_type::transform_type;
    using allocator_type = tensorwrapper::allocator::Eigen<FloatType>;
    allocator_type alloc(parallelzone::runtime::RuntimeView{});
    tensorwrapper::shape::Smooth shape{2};
    tensorwrapper::layout::Physical l(shape);
    auto e_buffer = alloc.allocate(l);
    e_buffer->set_elem({0}, -1.25330893);
    e_buffer->set_elem({1}, -0.47506974);
    tensor_type e(shape, std::move(e_buffer));
    return cmos_type(std::move(e), h2_aos(), h2_mos<FloatType>().transform());
}

template<typename FloatType>
inline auto he_cmos() {
    using cmos_type      = simde::type::cmos;
    using tensor_type    = typename cmos_type::transform_type;
    using allocator_type = tensorwrapper::allocator::Eigen<FloatType>;
    allocator_type alloc(parallelzone::runtime::RuntimeView{});
    tensorwrapper::shape::Smooth shape{1};
    tensorwrapper::layout::Physical l(shape);
    auto e_buffer = alloc.allocate(l);
    e_buffer->set_elem({0}, -0.876036);
    tensor_type e(shape, std::move(e_buffer));
    return cmos_type(std::move(e), he_aos(), he_mos<FloatType>().transform());
}

template<typename FloatType>
inline auto h2_density() {
    using density_type   = simde::type::decomposable_e_density;
    using tensor_type    = typename density_type::value_type;
    using allocator_type = tensorwrapper::allocator::Eigen<FloatType>;
    allocator_type alloc(parallelzone::runtime::RuntimeView{});
    tensorwrapper::shape::Smooth shape{2, 2};
    tensorwrapper::layout::Physical l(shape);
    auto pbuffer = alloc.construct(l, 0.31980835);
    tensor_type t(shape, std::move(pbuffer));
    return density_type(std::move(t), h2_mos<FloatType>());
}

template<typename FloatType>
inline auto he_density() {
    using density_type   = simde::type::decomposable_e_density;
    using tensor_type    = typename density_type::value_type;
    using allocator_type = tensorwrapper::allocator::Eigen<FloatType>;
    allocator_type alloc(parallelzone::runtime::RuntimeView{});
    tensorwrapper::shape::Smooth shape{1, 1};
    tensorwrapper::layout::Physical l(shape);
    auto pbuffer = alloc.construct(l, 1.000);
    tensor_type t(shape, std::move(pbuffer));
    return density_type(std::move(t), he_mos<FloatType>());
}

/// The Fock matrix consistent with h2_hamiltonian and h2_density
template<typename ElectronType, typename FloatType = double>
inline auto h2_fock() {
    ElectronType es;
    if constexpr(std::is_same_v<ElectronType, simde::type::many_electrons>) {
        es = simde::type::many_electrons(2);
    }

    auto h2  = make_h2<simde::type::nuclei>();
    auto rho = h2_density<FloatType>();
    simde::type::fock F;
    using namespace chemist::qm_operator;
    using t_type = Kinetic<ElectronType>;
    using v_type = Coulomb<ElectronType, chemist::Nuclei>;
    using j_type = Coulomb<ElectronType, simde::type::decomposable_e_density>;
    using k_type = Exchange<ElectronType, simde::type::decomposable_e_density>;

    F.emplace_back(1.0, std::make_unique<t_type>(es));
    F.emplace_back(1.0, std::make_unique<v_type>(es, h2));
    F.emplace_back(2.0, std::make_unique<j_type>(es, rho));
    F.emplace_back(-1.0, std::make_unique<k_type>(es, rho));
    return F;
}

} // namespace test_scf
