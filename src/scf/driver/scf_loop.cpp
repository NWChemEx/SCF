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

#include "driver.hpp"

namespace scf::driver {
namespace {

struct Kernel {
    explicit Kernel(parallelzone::runtime::RuntimeView rv) : m_rv(rv) {}
    template<typename FloatType>
    auto run(const tensorwrapper::buffer::BufferBase& a, double tol) {
        tensorwrapper::allocator::Eigen<FloatType> allocator(m_rv);
        const auto& eigen_a      = allocator.rebind(a);
        constexpr bool is_float  = std::is_same_v<FloatType, float>;
        constexpr bool is_double = std::is_same_v<FloatType, double>;
        if constexpr(is_float || is_double) {
            return std::fabs(eigen_a.at()) < tol;
        } else {
            return std::fabs(eigen_a.at().mean()) < tol;
        }
    }

    parallelzone::runtime::RuntimeView m_rv;
};

const auto desc = R"(
)";

} // namespace

using simde::type::electronic_hamiltonian;
using simde::type::hamiltonian;
using simde::type::op_base_type;

template<typename WfType>
using egy_pt = simde::eval_braket<WfType, hamiltonian, WfType>;

template<typename WfType>
using elec_egy_pt = simde::eval_braket<WfType, electronic_hamiltonian, WfType>;

template<typename WfType>
using pt = simde::Optimize<egy_pt<WfType>, WfType>;

template<typename WfType>
using update_pt = simde::UpdateGuess<WfType>;

using density_t = simde::type::decomposable_e_density;

using fock_pt = simde::FockOperator<density_t>;

using density_pt = simde::aos_rho_e_aos<simde::type::cmos>;

using v_nn_pt = simde::charge_charge_interaction;

using fock_matrix_pt = simde::aos_f_e_aos;
using s_pt           = simde::aos_s_e_aos;

struct GrabNuclear : chemist::qm_operator::OperatorVisitor {
    using V_nn_type = simde::type::V_nn_type;

    GrabNuclear() : chemist::qm_operator::OperatorVisitor(false) {}

    void run(const V_nn_type& V_nn) { m_pv = &V_nn; }

    const V_nn_type* m_pv;
};

MODULE_CTOR(SCFLoop) {
    using wf_type = simde::type::rscf_wf;
    description(desc);
    satisfies_property_type<pt<wf_type>>();

    const unsigned int max_itr = 20;
    add_input<unsigned int>("max iterations").set_default(max_itr);
    add_input<double>("energy tolerance").set_default(1.0E-6);
    add_input<double>("gradient tolerance").set_default(1.0E-6);

    add_submodule<elec_egy_pt<wf_type>>("Electronic energy");
    add_submodule<density_pt>("Density matrix");
    add_submodule<update_pt<wf_type>>("Guess update");
    add_submodule<fock_pt>("One-electron Fock operator");
    add_submodule<fock_pt>("Fock operator");
    add_submodule<fock_matrix_pt>("Fock matrix builder");
    add_submodule<v_nn_pt>("Charge-charge");
    add_submodule<s_pt>("Overlap matrix builder");
}

MODULE_RUN(SCFLoop) {
    using wf_type               = simde::type::rscf_wf;
    using density_op_type       = simde::type::rho_e<simde::type::cmos>;
    const auto&& [braket, psi0] = pt<wf_type>::unwrap_inputs(inputs);
    // TODO: Assert bra == ket == psi0
    const auto& H      = braket.op();
    const auto& H_elec = H.electronic_hamiltonian();
    const auto& H_core = H_elec.core_hamiltonian();
    const auto& aos    = psi0.orbitals().from_space();

    auto& egy_mod     = submods.at("Electronic energy");
    auto& density_mod = submods.at("Density matrix");
    auto& update_mod  = submods.at("Guess update");
    auto& fock_mod    = submods.at("One-electron Fock operator");
    auto& Fock_mod    = submods.at("Fock operator");
    auto& V_nn_mod    = submods.at("Charge-charge");

    // TODO: should be split off into orbital gradient module
    auto& F_mod = submods.at("Fock matrix builder");
    auto& S_mod = submods.at("Overlap matrix builder");

    // Step 1: Nuclear-nuclear repulsion
    GrabNuclear visitor;
    H.visit(visitor);
    // TODO: Clean up charges class to make this easier...
    const auto& V_nn       = *visitor.m_pv;
    const auto n_lhs       = V_nn.lhs_particle().as_nuclei();
    const auto qs_lhs_view = n_lhs.charges();
    const auto n_rhs       = V_nn.rhs_particle().as_nuclei();
    const auto qs_rhs_view = n_rhs.charges();
    simde::type::charges qs_lhs;
    simde::type::charges qs_rhs;
    for(const auto q_i : qs_lhs_view) {
        qs_lhs.push_back(q_i.as_point_charge());
    }
    for(const auto q_i : qs_rhs_view) {
        qs_rhs.push_back(q_i.as_point_charge());
    }
    auto e_nuclear = V_nn_mod.run_as<v_nn_pt>(qs_lhs, qs_rhs);

    // Compute S
    chemist::braket::BraKet s_mn(aos, simde::type::s_e_type{}, aos);
    const auto& S = S_mod.run_as<s_pt>(s_mn);

    wf_type psi_old = psi0;
    simde::type::tensor e_old;

    density_op_type rho_hat(psi_old.orbitals(), psi_old.occupations());
    chemist::braket::BraKet P_mn(aos, rho_hat, aos);
    const auto& P = density_mod.run_as<density_pt>(P_mn);
    density_t rho_old(P, psi_old.orbitals());

    const auto max_iter = inputs.at("max iterations").value<unsigned int>();
    const auto e_tol    = inputs.at("energy tolerance").value<double>();
    const auto g_tol    = inputs.at("gradient tolerance").value<double>();
    unsigned int iter   = 0;

    auto& logger = get_runtime().logger();

    while(iter < max_iter) {
        // Step 2: Old density is used to create the new Fock operator
        // TODO: Make easier to go from many-electron to one-electron
        // TODO: template fock_pt on Hamiltonian type and only pass H_elec
        const auto& f_new = fock_mod.run_as<fock_pt>(H, rho_old);
        const auto& F_new = Fock_mod.run_as<fock_pt>(H, rho_old);

        // Step 3: New Fock operator used to compute new wavefunction/density
        const auto& psi_new =
          update_mod.run_as<update_pt<wf_type>>(f_new, psi_old);

        density_op_type rho_hat_new(psi_new.orbitals(), psi_new.occupations());
        chemist::braket::BraKet P_mn_new(aos, rho_hat_new, aos);
        const auto& P_new = density_mod.run_as<density_pt>(P_mn_new);

        density_t rho_new(P_new, psi_new.orbitals());

        // Step 4: New electronic energy
        // Step 4a: New Fock operator to new electronic Hamiltonian
        // TODO: Should just be H_core + F;
        electronic_hamiltonian H_new;
        for(std::size_t i = 0; i < H_core.size(); ++i)
            H_new.emplace_back(H_core.coefficient(i),
                               H_core.get_operator(i).clone());
        for(std::size_t i = 0; i < F_new.size(); ++i)
            H_new.emplace_back(F_new.coefficient(i),
                               F_new.get_operator(i).clone());

        // Step 4b: New electronic hamiltonian to new electronic energy
        chemist::braket::BraKet H_00(psi_new, H_new, psi_new);
        auto e_new = egy_mod.run_as<elec_egy_pt<wf_type>>(H_00);

        logger.log("SCF iteration = " + std::to_string(iter) + ":");
        logger.log("  Electronic Energy = " + e_new.to_string());

        bool converged = false;
        // Step 5: Converged?
        if(iter > 0) {
            simde::type::tensor de;
            de("") = e_new("") - e_old("");

            // Orbital gradient: FPS-SPF
            // TODO: module satisfying BraKet(aos, Commutator(F,P), aos)
            chemist::braket::BraKet F_mn(aos, f_new, aos);
            const auto& F_matrix = F_mod.run_as<fock_matrix_pt>(F_mn);
            simde::type::tensor FPS;
            FPS("m,l") = F_matrix("m,n") * P_new("n,l");
            FPS("m,l") = FPS("m,n") * S("n,l");

            simde::type::tensor SPF;
            SPF("m,l") = P_new("m,n") * F_matrix("n,l");
            SPF("m,l") = S("m,n") * SPF("n,l");

            simde::type::tensor grad;
            simde::type::tensor grad_norm;
            grad("m,n")   = FPS("m,n") - SPF("m,n");
            grad_norm("") = grad("m,n") * grad("n,m");

            Kernel e_kernel(get_runtime());
            Kernel g_kernel(get_runtime());

            using tensorwrapper::utilities::floating_point_dispatch;
            auto e_conv = floating_point_dispatch(e_kernel, de.buffer(), e_tol);
            auto g_conv =
              floating_point_dispatch(g_kernel, grad_norm.buffer(), g_tol);

            logger.log("  dE = " + de.to_string());
            logger.log("  dG = " + grad_norm.to_string());

            if(e_conv && g_conv) converged = true;
        }

        // Step 6: Not converged so reset
        e_old   = e_new;
        psi_old = psi_new;
        rho_old = rho_new;
        if(converged) break;
        ++iter;
    }
    if(iter == max_iter) throw std::runtime_error("SCF failed to converge");

    simde::type::tensor e_total;

    // e_nuclear is a double. This hack converts it to udouble (if needed)
    tensorwrapper::allocator::Eigen<double> dalloc(get_runtime());
    using tensorwrapper::types::udouble;
    tensorwrapper::allocator::Eigen<udouble> ualloc(get_runtime());

    if(ualloc.can_rebind(e_old.buffer())) {
        simde::type::tensor temp(e_old);
        auto val = dalloc.rebind(e_nuclear.buffer()).at();
        ualloc.rebind(temp.buffer()).at() = val;
        e_nuclear                         = temp;
    }

    e_total("") = e_old("") + e_nuclear("");
    auto rv     = results();
    return pt<wf_type>::wrap_results(rv, e_total, psi_old);
}

} // namespace scf::driver