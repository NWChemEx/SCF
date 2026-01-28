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
#include <scf/driver/commutator.hpp>

namespace scf::driver {
namespace {

// Comparison between convergence metrics based on floating point type
struct Kernel {
    double m_tol;

    explicit Kernel(double tol) : m_tol(tol) {}

    template<typename FloatType>
    auto operator()(const std::span<FloatType>& a) {
        using tensorwrapper::types::fabs;
        return fabs(a[0]) < std::decay_t<FloatType>(m_tol);
    }
};

auto check_tolerance(const tensorwrapper::buffer::BufferBase& v, double tol) {
    Kernel kernel(tol);
    const auto& buffer = tensorwrapper::buffer::make_contiguous(v);
    return tensorwrapper::buffer::visit_contiguous_buffer(kernel, buffer);
}

// Pull out nuclear-nuclear interaction term, if there is one.
struct GrabNuclear : chemist::qm_operator::OperatorVisitor {
    using V_nn_type = simde::type::V_nn_type;

    GrabNuclear() : chemist::qm_operator::OperatorVisitor(false) {}

    void run(const V_nn_type& V_nn) { m_pv = &V_nn; }

    const V_nn_type* m_pv = nullptr;
};

const auto desc = R"(
)";

} // namespace

using tensor_t        = simde::type::tensor;
using cmos_t          = simde::type::cmos;
using diis_t          = tensorwrapper::diis::DIIS;
using density_t       = simde::type::decomposable_e_density;
using fock_pt         = simde::FockOperator<density_t>;
using density_pt      = simde::aos_rho_e_aos<cmos_t>;
using v_nn_pt         = simde::charge_charge_interaction;
using fock_matrix_pt  = simde::aos_f_e_aos;
using diagonalizer_pt = simde::GeneralizedEigenSolve;
using s_pt            = simde::aos_s_e_aos;
using simde::type::electronic_hamiltonian;

template<typename WfType>
using egy_pt = simde::eval_braket<WfType, simde::type::hamiltonian, WfType>;

template<typename WfType>
using elec_egy_pt = simde::eval_braket<WfType, electronic_hamiltonian, WfType>;

template<typename WfType>
using pt = simde::Optimize<egy_pt<WfType>, WfType>;

template<typename WfType>
using update_pt = simde::UpdateGuess<WfType>;

MODULE_CTOR(SCFLoop) {
    using wf_type = simde::type::rscf_wf;
    description(desc);
    satisfies_property_type<pt<wf_type>>();

    const unsigned int max_itr = 20;
    add_input<unsigned int>("max iterations").set_default(max_itr);
    add_input<double>("energy tolerance").set_default(1.0E-6);
    add_input<double>("density tolerance").set_default(1.0E-6);
    add_input<double>("gradient tolerance").set_default(1.0E-6);

    const std::size_t diis_sample_default = 5;
    add_input<bool>("DIIS").set_default(true);
    add_input<std::size_t>("DIIS max samples").set_default(diis_sample_default);

    add_submodule<elec_egy_pt<wf_type>>("Electronic energy");
    add_submodule<density_pt>("Density matrix");
    add_submodule<s_pt>("Overlap matrix builder");
    add_submodule<fock_matrix_pt>("Fock matrix builder");
    add_submodule<diagonalizer_pt>("Diagonalizer");
    add_submodule<fock_pt>("One-electron Fock operator");
    add_submodule<fock_pt>("Fock operator");
    add_submodule<v_nn_pt>("Charge-charge");
}

MODULE_RUN(SCFLoop) {
    using wf_type         = simde::type::rscf_wf;
    using density_op_type = simde::type::rho_e<cmos_t>;

    const auto&& [braket, psi0] = pt<wf_type>::unwrap_inputs(inputs);
    const auto max_iter = inputs.at("max iterations").value<unsigned int>();
    const auto e_tol    = inputs.at("energy tolerance").value<double>();
    const auto dp_tol   = inputs.at("density tolerance").value<double>();
    const auto g_tol    = inputs.at("gradient tolerance").value<double>();

    auto& egy_mod          = submods.at("Electronic energy");
    auto& density_mod      = submods.at("Density matrix");
    auto& diagonalizer_mod = submods.at("Diagonalizer");
    auto& fock_mod         = submods.at("One-electron Fock operator");
    auto& Fock_mod         = submods.at("Fock operator");
    auto& V_nn_mod         = submods.at("Charge-charge");
    auto& F_mod            = submods.at("Fock matrix builder");
    auto& S_mod            = submods.at("Overlap matrix builder");

    // Unwrap Braket and trial wavefunction
    // TODO: Assert bra == ket == psi0
    const auto& H      = braket.op();
    const auto& H_elec = H.electronic_hamiltonian();
    const auto& H_core = H_elec.core_hamiltonian();
    const auto& aos    = psi0.orbitals().from_space();

    // DIIS settings
    const auto diis_on = inputs.at("DIIS").value<bool>();
    const auto diis_max_samples =
      inputs.at("DIIS max samples").value<std::size_t>();
    diis_t diis(diis_max_samples);

    // Nuclear-nuclear repulsion
    GrabNuclear visitor;
    H.visit(visitor);
    bool has_nn = (visitor.m_pv != nullptr);

    // TODO: Clean up charges class to make this easier...
    tensor_t e_nuclear(0.0);
    if(has_nn) {
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
        e_nuclear = V_nn_mod.run_as<v_nn_pt>(qs_lhs, qs_rhs);
    }

    // Compute S
    chemist::braket::BraKet s_mn(aos, simde::type::s_e_type{}, aos);
    const auto& S = S_mod.run_as<s_pt>(s_mn);

    // For convergence checking
    wf_type psi_old;
    density_t rho_old;
    tensor_t e_old;
    tensor_t F_old;

    // Initialize loop
    unsigned int iter = 0;
    auto& logger      = get_runtime().logger();
    while(iter < max_iter) {
        // Step 1: Generate trial wavefunction
        wf_type psi;
        if(iter > 0) {
            // Diagonalize the Fock matrix
            const auto&& [evalues, evectors] =
              diagonalizer_mod.run_as<diagonalizer_pt>(F_old, S);

            // Construct new trial wavefunction
            cmos_t cmos(evalues, aos, evectors);
            psi = wf_type(psi_old.orbital_indices(), cmos);
        } else {
            // Use trial wavefunction provided from initial guess
            psi = psi0;
        }

        // Step 2: Construct electronic density
        density_op_type rho_hat(psi.orbitals(), psi.occupations());
        chemist::braket::BraKet P_mn(aos, rho_hat, aos);
        const auto& P = density_mod.run_as<density_pt>(P_mn);
        density_t rho(P, psi.orbitals());

        // Step 3: Construct new Fock matrix
        const auto& f_hat = fock_mod.run_as<fock_pt>(H, rho);
        chemist::braket::BraKet f_mn(aos, f_hat, aos);
        auto F = F_mod.run_as<fock_matrix_pt>(f_mn);

        // Step 4: New electronic energy
        // Step 4a: New Fock operator to new electronic Hamiltonian
        // TODO: Should just be H_core + F_hat;
        const auto& F_hat = Fock_mod.run_as<fock_pt>(H, rho);
        electronic_hamiltonian H_new;
        for(std::size_t i = 0; i < H_core.size(); ++i)
            H_new.emplace_back(H_core.coefficient(i),
                               H_core.get_operator(i).clone());
        for(std::size_t i = 0; i < F_hat.size(); ++i)
            H_new.emplace_back(F_hat.coefficient(i),
                               F_hat.get_operator(i).clone());

        // Step 4b: New electronic hamiltonian to new electronic energy
        chemist::braket::BraKet H_00(psi, H_new, psi);
        auto e = egy_mod.run_as<elec_egy_pt<wf_type>>(H_00);

        auto e_msg = "SCF iteration = " + std::to_string(iter) + ":";
        e_msg += "  Electronic Energy = " + e.to_string();
        logger.log(e_msg);

        // Step 5: Converged?
        bool converged = false;
        if(iter > 0) {
            // Change in the energy
            tensor_t de;
            de("") = e("") - e_old("");

            // Change in the density
            tensor_t dp;
            const auto& P_old = rho_old.value();
            dp("m,n")         = rho.value()("m,n") - P_old("m,n");
            auto dp_norm      = tensorwrapper::operations::infinity_norm(dp);

            // Orbital gradient: FPS-SPF
            auto grad = commutator(F, P, S);
            tensor_t grad_norm;
            grad_norm("") = grad("m,n") * grad("n,m");

            // Log convergence metrics
            logger.log("  dE = " + de.to_string());
            logger.log("  dP = " + dp_norm.to_string());
            logger.log("  dG = " + grad_norm.to_string());

            // Check for convergence
            auto e_conv  = check_tolerance(de.buffer(), e_tol);
            auto g_conv  = check_tolerance(grad_norm.buffer(), g_tol);
            auto dp_conv = check_tolerance(dp_norm.buffer(), dp_tol);
            if(e_conv && g_conv && dp_conv) converged = true;

            // If using DIIS and not converged, extrapolate new Fock matrix
            if(diis_on && !converged) { F = diis.extrapolate(F, grad); }
        } else if(diis_on) {
            // For DIIS, still need to the orbital gradient
            auto grad = commutator(F, P, S);

            // If using DIIS, extrapolate new Fock matrix
            F = diis.extrapolate(F, grad);
        }

        // Step 6: Not converged so reset
        e_old   = e;
        psi_old = psi;
        rho_old = rho;
        F_old   = F;
        if(converged) break;
        ++iter;
    }
    if(iter == max_iter) throw std::runtime_error("SCF failed to converge");

    tensor_t e_total;

    // e_nuclear is a double. This hack converts it to udouble (if needed)
    try {
        using tensorwrapper::buffer::make_contiguous;
        using tensorwrapper::types::udouble;
        const auto& e_contig = make_contiguous(e_nuclear.buffer());
        auto val = wtf::fp::float_cast<udouble>(e_contig.get_elem({}));
        std::vector<udouble> val_vector{val};
        tensorwrapper::shape::Smooth scalar{};
        tensorwrapper::buffer::Contiguous nuc_buffer(val_vector, scalar);
        e_nuclear = tensor_t(scalar, std::move(nuc_buffer));
    } catch(...) {
        // Nothing to do here
    }

    e_total("") = e_old("") + e_nuclear("");

    auto rv = results();
    return pt<wf_type>::wrap_results(rv, e_total, psi_old);
}

} // namespace scf::driver
