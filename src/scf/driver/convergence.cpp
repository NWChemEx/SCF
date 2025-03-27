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


MODULE_CTOR(ConvergenceMod) {
    using wf_type = simde::type::rscf_wf;
    description(desc);
    satisfies_property_type<pt<wf_type>>();

    const unsigned int max_itr = 20;
    add_input<unsigned int>("max iterations").set_default(max_itr);
    add_input<double>("energy tolerance").set_default(1.0E-6);
    add_input<double>("density tolerance").set_default(1.0E-6);
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

MODULE_RUN(ConvergenceMod) {
    using wf_type               = simde::type::rscf_wf;
    using density_op_type       = simde::type::rho_e<simde::type::cmos>;
    const auto&& [braket, psi0] = pt<wf_type>::unwrap_inputs(inputs);
        // Step 5: Converged?
            // Change in the energy
            simde::type::tensor de;
            de("") = e_new("") - e_old("");

            // Change in the density
            simde::type::tensor dp;
            dp("m,n")    = rho_new.value()("m,n") - rho_old.value()("m,n");
            auto dp_norm = tensorwrapper::operations::infinity_norm(dp);

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

            Kernel k(get_runtime());

            using tensorwrapper::utilities::floating_point_dispatch;
            auto e_conv = floating_point_dispatch(k, de.buffer(), e_tol);
            auto g_conv = floating_point_dispatch(k, grad_norm.buffer(), g_tol);
            auto dp_conv = floating_point_dispatch(k, dp_norm.buffer(), dp_tol);

            logger.log("  dE = " + de.to_string());
            logger.log("  dP = " + dp_norm.to_string());
            logger.log("  dG = " + grad_norm.to_string());

            if(e_conv && g_conv && dp_conv) converged = true;
        }

        // Step 6: Not converged so reset
        e_old   = e_new;
        psi_old = psi_new;
        rho_old = rho_new;
        if(converged) break;
        ++iter;
    }
