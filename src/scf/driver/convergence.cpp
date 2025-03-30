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
#include <driver/convergence.hpp>


namespace scf::driver {

template<typename KernelType>
TEMPLATED_MODULE_CTOR(ConvergenceMod, KernelType) {
    satisfies_property_type<scf::ConvergenceProp<KernelType>>();

    add_input<double>("energy tolerance").set_default(1.0E-6);
    add_input<double>("density tolerance").set_default(1.0E-6);
    add_input<double>("gradient tolerance").set_default(1.0E-6);
    add_submodule<simde::aos_f_e_aos>("Fock matrix builder");
}

template<typename KernelType>
TEMPLATED_MODULE_RUN(ConvergenceMod, KernelType) {
    const auto&& [e_new, e_old, rho_new, rho_old, P_new, S, f_new, aos, k, e_tol, dp_tol, g_tol] = ConvergenceProp<KernelType>::unwrap_inputs(inputs);

    bool converged = false;

    auto& F_mod = submods.at("Fock matrix builder");

    simde::type::tensor de;
    de("") = e_new("") - e_old("");

    // Change in the density
    simde::type::tensor dp;
    dp("m,n")    = rho_new.value()("m,n") - rho_old.value()("m,n");
    auto dp_norm = tensorwrapper::operations::infinity_norm(dp);

    // Orbital gradient: FPS-SPF
    // TODO: module satisfying BraKet(aos, Commutator(F,P), aos)
    chemist::braket::BraKet F_mn(aos, f_new, aos);
    const auto& F_matrix = F_mod.run_as<simde::aos_f_e_aos>(F_mn);
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

    using tensorwrapper::utilities::floating_point_dispatch;
    auto e_conv = floating_point_dispatch(k, de.buffer(), e_tol);
    auto g_conv = floating_point_dispatch(k, grad_norm.buffer(), g_tol);
    auto dp_conv = floating_point_dispatch(k, dp_norm.buffer(), dp_tol);

    // logger.log("  dE = " + de.to_string());
    // logger.log("  dP = " + dp_norm.to_string());
    // logger.log("  dG = " + grad_norm.to_string());

    if(e_conv && g_conv && dp_conv) converged = true;

    auto rv = results();
    return scf::ConvergenceProp<KernelType>::wrap_inputs(rv, converged);
}
}
    //     // Step 6: Not converged so reset
    //     e_old   = e_new;
    //     psi_old = psi_new;
    //     rho_old = rho_new;
    //     if(converged) break;
    //     ++iter;
    // }
