// #include "scf_modules.hpp"

// namespace scf {
// namespace {

// const auto desc = R"(
// )";

// }

// using simde::type::electronic_hamiltonian;
// using simde::type::op_base_type;

// template<typename WfType>
// using egy_pt = simde::type::eval_braket<WfType, electronic_hamiltonian,
// WfType>;

// template<typename WfType>
// using det_pt = simde::type::eval_braket<WfType, op_base_type, WfType>;

// template<typename WfType>
// using pt = simde::Optimize<egy_pt<WfType>, WfType>;

// template<typename WfType>
// using update_pt = simde::UpdateGuess<WfType>;

// using density_t = simde::type::decomposable_e_density

//   using fock_pt = simde::FockOperator<density_t>;

// using density_pt = simde::aos_rho_e_aos;

// MODULE_CTOR(SCFLoop) {
//     using wf_type = simde::type::rscf_wf;
//     description(desc);
//     satisfies_property_type<pt<wf_type>>();

//     add_submodule<egy_pt<wf_type>>("Determinant energy");
//     add_submodule<density_pt>("Density matrix");
//     add_submodule<update_pt<wf_type>>("Guess update");
//     add_submodule<fock_pt>("Fock operator");
// }

// MODULE_RUN(SCFLoop) {
//     using wf_type               = simde::type::rscf_wf;
//     const auto&& [braket, psi0] = pt<wf_type>::unwrap_inputs(inputs);
//     // TODO: Assert bra == ket == psi0
//     const auto&& H     = braket.op();
//     const auto& H_core = H.core_hamiltonian();

//     auto& density_mod = submods.at("Density matrix");
//     auto& egy_mod     = submods.at("Determinant energy");
//     auto& update_mod  = submods.at("Guess update");
//     auto& fock_mod    = submods.at("Fock operator");

//     double e_old                = 0;
//     const unsigned int max_iter = 100;
//     unsigned int iter           = 0;
//     while(iter < max_iter) {
//         // Step 1: Build old density
//         simde::type::rho_e rho_hat(psi0.orbitals(), psi0.occupations());

//         auto P_mn           = chemist::simde::BraKet(aos, rho_hat, aos);
//         const auto& rho_old = density_mod.run_as<density_pt>(P_mn);

//         // Step 2: Old density becomes new Fock operator
//         const auto& f_new = fock_mod.run_as<fock_pt>(H, rho_old);
//         const auto& F_new =
//           Fock_mod.run_as<fock_pt>(H, rho_old); // This is a hack

//         // Step 3: New Fock operator becomes the new guess
//         auto new_psi0 = update_mod.run_as<update_pt<wf_type>>(f_new, psi0);

//         // Step 4: New Fock operator to new Hamiltonian
//         electronic_hamiltonian H_new;
//         for(std::size_t i = 0; i < H_core.size(); ++i)
//             H_new.insert(H_core.coefficient(i), H_core.get_op(i));
//         for(std::size_t i = 0; i < F_new.size(); ++i)
//             H_new.insert(F_new.coefficient(i), F_new.get_op(i));

//         // Step 5: New hamiltonian to new energy
//         chemist::braket::BraKet H_00(new_psi0, H_new, new_psi0);
//         auto e_new = egy_mod.run_as<egy_pt<wf_type>>(H_00);

//         // Step 6: Converged?

//         // Step 7: Not converged reset
//         e_old = e_new;
//         psi0  = new_psi0;
//         ++iter;
//     }
//     auto rv = results();
//     return pt<wf_type>::wrap_results(rv, e_old, psi0);
// }

// } // namespace scf