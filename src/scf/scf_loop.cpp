// #include "scf_modules.hpp"

// namespace scf {
// namespace {
// const auto desc = R"(
// )";

// }

// using simde::type::hamiltonian;

// template<typename WfType>
// using braket_t = chemist::braket::BraKet<WfType, hamiltonian, WfType>;

// template<typename WfType>
// using egy_pt = simde::EvaluateBraKet<braket_t<WfType>>;

// template<typename WfType>
// using pt = simde::Optimize<egy_pt<WfType>, WfType>;

// template<typename WfType>
// using update_pt = simde::UpdateGuess<WfType>;

// MODULE_CTOR(SCFLoop) {
//     using wf_type = simde::type::rscf_wf;
//     description(desc);
//     satisfies_property_type<pt<wf_type>>();

//     add_submodule<egy_pt<wf_type>>("Determinant integrals");
//     add_submodule<update_pt<wf_type>>("Guess update");
// }

// MODULE_RUN(SCFLoop) {
//     using wf_type               = simde::type::rscf_wf;
//     const auto&& [braket, psi0] = pt<wf_type>::unwrap_inputs(inputs);
//     // TODO: Assert bra == ket == psi0
//     const auto&& H = braket.op();

//     auto& egy_mod    = submods.at("Determinant integrals");
//     auto& update_mod = submods.at("Guess update");

//     auto e_old     = egy_mod.run_as<egy_pt<wf_type>>(braket);
//     auto f_old     = fock_mod.run_as<>(H, rho);
//     auto new_guess = update_mod.run_as<update_pt<wf_type>>(f_old, psi0);
//     braket_t<wf_type> new_braket(new_guess, H, new_guess);
//     auto e_new = egy_mod.run_as<egy_pt<wf_type>>()

//                    auto rv = results();
//     return pt<wf_type>::wrap_results(rv, e, psi);
// }

// } // namespace scf