#include "guess.hpp"

namespace scf::guess {

using rscf_wf    = simde::type::rscf_wf;
using density_t  = simde::type::decomposable_e_density;
using pt         = simde::InitialGuess<rscf_wf>;
using fock_op_pt = simde::FockOperator<density_t>;
using update_pt  = simde::UpdateGuess<rscf_wf>;

const auto desc = R"(
Core Guess
----------

TODO: Write me!!!
)";

MODULE_CTOR(Core) {
    description(desc);
    satisfies_property_type<pt>();
    add_submodule<fock_op_pt>("Build Fock operator");
    add_submodule<update_pt>("Guess updater");
}

MODULE_RUN(Core) {
    // const auto&& [H, aos] = pt::unwrap_inputs(inputs);

    // // Step 1: Build Fock Operator with zero density
    // density_t rho;
    // auto& fock_op_mod = submods.at("Build Fock operator");
    // const auto& F     = fock_op_mod.run_as<fock_op_pt>(rho);

    // // Step 2: Update guess Call module to update guess
    // simde::type::cmos cmos(zero_vector, aos, Tensor{});
    // typename rscf_wf::orbital_index_set_type occupations(aos.size(), 0);
    // // Fill in occupations

    // rscf_wf zero_guess(occupations aos, cmos);
    // auto& update_mod = submods.at("Guess updater");
    // const auto& Psi0 = update_mod.run_as<update_pt>(F, zero_guess);

    // auto rv = results();
    // return pt::wrap_results(rv, Psi0);
}

} // namespace scf::guess