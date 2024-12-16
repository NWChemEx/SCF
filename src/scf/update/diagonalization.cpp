#include "update.hpp"

namespace scf::update {

using rscf_wf         = simde::type::rscf_wf;
using fock_matrix_pt  = simde::aos_f_e_aos;
using pt              = simde::UpdateGuess<rscf_wf>;
using diagonalizer_pt = simde::GeneralizedEigenSolve;

const auto desc = R"(
)";

MODULE_CTOR(Diagonalization) {
    description(desc);
    satisfies_property_type<pt>();
    add_submodule<fock_matrix_pt>("Fock matrix builder");
    add_submodule<diagonalizer_pt>("Diagonalizer");
}

MODULE_RUN(Diagonalization) {
    // const auto&& [F, old_guess] = pt::unwrap_inputs(inputs);

    // // To one-electron F
    // const auto& aos       = old_guess.orbitals().from_space();
    // auto& fock_matrix_mod = submods.at("Fock matrix builder");
    // const auto& f_matrix  = fock_matrix_mod.run_as<fock_matrix_pt>(aos, f,
    // aos);

    // // Compute S
    // Tensor s;

    // // Diagonalize
    // const auto& diagonalizer_mod = submods.at("Diagonalizer");
    // const auto&& [evalues, evectors] =
    //   diagonalizer_mod.run_as<diagonalizer_pt>(f_matrix, s);

    // // Create new guess
    // simde::type::cmos cmos(evalues, aos, evectors);
    // rscf_wf new_guess(old_guess.orbital_indices(), cmos);

    // auto rv = results();
    // return pt::wrap_results(rv, new_guess);
}

} // namespace scf::update