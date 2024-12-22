#include "driver.hpp"

namespace scf::driver {
namespace {

const auto desc = R"(
)";

}

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

    add_submodule<elec_egy_pt<wf_type>>("Electronic energy");
    add_submodule<density_pt>("Density matrix");
    add_submodule<update_pt<wf_type>>("Guess update");
    add_submodule<fock_pt>("One-electron Fock operator");
    add_submodule<fock_pt>("Fock operator");
    add_submodule<v_nn_pt>("Charge-charge");
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

    wf_type psi_old             = psi0;
    double e_old                = 0;
    const unsigned int max_iter = 3;
    unsigned int iter           = 0;

    while(iter < max_iter) {
        // Step 2: Build old density
        density_op_type rho_hat(psi_old.orbitals(), psi_old.occupations());
        chemist::braket::BraKet P_mn(aos, rho_hat, aos);
        const auto& P = density_mod.run_as<density_pt>(P_mn);
        density_t rho_old(P, psi_old.orbitals());

        // Step 3: Old density is used to create the new Fock operator
        // TODO: Make easier to go from many-electron to one-electron
        // TODO: template fock_pt on Hamiltonian type and only pass H_elec
        const auto& f_new = fock_mod.run_as<fock_pt>(H, rho_old);
        const auto& F_new = Fock_mod.run_as<fock_pt>(H, rho_old);

        // Step 4: New Fock operator is used to compute the new wavefunction
        auto new_psi = update_mod.run_as<update_pt<wf_type>>(f_new, psi_old);
        const auto& new_cmos  = new_psi.orbitals();
        const auto& new_evals = new_cmos.diagonalized_matrix();
        const auto& new_c     = new_cmos.transform();

        // Step 5: New electronic energy
        // Step 5a: New Fock operator to new electronic Hamiltonian
        // TODO: Should just be H_core + F;
        electronic_hamiltonian H_new;
        for(std::size_t i = 0; i < H_core.size(); ++i)
            H_new.emplace_back(H_core.coefficient(i),
                               H_core.get_operator(i).clone());
        for(std::size_t i = 0; i < F_new.size(); ++i)
            H_new.emplace_back(F_new.coefficient(i),
                               F_new.get_operator(i).clone());

        // Step 5b: New electronic hamiltonian to new electronic energy
        chemist::braket::BraKet H_00(new_psi, H_new, new_psi);
        auto e_new = egy_mod.run_as<elec_egy_pt<wf_type>>(H_00);

        // Step 6: Converged?
        std::cout << "E_new: " << e_new + e_nuclear << std::endl;
        std::cout << "E_old: " << e_old + e_nuclear << std::endl;
        std::cout << e_new - e_old << std::endl;

        // Step 7: Not converged so reset
        e_old   = e_new;
        psi_old = new_psi;
        ++iter;
    }
    auto rv = results();
    return pt<wf_type>::wrap_results(rv, e_old + e_nuclear, psi_old);
}

} // namespace scf::driver