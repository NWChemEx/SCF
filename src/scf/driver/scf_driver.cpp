#include "driver.hpp"

namespace scf::driver {
namespace {

const auto desc = R"(
)";

}

using simde::type::hamiltonian;

using pt = simde::AOEnergy;

using ham_pt = simde::Convert<hamiltonian, simde::type::chemical_system>;

template<typename WfType>
using guess_pt = simde::InitialGuess<WfType>;

template<typename WfType>
using braket_t = chemist::braket::BraKet<WfType, hamiltonian, WfType>;

template<typename WfType>
using egy_pt = simde::EvaluateBraKet<braket_t<WfType>>;

template<typename WfType>
using opt_pt = simde::Optimize<egy_pt<WfType>, WfType>;

MODULE_CTOR(SCFDriver) {
    using wf_type = simde::type::rscf_wf;

    description(desc);
    satisfies_property_type<pt>();
    add_submodule<ham_pt>("Hamiltonian");
    add_submodule<guess_pt<wf_type>>("Guess");
    add_submodule<opt_pt<wf_type>>("Optimizer");
}

MODULE_RUN(SCFDriver) {
    using wf_type = simde::type::rscf_wf;

    const auto& [ao_params, sys] = pt::unwrap_inputs(inputs);
    simde::type::aos aos(ao_params);

    auto& ham_mod = submods.at("Hamiltonian");
    const auto& H = ham_mod.run_as<ham_pt>(sys);

    auto& guess_mod  = submods.at("Guess");
    const auto& Psi0 = guess_mod.run_as<guess_pt<wf_type>>(H, aos);

    auto& opt_mod = submods.at("Optimizer");
    braket_t<wf_type> H_00(Psi0, H, Psi0);
    const auto&& [e, Psi] = opt_mod.run_as<opt_pt<wf_type>>(H_00, Psi0);

    auto rv = results();
    return pt::wrap_results(rv, e);
}

} // namespace scf::driver