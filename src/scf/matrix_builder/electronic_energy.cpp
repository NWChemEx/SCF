#include "matrix_builder.hpp"

namespace scf::matrix_builder {

namespace {

const auto desc = R"(
)";

}

using simde::type::electronic_hamiltonian;

template<typename WFType>
using pt = simde::eval_braket<WFType, electronic_hamiltonian, WFType>;

template<typename WFType>
using det_pt = simde::eval_braket<WFType, simde::type::op_base_type, WFType>;

MODULE_CTOR(ElectronicEnergy) {
    using wf_type = simde::type::rscf_wf;
    description(desc);
    satisfies_property_type<pt<wf_type>>();
    add_submodule<det_pt<wf_type>>("determinant driver");
}

MODULE_RUN(ElectronicEnergy) {
    using wf_type       = simde::type::rscf_wf;
    const auto&& [H_ij] = pt<wf_type>::unwrap_inputs(inputs);
    const auto& bra     = H_ij.bra();
    const auto& H       = H_ij.op();
    const auto& ket     = H_ij.ket();

    auto& mod = submods.at("determinant driver");
    double e  = 0.0;

    auto n_ops = H.size();
    for(decltype(n_ops) i = 0; i < n_ops; ++i) {
        const auto& ci  = H.coefficient(i);
        const auto& O_i = H.get_operator(i);

        chemist::braket::BraKet O_ij(bra, O_i, ket);
        e += ci * mod.run_as<det_pt<wf_type>>(O_ij);
    }

    auto rv = results();
    return pt<wf_type>::wrap_results(rv, e);
}

} // namespace scf::matrix_builder