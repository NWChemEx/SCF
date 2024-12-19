#include "matrix_builder.hpp"

namespace scf::matrix_builder {
namespace {

const auto desc = R"(
)";

}

using pt    = simde::aos_f_e_aos;
using ao_pt = simde::aos_op_base_aos;

MODULE_CTOR(Fock) {
    description(desc);
    satisfies_property_type<pt>();
    add_submodule<ao_pt>("Two center evaluator");
}

MODULE_RUN(Fock) {
    const auto&& [braket] = pt::unwrap_inputs(inputs);
    const auto& bra_aos   = braket.bra();
    const auto& f         = braket.op();
    const auto& ket_aos   = braket.ket();
    auto& ao_dispatcher   = submods.at("Two center evaluator");

    simde::type::tensor F;
    auto n_terms = f.size();
    if(n_terms > 0) {
        // Initialize F to zero-th term
        const auto c0    = f.coefficient(0);
        const auto& op_0 = f.get_operator(0);
        chemist::braket::BraKet term0(bra_aos, op_0, ket_aos);
        F = ao_dispatcher.run_as<ao_pt>(term0);

        using allocator_type = tensorwrapper::allocator::Eigen<double, 2>;
        auto& f_buffer       = allocator_type::rebind(F.buffer());
        f_buffer.value()     = f_buffer.value() * c0;

        for(decltype(n_terms) i = 1; i < n_terms; ++i) {
            const auto ci    = f.coefficient(i);
            const auto& op_i = f.get_operator(i);
            chemist::braket::BraKet termi(bra_aos, op_i, ket_aos);
            auto matrix = ao_dispatcher.run_as<ao_pt>(termi);

            auto& matrix_buffer = allocator_type::rebind(matrix.buffer());
            f_buffer.value() += ci * matrix_buffer.value();
        }
    }

    auto rv = results();
    return pt::wrap_results(rv, F);
}

} // namespace scf::matrix_builder