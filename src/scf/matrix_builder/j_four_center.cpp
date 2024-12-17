#include "matrix_builder.hpp"

namespace scf::matrix_builder {

using pt    = simde::aos_j_e_aos;
using pt_4c = simde::ERI4;

auto desc = R"(
Four-Center J Builder
---------------------
)";

MODULE_CTOR(JFourCenter) {
    description(desc);
    satisfies_property_type<pt>();
    add_submodule<pt_4c>("Four-center ERI");
}

MODULE_RUN(JFourCenter) {
    const auto&& [braket] = pt::unwrap_inputs(inputs);
    // TODO: avoid copying AOs
    simde::type::aos bra_e0 = braket.bra();
    const auto& j_hat       = braket.op();
    simde::type::aos ket_e0 = braket.ket();
    const auto& rho         = j_hat.rhs_particle();
    simde::type::aos aos_e1 = rho.basis_set();
    const auto& p           = rho.value();
    auto& eri_mod           = submods.at("Four-center ERI");

    // auto aos2_v_aos2 = (aos_e1 * aos_e1 | v_ee | bra_e0 * ket_e0);
    simde::type::v_ee_type v_ee;
    simde::type::aos_squared e1_pair(aos_e1, aos_e1);
    simde::type::aos_squared e2_pair(bra_e0, ket_e0);
    chemist::braket::BraKet aos2_v_aos2(e1_pair, v_ee, e2_pair);
    const auto& I = eri_mod.run_as<pt_4c>(std::move(aos2_v_aos2));

    // This goes away when j("m,n") = p("l,s")*I("l,s,m,n") works
    // {
    using eri_alloc  = tensorwrapper::allocator::Eigen<double, 4>;
    using rho_alloc  = tensorwrapper::allocator::Eigen<double, 2>;
    using index_pair = Eigen::IndexPair<int>;
    using array_t    = Eigen::array<index_pair, 2>;

    const auto& I_eigen  = eri_alloc::rebind(I.buffer()).value();
    const auto& p_buffer = rho_alloc::rebind(p.buffer());
    const auto& p_eigen  = p_buffer.value();

    array_t contract_modes{index_pair{0, 0}, index_pair(1, 1)};
    using eigen_buffer_type = typename rho_alloc::eigen_buffer_type;
    using eigen_tensor_type = typename eigen_buffer_type::data_type;

    eigen_tensor_type j_eigen = I_eigen.contract(p_eigen, contract_modes);
    auto j_buffer = std::make_unique<eigen_buffer_type>(std::move(j_eigen),
                                                        p_buffer.layout());

    using logical_type = tensorwrapper::layout::Logical;
    auto playout       = p.logical_layout().clone_as<logical_type>();
    simde::type::tensor j(std::move(playout), std::move(j_buffer));
    //}
    auto rv = results();
    return pt::wrap_results(rv, std::move(j));
}

} // namespace scf::matrix_builder