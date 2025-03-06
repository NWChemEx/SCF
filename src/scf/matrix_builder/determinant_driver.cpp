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

#include "matrix_builder.hpp"

namespace scf::matrix_builder {

namespace {

const auto desc = R"(
)";

using ao_pt   = simde::aos_op_base_aos;
using fock_pt = simde::aos_f_e_aos;

template<typename WFType>
using pt = simde::eval_braket<WFType, simde::type::op_base_type, WFType>;

template<typename WFType>
using density_op_t = simde::type::rho_e<typename WFType::orbital_space_type>;

using ao_type     = simde::type::aos;
using op_type     = simde::type::op_base_type;
using braket_base = simde::type::braket<ao_type, op_type, ao_type>;

template<typename WFType>
class DeterminantDispatcher : public chemist::qm_operator::OperatorVisitor {
public:
    using T_e_type  = simde::type::T_e_type;
    using V_en_type = simde::type::V_en_type;
    using J_e_type  = simde::type::J_e_type;
    using K_e_type  = simde::type::K_e_type;

    using submodule_map = pluginplay::type::submodule_map;

    DeterminantDispatcher(const WFType& bra, const WFType& ket,
                          submodule_map& submods) :
      m_pbra_(&bra), m_pket_(&ket), m_psubmods_(&submods) {}

    void run(const T_e_type& T_e) {
        simde::type::electron e;
        simde::type::t_e_type t_e(e);
        run_(t_e);
    }

    void run(const V_en_type& V_en) {
        simde::type::electron e;
        simde::type::v_en_type v_en(e, V_en.rhs_particle().as_nuclei());
        run_(v_en);
    }

    void run(const J_e_type& J_e) {
        simde::type::electron e;
        simde::type::j_e_type j_e(e, J_e.rhs_particle());
        run_(j_e);
    }

    void run(const K_e_type& K_e) {
        simde::type::electron e;
        simde::type::k_e_type k_e(e, K_e.rhs_particle());
        run_(k_e);
    }

    template<typename OpType>
    void run_(const OpType& o) {
        const auto& aos = m_pbra_->orbitals().from_space();
        braket_base o_mn(aos, o, aos);
        auto& mod     = m_psubmods_->at("Two center evaluator");
        const auto& t = mod.run_as<ao_pt>(o_mn);
        m_pt          = t;
    }

    simde::type::tensor m_pt;

private:
    const WFType* m_pbra_;
    const WFType* m_pket_;
    submodule_map* m_psubmods_;
};

} // namespace

MODULE_CTOR(DeterminantDriver) {
    using wf_type = simde::type::rscf_wf;
    description(desc);
    satisfies_property_type<pt<wf_type>>();
    add_submodule<ao_pt>("Two center evaluator");
    // For now Fock matrix is special
    add_submodule<fock_pt>("Fock matrix");
}

MODULE_RUN(DeterminantDriver) {
    using wf_type         = simde::type::rscf_wf;
    const auto&& [braket] = pt<wf_type>::unwrap_inputs(inputs);
    const auto& bra       = braket.bra();
    const auto& op        = braket.op();
    const auto& ket       = braket.ket();

    // Step 1: Build density matrix
    const auto& mos = bra.orbitals();
    const auto& aos = mos.from_space();
    const auto& occ = bra.occupations();

    density_op_t<wf_type> rho_hat(mos, occ);
    braket_base p_mn(aos, rho_hat, aos);

    auto& mod       = submods.at("Two center evaluator");
    const auto& rho = mod.run_as<ao_pt>(p_mn);

    // Step 2: Build the other matrix
    DeterminantDispatcher<wf_type> visitor(bra, ket, submods);
    op.visit(visitor);
    const auto& t = visitor.m_pt;

    // Step 3: Contract
    tensorwrapper::Tensor x;
    x("") = rho("i,j") * t("i,j");

    auto rv = results();
    return pt<wf_type>::wrap_results(rv, x);
}

} // namespace scf::matrix_builder
