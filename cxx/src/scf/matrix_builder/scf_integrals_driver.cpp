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
AO Integrals Driver
-------------------
)";

}

using pluginplay::type::submodule_map;
using simde::type::aos;
using simde::type::tensor;

using pt       = simde::aos_op_base_aos;
using f_e_pt   = simde::aos_f_e_aos;
using rho_e_pt = simde::aos_rho_e_aos<simde::type::cmos>;
using xc_e_pt  = simde::aos_xc_e_aos;

class AODispatcher : public chemist::qm_operator::OperatorVisitor {
private:
    using base_type = chemist::qm_operator::OperatorVisitor;

public:
    using f_e_type     = simde::type::fock;
    using rho_e_type   = simde::type::rho_e<simde::type::cmos>;
    using xc_e_type    = simde::type::xc_e_type;
    using submods_type = pluginplay::type::submodule_map;

    AODispatcher(const aos& bra, const aos& ket, submodule_map& submods,
                 tensor& t) :
      base_type(false),
      m_pbra_(&bra),
      m_pket_(&ket),
      m_evaluated_(false),
      m_psubmods_(&submods),
      m_ptensor_(&t) {}

    // void run(const f_e_type& f_e) {
    //     chemist::braket::BraKet input(*m_pbra_, f_e, *m_pket_);
    //     const auto key = "Fock matrix";
    //     *m_ptensor_    = m_psubmods_->at(key).run_as<f_e_pt>(input);
    //     m_evaluated_   = true;
    // }

    void run(const rho_e_type& rho_e) {
        chemist::braket::BraKet input(*m_pbra_, rho_e, *m_pket_);
        const auto key = "Density matrix";
        *m_ptensor_    = m_psubmods_->at(key).run_as<rho_e_pt>(input);
        m_evaluated_   = true;
    }

    void run(const xc_e_type& xc_e) {
        chemist::braket::BraKet input(*m_pbra_, xc_e, *m_pket_);
        const auto key = "XC Potential";
        *m_ptensor_    = m_psubmods_->at(key).run_as<xc_e_pt>(input);
        m_evaluated_   = true;
    }

    bool evaluated() const noexcept { return m_evaluated_; }

private:
    const aos* m_pbra_;
    const aos* m_pket_;
    bool m_evaluated_ = false;
    submodule_map* m_psubmods_;
    tensor* m_ptensor_;
};

MODULE_CTOR(SCFIntegralsDriver) {
    description(desc);
    satisfies_property_type<pt>();
    add_submodule<pt>("Fundamental matrices");
    add_submodule<rho_e_pt>("Density matrix");
    add_submodule<xc_e_pt>("XC Potential");
}

MODULE_RUN(SCFIntegralsDriver) {
    const auto&& [braket] = pt::unwrap_inputs(inputs);
    const auto& bra       = braket.bra();
    const auto& op        = braket.op();
    const auto& ket       = braket.ket();

    tensor t;
    AODispatcher visitor(bra, ket, submods, t);
    op.visit(visitor);
    if(!visitor.evaluated()) {
        t = submods.at("Fundamental matrices").run_as<pt>(braket);
    }

    auto rv = results();
    return pt::wrap_results(rv, std::move(t));
}

} // namespace scf::matrix_builder
