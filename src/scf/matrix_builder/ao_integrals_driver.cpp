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
using t_e_pt   = simde::aos_t_e_aos;
using v_en_pt  = simde::aos_v_en_aos;
using j_e_pt   = simde::aos_j_e_aos;
using k_e_pt   = simde::aos_k_e_aos;
using f_e_pt   = simde::aos_f_e_aos;
using rho_e_pt = simde::aos_rho_e_aos<simde::type::cmos>;

class AODispatcher : public chemist::qm_operator::OperatorVisitor {
public:
    using t_e_type   = simde::type::t_e_type;
    using v_en_type  = simde::type::v_en_type;
    using j_e_type   = simde::type::j_e_type;
    using k_e_type   = simde::type::k_e_type;
    using f_e_type   = simde::type::fock;
    using rho_e_type = simde::type::rho_e<simde::type::cmos>;

    using submods_type = pluginplay::type::submodule_map;

    AODispatcher(const aos& bra, const aos& ket, submodule_map& submods,
                 tensor& t) :
      m_pbra_(&bra), m_pket_(&ket), m_psubmods_(&submods), m_ptensor_(&t) {}

    void run(const t_e_type& t_e) {
        chemist::braket::BraKet input(*m_pbra_, t_e, *m_pket_);
        *m_ptensor_ = m_psubmods_->at("Kinetic").run_as<t_e_pt>(input);
    }

    void run(const v_en_type& v_en) {
        chemist::braket::BraKet input(*m_pbra_, v_en, *m_pket_);
        const auto key = "Electron-Nuclear attraction";
        *m_ptensor_    = m_psubmods_->at(key).run_as<v_en_pt>(input);
    }

    void run(const j_e_type& j_e) {
        chemist::braket::BraKet input(*m_pbra_, j_e, *m_pket_);
        const auto key = "Coulomb matrix";
        *m_ptensor_    = m_psubmods_->at(key).run_as<j_e_pt>(input);
    }

    void run(const k_e_type& k_e) {
        chemist::braket::BraKet input(*m_pbra_, k_e, *m_pket_);
        const auto key = "Exchange matrix";
        *m_ptensor_    = m_psubmods_->at(key).run_as<k_e_pt>(input);
    }

    // void run(const f_e_type& f_e) {
    //     chemist::braket::BraKet input(*m_pbra_, f_e, *m_pket_);
    //     const auto key = "Fock matrix";
    //     *m_ptensor_    = m_psubmods_->at(key).run_as<f_e_pt>(input);
    // }

    void run(const rho_e_type& rho_e) {
        chemist::braket::BraKet input(*m_pbra_, rho_e, *m_pket_);
        const auto key = "Density matrix";
        *m_ptensor_    = m_psubmods_->at(key).run_as<rho_e_pt>(input);
    }

private:
    const aos* m_pbra_;
    const aos* m_pket_;
    submodule_map* m_psubmods_;
    tensor* m_ptensor_;
};

MODULE_CTOR(AOIntegralsDriver) {
    description(desc);
    satisfies_property_type<pt>();
    add_submodule<t_e_pt>("Kinetic");
    add_submodule<v_en_pt>("Electron-Nuclear attraction");
    add_submodule<j_e_pt>("Coulomb matrix");
    add_submodule<k_e_pt>("Exchange matrix");
    // add_submodule<f_e_pt>("Fock matrix");
    add_submodule<rho_e_pt>("Density matrix");
}

MODULE_RUN(AOIntegralsDriver) {
    const auto&& [braket] = pt::unwrap_inputs(inputs);
    const auto& bra       = braket.bra();
    const auto& op        = braket.op();
    const auto& ket       = braket.ket();

    tensor t;
    AODispatcher visitor(bra, ket, submods, t);
    op.visit(visitor);

    auto rv = results();
    return pt::wrap_results(rv, std::move(t));
}

} // namespace scf::matrix_builder