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

#include "fock_operator.hpp"

namespace scf::fock_operator {

template<typename DensityType>
class RestrictedVisitor : public chemist::qm_operator::OperatorVisitor {
public:
    using many_electrons = simde::type::many_electrons;
    using J_type = chemist::qm_operator::Coulomb<many_electrons, DensityType>;
    using K_type = chemist::qm_operator::Exchange<many_electrons, DensityType>;

    using T_e_term  = simde::type::T_e_type;
    using V_en_term = simde::type::V_en_type;
    using V_ee_term = simde::type::V_ee_type;
    using V_nn_term = simde::type::V_nn_type;

    RestrictedVisitor(simde::type::fock& F, const DensityType& rho) :
      m_pF_(&F), m_prho_(&rho) {}

    void run(const T_e_term& T_e) { m_pF_->emplace_back(1.0, T_e.clone()); }

    void run(const V_en_term& V_en) { m_pF_->emplace_back(1.0, V_en.clone()); }

    void run(const V_ee_term& V_ee) {
        auto es = V_ee.at<0>();
        m_pF_->emplace_back(2.0, std::make_unique<J_type>(es, *m_prho_));
        m_pF_->emplace_back(-1.0, std::make_unique<K_type>(es, *m_prho_));
    }

private:
    simde::type::fock* m_pF_;
    const DensityType* m_prho_;
};

const auto desc = R"(
Restricted Fock Operator Builder
--------------------------------

This module builds the Fock operator given a restricted density. Under this 
assumption, the Fock operator is created by mapping each term in the 
electronic Hamiltonian to itself with the exception of the electron-electron
repulsion. The electron-electron repulsion is mapped to 2J-K where J is the
Coulomb operator for the electrons interacting with a density and K is the 
exchange operator for the same density.
)";

template<typename DensityType>
TEMPLATED_MODULE_CTOR(Restricted, DensityType) {
    using pt = simde::FockOperator<DensityType>;
    satisfies_property_type<pt>();
    description(desc);
}

template<typename DensityType>
TEMPLATED_MODULE_RUN(Restricted, DensityType) {
    using pt = simde::FockOperator<DensityType>;
    using simde::type::many_electrons;

    const auto& [H, rho] = pt::unwrap_inputs(inputs);
    auto H_elec          = H.electronic_hamiltonian();

    simde::type::fock F;
    RestrictedVisitor<DensityType> visitor(F, rho);
    H_elec.visit(visitor);

    auto rv = results();
    return pt::wrap_results(rv, F);
}

template class Restricted<simde::type::e_density>;
template class Restricted<simde::type::decomposable_e_density>;

} // namespace scf::fock_operator