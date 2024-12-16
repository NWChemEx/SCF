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

using namespace chemist::qm_operator;

template<typename DensityType, typename ElectronType>
class RestrictedVisitor : public chemist::qm_operator::OperatorVisitor {
public:
    using T_e_term  = simde::type::T_e_type;
    using V_en_term = simde::type::V_en_type;
    using V_ee_term = simde::type::V_ee_type;
    using V_nn_term = simde::type::V_nn_type;

    RestrictedVisitor(simde::type::fock& F, const DensityType& rho) :
      m_pF_(&F), m_prho_(&rho) {}

    void run(const T_e_term& T_e) {
        auto t = std::make_unique<Kinetic<ElectronType>>(get_e_(T_e));
        m_pF_->emplace_back(1.0, std::move(t));
    }

    void run(const V_en_term& V_en) {
        auto rhs     = V_en.rhs_particle().as_nuclei();
        using v_type = Coulomb<ElectronType, simde::type::nuclei>;
        auto v       = std::make_unique<v_type>(get_e_(V_en), rhs);
        m_pF_->emplace_back(1.0, std::move(v));
    }

    void run(const V_ee_term& V_ee) {
        if(*m_prho_ == DensityType{}) return;
        using j_type = Coulomb<ElectronType, DensityType>;
        using k_type = Exchange<ElectronType, DensityType>;

        auto j = std::make_unique<j_type>(get_e_(V_ee), *m_prho_);
        auto k = std::make_unique<k_type>(get_e_(V_ee), *m_prho_);
        m_pF_->emplace_back(2.0, std::move(j));
        m_pF_->emplace_back(-1.0, std::move(k));
    }

private:
    template<typename T>
    auto get_e_(T&& op) {
        if constexpr(std::is_same_v<ElectronType, simde::type::electron>) {
            return simde::type::electron{};
        } else {
            return op.template at<0>();
        }
    }

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

N.b. Empty densities are treated as zero densities and this module will skip
adding 2J-K if provided an empty density.
)";

template<typename DensityType, typename ElectronType>
TEMPLATED_MODULE_CTOR(Restricted, DensityType, ElectronType) {
    using pt = simde::FockOperator<DensityType>;
    satisfies_property_type<pt>();
    description(desc);
}

template<typename DensityType, typename ElectronType>
TEMPLATED_MODULE_RUN(Restricted, DensityType, ElectronType) {
    using pt = simde::FockOperator<DensityType>;

    const auto& [H, rho] = pt::unwrap_inputs(inputs);
    auto H_elec          = H.electronic_hamiltonian();

    simde::type::fock F;
    RestrictedVisitor<DensityType, ElectronType> visitor(F, rho);
    H_elec.visit(visitor);

    auto rv = results();
    return pt::wrap_results(rv, F);
}

#define RESTRICTED(density)                                          \
    template class Restricted<density, simde::type::many_electrons>; \
    template class Restricted<density, simde::type::electron>

RESTRICTED(simde::type::e_density);
RESTRICTED(simde::type::decomposable_e_density);

#undef RESTRICTED

} // namespace scf::fock_operator