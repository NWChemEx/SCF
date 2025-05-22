/*
 * Copyright 2025 NWChemEx-Project
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
#include "fock_operator_common.hpp"

namespace scf::fock_operator {
using xc_functional = chemist::qm_operator::xc_functional;

template<typename DensityType, typename ElectronType>
class RKohnShamVisitor : public FockOperatorCommon<ElectronType> {
private:
    using base_type = FockOperatorCommon<ElectronType>;

public:
    using V_ee_term = simde::type::V_ee_type;

    RKohnShamVisitor(xc_functional func, simde::type::fock& F,
                     const DensityType& rho) :
      base_type(F), m_func_(func), m_prho_(&rho) {}

    using base_type::run;
    void run(const V_ee_term& V_ee) override {
        if(*m_prho_ == DensityType{}) return;
        using j_type  = Coulomb<ElectronType, DensityType>;
        using xc_type = ExchangeCorrelation<ElectronType, DensityType>;

        const auto& electrons = this->get_e_(V_ee);
        auto j                = std::make_unique<j_type>(electrons, *m_prho_);
        auto xc = std::make_unique<xc_type>(m_func_, electrons, *m_prho_);
        this->m_pF_->emplace_back(2.0, std::move(j));
        this->m_pF_->emplace_back(1.0, std::move(xc));
    }

private:
    xc_functional m_func_;
    const DensityType* m_prho_;
};

const auto desc = R"(
Restricted Kohn-Sham Operator Builder
-------------------------------------

This module builds the Kohn-Sham operator given a restricted density. Under this 
assumption, the Kohn-Sham operator is created by mapping each term in the 
electronic Hamiltonian to itself with the exception of the electron-electron
repulsion. The electron-electron repulsion is mapped to 2J + X where J is the
Coulomb operator for the electrons interacting with a density and X is the 
exchange-correlation potential for the same density.

N.b. Empty densities are treated as zero densities and this module will skip
adding 2J + X if provided an empty density.
)";

template<typename DensityType, typename ElectronType>
TEMPLATED_MODULE_CTOR(RKohnSham, DensityType, ElectronType) {
    using pt = simde::FockOperator<DensityType>;
    satisfies_property_type<pt>();
    description(desc);

    add_input<xc_functional>("XC Potential");
}

template<typename DensityType, typename ElectronType>
TEMPLATED_MODULE_RUN(RKohnSham, DensityType, ElectronType) {
    using pt = simde::FockOperator<DensityType>;

    const auto& [H, rho] = pt::unwrap_inputs(inputs);
    auto H_elec          = H.electronic_hamiltonian();
    auto func            = inputs.at("XC Potential").value<xc_functional>();
    simde::type::fock F;
    RKohnShamVisitor<DensityType, ElectronType> visitor(func, F, rho);
    H_elec.visit(visitor);

    auto rv = results();
    return pt::wrap_results(rv, F);
}

#define RESTRICTED(density)                                          \
    template struct RKohnSham<density, simde::type::many_electrons>; \
    template struct RKohnSham<density, simde::type::electron>

RESTRICTED(simde::type::e_density);
RESTRICTED(simde::type::decomposable_e_density);

#undef RESTRICTED

} // namespace scf::fock_operator