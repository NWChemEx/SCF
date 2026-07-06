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
#include "fock_operator_common.hpp"
namespace scf::fock_operator {

template<typename DensityType, typename ElectronType>
class RestrictedVisitor : public FockOperatorCommon<ElectronType> {
private:
    using base_type = FockOperatorCommon<ElectronType>;

public:
    using V_ee_term = simde::type::V_ee_type;

    RestrictedVisitor(simde::type::fock& F, const DensityType& rho) :
      base_type(F), m_prho_(&rho) {}

    using base_type::run;
    void run(const V_ee_term& V_ee) override {
        if(*m_prho_ == DensityType{}) return;
        using j_type = Coulomb<ElectronType, DensityType>;
        using k_type = Exchange<ElectronType, DensityType>;

        const auto& electrons = this->get_e_(V_ee);
        auto j                = std::make_unique<j_type>(electrons, *m_prho_);
        auto k                = std::make_unique<k_type>(electrons, *m_prho_);
        this->m_pF_->emplace_back(2.0, std::move(j));
        this->m_pF_->emplace_back(-1.0, std::move(k));
    }

private:
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

#define RESTRICTED(density)                                           \
    template struct Restricted<density, simde::type::many_electrons>; \
    template struct Restricted<density, simde::type::electron>

RESTRICTED(simde::type::e_density);
RESTRICTED(simde::type::decomposable_e_density);

#undef RESTRICTED

} // namespace scf::fock_operator
