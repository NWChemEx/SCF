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

#include <simde/simde.hpp>

namespace scf::fock_operator {

using namespace chemist::qm_operator;

template<typename ElectronType>
class FockOperatorCommon : public chemist::qm_operator::OperatorVisitor {
public:
    using T_e_term  = simde::type::T_e_type;
    using V_en_term = simde::type::V_en_type;
    using V_nn_term = simde::type::V_nn_type;

    explicit FockOperatorCommon(simde::type::fock& F) : m_pF_(&F) {}

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

protected:
    template<typename T>
    auto get_e_(T&& op) {
        if constexpr(std::is_same_v<ElectronType, simde::type::electron>) {
            return simde::type::electron{};
        } else {
            return op.template at<0>();
        }
    }

    simde::type::fock* m_pF_;
};
} // namespace scf::fock_operator