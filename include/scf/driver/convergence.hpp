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

#pragma once
#include <pluginplay/pluginplay.hpp>
#include <simde/types.hpp>


namespace scf {

template<typename KernelType>
DECLARE_TEMPLATED_PROPERTY_TYPE(ConvergenceProp, KernelType);

template<typename KernelType>
TEMPLATED_PROPERTY_TYPE_INPUTS(ConvergenceProp, KernelType) {
    auto rv     = pluginplay::declare_input()
                    .add_field<simde::type::tensor>("New Energy")
                    .template add_field<simde::type::tensor>("Old Energy")
                    .template add_field<chemist::DecomposableDensity<chemist::Electron>>("New Rho")
                    .template add_field<chemist::DecomposableDensity<chemist::Electron>>("Old Rho")
                    .template add_field<simde::type::tensor>("P_new")
                    .template add_field<simde::type::tensor>("S")
                    .template add_field<simde::type::fock>("Fock Operator")
                    .template add_field<chemist::wavefunction::AOs>("Wave Function")
                    .template add_field<KernelType>("K")
                    .template add_field<double>("Energy Tolerance")
                    .template add_field<double>("Density Tolerance")
                    .template add_field<double>("Gradient Tolerance");
    return rv;
}

template<typename KernelType>
TEMPLATED_PROPERTY_TYPE_RESULTS(ConvergenceProp, KernelType) {
    auto rv     = pluginplay::declare_result().add_field<bool>("Convergence Status");
    rv.at("Convergence Status")
      .set_description("The molecule corresponding to the input string");
    return rv;
}

} // namespace simde

