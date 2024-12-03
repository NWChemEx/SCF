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

#pragma once
#include <simde/simde.hpp>

namespace scf::fock_operator {

template<typename DensityType>
DECLARE_MODULE(Restricted);

inline void load_modules(pluginplay::ModuleManager& mm) {
    using simde::type::decomposable_e_density;
    using simde::type::e_density;
    mm.add_module<Restricted<e_density>>("AO Restricted Fock Op");
    mm.add_module<Restricted<decomposable_e_density>>("Restricted Fock Op");
}

extern template class Restricted<simde::type::e_density>;
extern template class Restricted<simde::type::decomposable_e_density>;

} // namespace scf::fock_operator