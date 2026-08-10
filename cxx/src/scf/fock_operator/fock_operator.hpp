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

template<typename DensityType, typename ElectronType>
DECLARE_MODULE(Restricted);

template<typename DensityType, typename ElectronType>
DECLARE_MODULE(RKohnSham);

inline void load_modules(pluginplay::ModuleManager& mm) {
    using simde::type::decomposable_e_density;
    using simde::type::e_density;
    using simde::type::electron;
    using simde::type::many_electrons;
    using e_many    = Restricted<e_density, many_electrons>;
    using e_e       = Restricted<e_density, electron>;
    using d_many    = Restricted<decomposable_e_density, many_electrons>;
    using d_e       = Restricted<decomposable_e_density, electron>;
    using e_ks_many = RKohnSham<e_density, many_electrons>;
    using e_ks_e    = RKohnSham<e_density, electron>;
    using d_ks_many = RKohnSham<decomposable_e_density, many_electrons>;
    using d_ks_e    = RKohnSham<decomposable_e_density, electron>;
    mm.add_module<e_many>("AO Restricted Fock Op");
    mm.add_module<d_many>("Restricted Fock Op");
    mm.add_module<e_e>("AO Restricted One-Electron Fock Op");
    mm.add_module<d_e>("Restricted One-Electron Fock Op");
    mm.add_module<e_ks_many>("AO Restricted Kohn-Sham Op");
    mm.add_module<d_ks_many>("Restricted Kohn-Sham Op");
    mm.add_module<e_ks_e>("AO Restricted One-Electron Kohn-Sham Op");
    mm.add_module<d_ks_e>("Restricted One-Electron Kohn-Sham Op");
}

#define EXTERN_RESTRICTED(density)                                           \
    extern template struct Restricted<density, simde::type::many_electrons>; \
    extern template struct Restricted<density, simde::type::electron>

#define EXTERN_RKOHN_SHAM(density)                                          \
    extern template struct RKohnSham<density, simde::type::many_electrons>; \
    extern template struct RKohnSham<density, simde::type::electron>

EXTERN_RESTRICTED(simde::type::e_density);
EXTERN_RESTRICTED(simde::type::decomposable_e_density);
EXTERN_RKOHN_SHAM(simde::type::e_density);
EXTERN_RKOHN_SHAM(simde::type::decomposable_e_density);

#undef EXTERN_RKOHN_SHAM
#undef EXTERN_RESTRICTED
} // namespace scf::fock_operator
