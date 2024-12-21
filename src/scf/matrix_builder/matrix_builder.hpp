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

namespace scf::matrix_builder {

DECLARE_MODULE(AOIntegralsDriver);
DECLARE_MODULE(DensityMatrix);
DECLARE_MODULE(DeterminantDriver);
DECLARE_MODULE(ElectronicEnergy);
DECLARE_MODULE(Fock);
DECLARE_MODULE(JFourCenter);
DECLARE_MODULE(KFourCenter);

inline void load_modules(pluginplay::ModuleManager& mm) {
    mm.add_module<AOIntegralsDriver>("AO integral driver");
    mm.add_module<DensityMatrix>("Density matrix builder");
    mm.add_module<DeterminantDriver>("Determinant driver");
    mm.add_module<ElectronicEnergy>("Electronic energy");
    mm.add_module<Fock>("Fock matrix builder");
    mm.add_module<JFourCenter>("Four center J builder");
    mm.add_module<KFourCenter>("Four center K builder");
}

inline void set_defaults(pluginplay::ModuleManager& mm) {
    const auto ao_driver = "AO integral driver";
    mm.change_submod(ao_driver, "Coulomb matrix", "Four center J builder");
    mm.change_submod(ao_driver, "Exchange matrix", "Four center K builder");
    // TODO: Re-enable when PluginPlay doesn't choke on loops in modules
    // mm.change_submod(ao_driver, "Fock matrix", "Fock Matrix Builder");
    mm.change_submod(ao_driver, "Density matrix", "Density matrix builder");

    mm.change_submod("Fock matrix builder", "Two center evaluator", ao_driver);

    const auto det_driver = "Determinant driver";
    mm.change_submod(det_driver, "Two center evaluator", ao_driver);
    mm.change_submod(det_driver, "Fock matrix", "Fock matrix builder");

    mm.change_submod("Electronic energy", "determinant driver", det_driver);
}

} // namespace scf::matrix_builder