#pragma once
#include <simde/simde.hpp>

namespace scf::matrix_builder {

DECLARE_MODULE(AOIntegralsDriver);
DECLARE_MODULE(Fock);
DECLARE_MODULE(JFourCenter);
DECLARE_MODULE(KFourCenter);

inline void load_modules(pluginplay::ModuleManager& mm) {
    mm.add_module<AOIntegralsDriver>("AO integral driver");
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

    mm.change_submod("Fock matrix builder", "Two center evaluator", ao_driver);
}

} // namespace scf::matrix_builder