#pragma once
#include <simde/simde.hpp>

namespace scf::matrix_builder {

// DECLARE_MODULE(Core);
// DECLARE_MODULE(Fock);
DECLARE_MODULE(JFourCenter);

inline void load_modules(pluginplay::ModuleManager& mm) {
    // mm.add_module<Core>("Core Matrix Builder");
    // mm.add_module<Fock>("Fock matrix builder");
    mm.add_module<JFourCenter>("Four center J builder");
}

} // namespace scf::matrix_builder