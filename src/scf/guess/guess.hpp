#pragma once
#include <simde/simde.hpp>

namespace scf::guess {

DECLARE_MODULE(Core);

inline void load_modules(pluginplay::ModuleManager& mm) {
    mm.add_module<Core>("Core guess");
}

inline void set_defaults(pluginplay::ModuleManager& mm) {
    mm.change_submod("Core guess", "Build Fock operator",
                     "Restricted One-Electron Fock Op");
    mm.change_submod("Core guess", "Guess updater",
                     "Diagonalization Fock update");
}

} // namespace scf::guess