#pragma once
#include <simde/simde.hpp>

namespace scf::driver {

DECLARE_MODULE(SCFDriver);
DECLARE_MODULE(SCFLoop);

inline void load_modules(pluginplay::ModuleManager& mm) {
    mm.add_module<SCFDriver>("SCF Driver");
    mm.add_module<SCFLoop>("Loop");
}

inline void set_defaults(pluginplay::ModuleManager& mm) {
    mm.change_submod("Loop", "Electronic energy", "Electronic energy");
    mm.change_submod("Loop", "Density matrix", "Density matrix builder");
    mm.change_submod("Loop", "Guess update", "Diagonalization Fock update");
    mm.change_submod("Loop", "One-electron Fock operator",
                     "Restricted One-Electron Fock op");
    mm.change_submod("Loop", "Fock operator", "Restricted Fock Op");
    mm.change_submod("Loop", "Charge-charge", "Coulomb's Law");

    mm.change_submod("SCF Driver", "Guess", "Core guess");
    mm.change_submod("SCF Driver", "Optimizer", "Loop");
}

} // namespace scf::driver