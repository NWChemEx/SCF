#pragma once
#include <simde/simde.hpp>

namespace scf::update {

DECLARE_MODULE(Diagonalization);

inline void load_modules(pluginplay::ModuleManager& mm) {
    mm.add_module<Diagonalization>("Diagonalization Fock update.");
}

inline void set_defaults(pluginplay::ModuleManager& mm) {
    mm.change_submod("Diagonalization Fock update", "Diagonalizer",
                     "Generalized eigensolve via Eigen");
}

} // namespace scf::update