#pragma once
#include "../test_scf.hpp"
#include <integrals/integrals.hpp>
#include <scf/scf.hpp>
#include <simde/simde.hpp>

namespace test_scf {

/// Factors out setting submodules for SCF plugin from other plugins
inline auto load_modules() {
    pluginplay::ModuleManager mm;
    scf::load_modules(mm);
    integrals::load_modules(mm);

    mm.change_submod("Four center J builder", "Four-center ERI", "ERI4");
    mm.change_submod("Four center K builder", "Four-center ERI", "ERI4");

    const auto ao_driver_key = "AO integral driver";
    mm.change_submod(ao_driver_key, "Kinetic", "Kinetic");
    mm.change_submod(ao_driver_key, "Electron-Nuclear attraction", "Nuclear");

    mm.change_submod("Diagonalization Fock update", "Overlap matrix builder",
                     "Overlap");

    return mm;
}

} // namespace test_scf