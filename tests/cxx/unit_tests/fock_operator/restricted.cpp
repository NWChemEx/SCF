#include "../../test_scf.hpp"
#include <scf/scf.hpp>
#include <simde/simde.hpp>

TEST_CASE("Restricted") {
    pluginplay::ModuleManager mm;
    scf::load_modules(mm);

    using density_type = simde::type::decomposable_e_density;
    using pt           = simde::FockOperator<density_type>;

    auto& mod = mm.at("Restricted Fock Op");

    SECTION("H2 Molecule") {
        auto H   = test_scf::h2_hamiltonian();
        auto rho = test_scf::h2_density();

        auto F      = mod.run_as<pt>(H, rho);
        auto F_corr = test_scf::h2_fock();
        REQUIRE(F == F_corr);
    }
}