#include <scf/scf.hpp>
#include <simde/simde.hpp>

TEST_CASE("Restricted") {
    pluginplay::ModuleManager mm;
    scf::load_modules(mm);

    using density_type = simde::type::e_density;
    using pt           = simde::FockOperator<density_type>;

    auto& mod = mm.at("AO Restricted Fock Op");
}