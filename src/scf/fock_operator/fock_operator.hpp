#pragma once
#include <simde/simde.hpp>

namespace scf::fock_operator {

template<typename DensityType>
DECLARE_MODULE(Restricted);

inline void load_modules(pluginplay::ModuleManager& mm) {
    using simde::type::decomposable_e_density;
    using simde::type::e_density;
    mm.add_module<Restricted<e_density>>("AO Restricted Fock Op");
    mm.add_module<Restricted<decomposable_e_density>>("Restricted Fock Op");
}

extern template class Restricted<simde::type::e_density>;
extern template class Restricted<simde::type::decomposable_e_density>;

} // namespace scf::fock_operator