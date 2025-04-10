#pragma once
#include <simde/simde.hpp>

namespace scf::driver {

simde::type::tensor commutator(simde::type::tensor Fock_Matrix,
                               simde::type::tensor Density_Matrix,
                               simde::type::tensor Overlap_Matrix);
} // namespace scf::driver
