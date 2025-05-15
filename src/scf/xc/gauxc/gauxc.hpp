#pragma once
#include "gauxc_property_types.hpp"
#include <simde/simde.hpp>

namespace scf::xc::gauxc {

// DECLARE_MODULE(XC);

// DECLARE_MODULE(snLinK);

// Module to convert Chemist Basis -> GauXC Basis
DECLARE_MODULE(BasisConversion);

// Module to convert Chemist AO Basis -> GauXC Molecule
DECLARE_MODULE(MoleculeConversion);

// Modules to generate XC Quadrature Batches
// DECLARE_MODULE(QuadratureBatches);

void set_defaults(pluginplay::ModuleManager& mm);
void load_modules(pluginplay::ModuleManager& mm);

} // namespace scf::xc::gauxc