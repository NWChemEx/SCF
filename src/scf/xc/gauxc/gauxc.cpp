#include "gauxc.hpp"

namespace scf::xc::gauxc {

void load_modules(pluginplay::ModuleManager& mm) {
    // Conversion Utilities
    mm.add_module<BasisConversion>("GauXC Basis Converter");
    mm.add_module<MoleculeConversion>("GauXC Molecule Converter");

    // Data geneation utilities
    mm.add_module<QuadratureBatches>("GauXC Quadrature Batches");

    // XC Integration
    mm.add_module<GauXCDriver>("GauXC Driver");
    mm.add_module<XCPotential>("GauXC XC Potential");
    mm.add_module<XCEnergy>("GauXC XC Energy");

    // sn-LinK Integration
    // mm.add_module<snLinK>("snLinK");
}

void set_defaults(pluginplay::ModuleManager& mm) {
    // Data geneation utilities
    mm.change_submod("GauXC Quadrature Batches", "GauXC Basis Converter",
                     "GauXC Basis Converter");
    mm.change_submod("GauXC Quadrature Batches", "GauXC Molecule Converter",
                     "GauXC Molecule Converter");

    // XC Integration
    mm.change_submod("GauXC Driver", "Quadrature Batches",
                     "GauXC Quadrature Batches");

    mm.change_submod("GauXC XC Potential", "XC Driver", "GauXC Driver");
    mm.change_submod("GauXC XC Energy", "XC Driver", "GauXC Driver");

    // sn-LinK Integration
    // mm.change_submod("snLinK", "Quadrature Batches",
    //                 "GauXC Quadrature Batches");
}

} // namespace scf::xc::gauxc