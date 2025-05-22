/*
 * Copyright 2022 NWChemEx-Project
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "gauxc.hpp"

namespace scf::xc::gauxc {
namespace {

const auto desc = R"(
Chemist Molecule to GauXC Molecule
----------------------------------
TODO: Write me!!!
)";

}

using pt = basis_to_gauxc_molecule_conversion_t;

// Module to convert NWX AO Basis -> GauXC Molecule
MODULE_CTOR(MoleculeConversion) {
    satisfies_property_type<pt>();
    description(desc);
}

MODULE_RUN(MoleculeConversion) {
    const auto& [nwx_basis] = pt::unwrap_inputs(inputs);

    GauXC::Molecule gauxc_mol;
    gauxc_mol.reserve(nwx_basis.size());
    for(const auto& center : nwx_basis) {
        auto Z_optional = center.atomic_number();
        auto Z          = Z_optional.has_value() ? Z_optional.value() : 0;
        auto x          = center.center().coord(0);
        auto y          = center.center().coord(1);
        auto z          = center.center().coord(2);
        gauxc_mol.emplace_back(GauXC::AtomicNumber(Z), x, y, z);
    }
    auto rv = results();
    return pt::wrap_results(rv, gauxc_mol);
}

} // namespace scf::xc::gauxc
