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

#ifdef BUILD_GAUXC
#include "gauxc/gauxc.hpp"
#endif
#include "libxc/libxc.hpp"
#include "native/native.hpp"
#include "xc.hpp"

namespace scf::xc {
namespace {

/// The angular size of "Reference IntegratorXX Grid": the largest
/// Lebedev-Laikov grid IntegratorXX has (algebraic order 131)
constexpr std::size_t reference_angular_size = 5810;

/// The radial size of "Reference IntegratorXX Grid". For SVWN3 water/STO-3G
/// with 5810 angular points, E_xc and V_xc change by less than 1.1E-11 from
/// 150 to 250 radial points (see xc_uq.rst)
constexpr std::size_t reference_radial_size = 150;

} // namespace

void load_modules(pluginplay::ModuleMap& modules) {
#ifdef BUILD_GAUXC
    gauxc::load_modules(modules);
#endif
    libxc::load_modules(modules);
    native::load_modules(modules);
    modules.add_module<AOsOnGrid>("AOs on a grid");
    modules.add_module<Gau2Grid>("Gau2Grid");
    modules.add_module<GridFromFile>("Grid From File");
    modules.add_module<GridFromIntegratorXX>("Grid From IntegratorXX");
    modules.add_module<Density2Grid>("Density2Grid");
    modules.add_module<UQXCEnergy>("UQ XC Energy");
    modules.add_module<UQXCPotential>("UQ XC Potential");

    using scf::property_types::xc_energy_quadrature;
    using scf::property_types::xc_potential_quadrature;
    modules.add_module<ReferenceGridError<xc_energy_quadrature>>(
      "Reference Grid XC Energy Error");
    modules.add_module<ReferenceGridError<xc_potential_quadrature>>(
      "Reference Grid XC Potential Error");

    // The reference grid of the reference grid error models
    const auto ref_grid_key = "Reference IntegratorXX Grid";
    modules.copy_module("Grid From IntegratorXX", ref_grid_key);
    modules.change_input(ref_grid_key, "Angular Size",
                         std::size_t{reference_angular_size});
    modules.change_input(ref_grid_key, "Radial Size",
                         std::size_t{reference_radial_size});
}

void set_defaults(pluginplay::ModuleMap& modules) {
#ifdef BUILD_GAUXC
    gauxc::set_defaults(modules);
#endif
    libxc::set_defaults(modules);
    native::set_defaults(modules);
    modules.change_submod("Density2Grid", "AOs on a grid", "AOs on a grid");

    // The UQ XC modules default to the native SVWN3 quadratures, since they
    // support the UQ types, and to the reference grid error models. The
    // "UQ Initializer" of the error models comes from Integrals (see
    // scf_mm.cpp).
    const auto grid_key  = "Grid From IntegratorXX";
    const auto e_quad    = "SVWN3 Energy Quadrature";
    const auto v_quad    = "SVWN3 Potential Quadrature";
    const auto e_error   = "Reference Grid XC Energy Error";
    const auto v_error   = "Reference Grid XC Potential Error";
    const auto ref_grid  = "Reference IntegratorXX Grid";
    const auto error_key = "Incompleteness Error Model";
    modules.change_submod("UQ XC Energy", "Integration grid", grid_key);
    modules.change_submod("UQ XC Potential", "Integration grid", grid_key);
    modules.change_submod("UQ XC Energy", "XC Energy", e_quad);
    modules.change_submod("UQ XC Potential", "XC Potential", v_quad);
    modules.change_submod("UQ XC Energy", error_key, e_error);
    modules.change_submod("UQ XC Potential", error_key, v_error);
    modules.change_submod(e_error, "Reference grid", ref_grid);
    modules.change_submod(v_error, "Reference grid", ref_grid);
    modules.change_submod(e_error, "Reference XC Module", e_quad);
    modules.change_submod(v_error, "Reference XC Module", v_quad);
}

} // namespace scf::xc
