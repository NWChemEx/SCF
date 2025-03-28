/*
 * Copyright 2024 NWChemEx-Project
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

#ifdef BUILD_TAMM_SCF
#include "exachem/common/chemenv.hpp"
#include "exachem/common/initialize_system_data.hpp"
#include "exachem/scf/scf_main.hpp"
#include "scf_modules.hpp"
#include <libint2.hpp>
#include <simde/simde.hpp>

namespace scf {

using energy_pt = simde::AOEnergy;

inline libint2::BasisSet
make_libint_basis(const simde::type::ao_basis_set &bs) {
  /// Typedefs for everything
  using atom_t = libint2::Atom;
  using shell_t = libint2::Shell;
  using basis_t = libint2::BasisSet;
  using cont_t = libint2::Shell::Contraction;
  using svec_d_t = libint2::svector<double>;
  using conts_t = libint2::svector<cont_t>;
  using centers_t = std::vector<atom_t>;
  using atom_bases_t = std::vector<shell_t>;
  using element_bases_t = std::vector<atom_bases_t>;

  /// Inputs for BasisSet constructor
  centers_t centers{};
  element_bases_t element_bases{};

  /// Atom doesn't have a value ctor, so here's a stand in
  auto atom_ctor = [](int Z, double x, double y, double z) {
    atom_t atom{};
    atom.atomic_number = Z;
    atom.x = x;
    atom.y = y;
    atom.z = z;
    return atom;
  };

  /// Origin for shell construction
  std::array<double, 3> origin = {0.0, 0.0, 0.0};

  /// Convert centers and their shells to libint equivalents.
  for (decltype(bs.size()) abs_i = 0; abs_i < bs.size(); ++abs_i) {
    /// Add current center to atoms list
    const auto &abs = bs[abs_i];
    centers.push_back(
        atom_ctor(abs_i, abs.center().x(), abs.center().y(), abs.center().z()));

    /// Gather shells for this center and add them to element_bases
    atom_bases_t atom_bases{};
    for (const auto &&shelli : abs) {
      const auto nprims = shelli.n_primitives();
      const auto prim0 = shelli.primitive(0);
      const auto primN = shelli.primitive(nprims - 1);
      const bool pure = shelli.pure() == chemist::ShellType::pure;
      const int l = shelli.l();

      svec_d_t alphas(&prim0.exponent(), &primN.exponent() + 1);
      svec_d_t coefs(&prim0.coefficient(), &primN.coefficient() + 1);
      conts_t conts{cont_t{l, pure, coefs}};
      /// Use origin for position, because BasisSet moves shells to center
      atom_bases.push_back(shell_t(alphas, conts, origin));
    }
    element_bases.push_back(atom_bases);
  }

  /// Return the new basis set
  return basis_t(centers, element_bases);
}

MODULE_CTOR(TAMMEnergy) {
  satisfies_property_type<energy_pt>();

  add_input<std::string>("molecule_name")
      .set_description("The name of the molecule");
  add_input<std::string>("units")
      .set_default("angstrom")
      .set_description("Specifies the units as bohr or angstrom");

  add_input<int>("charge").set_default(0).set_description("Charge");
  add_input<int>("multiplicity")
      .set_default(1)
      .set_description(
          "number of singly occupied orbitals for a particular calculation");
  add_input<double>("lshift").set_default(0.0).set_description(
      "level shift factor denoting the amount of shift applied to the diagonal "
      "elements of the unoccupied block of the Fock matrix");

  add_input<double>("tol_int").set_default(1e-22).set_description(
      "integral primitive screening threshold");

  add_input<double>("tol_sch").set_default(1.0e-10).set_description(
      "The Schwarz inequality is used to screen the product of integrals and "
      "density matrices");

  add_input<double>("tol_lindep")
      .set_default(1e-5)
      .set_description(
          "Tolerance for detecting the linear dependence of basis set");

  add_input<double>("conve").set_default(1e-8).set_description(
      "Specifies the energy convergence threshold");
  add_input<double>("convd").set_default(1e-6).set_description(
      "Specifies the density convergence threshold");

  add_input<int>("diis_hist")
      .set_default(10)
      .set_description("number of DIIS history entries to store for the fock "
                       "and error matrices");

  add_input<int>("damp").set_default(100).set_description(
      "percentage of the current iterations density mixed with the previous "
      "iterations density. 100\% indicates no damping");

  add_input<int>("writem").set_default(1).set_description(
      "An integer specifying the frequency (as number of iterations) after "
      "which the movecs and density matrices are written to disk for "
      "restarting the calculation");

  add_input<bool>("debug").set_default(false).set_description(
      "enable verbose printing for debugging a calculation");

  add_input<bool>("restart").set_default(false).set_description(
      "indicates the calculation be restarted");
  add_input<bool>("noscf").set_default(false).set_description(
      "computes only the SCF energy upon restart");

  add_input<std::string>("scf_type")
      .set_default("restricted")
      .set_description("options supported are restricted and unrestricted");
  add_input<bool>("direct_df")
      .set_default(false)
      .set_description("Requests the direct computation of the density-fitted "
                       "Coulomb contribution. Works only for pure Kohn-Sham "
                       "fnctionals (no exact exchange) ");

  // DFT
  add_input<bool>("snK").set_default(false).set_description(
      "Computes the exact exchange contribution using the seminumerical "
      "approach implemented in GauXC");
  add_input<std::vector<std::string>>("xc_type")
      .set_default(std::vector<std::string>{})
      .set_description("A list of strings specifying the exchange and "
                       "correlation functionals for DFT calculations");

  add_input<std::string>("xc_grid_type")
      .set_default("UltraFine")
      .set_description(
          "Specifies the quality of the numerical integration grid");

  add_input<std::string>("xc_pruning_scheme")
      .set_default("Robust")
      .set_description("GauXC pruning scheme. Options supported are Robust, "
                       "Treutler, Unpruned");
  add_input<std::string>("xc_rad_quad")
      .set_default("MK")
      .set_description(
          "Specifies the GauXC radial quadrature. Options are MK,TA,MHL");
  add_input<std::string>("xc_weight_scheme")
      .set_default("SSF")
      .set_description(
          "Specifies the GauXC partitioning scheme. Can be SSF, Becke, LKO");
  add_input<std::string>("xc_exec_space")
      .set_default("Host")
      .set_description("Specifies the GauXC execution space (Host or Device) "
                       "for the load balancer and integrator");

  add_input<double>("xc_basis_tol")
      .set_default(1e-10)
      .set_description("Specifies the GauXC basis tolerance");
  add_input<int>("xc_batch_size")
      .set_default(2048)
      .set_description("Specifies the GauXC batch size");

  add_input<double>("xc_snK_etol")
      .set_default(1e-10)
      .set_description(
          "snK energy tolerance. If conve < xc_snK_etol, this tolerance will "
          "be automatically set to the conve value");
  add_input<double>("xc_snK_ktol")
      .set_default(1e-10)
      .set_description(
          "K matrix tolerance. If conve * 1e-2 < xc_snK_ktol, this tolerance "
          "will be automatically set to conve * 1e-2");

  add_input<std::string>("xc_lb_kernel")
      .set_default("Default")
      .set_description("Specifies the GauXC Load Balancer Kernel");
  add_input<std::string>("xc_mw_kernel")
      .set_default("Default")
      .set_description("Specifies the GauXC Molecular Weights Kernel");
  add_input<std::string>("xc_int_kernel")
      .set_default("Default")
      .set_description("Specifies the GauXC Integrator Kernel");
  add_input<std::string>("xc_red_kernel")
      .set_default("Default")
      .set_description("Specifies the GauXC Reduction Kernel");
  add_input<std::string>("xc_lwd_kernel")
      .set_default("Default")
      .set_description("Specifies the GauXC Local Work Driver Kernel");
}

MODULE_RUN(TAMMEnergy) {

  const auto rank = ProcGroup::world_rank();
  ProcGroup pg = ProcGroup::create_world_coll();
  ExecutionContext ec{pg, DistributionKind::nw, MemoryManagerKind::ga};

  const auto &[aos, cs] = energy_pt::unwrap_inputs(inputs);

  const double angstrom_to_bohr = 1.8897259878858;

  libint2::BasisSet li_shells = make_libint_basis(aos);

  ChemEnv chem_env;
  chem_env.input_file = inputs.at("molecule_name").value<std::string>();
  auto mol = cs.molecule();
  chem_env.ec_atoms.resize(mol.size());
  chem_env.atoms.resize(mol.size());

  CommonOptions &coptions = chem_env.ioptions.common_options;
  // parse common options
  coptions.geom_units = inputs.at("units").value<std::string>();
  const double convert_units =
      (coptions.geom_units == "angstrom") ? angstrom_to_bohr : 1.0;
  coptions.basis = aos[0].basis_set_name().value_or("sto-3g");

  std::cout << std::endl;
  for (decltype(mol.size()) i = 0; i < mol.size(); i++) {
    auto atom_i = mol[i];
    chem_env.atoms[i] = {(int)atom_i.Z(), atom_i.x() * convert_units,
                         atom_i.y() * convert_units,
                         atom_i.z() * convert_units};
    chem_env.ec_atoms[i].atom = chem_env.atoms[i];
    chem_env.ec_atoms[i].esymbol = atom_i.name();
    chem_env.ec_atoms[i].basis = coptions.basis;

    std::cout << std::setw(3) << std::left << chem_env.ec_atoms[i].esymbol
              << " " << std::right << std::setw(14) << std::fixed
              << std::setprecision(10) << chem_env.atoms[i].x << " "
              << std::right << std::setw(14) << std::fixed
              << std::setprecision(10) << chem_env.atoms[i].y << " "
              << std::right << std::setw(14) << std::fixed
              << std::setprecision(10) << chem_env.atoms[i].z << "\n";
  }

  chem_env.sys_data.input_molecule = chem_env.input_file;

  if (chem_env.ioptions.common_options.file_prefix.empty()) {
    chem_env.ioptions.common_options.file_prefix =
        chem_env.sys_data.input_molecule;
  }

  chem_env.sys_data.output_file_prefix =
      chem_env.ioptions.common_options.file_prefix + "." +
      chem_env.ioptions.common_options.basis;
  chem_env.workspace_dir = chem_env.sys_data.output_file_prefix + "_files/";

  // Set SCF options
  SCFOptions &scf = chem_env.ioptions.scf_options;
  scf.charge = inputs.at("charge").value<int>();
  scf.multiplicity = inputs.at("multiplicity").value<int>();
  scf.lshift = inputs.at("lshift").value<double>();
  scf.tol_int = inputs.at("tol_int").value<double>();
  scf.tol_sch = inputs.at("tol_sch").value<double>();
  scf.tol_lindep = inputs.at("tol_lindep").value<double>();
  scf.conve = inputs.at("conve").value<double>();
  scf.conve = inputs.at("convd").value<double>();
  scf.diis_hist = inputs.at("diis_hist").value<int>();
  scf.damp = inputs.at("damp").value<int>();
  scf.writem = inputs.at("writem").value<int>();
  scf.debug = inputs.at("debug").value<bool>();
  scf.restart = inputs.at("restart").value<bool>();
  scf.noscf = inputs.at("noscf").value<bool>();
  scf.scf_type = inputs.at("scf_type").value<std::string>();
  scf.direct_df = inputs.at("direct_df").value<bool>();

  // DFT
  scf.snK = inputs.at("snK").value<bool>();
  scf.xc_type = inputs.at("xc_type").value<std::vector<std::string>>();

  scf.xc_grid_type = inputs.at("xc_grid_type").value<std::string>();
  scf.xc_pruning_scheme = inputs.at("xc_pruning_scheme").value<std::string>();
  scf.xc_rad_quad = inputs.at("xc_rad_quad").value<std::string>();
  scf.xc_weight_scheme = inputs.at("xc_weight_scheme").value<std::string>();
  scf.xc_exec_space = inputs.at("xc_exec_space").value<std::string>();

  scf.xc_basis_tol = inputs.at("xc_basis_tol").value<double>();
  scf.xc_batch_size = inputs.at("xc_batch_size").value<int>();
  scf.xc_snK_etol = inputs.at("xc_snK_etol").value<double>();
  scf.xc_snK_ktol = inputs.at("xc_snK_ktol").value<double>();

  scf.xc_lb_kernel = inputs.at("xc_lb_kernel").value<std::string>();
  scf.xc_mw_kernel = inputs.at("xc_mw_kernel").value<std::string>();
  scf.xc_int_kernel = inputs.at("xc_int_kernel").value<std::string>();
  scf.xc_red_kernel = inputs.at("xc_red_kernel").value<std::string>();
  scf.xc_lwd_kernel = inputs.at("xc_lwd_kernel").value<std::string>();

  IniSystemData ini_sys_data(chem_env);
  SCFOptions &scf_options = chem_env.ioptions.scf_options;
  chem_env.ec_basis =
      ECBasis(ec, scf_options.basis, scf_options.basisfile,
              scf_options.gaussian_type, chem_env.atoms, chem_env.ec_atoms);
  chem_env.shells = chem_env.ec_basis.shells;
  chem_env.sys_data.has_ecp = chem_env.ec_basis.has_ecp;

  exachem::scf::scf(ec, chem_env);

  // This is a total energy in Hartree
  double E0 = chem_env.scf_context.hf_energy;
  auto rv = results();
  return energy_pt::wrap_results(rv, E0);
}

} // namespace scf
#endif