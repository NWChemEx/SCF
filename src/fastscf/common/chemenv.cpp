/*
 * ExaChem: Open Source Exascale Computational Chemistry Software.
 *
 * Copyright 2024 Pacific Northwest National Laboratory, Battelle Memorial Institute.
 *
 * See LICENSE.txt for details
 */

#include "chemenv.hpp"

void ChemEnv::write_sinfo() {
  std::string basis = ioptions.scf_options.basis;

  std::string out_fp       = workspace_dir;
  std::string files_dir    = out_fp + ioptions.scf_options.scf_type;
  std::string files_prefix = /*out_fp;*/ files_dir + "/" + sys_data.output_file_prefix;
  if(!fs::exists(files_dir)) fs::create_directories(files_dir);

  json results;

  results["molecule"]["name"]         = sys_data.input_molecule;
  results["molecule"]["basis"]["all"] = basis;

  for(size_t i = 0; i < atoms.size(); i++) {
    ECAtom& iatom = ec_atoms[i];
    if(iatom.basis != basis) results["molecule"]["basis"][iatom.esymbol] = iatom.basis;
  }

  results["molecule"]["nbf"]              = shells.nbf();
  results["molecule"]["nshells"]          = shells.size();
  results["molecule"]["nelectrons"]       = sys_data.nelectrons;
  results["molecule"]["nelectrons_alpha"] = sys_data.nelectrons_alpha;
  results["molecule"]["nelectrons_beta"]  = sys_data.nelectrons_beta;

  std::string json_file   = files_prefix + ".sinfo.json";
  bool        json_exists = std::filesystem::exists(json_file);
  if(json_exists) std::filesystem::remove(json_file);

  std::ofstream res_file(json_file);
  res_file << std::setw(2) << results << std::endl;
}

void ChemEnv::write_json_data(const std::string cmodule) {
  auto  scf     = ioptions.scf_options;
  json& results = sys_data.results;

  auto str_bool = [=](const bool val) {
    if(val) return "true";
    return "false";
  };

  results["input"]["molecule"]["name"]     = sys_data.input_molecule;
  results["input"]["molecule"]["basisset"] = scf.basis;
  // results["input"]["molecule"]["gaussian_type"]  = scf.gaussian_type;
  results["input"]["molecule"]["geometry_units"] = scf.geom_units;
  // SCF options
  results["input"]["SCF"]["tol_int"]        = scf.tol_int;
  results["input"]["SCF"]["tol_sch"]        = scf.tol_sch;
  results["input"]["SCF"]["tol_lindep"]     = scf.tol_lindep;
  results["input"]["SCF"]["conve"]          = scf.conve;
  results["input"]["SCF"]["convd"]          = scf.convd;
  results["input"]["SCF"]["diis_hist"]      = scf.diis_hist;
  results["input"]["SCF"]["AO_tilesize"]    = scf.AO_tilesize;
  results["input"]["SCF"]["force_tilesize"] = str_bool(scf.force_tilesize);
  results["input"]["SCF"]["scf_type"]       = scf.scf_type;
  results["input"]["SCF"]["multiplicity"]   = scf.multiplicity;

  std::string l_module = cmodule;
  txt_utils::to_lower(l_module);

  std::string files_dir =
    sys_data.output_file_prefix + "_files/" + ioptions.scf_options.scf_type + "/json";
  if(!fs::exists(files_dir)) fs::create_directories(files_dir);
  std::string files_prefix = files_dir + "/" + sys_data.output_file_prefix;
  std::string json_file    = files_prefix + "." + l_module + ".json";
  bool        json_exists  = std::filesystem::exists(json_file);
  if(json_exists) {
    // std::ifstream jread(json_file);
    // jread >> results;
    std::filesystem::remove(json_file);
  }

  // std::cout << std::endl << std::endl << results.dump() << std::endl;
  std::ofstream res_file(json_file);
  res_file << std::setw(2) << results << std::endl;
}

void ChemEnv::sinfo() {
  SCFOptions& scf_options = ioptions.scf_options;

  // TODO: dont create ec
  ProcGroup        pg = ProcGroup::create_world_coll();
  ExecutionContext ec{pg, DistributionKind::nw, MemoryManagerKind::ga};
  auto             rank   = ec.pg().rank();
  std::string      basis  = scf_options.basis;
  int              charge = scf_options.charge;

  libint2::initialize(false);

  std::string basis_set_file = std::string(DATADIR) + "/basis/" + basis + ".g94";

  int basis_file_exists = 0;
  if(rank == 0) basis_file_exists = std::filesystem::exists(basis_set_file);
  ec.pg().broadcast(&basis_file_exists, 0);

  if(!basis_file_exists)
    tamm_terminate("ERROR: basis set file " + basis_set_file + " does not exist");

  libint2::BasisSet shells;
  {
    std::vector<libint2::Shell> bset_vec;
    for(size_t i = 0; i < atoms.size(); i++) {
      // const auto        Z = atoms[i].atomic_number;
      libint2::BasisSet ashells(ec_atoms[i].basis, {atoms[i]});
      bset_vec.insert(bset_vec.end(), ashells.begin(), ashells.end());
    }
    libint2::BasisSet bset(bset_vec);
    shells = std::move(bset);
  }

  shells.set_pure(true);

  const int N = shells.nbf();

  auto nelectrons = 0;
  for(size_t i = 0; i < atoms.size(); ++i) nelectrons += atoms[i].atomic_number;
  nelectrons -= charge;

  // sys_data.nelectrons = nelectrons;
  EXPECTS((nelectrons + scf_options.multiplicity - 1) % 2 == 0);

  int nelectrons_alpha = (nelectrons + scf_options.multiplicity - 1) / 2;
  int nelectrons_beta  = nelectrons - nelectrons_alpha;

  if(rank == 0) {
    std::cout << std::endl << "Number of basis functions = " << N << std::endl;
    std::cout << std::endl << "Total number of shells = " << shells.size() << std::endl;
    std::cout << std::endl << "Total number of electrons = " << nelectrons << std::endl;
    std::cout << "  # of alpha electrons    = " << nelectrons_alpha << std::endl;
    std::cout << "  # of beta electons      = " << nelectrons_beta << std::endl << std::flush;
  }

  if(rank == 0) write_sinfo();

  ec.flush_and_sync();
}