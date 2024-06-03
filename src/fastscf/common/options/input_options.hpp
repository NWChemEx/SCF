/*
 * ExaChem: Open Source Exascale Computational Chemistry Software.
 *
 * Copyright 2024 Pacific Northwest National Laboratory, Battelle Memorial Institute.
 *
 * See LICENSE.txt for details
 */

#pragma once

#include "common/ecatom.hpp"
#include "common/txt_utils.hpp"
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <vector>

class CommonOptions {
public:
  bool        debug{false};
  int         maxiter{50};
  std::string basis{"sto-3g"};
  std::string dfbasis{};
  std::string basisfile{}; // supports only ECPs for now
  std::string gaussian_type{"spherical"};
  std::string geom_units{"angstrom"};
  std::string file_prefix{};
  std::string ext_data_path{};
  void        print();
};

class SCFOptions: public CommonOptions {
public:
  int      charge{0};
  int      multiplicity{1};
  double   lshift{0};        // level shift factor, +ve value b/w 0 and 1
  double   tol_int{1e-22};   // tolerance for integral primitive screening
  double   tol_sch{1e-10};   // tolerance for schwarz screening
  double   tol_lindep{1e-5}; // tolerance for linear dependencies
  double   conve{1e-8};      // energy convergence
  double   convd{1e-7};      // density convergence
  int      diis_hist{10};    // number of diis history entries
  uint32_t AO_tilesize{30};
  uint32_t dfAO_tilesize{50};
  bool     restart{false}; // Read movecs from disk
  bool     noscf{false};   // only recompute energy from movecs
  bool     sad{true};
  bool     force_tilesize{false};
  bool     direct_df{false};
  bool     snK{false};
  int restart_size{2000}; // read/write orthogonalizer, schwarz, etc matrices when N>=restart_size
  int scalapack_nb{256};
  int nnodes{1};
  int scalapack_np_row{0};
  int scalapack_np_col{0};
  std::string         moldenfile{""};
  int                 n_lindep{0};
  int                 writem{1};
  int                 damp{100}; // density mixing parameter
  std::string         scf_type{"restricted"};
  std::string         xc_grid_type{"UltraFine"};
  std::string         xc_pruning_scheme{"Robust"};
  std::string         xc_rad_quad{"MK"};
  std::string         xc_weight_scheme{"SSF"};
  std::string         xc_exec_space{"Host"};
  std::string         xc_lb_kernel{"Default"};
  std::string         xc_mw_kernel{"Default"};
  std::string         xc_int_kernel{"Default"};
  std::string         xc_red_kernel{"Default"};
  std::string         xc_lwd_kernel{"Default"};
  std::pair<int, int> xc_radang_size{0, 0};
  int                 xc_batch_size{2048};
  double              xc_basis_tol{1e-8};
  double              xc_snK_etol{1e-10};
  double              xc_snK_ktol{1e-10};

  std::map<std::string, std::tuple<int, int>> guess_atom_options;

  std::vector<std::string> xc_type;
  // mos_txt: write lcao, mo transformed core H, fock, and v2 to disk as text files.
  bool                             mos_txt{false};
  bool                             mulliken_analysis{false};
  std::pair<bool, double>          mo_vectors_analysis{false, 0.15};
  std::vector<double>              qed_omegas{};
  std::vector<double>              qed_volumes{};
  std::vector<double>              qed_lambdas{};
  std::vector<std::vector<double>> qed_polvecs{};
  void                             print();
};

class ECOptions {
public:
  ECOptions() = default;

  CommonOptions common_options;
  SCFOptions    scf_options;
};
