/*
 * ExaChem: Open Source Exascale Computational Chemistry Software.
 *
 * Copyright 2024 Pacific Northwest National Laboratory, Battelle Memorial Institute.
 *
 * See LICENSE.txt for details
 */

#include "scf_iter.hpp"

template<typename TensorType>
std::tuple<TensorType, TensorType> exachem::scf::SCFIter::scf_iter_body(
  ExecutionContext& ec, ChemEnv& chem_env, ScalapackInfo& scalapack_info, const int& iter,
  SCFVars& scf_vars, TAMMTensors& ttensors, EigenTensors& etensors
#if defined(USE_GAUXC)
  ,
  GauXC::XCIntegrator<Matrix>& gauxc_integrator
#endif
) {

  SystemData& sys_data    = chem_env.sys_data;
  SCFOptions& scf_options = chem_env.ioptions.scf_options;

  const bool   is_uhf = sys_data.is_unrestricted;
  const bool   is_rhf = sys_data.is_restricted;
  const bool   do_snK = sys_data.do_snK;
  const double lshift = scf_vars.lshift;

  Tensor<TensorType>& H1       = ttensors.H1;
  Tensor<TensorType>& S1       = ttensors.S1;
  Tensor<TensorType>& D_diff   = ttensors.D_diff;
  Tensor<TensorType>& ehf_tmp  = ttensors.ehf_tmp;
  Tensor<TensorType>& ehf_tamm = ttensors.ehf_tamm;

  Matrix&             D_alpha           = etensors.D_alpha;
  Tensor<TensorType>& F_alpha           = ttensors.F_alpha;
  Tensor<TensorType>& FD_alpha_tamm     = ttensors.FD_alpha;
  Tensor<TensorType>& FDS_alpha_tamm    = ttensors.FDS_alpha;
  Tensor<TensorType>& D_alpha_tamm      = ttensors.D_alpha;
  Tensor<TensorType>& D_last_alpha_tamm = ttensors.D_last_alpha;

  Matrix&             D_beta           = etensors.D_beta;
  Tensor<TensorType>& F_beta           = ttensors.F_beta;
  Tensor<TensorType>& FD_beta_tamm     = ttensors.FD_beta;
  Tensor<TensorType>& FDS_beta_tamm    = ttensors.FDS_beta;
  Tensor<TensorType>& D_beta_tamm      = ttensors.D_beta;
  Tensor<TensorType>& D_last_beta_tamm = ttensors.D_last_beta;

  Scheduler sch{ec};

  const int64_t          N   = sys_data.nbf_orig;
  const TiledIndexSpace& tAO = scf_vars.tAO;

  auto rank          = ec.pg().rank();
  auto debug         = scf_options.debug;
  auto [mu, nu, ku]  = tAO.labels<3>("all");
  const int max_hist = scf_options.diis_hist;

  double ehf = 0.0;

#if defined(USE_GAUXC)
  if(do_snK) {
    const auto snK_start = std::chrono::high_resolution_clock::now();
    scf::gauxc::compute_exx<TensorType>(ec, chem_env, scf_vars, ttensors, etensors,
                                        gauxc_integrator);
    const auto snK_stop = std::chrono::high_resolution_clock::now();
    const auto snK_time =
      std::chrono::duration_cast<std::chrono::duration<double>>((snK_stop - snK_start)).count();
    if(rank == 0 && debug)
      std::cout << std::fixed << std::setprecision(2) << "snK: " << snK_time << "s, ";
  }
#endif

  if(is_rhf) {
    // clang-format off
    sch
      (ehf_tmp(mu,nu)  = H1(mu,nu))
      (ehf_tmp(mu,nu) += F_alpha(mu,nu))
      (ehf_tamm()      = 0.5 * D_last_alpha_tamm() * ehf_tmp())
      .execute();
    // clang-format on
  }

  if(is_uhf) {
    // clang-format off
    sch
      (ehf_tmp(mu,nu)  = H1(mu,nu))
      (ehf_tmp(mu,nu) += F_alpha(mu,nu))
      (ehf_tamm()      = 0.5 * D_last_alpha_tamm() * ehf_tmp())
      (ehf_tmp(mu,nu)  = H1(mu,nu))
      (ehf_tmp(mu,nu) += F_beta(mu,nu))
      (ehf_tamm()     += 0.5 * D_last_beta_tamm()  * ehf_tmp())
      .execute();
    // clang-format on
  }

  ehf = get_scalar(ehf_tamm);

#if defined(USE_GAUXC)
  const bool is_ks     = sys_data.is_ks;
  double     gauxc_exc = 0;
  if(is_ks) {
    const auto xcf_start = std::chrono::high_resolution_clock::now();
    gauxc_exc =
      scf::gauxc::compute_xcf<TensorType>(ec, chem_env, ttensors, etensors, gauxc_integrator);

    const auto xcf_stop = std::chrono::high_resolution_clock::now();
    const auto xcf_time =
      std::chrono::duration_cast<std::chrono::duration<double>>((xcf_stop - xcf_start)).count();
    if(rank == 0 && debug)
      std::cout << std::fixed << std::setprecision(2) << "xcf: " << xcf_time << "s, ";
  }

  ehf += gauxc_exc;
  scf_vars.exc = gauxc_exc;

  if(is_ks) {
    sch(F_alpha() += ttensors.VXC_alpha());
    if(is_uhf) {
      // clang-format off
      sch
        (F_alpha() += ttensors.VXC_beta())
        (F_beta()  += ttensors.VXC_alpha())
        (F_beta()  += -1.0 * ttensors.VXC_beta());
      // clang-format on
    }
    sch.execute();
  }
#endif

  Tensor<TensorType> err_mat_alpha_tamm{tAO, tAO};
  Tensor<TensorType> err_mat_beta_tamm{tAO, tAO};
  Tensor<TensorType>::allocate(&ec, err_mat_alpha_tamm);
  if(is_uhf) Tensor<TensorType>::allocate(&ec, err_mat_beta_tamm);

  // clang-format off
  sch
    (FD_alpha_tamm(mu,nu)       = F_alpha(mu,ku)      * D_last_alpha_tamm(ku,nu))
    (FDS_alpha_tamm(mu,nu)      = FD_alpha_tamm(mu,ku) * S1(ku,nu))
    (err_mat_alpha_tamm(mu,nu)  = FDS_alpha_tamm(mu,nu))
    (err_mat_alpha_tamm(mu,nu) -= FDS_alpha_tamm(nu,mu))
    .execute();
  // clang-format on

  if(is_uhf) {
    // clang-format off
    sch
      (FD_beta_tamm(mu,nu)        = F_beta(mu,ku)       * D_last_beta_tamm(ku,nu))
      (FDS_beta_tamm(mu,nu)       = FD_beta_tamm(mu,ku)  * S1(ku,nu))
      (err_mat_beta_tamm(mu,nu)   = FDS_beta_tamm(mu,nu))
      (err_mat_beta_tamm(mu,nu)  -= FDS_beta_tamm(nu,mu))
      .execute();
    // clang-format on
  }

  if(iter >= 1) {
    auto do_t1 = std::chrono::high_resolution_clock::now();
    if(is_rhf) {
      ++scf_vars.idiis;
      scf_diis(ec, chem_env, tAO, F_alpha, F_alpha, err_mat_alpha_tamm, err_mat_alpha_tamm, iter,
               max_hist, scf_vars, sys_data.n_lindep, ttensors.diis_hist, ttensors.diis_hist,
               ttensors.fock_hist, ttensors.fock_hist);
    }
    if(is_uhf) {
      ++scf_vars.idiis;
      scf_diis(ec, chem_env, tAO, F_alpha, F_beta, err_mat_alpha_tamm, err_mat_beta_tamm, iter,
               max_hist, scf_vars, sys_data.n_lindep, ttensors.diis_hist, ttensors.diis_beta_hist,
               ttensors.fock_hist, ttensors.fock_beta_hist);
      // scf_diis(ec, chem_env, tAO, D_beta_tamm, F_beta, err_mat_beta_tamm, iter, max_hist,
      // scf_vars,
      //         sys_data.n_lindep, ttensors.diis_beta_hist, ttensors.fock_beta_hist);
    }
    auto do_t2 = std::chrono::high_resolution_clock::now();
    auto do_time =
      std::chrono::duration_cast<std::chrono::duration<double>>((do_t2 - do_t1)).count();

    if(rank == 0 && debug)
      std::cout << std::fixed << std::setprecision(2) << "diis: " << do_time << "s, ";
  }

  if(lshift > 0) {
    double lval = is_rhf ? 0.5 * lshift : lshift;
    // clang-format off
    sch
    (ehf_tmp(mu,ku) = S1(mu,nu) * D_last_alpha_tamm(nu,ku))
    (F_alpha(mu,ku) -= lval * ehf_tmp(mu,nu) * S1(nu,ku))
    .execute();
    // clang-format on

    if(is_uhf) {
      // clang-format off
      sch
      (ehf_tmp(mu,ku) = S1(mu,nu) * D_last_beta_tamm(nu,ku))
      (F_beta(mu,ku) -= lval * ehf_tmp(mu,nu) * S1(nu,ku))
      .execute();
      // clang-format on
    }
  }

  auto     do_t1 = std::chrono::high_resolution_clock::now();
  SCFGuess scf_guess;
  scf_guess.scf_diagonalize<TensorType>(sch, chem_env, scf_vars, scalapack_info, ttensors,
                                        etensors);

  auto do_t2   = std::chrono::high_resolution_clock::now();
  auto do_time = std::chrono::duration_cast<std::chrono::duration<double>>((do_t2 - do_t1)).count();

  if(rank == 0 && debug)
    std::cout << std::fixed << std::setprecision(2) << "diagonalize: " << do_time << "s, ";

  compute_density<TensorType>(ec, chem_env, scf_vars, scalapack_info, ttensors, etensors);

  double rmsd = 0.0;
  // clang-format off
  sch
    (D_diff()  = D_alpha_tamm())
    (D_diff() -= D_last_alpha_tamm())
    .execute();
  // clang-format on
  rmsd = norm(D_diff) / (double) (1.0 * N);

  if(is_uhf) {
    // clang-format off
    sch(D_diff() = D_beta_tamm())(D_diff() -= D_last_beta_tamm()).execute();
    // clang-format on
    rmsd += norm(D_diff) / (double) (1.0 * N);
  }

  const auto   damp  = scf_options.damp;
  const double alpha = damp / 100.0;
  // if(rmsd < 1e-6) { scf_vars.switch_diis = true; } // rmsd check

  if(damp < 100) {
    // D = alpha*D + (1.0-alpha)*D_last;
    if(is_rhf) {
      tamm::scale_ip(D_alpha_tamm(), alpha);
      sch(D_alpha_tamm() += (1.0 - alpha) * D_last_alpha_tamm()).execute();
      tamm_to_eigen_tensor(D_alpha_tamm, D_alpha);
    }
    if(is_uhf) {
      tamm::scale_ip(D_alpha_tamm(), alpha);
      sch(D_alpha_tamm() += (1.0 - alpha) * D_last_alpha_tamm()).execute();
      tamm_to_eigen_tensor(D_alpha_tamm, D_alpha);
      tamm::scale_ip(D_beta_tamm(), alpha);
      sch(D_beta_tamm() += (1.0 - alpha) * D_last_beta_tamm()).execute();
      tamm_to_eigen_tensor(D_beta_tamm, D_beta);
    }
  }

  return std::make_tuple(ehf, rmsd);
}

template<typename TensorType>
std::tuple<std::vector<int>, std::vector<int>, std::vector<int>>
exachem::scf::SCFIter::compute_2bf_taskinfo(ExecutionContext& ec, ChemEnv& chem_env,
                                            const SCFVars& scf_vars, const bool do_schwarz_screen,
                                            const std::vector<size_t>& shell2bf,
                                            const Matrix& SchwarzK, const size_t& max_nprim4,
                                            TAMMTensors& ttensors, EigenTensors& etensors,
                                            const bool cs1s2) {
  Matrix&             D       = etensors.D_alpha;
  Matrix&             D_beta  = etensors.D_beta;
  Tensor<TensorType>& F_dummy = ttensors.F_dummy;

  SystemData&              sys_data    = chem_env.sys_data;
  SCFOptions&              scf_options = chem_env.ioptions.scf_options;
  const libint2::BasisSet& obs         = chem_env.shells;

  double fock_precision = std::min(scf_options.tol_sch, 1e-2 * scf_options.conve);
  // auto       rank           = ec.pg().rank();
  const bool is_uhf = sys_data.is_unrestricted;

  Matrix D_shblk_norm = compute_shellblock_norm(obs, D); // matrix of infty-norms of shell blocks
  if(is_uhf) D_shblk_norm += compute_shellblock_norm(obs, D_beta);

  std::vector<int> s1vec;
  std::vector<int> s2vec;
  std::vector<int> ntask_vec;

  auto comp_2bf_lambda = [&](IndexVector blockid) {
    auto s1        = blockid[0];
    auto sp12_iter = scf_vars.obs_shellpair_data.at(s1).begin();

    auto s2     = blockid[1];
    auto s2spl  = scf_vars.obs_shellpair_list.at(s1);
    auto s2_itr = std::find(s2spl.begin(), s2spl.end(), s2);
    if(s2_itr == s2spl.end()) return;
    auto s2_pos = std::distance(s2spl.begin(), s2_itr);

    std::advance(sp12_iter, s2_pos);
    // const auto* sp12 = sp12_iter->get();

    const auto Dnorm12 = do_schwarz_screen ? D_shblk_norm(s1, s2) : 0.;

    size_t taskid = 0;

    for(decltype(s1) s3 = 0; s3 <= s1; ++s3) {
      const auto Dnorm123 =
        do_schwarz_screen ? std::max(D_shblk_norm(s1, s3), std::max(D_shblk_norm(s2, s3), Dnorm12))
                          : 0.;

      auto sp34_iter = scf_vars.obs_shellpair_data.at(s3).begin();

      const auto s4_max = (s1 == s3) ? s2 : s3;
      for(const auto& s4: scf_vars.obs_shellpair_list.at(s3)) {
        if(s4 > s4_max)
          break; // for each s3, s4 are stored in monotonically increasing
                 // order

        // must update the iter even if going to skip s4
        // const auto* sp34 = sp34_iter->get();
        ++sp34_iter;

        const auto Dnorm1234 =
          do_schwarz_screen
            ? std::max(D_shblk_norm(s1, s4),
                       std::max(D_shblk_norm(s2, s4), std::max(D_shblk_norm(s3, s4), Dnorm123)))
            : 0.;

        if(do_schwarz_screen && Dnorm1234 * SchwarzK(s1, s2) * SchwarzK(s3, s4) < fock_precision)
          continue;

        taskid++;
      }
    }
    // taskinfo << rank.value() << ", " << s1 << ", " << s2 << ", " << taskid << "\n";
    s1vec.push_back(s1);
    s2vec.push_back(s2);
    ntask_vec.push_back(taskid);
  };

  block_for(ec, F_dummy(), comp_2bf_lambda);
  return std::make_tuple(s1vec, s2vec, ntask_vec);
}

template<typename TensorType>
void exachem::scf::SCFIter::compute_3c_ints(ExecutionContext& ec, ChemEnv& chem_env,
                                            const SCFVars& scf_vars, Tensor<TensorType>& xyZ) {
  using libint2::BraKet;
  using libint2::Engine;
  using libint2::Operator;

  SCFOptions&              scf_options = chem_env.ioptions.scf_options;
  const libint2::BasisSet& obs         = chem_env.shells;

  auto                       rank           = ec.pg().rank();
  auto                       debug          = scf_options.debug;
  const std::vector<Tile>&   AO_tiles       = scf_vars.AO_tiles;
  const std::vector<size_t>& shell_tile_map = scf_vars.shell_tile_map;

  const libint2::BasisSet&   dfbs              = scf_vars.dfbs;
  const std::vector<Tile>&   dfAO_tiles        = scf_vars.dfAO_tiles;
  const std::vector<size_t>& df_shell_tile_map = scf_vars.df_shell_tile_map;

  Scheduler sch{ec};

  double      engine_precision = scf_options.tol_int; // default: 1e-22
  const auto& unitshell        = libint2::Shell::unit();
  auto        engine           = libint2::Engine(libint2::Operator::coulomb,
                                                 std::max(obs.max_nprim(), dfbs.max_nprim()),
                                                 std::max(obs.max_l(), dfbs.max_l()), 0);
  engine.set(libint2::BraKet::xs_xx);
  engine.set_precision(engine_precision);

  auto        shell2bf    = obs.shell2bf();
  auto        shell2bf_df = dfbs.shell2bf();
  const auto& results     = engine.results();

  auto compute_2body_fock_dfC_lambda = [&](const IndexVector& blockid) {
    auto bi0 = blockid[0];
    auto bi1 = blockid[1];
    auto bi2 = blockid[2];

    const TAMM_SIZE         size       = xyZ.block_size(blockid);
    auto                    block_dims = xyZ.block_dims(blockid);
    std::vector<TensorType> dbuf(size);

    auto bd1 = block_dims[1];
    auto bd2 = block_dims[2];

    auto s0range_end = df_shell_tile_map[bi2];
    // auto n0 = dfbs[s0].size();
    auto                  s1range_end   = shell_tile_map[bi0];
    decltype(s1range_end) s1range_start = 0l;
    if(bi0 > 0) s1range_start = shell_tile_map[bi0 - 1] + 1;

    tamm::Tile curshelloffset_i = 0U;
    for(auto s1 = s1range_start; s1 <= s1range_end; ++s1) {
      // auto n1 = shells[s1].size();
      auto                  dimi          = curshelloffset_i + AO_tiles[s1];
      auto                  s2range_end   = shell_tile_map[bi1];
      decltype(s2range_end) s2range_start = 0l;

      if(bi1 > 0) s2range_start = shell_tile_map[bi1 - 1] + 1;

      tamm::Tile curshelloffset_j = 0U;
      for(auto s2 = s2range_start; s2 <= s2range_end; ++s2) {
        //// if (s2>s1) continue;
        // auto n2 = shells[s2].size();
        // auto n123 = n0*n1*n2;
        // std::vector<TensorType> tbuf(n123);
        auto                  dimj          = curshelloffset_j + AO_tiles[s2];
        decltype(s0range_end) s0range_start = 0l;
        if(bi2 > 0) s0range_start = df_shell_tile_map[bi2 - 1] + 1;

        tamm::Tile curshelloffset_k = 0U;
        for(auto s0 = s0range_start; s0 <= s0range_end; ++s0) {
          engine.compute2<Operator::coulomb, BraKet::xs_xx, 0>(dfbs[s0], unitshell, obs[s1],
                                                               obs[s2]);
          const auto* buf = results[0];
          if(buf == nullptr) continue;
          // std::copy(buf, buf + n123, tbuf.begin());

          size_t c    = 0;
          auto   dimk = curshelloffset_k + dfAO_tiles[s0];

          for(auto k = curshelloffset_k; k < dimk; k++)
            for(auto i = curshelloffset_i; i < dimi; i++)
              for(auto j = curshelloffset_j; j < dimj; j++, c++)
                dbuf[(i * bd1 + j) * bd2 + k] = buf[c];

          curshelloffset_k += dfAO_tiles[s0];
        } // s2
        curshelloffset_j += AO_tiles[s2];
      } // s1
      curshelloffset_i += AO_tiles[s1];
    } // s0
    xyZ.put(blockid, dbuf);
  };

  auto do_t1 = std::chrono::high_resolution_clock::now();
  block_for(ec, xyZ(), compute_2body_fock_dfC_lambda);

  auto   do_t2 = std::chrono::high_resolution_clock::now();
  double do_time =
    std::chrono::duration_cast<std::chrono::duration<double>>((do_t2 - do_t1)).count();
  if(rank == 0 && debug) std::cout << "2BF-DFC: " << do_time << "s, ";
}

template<typename TensorType>
void exachem::scf::SCFIter::compute_2c_ints(ExecutionContext& ec, ChemEnv& chem_env,
                                            EigenTensors& etensors, const SCFVars& scf_vars,
                                            TAMMTensors& ttensors) {
  using libint2::BraKet;
  using libint2::Engine;
  using libint2::Operator;

  SCFOptions& scf_options = chem_env.ioptions.scf_options;

  auto rank = ec.pg().rank();

  const libint2::BasisSet&   dfbs              = scf_vars.dfbs;
  const std::vector<Tile>&   dfAO_tiles        = scf_vars.dfAO_tiles;
  const std::vector<size_t>& df_shell_tile_map = scf_vars.df_shell_tile_map;

  Tensor<TensorType>& Vm1 = ttensors.Vm1;

  auto& dfNorm = etensors.dfNorm;

  auto d_mu = scf_vars.d_mu, d_nu = scf_vars.d_nu, d_ku = scf_vars.d_ku;

  Scheduler sch{ec};

  double engine_precision = scf_options.tol_int; // default: 1e-22
  auto   engine = libint2::Engine(libint2::Operator::coulomb, dfbs.max_nprim(), dfbs.max_l(), 0);
  engine.set(libint2::BraKet::xs_xs);
  engine.set_precision(engine_precision);

  auto        shell2bf_df = dfbs.shell2bf();
  const auto& results     = engine.results();

  dfNorm.resize(dfbs.size());
  dfNorm.setZero();

  auto compute_2body_2index_ints_lambda = [&](const IndexVector& blockid) {
    auto bi0 = blockid[0];
    auto bi1 = blockid[1];

    const TAMM_SIZE         size       = Vm1.block_size(blockid);
    auto                    block_dims = Vm1.block_dims(blockid);
    std::vector<TensorType> dbuf(size);

    auto                  bd1           = block_dims[1];
    auto                  s1range_end   = df_shell_tile_map[bi0];
    decltype(s1range_end) s1range_start = 0l;

    if(bi0 > 0) s1range_start = df_shell_tile_map[bi0 - 1] + 1;

    for(auto s1 = s1range_start; s1 <= s1range_end; ++s1) {
      auto n1 = dfbs[s1].size();

      auto                  s2range_end   = df_shell_tile_map[bi1];
      decltype(s2range_end) s2range_start = 0l;
      if(bi1 > 0) s2range_start = df_shell_tile_map[bi1 - 1] + 1;

      for(auto s2 = s2range_start; s2 <= s2range_end; ++s2) {
        // if (s2>s1) continue;
        // if(s2>s1){ TODO: screening doesnt work - revisit
        //   auto s2spl = scf_vars.dfbs_shellpair_list.at(s2);
        //   if(std::find(s2spl.begin(),s2spl.end(),s1) == s2spl.end()) continue;
        // }
        // else{
        //   auto s2spl = scf_vars.dfbs_shellpair_list.at(s1);
        //   if(std::find(s2spl.begin(),s2spl.end(),s2) == s2spl.end()) continue;
        // }

        auto n2 = dfbs[s2].size();

        std::vector<TensorType> tbuf(n1 * n2);
        engine.compute(dfbs[s1], dfbs[s2]);
        if(results[0] == nullptr) continue;

        Eigen::Map<const Matrix> buf_mat(results[0], n1, n2);
        Eigen::Map<Matrix>(&tbuf[0], n1, n2) = buf_mat;

        if(s1 == s2) dfNorm(s1) = std::sqrt(buf_mat.maxCoeff());

        tamm::Tile curshelloffset_i = 0U;
        tamm::Tile curshelloffset_j = 0U;
        for(decltype(s1) x = s1range_start; x < s1; x++) curshelloffset_i += dfAO_tiles[x];
        for(decltype(s2) x = s2range_start; x < s2; x++) curshelloffset_j += dfAO_tiles[x];

        size_t c    = 0;
        auto   dimi = curshelloffset_i + dfAO_tiles[s1];
        auto   dimj = curshelloffset_j + dfAO_tiles[s2];
        for(auto i = curshelloffset_i; i < dimi; i++)
          for(auto j = curshelloffset_j; j < dimj; j++, c++) dbuf[i * bd1 + j] = tbuf[c];

        // TODO: not needed if screening works
        //  if (s1 != s2)  // if s1 >= s2, copy {s1,s2} to the corresponding {s2,s1}
        //               // block, note the transpose!
        //  result.block(bf2, bf1, n2, n1) = buf_mat.transpose();
      }
    }
    Vm1.put(blockid, dbuf);
  };
  block_for(ec, Vm1(), compute_2body_2index_ints_lambda);

  std::vector<double> norms(dfbs.size());
  ec.pg().allreduce(&dfNorm(0), norms.data(), dfbs.size(), tamm::ReduceOp::sum);
  for(size_t s1 = 0; s1 < dfbs.size(); ++s1) dfNorm(s1) = norms[s1];
}

template<typename TensorType>
void exachem::scf::SCFIter::init_ri(ExecutionContext& ec, ChemEnv& chem_env,
                                    ScalapackInfo& scalapack_info, const SCFVars& scf_vars,
                                    EigenTensors& etensors, TAMMTensors& ttensors) {
  SystemData& sys_data    = chem_env.sys_data;
  SCFOptions& scf_options = chem_env.ioptions.scf_options;

  const auto ndf    = sys_data.ndf;
  const bool direct = scf_vars.direct_df;
  auto       rank   = ec.pg().rank();
  auto       debug  = scf_options.debug;

  auto mu = scf_vars.mu, nu = scf_vars.nu, ku = scf_vars.ku;
  auto d_mu = scf_vars.d_mu, d_nu = scf_vars.d_nu, d_ku = scf_vars.d_ku;

  Scheduler   sch{ec};
  ExecutionHW exhw = ec.exhw();

  Tensor<TensorType>& xyK = ttensors.xyK;
  Tensor<TensorType>& xyZ = ttensors.xyZ;
  Tensor<TensorType>& Vm1 = ttensors.Vm1;

  // Compute 2c integrals
  Tensor<TensorType>::allocate(&ec, Vm1);
  sch(Vm1() = 0.0).execute();
  ec.pg().barrier();

  compute_2c_ints<TensorType>(ec, chem_env, etensors, scf_vars, ttensors);

  // Obtain inverse (square root) of V
  Matrix                  V;
  std::vector<TensorType> eps(ndf);

  auto ig1 = std::chrono::high_resolution_clock::now();

  Tensor<TensorType> v_tmp{scf_vars.tdfAO, scf_vars.tdfAO};
  Tensor<TensorType> eps_tamm{scf_vars.tdfAO};
  Tensor<TensorType>::allocate(&ec, v_tmp, eps_tamm);

#if defined(USE_SCALAPACK)
  Tensor<TensorType> V_sca;
  if(scalapack_info.comm != MPI_COMM_NULL) {
    blacspp::Grid*                  blacs_grid       = scalapack_info.blacs_grid.get();
    const auto&                     grid             = *blacs_grid;
    scalapackpp::BlockCyclicDist2D* blockcyclic_dist = scalapack_info.blockcyclic_dist.get();
    const tamm::Tile                mb               = blockcyclic_dist->mb();

    TiledIndexSpace    tN_bc{IndexSpace{range(ndf)}, mb};
    Tensor<TensorType> S_BC{tN_bc, tN_bc};
    V_sca = {tN_bc, tN_bc};
    S_BC.set_block_cyclic({scalapack_info.npr, scalapack_info.npc});
    V_sca.set_block_cyclic({scalapack_info.npr, scalapack_info.npc});
    Tensor<TensorType>::allocate(&scalapack_info.ec, S_BC, V_sca);

    tamm::to_block_cyclic_tensor(Vm1, S_BC);

    auto desc_lambda = [&](const int64_t M, const int64_t N) {
      auto [M_loc, N_loc] = (*blockcyclic_dist).get_local_dims(M, N);
      return (*blockcyclic_dist).descinit_noerror(M, N, M_loc);
    };

    if(grid.ipr() >= 0 and grid.ipc() >= 0) {
      auto desc_S = desc_lambda(ndf, ndf);
      auto desc_V = desc_lambda(ndf, ndf);

      /*info=*/scalapackpp::hereig(scalapackpp::Job::Vec, scalapackpp::Uplo::Lower, desc_S[2],
                                   S_BC.access_local_buf(), 1, 1, desc_S, eps.data(),
                                   V_sca.access_local_buf(), 1, 1, desc_V);
    }

    Tensor<TensorType>::deallocate(S_BC);
    tamm::from_block_cyclic_tensor(V_sca, v_tmp);
  }
#else
  if(rank == 0) {
    V.setZero(ndf, ndf);
    tamm_to_eigen_tensor(Vm1, V);

    lapack::syevd(lapack::Job::Vec, lapack::Uplo::Lower, ndf, V.data(), ndf, eps.data());
    tamm::eigen_to_tamm_tensor(v_tmp, V);
  }
#endif

  if(rank == 0) {
    std::transform(eps.begin(), eps.end(), eps.begin(), [](auto& c) {
      if(c < 1.0e-6) {
        std::cout << "WARNING: small eigenvalue in V: " << std::scientific << c << std::fixed
                  << std::endl;
        return 0.0;
      }
      return 1.0 / std::sqrt(c);
    });
    tamm::vector_to_tamm_tensor(eps_tamm, eps);
  }
  ec.pg().barrier();

  sch(Vm1(d_mu, d_nu) = v_tmp(d_nu, d_mu) * eps_tamm(d_nu)).deallocate(v_tmp, eps_tamm).execute();

  auto ig2    = std::chrono::high_resolution_clock::now();
  auto igtime = std::chrono::duration_cast<std::chrono::duration<double>>((ig2 - ig1)).count();

  if(rank == 0 && debug)
    std::cout << std::fixed << std::setprecision(2) << "V^-1/2: " << igtime << "s, ";

  if(!direct) {
    // Compute 3c ints
    Tensor<TensorType>::allocate(&ec, xyZ);
    compute_3c_ints<TensorType>(ec, chem_env, scf_vars, xyZ);

    // Orthonormalize DF basis
    sch(xyK(mu, nu, d_nu) = xyZ(mu, nu, d_mu) * Vm1(d_mu, d_nu))
      .deallocate(Vm1, xyZ) // release memory Zxy
      .execute(exhw);
  }
}

template<typename TensorType>
void exachem::scf::SCFIter::compute_2bf_ri_direct(ExecutionContext& ec, ChemEnv& chem_env,
                                                  const SCFVars&             scf_vars,
                                                  const std::vector<size_t>& shell2bf,
                                                  TAMMTensors& ttensors, EigenTensors& etensors,
                                                  const Matrix& SchwarzK) {
  using libint2::BraKet;
  using libint2::Engine;
  using libint2::Operator;

  SystemData& sys_data    = chem_env.sys_data;
  SCFOptions& scf_options = chem_env.ioptions.scf_options;

  Scheduler sch{ec};
  // ExecutionHW exhw = ec.exhw();

  const libint2::BasisSet& obs = chem_env.shells;

  const bool debug  = scf_options.debug;
  const bool is_uhf = sys_data.is_unrestricted;
  // const bool is_spherical = (scf_options.gaussian_type == "spherical");

  auto rank = ec.pg().rank();

  const auto               ndf         = sys_data.ndf;
  const libint2::BasisSet& dfbs        = scf_vars.dfbs;
  auto                     shell2bf_df = dfbs.shell2bf();

  auto mu = scf_vars.mu, nu = scf_vars.nu, ku = scf_vars.ku;
  auto d_mu = scf_vars.d_mu, d_nu = scf_vars.d_nu, d_ku = scf_vars.d_ku;

  double engine_precision = scf_options.tol_int; // default: 1e-22
  double fock_precision   = std::min(scf_options.tol_sch, 1e-2 * scf_options.conve);

  const auto& unitshell = libint2::Shell::unit();
  auto        engine    = libint2::Engine(libint2::Operator::coulomb,
                                          std::max(dfbs.max_nprim(), obs.max_nprim()),
                                          std::max(obs.max_l(), dfbs.max_l()), 0);
  engine.set(libint2::BraKet::xs_xx);
  engine.set_precision(engine_precision);
  const auto& buf = engine.results();

  Matrix& D      = etensors.D_alpha;
  Matrix& D_beta = etensors.D_beta;
  Matrix& G      = etensors.G_alpha;
  auto&   dfNorm = etensors.dfNorm;

  // auto   shblk = is_spherical ? 2 * obs.max_l() + 1 : ((obs.max_l() + 1) * (obs.max_l() + 2)) /
  // 2; Matrix Dblk(shblk, shblk);

  Matrix              Jtmp    = Matrix::Zero(1, ndf);
  Tensor<TensorType>& F_dummy = ttensors.F_dummy;
  Tensor<TensorType>& Vm1     = ttensors.Vm1;
  IndexSpace          dummy{range(1)};
  TiledIndexSpace     tdummy{dummy};
  Tensor<TensorType>  Jtmp_tamm{tdummy, scf_vars.tdfAO}, Xtmp_tamm{tdummy, scf_vars.tdfAO};

  const auto buildJ_start = std::chrono::high_resolution_clock::now();

  sch.allocate(Jtmp_tamm, Xtmp_tamm).execute();
  sch(Jtmp_tamm() = 0.0).execute();

  Matrix D_shblk_norm;
  D_shblk_norm = compute_shellblock_norm(obs, D); // matrix of infty-norms of shell blocks
  if(is_uhf) D_shblk_norm += compute_shellblock_norm(obs, D_beta);

  const auto       max_engine_precision    = std::numeric_limits<double>::epsilon() / 1e10;
  const auto       ln_max_engine_precision = std::log(max_engine_precision);
  shellpair_data_t spdata(1);
  for(size_t s1 = 0l; s1 != dfbs.size(); ++s1) {
    spdata[0].emplace_back(
      std::make_shared<libint2::ShellPair>(dfbs[s1], unitshell, ln_max_engine_precision));
  }

  auto comp_J_lambda = [&](IndexVector blockid) {
    auto s1        = blockid[0];
    auto bf1_first = shell2bf[s1];
    auto n1        = obs[s1].size();
    auto sp12_iter = scf_vars.obs_shellpair_data.at(s1).begin();

    auto s2     = blockid[1];
    auto s2spl  = scf_vars.obs_shellpair_list.at(s1);
    auto s2_itr = std::find(s2spl.begin(), s2spl.end(), s2);
    if(s2_itr == s2spl.end()) return;
    auto s2_pos    = std::distance(s2spl.begin(), s2_itr);
    auto bf2_first = shell2bf[s2];
    auto n2        = obs[s2].size();

    std::advance(sp12_iter, s2_pos);
    const auto* sp12 = sp12_iter->get();

    const auto Norm12  = D_shblk_norm(s1, s2) * SchwarzK(s1, s2);
    const auto Jfactor = (s1 == s2) ? 1.0 : 2.0;

    if(Norm12 * dfNorm.maxCoeff() < fock_precision) return;

    Matrix Dblk = D.block(bf1_first, bf2_first, n1, n2);
    if(is_uhf) Dblk.block(0, 0, n1, n2) += D_beta.block(bf1_first, bf2_first, n1, n2);

    auto sp_iter = spdata.at(0).begin();
    for(decltype(s1) s3 = 0; s3 < dfbs.size(); ++s3) {
      auto        bf3_first = shell2bf_df[s3];
      auto        n3        = dfbs[s3].size();
      const auto* sp        = sp_iter->get();
      ++sp_iter;

      if(Norm12 * dfNorm(s3) < fock_precision) continue;

      engine.prescale_by(Jfactor);
      engine.compute2<Operator::coulomb, BraKet::xs_xx, 0>(dfbs[s3], unitshell, obs[s1], obs[s2],
                                                           sp, sp12);
      const auto* buf_312 = buf[0];
      if(buf_312 == nullptr) continue;

      for(decltype(n1) bf3 = bf3_first, f312 = 0; bf3 != bf3_first + n3; ++bf3)
        for(decltype(n1) f1 = 0; f1 != n1; ++f1)
          for(decltype(n2) f2 = 0; f2 != n2; ++f2, ++f312)
            Jtmp(0, bf3) += buf_312[f312] * Dblk(f1, f2);
    }
  };
  block_for(ec, F_dummy(), comp_J_lambda);

  eigen_to_tamm_tensor_acc(Jtmp_tamm, Jtmp);
  ec.pg().barrier();

  // const auto buildJ_stop = std::chrono::high_resolution_clock::now();
  // const auto buildJ_time =
  //   std::chrono::duration_cast<std::chrono::duration<double>>((buildJ_stop -
  //   buildJ_start)).count();
  // if(rank == 0 && debug) std::cout << "buildJ: " << buildJ_time << "s, ";

  // const auto buildX_start = std::chrono::high_resolution_clock::now();
  sch(Xtmp_tamm("i", "j") = Jtmp_tamm("i", "k") * Vm1("k", "j"))(
    Jtmp_tamm("i", "j") = Xtmp_tamm("i", "k") * Vm1("j", "k"))
    .deallocate(Xtmp_tamm)
    .execute();

  Jtmp.setZero();
  tamm_to_eigen_tensor(Jtmp_tamm, Jtmp);
  ec.pg().barrier();
  Tensor<TensorType>::deallocate(Jtmp_tamm);
  // const auto buildX_stop = std::chrono::high_resolution_clock::now();
  // const auto buildX_time =
  //   std::chrono::duration_cast<std::chrono::duration<double>>((buildX_stop -
  //   buildX_start)).count();
  // if(rank == 0 && debug) std::cout << "buildX: " << buildX_time << "s, ";

  // const auto                            buildF_start = std::chrono::high_resolution_clock::now();
  Eigen::VectorXd                       Jvec(Eigen::Map<Eigen::VectorXd>(Jtmp.data(), Jtmp.cols()));
  Eigen::Vector<double, Eigen::Dynamic> Jnorm;
  Jnorm.resize(dfbs.size());
  for(size_t s3 = 0; s3 < dfbs.size(); ++s3) {
    Jnorm(s3) =
      Jvec.segment(shell2bf_df[s3], dfbs[s3].size()).lpNorm<Eigen::Infinity>() * dfNorm(s3);
  }

  // Scale JVEC to account for permutational symmetry
  Jvec *= 2.0;

  auto comp_df_lambda = [&](IndexVector blockid) {
    auto s1        = blockid[0];
    auto bf1_first = shell2bf[s1];
    auto n1        = obs[s1].size();
    auto sp12_iter = scf_vars.obs_shellpair_data.at(s1).begin();

    auto s2     = blockid[1];
    auto s2spl  = scf_vars.obs_shellpair_list.at(s1);
    auto s2_itr = std::find(s2spl.begin(), s2spl.end(), s2);
    if(s2_itr == s2spl.end()) return;
    auto s2_pos    = std::distance(s2spl.begin(), s2_itr);
    auto bf2_first = shell2bf[s2];
    auto n2        = obs[s2].size();

    std::advance(sp12_iter, s2_pos);
    const auto* sp12 = sp12_iter->get();

    const auto Norm12  = SchwarzK(s1, s2);
    const auto Jfactor = (s1 == s2) ? 0.5 : 1.0;
    if(Norm12 * Jnorm.maxCoeff() < fock_precision) return;

    Matrix J12 = Matrix::Zero(n1, n2);

    auto sp_iter = spdata.at(0).begin();
    for(decltype(s1) s3 = 0; s3 < dfbs.size(); ++s3) {
      auto        bf3_first = shell2bf_df[s3];
      auto        n3        = dfbs[s3].size();
      const auto* sp        = sp_iter->get();
      ++sp_iter;
      if(Norm12 * Jnorm(s3) < fock_precision) continue;

      engine.prescale_by(Jfactor);
      engine.compute2<Operator::coulomb, BraKet::xs_xx, 0>(dfbs[s3], unitshell, obs[s1], obs[s2],
                                                           sp, sp12);
      const auto* buf_312 = buf[0];
      if(buf_312 == nullptr) continue;

      for(decltype(n1) bf3 = bf3_first, f312 = 0; bf3 != bf3_first + n3; ++bf3)
        for(decltype(n1) f1 = 0; f1 != n1; ++f1)
          for(decltype(n2) f2 = 0; f2 != n2; ++f2, ++f312) J12(f1, f2) += buf_312[f312] * Jvec(bf3);
    }
    G.block(bf1_first, bf2_first, n1, n2) = J12;
  };
  block_for(ec, F_dummy(), comp_df_lambda);

  eigen_to_tamm_tensor_acc(ttensors.F_alpha_tmp, G);
  ec.pg().barrier();

  if(is_uhf) sch(ttensors.F_beta_tmp() = ttensors.F_alpha_tmp()).execute();

  const auto buildF_stop = std::chrono::high_resolution_clock::now();
  // const auto buildF_time =
  //   std::chrono::duration_cast<std::chrono::duration<double>>((buildF_stop -
  //   buildF_start)).count();
  // if(rank == 0 && debug) std::cout << "buildF: " << buildF_time << "s, ";

  const auto total_time =
    std::chrono::duration_cast<std::chrono::duration<double>>((buildF_stop - buildJ_start)).count();
  if(rank == 0 && debug)
    std::cout << std::fixed << std::setprecision(2) << "DF-J: " << total_time << "s, ";
};

template<typename TensorType>
void exachem::scf::SCFIter::compute_2bf_ri(ExecutionContext& ec, ChemEnv& chem_env,
                                           ScalapackInfo& scalapack_info, const SCFVars& scf_vars,
                                           const std::vector<size_t>& shell2bf,
                                           TAMMTensors& ttensors, EigenTensors& etensors,
                                           bool& is_3c_init, double xHF) {
  SystemData& sys_data    = chem_env.sys_data;
  SCFOptions& scf_options = chem_env.ioptions.scf_options;

  const bool is_uhf = sys_data.is_unrestricted;
  const bool is_rhf = sys_data.is_restricted;
  const bool do_snK = sys_data.do_snK;

  auto rank  = ec.pg().rank();
  auto debug = scf_options.debug;

  auto mu = scf_vars.mu, nu = scf_vars.nu, ku = scf_vars.ku;
  auto d_mu = scf_vars.d_mu, d_nu = scf_vars.d_nu, d_ku = scf_vars.d_ku;
  // const tamm::TiledIndexLabel& dCocc_til = scf_vars.dCocc_til;

  Scheduler   sch{ec};
  ExecutionHW exhw = ec.exhw();

  Tensor<TensorType>& xyK         = ttensors.xyK;
  Tensor<TensorType>& F_alpha_tmp = ttensors.F_alpha_tmp;
  Tensor<TensorType>& F_beta_tmp  = ttensors.F_beta_tmp;

  Tensor<TensorType> Jtmp_tamm{scf_vars.tdfAO};                          // ndf
  Tensor<TensorType> tmp_df{scf_vars.tAO, scf_vars.tAO, scf_vars.tdfAO}; // n, n, ndf

  auto ig1  = std::chrono::high_resolution_clock::now();
  auto tig1 = ig1;

  // contract(1.0, xyK, {1, 2, 3}, Co, {2, 4}, 0.0, xiK, {1, 4, 3});

  sch.allocate(Jtmp_tamm).execute();

  ig1 = std::chrono::high_resolution_clock::now();
  // compute Coulomb
  // clang-format off
    if (is_uhf) sch(ttensors.D_alpha() += ttensors.D_beta());

    sch
      (Jtmp_tamm(d_mu) = xyK(mu,nu,d_mu) * ttensors.D_alpha(mu,nu))
      (F_alpha_tmp(mu, nu) = xyK(mu, nu, d_mu) * Jtmp_tamm(d_mu));
    
    if (is_uhf) {
      sch
        (F_beta_tmp() = F_alpha_tmp())
        (ttensors.D_alpha() -= ttensors.D_beta());
    }
    sch.deallocate(Jtmp_tamm).execute();
  // clang-format on

  auto ig2    = std::chrono::high_resolution_clock::now();
  auto igtime = std::chrono::duration_cast<std::chrono::duration<double>>((ig2 - ig1)).count();
  if(rank == 0 && debug)
    std::cout << " J: " << std::fixed << std::setprecision(2) << igtime << "s, ";

  if(xHF > 0.0 && !do_snK) {
    ig1 = std::chrono::high_resolution_clock::now();

    TensorType factor = is_rhf ? 0.5 : 1.0;
    // clang-format off
      sch.allocate(tmp_df);
      sch
        (tmp_df(mu, ku, d_mu) = xyK(mu, nu, d_mu) * ttensors.D_alpha(nu, ku))
        (F_alpha_tmp(mu, nu) += -1.0 * factor * xHF * tmp_df(mu, ku, d_mu) * xyK(nu, ku, d_mu));
    
      if (is_uhf) {
        sch
          (tmp_df(mu, ku, d_mu) = xyK(mu, nu, d_mu) * ttensors.D_beta(nu, ku))
          (F_beta_tmp(mu, nu) += -1.0 * xHF * tmp_df(mu, ku, d_mu) * xyK(nu, ku, d_mu));
      }
      sch.deallocate(tmp_df).execute(exhw);
    // clang-format on

    ig2    = std::chrono::high_resolution_clock::now();
    igtime = std::chrono::duration_cast<std::chrono::duration<double>>((ig2 - ig1)).count();
    if(rank == 0 && debug)
      std::cout << " K: " << std::fixed << std::setprecision(2) << igtime << "s, ";
  }

  auto tig2    = std::chrono::high_resolution_clock::now();
  auto tigtime = std::chrono::duration_cast<std::chrono::duration<double>>((tig2 - tig1)).count();

  if(rank == 0 && debug)
    std::cout << "3c contractions: " << std::fixed << std::setprecision(2) << tigtime << "s, ";
}

template<typename TensorType>
void exachem::scf::SCFIter::compute_2bf(ExecutionContext& ec, ChemEnv& chem_env,
                                        ScalapackInfo& scalapack_info, const SCFVars& scf_vars,
                                        const bool                 do_schwarz_screen,
                                        const std::vector<size_t>& shell2bf, const Matrix& SchwarzK,
                                        const size_t& max_nprim4, TAMMTensors& ttensors,
                                        EigenTensors& etensors, bool& is_3c_init,
                                        const bool do_density_fitting, double xHF) {
  using libint2::Operator;

  SystemData&              sys_data    = chem_env.sys_data;
  SCFOptions&              scf_options = chem_env.ioptions.scf_options;
  Scheduler                sch{ec};
  const libint2::BasisSet& obs = chem_env.shells;

  const bool is_uhf       = sys_data.is_unrestricted;
  const bool is_rhf       = sys_data.is_restricted;
  const bool do_snK       = sys_data.do_snK;
  const bool doK          = xHF != 0.0 && !do_snK;
  const bool is_spherical = (scf_options.gaussian_type == "spherical");

  Matrix& G      = etensors.G_alpha;
  Matrix& D      = etensors.D_alpha;
  Matrix& G_beta = etensors.G_beta;
  Matrix& D_beta = etensors.D_beta;

  Tensor<TensorType>& F_dummy     = ttensors.F_dummy;
  Tensor<TensorType>& F_alpha_tmp = ttensors.F_alpha_tmp;
  Tensor<TensorType>& F_beta_tmp  = ttensors.F_beta_tmp;

  double fock_precision = std::min(scf_options.tol_sch, 1e-2 * scf_options.conve);
  auto   rank           = ec.pg().rank();
  auto   N              = sys_data.nbf_orig;
  auto   debug          = scf_options.debug;

  auto   do_t1 = std::chrono::high_resolution_clock::now();
  Matrix D_shblk_norm;

  double Kfactor = is_rhf ? -0.25 * xHF : -0.5 * xHF;

  double engine_precision = scf_options.tol_int; // default: 1e-22

  // construct the 2-electron repulsion integrals engine pool
  using libint2::Engine;
  Engine engine(Operator::coulomb, obs.max_nprim(), obs.max_l(), 0);

  engine.set_precision(engine_precision);
  const auto& buf = engine.results();

  auto shblk = is_spherical ? 2 * obs.max_l() + 1 : ((obs.max_l() + 1) * (obs.max_l() + 2)) / 2;

  // To avoid skipping over the whole Fock matrix
  Matrix J12(shblk, shblk), J34(shblk, shblk);
  Matrix D12(shblk, shblk), D34(shblk, shblk);
  Matrix K1_alpha(shblk, N), K1_beta(shblk, N);
  Matrix K2_alpha(shblk, N), K2_beta(shblk, N);

  auto comp_2bf_lambda = [&](IndexVector blockid) {
    auto s1        = blockid[0];
    auto bf1_first = shell2bf[s1];
    auto n1        = obs[s1].size();
    auto sp12_iter = scf_vars.obs_shellpair_data.at(s1).begin();

    auto s2     = blockid[1];
    auto s2spl  = scf_vars.obs_shellpair_list.at(s1);
    auto s2_itr = std::find(s2spl.begin(), s2spl.end(), s2);
    if(s2_itr == s2spl.end()) return;
    auto s2_pos    = std::distance(s2spl.begin(), s2_itr);
    auto bf2_first = shell2bf[s2];
    auto n2        = obs[s2].size();

    std::advance(sp12_iter, s2_pos);
    const auto* sp12 = sp12_iter->get();

    const auto Dnorm12 = do_schwarz_screen ? 2 * D_shblk_norm(s1, s2) : 0.;

    // To avoid skipping over the whole Fock matrix
    J12.setZero();

    if(doK) {
      K1_alpha.setZero();
      K2_alpha.setZero();
      if(is_uhf) {
        K1_beta.setZero();
        K2_beta.setZero();
      }
    }

    // For Coulomb part
    D12.block(0, 0, n1, n2) = D.block(bf1_first, bf2_first, n1, n2);
    if(is_uhf) D12.block(0, 0, n1, n2) += D_beta.block(bf1_first, bf2_first, n1, n2);

    for(decltype(s1) s3 = 0; s3 <= s1; ++s3) {
      auto bf3_first = shell2bf[s3];
      auto n3        = obs[s3].size();

      const auto Dnorm123 =
        do_schwarz_screen ? std::max(D_shblk_norm(s1, s3), std::max(D_shblk_norm(s2, s3), Dnorm12))
                          : 0.;

      auto sp34_iter = scf_vars.obs_shellpair_data.at(s3).begin();

      const auto s4_max = (s1 == s3) ? s2 : s3;
      for(const auto& s4: scf_vars.obs_shellpair_list.at(s3)) {
        if(s4 > s4_max)
          break; // for each s3, s4 are stored in monotonically increasing
                 // order

        // must update the iter even if going to skip s4
        const auto* sp34 = sp34_iter->get();
        ++sp34_iter;

        const auto Dnorm34 = do_schwarz_screen ? 2 * D_shblk_norm(s3, s4) : 0.0;

        const auto Dnorm1234 =
          do_schwarz_screen
            ? std::max(D_shblk_norm(s1, s4),
                       std::max(D_shblk_norm(s2, s4), std::max(D_shblk_norm(s3, s4), Dnorm123)))
            : 0.;

        if(do_schwarz_screen && Dnorm1234 * SchwarzK(s1, s2) * SchwarzK(s3, s4) < fock_precision)
          continue;

        if(!doK) {
          if(do_schwarz_screen && Dnorm12 * SchwarzK(s1, s2) * SchwarzK(s3, s4) < fock_precision &&
             Dnorm34 * SchwarzK(s1, s2) * SchwarzK(s3, s4) < fock_precision)
            continue;
        }

        auto bf4_first = shell2bf[s4];
        auto n4        = obs[s4].size();

        // For Coulomb part
        D34.block(0, 0, n3, n4) = D.block(bf3_first, bf4_first, n3, n4);
        if(is_uhf) D34.block(0, 0, n3, n4) += D_beta.block(bf3_first, bf4_first, n3, n4);

        // compute the permutational degeneracy (i.e. # of equivalents) of
        // the given shell set
        auto s12_deg    = (s1 == s2) ? 1 : 2;
        auto s34_deg    = (s3 == s4) ? 1 : 2;
        auto s12_34_deg = (s1 == s3) ? (s2 == s4 ? 1 : 2) : 2;
        auto s1234_deg  = s12_deg * s34_deg * s12_34_deg;

        // prescale the integrals inside Libint
        engine.prescale_by(0.5 * s1234_deg);

        // compute integrals
        engine.compute2<Operator::coulomb, libint2::BraKet::xx_xx, 0>(obs[s1], obs[s2], obs[s3],
                                                                      obs[s4], sp12, sp34);

        const auto* buf_1234 = buf[0];
        if(buf_1234 == nullptr) continue; // if all integrals screened out, skip to next quartet

        // 1) each shell set of integrals contributes up to 6 shell sets of
        // the Fock matrix:
        //    F(a,b) += 1/2 * (ab|cd) * D(c,d)
        //    F(c,d) += 1/2 * (ab|cd) * D(a,b)
        //    F(b,d) -= 1/8 * (ab|cd) * D(a,c)
        //    F(b,c) -= 1/8 * (ab|cd) * D(a,d)
        //    F(a,c) -= 1/8 * (ab|cd) * D(b,d)
        //    F(a,d) -= 1/8 * (ab|cd) * D(b,c)
        // 2) each permutationally-unique integral (shell set) must be
        // scaled by its degeneracy,
        //    i.e. the number of the integrals/sets equivalent to it
        // 3) the end result must be symmetrized

        // J and K for closed-shell calculations
        if(doK && is_rhf) {
          for(decltype(n1) f1 = 0, f1234 = 0; f1 != n1; ++f1) {
            const auto bf1 = f1 + bf1_first;
            for(decltype(n2) f2 = 0; f2 != n2; ++f2) {
              const auto bf2 = f2 + bf2_first;
              for(decltype(n3) f3 = 0; f3 != n3; ++f3) {
                const auto bf3 = f3 + bf3_first;
                for(decltype(n4) f4 = 0; f4 != n4; ++f4, ++f1234) {
                  const auto bf4               = f4 + bf4_first;
                  auto       value_scal_by_deg = buf_1234[f1234];
                  J12(f1, f2) += D34(f3, f4) * value_scal_by_deg;
                  G(bf3, bf4) += D12(f1, f2) * value_scal_by_deg;

                  value_scal_by_deg *= Kfactor;
                  K1_alpha(f1, bf3) += D(bf2, bf4) * value_scal_by_deg;
                  K1_alpha(f1, bf4) += D(bf2, bf3) * value_scal_by_deg;
                  K2_alpha(f2, bf3) += D(bf1, bf4) * value_scal_by_deg;
                  K2_alpha(f2, bf4) += D(bf1, bf3) * value_scal_by_deg;
                }
              }
            }
          }
        }
        else if(doK && is_uhf) {
          for(decltype(n1) f1 = 0, f1234 = 0; f1 != n1; ++f1) {
            const auto bf1 = f1 + bf1_first;
            for(decltype(n2) f2 = 0; f2 != n2; ++f2) {
              const auto bf2 = f2 + bf2_first;
              for(decltype(n3) f3 = 0; f3 != n3; ++f3) {
                const auto bf3 = f3 + bf3_first;
                for(decltype(n4) f4 = 0; f4 != n4; ++f4, ++f1234) {
                  const auto bf4               = f4 + bf4_first;
                  auto       value_scal_by_deg = buf_1234[f1234];
                  auto       J34               = D12(f1, f2) * value_scal_by_deg;
                  J12(f1, f2) += D34(f3, f4) * value_scal_by_deg;
                  G(bf3, bf4) += J34;
                  G_beta(bf3, bf4) += J34;

                  value_scal_by_deg *= Kfactor;
                  K1_alpha(f1, bf3) += D(bf2, bf4) * value_scal_by_deg;
                  K1_alpha(f1, bf4) += D(bf2, bf3) * value_scal_by_deg;
                  K2_alpha(f2, bf3) += D(bf1, bf4) * value_scal_by_deg;
                  K2_alpha(f2, bf4) += D(bf1, bf3) * value_scal_by_deg;

                  K1_beta(f1, bf3) += D_beta(bf2, bf4) * value_scal_by_deg;
                  K1_beta(f1, bf4) += D_beta(bf2, bf3) * value_scal_by_deg;
                  K2_beta(f2, bf3) += D_beta(bf1, bf4) * value_scal_by_deg;
                  K2_beta(f2, bf4) += D_beta(bf1, bf3) * value_scal_by_deg;
                }
              }
            }
          }
        }
        else if(is_rhf) {
          for(decltype(n1) f1 = 0, f1234 = 0; f1 != n1; ++f1) {
            for(decltype(n2) f2 = 0; f2 != n2; ++f2) {
              for(decltype(n3) f3 = 0; f3 != n3; ++f3) {
                const auto bf3 = f3 + bf3_first;
                for(decltype(n4) f4 = 0; f4 != n4; ++f4, ++f1234) {
                  const auto bf4               = f4 + bf4_first;
                  auto       value_scal_by_deg = buf_1234[f1234];
                  J12(f1, f2) += D34(f3, f4) * value_scal_by_deg;
                  G(bf3, bf4) += D12(f1, f2) * value_scal_by_deg;
                }
              }
            }
          }
        }
        else if(is_uhf) {
          for(decltype(n1) f1 = 0, f1234 = 0; f1 != n1; ++f1) {
            for(decltype(n2) f2 = 0; f2 != n2; ++f2) {
              for(decltype(n3) f3 = 0; f3 != n3; ++f3) {
                const auto bf3 = f3 + bf3_first;
                for(decltype(n4) f4 = 0; f4 != n4; ++f4, ++f1234) {
                  const auto bf4               = f4 + bf4_first;
                  auto       value_scal_by_deg = buf_1234[f1234];
                  auto       J34               = D12(f1, f2) * value_scal_by_deg;
                  J12(f1, f2) += D34(f3, f4) * value_scal_by_deg;
                  G(bf3, bf4) += J34;
                  G_beta(bf3, bf4) += J34;
                }
              }
            }
          }
        }
      }
    }
    // Add contributions to (s1,s2) block
    G.block(bf1_first, bf2_first, n1, n2) += J12.block(0, 0, n1, n2);
    if(is_uhf) G_beta.block(bf1_first, bf2_first, n1, n2) += J12.block(0, 0, n1, n2);

    // Add contributions to (s1,N) and (s2,N) blocks
    if(doK) {
      G.block(bf1_first, 0, n1, N) += K1_alpha.block(0, 0, n1, N);
      G.block(bf2_first, 0, n2, N) += K2_alpha.block(0, 0, n2, N);
      if(is_uhf) {
        G_beta.block(bf1_first, 0, n1, N) += K1_beta.block(0, 0, n1, N);
        G_beta.block(bf2_first, 0, n2, N) += K2_beta.block(0, 0, n2, N);
      }
    }
  };

  decltype(do_t1) do_t2;
  double          do_time;

  if(!do_density_fitting) {
    D_shblk_norm = compute_shellblock_norm(obs, D); // matrix of infty-norms of shell blocks
    if(is_uhf) D_shblk_norm = D_shblk_norm.cwiseMax(compute_shellblock_norm(obs, D_beta));

    G.setZero(N, N);
    if(is_uhf) G_beta.setZero(N, N);
    if(!scf_vars.do_load_bal) block_for(ec, F_dummy(), comp_2bf_lambda);
    else {
      for(Eigen::Index i1 = 0; i1 < etensors.taskmap.rows(); i1++)
        for(Eigen::Index j1 = 0; j1 < etensors.taskmap.cols(); j1++) {
          if(etensors.taskmap(i1, j1) == -1 || etensors.taskmap(i1, j1) != rank) continue;
          IndexVector blockid{(tamm::Index) i1, (tamm::Index) j1};
          comp_2bf_lambda(blockid);
        }
      ec.pg().barrier();
    }

    // Matrix Gt = 0.5 * (G + G.transpose()); G=Gt
    // Gt     = 0.5 * (G_beta + G_beta.transpose());
    eigen_to_tamm_tensor_acc(F_alpha_tmp, G);
    if(is_uhf) eigen_to_tamm_tensor_acc(F_beta_tmp, G_beta);

    do_t2   = std::chrono::high_resolution_clock::now();
    do_time = std::chrono::duration_cast<std::chrono::duration<double>>((do_t2 - do_t1)).count();

    if(rank == 0 && debug)
      std::cout << std::fixed << std::setprecision(2) << "Fock build: " << do_time << "s, ";

    // ec.pg().barrier();
  }
  else {
    if(scf_vars.direct_df) {
      G.setZero(N, N);
      compute_2bf_ri_direct<TensorType>(ec, chem_env, scf_vars, shell2bf, ttensors, etensors,
                                        SchwarzK);
    }
    else {
      compute_2bf_ri<TensorType>(ec, chem_env, scalapack_info, scf_vars, shell2bf, ttensors,
                                 etensors, is_3c_init, xHF);
    }
  } // end density fitting

  auto [mu, nu]               = scf_vars.tAO.labels<2>("all");
  Tensor<TensorType>& H1      = ttensors.H1;
  Tensor<TensorType>& F_alpha = ttensors.F_alpha;
  Tensor<TensorType>& F_beta  = ttensors.F_beta;

  // symmetrize the result
  // clang-format off
  sch 
      (F_alpha()       = 0.5 * F_alpha_tmp())
      (F_alpha(mu,nu) += 0.5 * F_alpha_tmp(nu,mu))
      (F_alpha()      += H1())
      .execute();
  // clang-format on

  if(is_uhf) {
    // clang-format off
    sch 
      (F_beta()       = 0.5 * F_beta_tmp())
      (F_beta(mu,nu) += 0.5 * F_beta_tmp(nu,mu))
      (F_beta()      += H1())
      .execute();
    // clang-format on
  }
}

template<typename TensorType>
void exachem::scf::SCFIter::scf_diis(ExecutionContext& ec, ChemEnv& chem_env,
                                     const TiledIndexSpace& tAO, Tensor<TensorType> F_alpha,
                                     Tensor<TensorType> F_beta, Tensor<TensorType> err_mat_alpha,
                                     Tensor<TensorType> err_mat_beta, int iter, int max_hist,
                                     const SCFVars& scf_vars, const int n_lindep,
                                     std::vector<Tensor<TensorType>>& diis_hist_alpha,
                                     std::vector<Tensor<TensorType>>& diis_hist_beta,
                                     std::vector<Tensor<TensorType>>& fock_hist_alpha,
                                     std::vector<Tensor<TensorType>>& fock_hist_beta) {
  using Vector = Eigen::Matrix<TensorType, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

  tamm::Scheduler sch{ec};
  SystemData&     sys_data = chem_env.sys_data;
  bool const      is_uhf   = sys_data.is_unrestricted;

  auto rank  = ec.pg().rank().value();
  auto ndiis = scf_vars.idiis;

  if(ndiis > max_hist) {
    auto maxe = 0;
    if(!scf_vars.switch_diis) {
      std::vector<TensorType> max_err(diis_hist_alpha.size());
      max_err[diis_hist_alpha.size() - 1] = 0.0;
      for(size_t i = 0; i < diis_hist_alpha.size() - 1; i++) {
        max_err[i] = tamm::norm(diis_hist_alpha[i]);
      }
      if(is_uhf) {
        for(size_t i = 0; i < diis_hist_beta.size() - 1; i++) {
          max_err[i] = pow(max_err[i], 2) + pow(tamm::norm(diis_hist_beta[i]), 2);
        }
      }
      maxe = std::distance(max_err.begin(), std::max_element(max_err.begin(), max_err.end()));
    }
    Tensor<TensorType>::deallocate(diis_hist_alpha[maxe]);
    Tensor<TensorType>::deallocate(fock_hist_alpha[maxe]);
    diis_hist_alpha.erase(diis_hist_alpha.begin() + maxe);
    fock_hist_alpha.erase(fock_hist_alpha.begin() + maxe);
    if(is_uhf) {
      Tensor<TensorType>::deallocate(diis_hist_beta[maxe]);
      Tensor<TensorType>::deallocate(fock_hist_beta[maxe]);
      diis_hist_beta.erase(diis_hist_beta.begin() + maxe);
      fock_hist_beta.erase(fock_hist_beta.begin() + maxe);
    }
  }
  else {
    if(ndiis == (int) (max_hist / 2) && n_lindep > 1) {
      for(auto x: diis_hist_alpha) Tensor<TensorType>::deallocate(x);
      for(auto x: fock_hist_alpha) Tensor<TensorType>::deallocate(x);
      diis_hist_alpha.clear();
      fock_hist_alpha.clear();
      if(is_uhf) {
        for(auto x: diis_hist_beta) Tensor<TensorType>::deallocate(x);
        for(auto x: fock_hist_beta) Tensor<TensorType>::deallocate(x);
        diis_hist_beta.clear();
        fock_hist_beta.clear();
      }
    }
  }

  Tensor<TensorType> Fcopy{tAO, tAO};
  Tensor<TensorType> Fcopy_beta{tAO, tAO};
  Tensor<TensorType>::allocate(&ec, Fcopy);
  sch(Fcopy() = F_alpha()).execute();

  diis_hist_alpha.push_back(err_mat_alpha);
  fock_hist_alpha.push_back(Fcopy);

  if(is_uhf) {
    Tensor<TensorType>::allocate(&ec, Fcopy_beta);
    sch(Fcopy_beta() = F_beta()).execute();
    diis_hist_beta.push_back(err_mat_beta);
    fock_hist_beta.push_back(Fcopy_beta);
  }

  Matrix                  A;
  Vector                  b;
  int                     idim = std::min((int) diis_hist_alpha.size(), max_hist);
  int64_t                 info = -1;
  int64_t                 N    = idim + 1;
  std::vector<TensorType> X;

  // Tensor<TensorType> dhi_trans{tAO, tAO};
  Tensor<TensorType> dhi_trace{};
  auto [mu, nu] = tAO.labels<2>("all");
  Tensor<TensorType>::allocate(&ec, dhi_trace); // dhi_trans

  // ----- Construct Pulay matrix -----
  A = Matrix::Zero(idim + 1, idim + 1);

  for(int i = 0; i < idim; i++) {
    for(int j = i; j < idim; j++) {
      // A(i, j) = (diis_hist[i].transpose() * diis_hist[j]).trace();
      sch(dhi_trace() = diis_hist_alpha[i](nu, mu) * diis_hist_alpha[j](nu, mu)).execute();
      TensorType dhi  = get_scalar(dhi_trace); // dhi_trace.trace();
      A(i + 1, j + 1) = dhi;
    }
  }

  if(is_uhf) {
    for(int i = 0; i < idim; i++) {
      for(int j = i; j < idim; j++) {
        // A(i, j) = (diis_hist[i].transpose() * diis_hist[j]).trace();
        sch(dhi_trace() = diis_hist_beta[i](nu, mu) * diis_hist_beta[j](nu, mu)).execute();
        TensorType dhi = get_scalar(dhi_trace); // dhi_trace.trace();
        A(i + 1, j + 1) += dhi;
      }
    }
  }

  for(int i = 1; i <= idim; i++) {
    for(int j = i; j <= idim; j++) { A(j, i) = A(i, j); }
  }

  for(int i = 1; i <= idim; i++) {
    A(i, 0) = -1.0;
    A(0, i) = -1.0;
  }

  while(info != 0) {
    if(idim == 1) {
      Tensor<TensorType>::deallocate(dhi_trace);
      return;
    }

    N = idim + 1;
    std::vector<TensorType> AC(N * (N + 1) / 2);
    std::vector<TensorType> AF(AC.size());
    std::vector<TensorType> B(N, 0.);
    B.front() = -1.0;
    X.resize(N);

    int ac_i = 0;
    for(int i = 0; i <= idim; i++) {
      for(int j = i; j <= idim; j++) {
        AC[ac_i] = A(j, i);
        ac_i++;
      }
    }

    TensorType              RCOND;
    std::vector<int64_t>    IPIV(N);
    std::vector<TensorType> BERR(1), FERR(1); // NRHS
    info = lapack::spsvx(lapack::Factored::NotFactored, lapack::Uplo::Lower, N, 1, AC.data(),
                         AF.data(), IPIV.data(), B.data(), N, X.data(), N, &RCOND, FERR.data(),
                         BERR.data());

    if(info != 0) {
      if(rank == 0)
        cout
          << "<DIIS> Singularity in Pulay matrix detected." /*Error and Fock matrices removed." */
          << endl;
      Tensor<TensorType>::deallocate(diis_hist_alpha[0]);
      Tensor<TensorType>::deallocate(fock_hist_alpha[0]);
      diis_hist_alpha.erase(diis_hist_alpha.begin());
      fock_hist_alpha.erase(fock_hist_alpha.begin());

      if(is_uhf) {
        Tensor<TensorType>::deallocate(diis_hist_beta[0]);
        Tensor<TensorType>::deallocate(fock_hist_beta[0]);
        diis_hist_beta.erase(diis_hist_beta.begin());
        fock_hist_beta.erase(fock_hist_beta.begin());
      }

      idim--;
      if(idim == 1) {
        Tensor<TensorType>::deallocate(dhi_trace);
        return;
      }
      Matrix A_pl  = A.block(2, 2, idim, idim);
      Matrix A_new = Matrix::Zero(idim + 1, idim + 1);

      for(int i = 1; i <= idim; i++) {
        A_new(i, 0) = -1.0;
        A_new(0, i) = -1.0;
      }
      A_new.block(1, 1, idim, idim) = A_pl;
      A                             = A_new;
    }

  } // while

  Tensor<TensorType>::deallocate(dhi_trace); // dhi_trans

  std::vector<TensorType> X_final(N);
  // Reordering [0...N] -> [1...N,0] to match eigen's lu solve
  X_final.back() = X.front();
  std::copy(X.begin() + 1, X.end(), X_final.begin());

  // if(rank==0 && scf_options.debug)
  //  std::cout << "diis weights sum, vector: " << std::setprecision(5) <<
  //  std::accumulate(X.begin(),X.end(),0.0d) << std::endl << X << std::endl << X_final <<
  //  std::endl;

  sch(F_alpha() = 0);
  for(int j = 0; j < idim; j++) { sch(F_alpha() += X_final[j] * fock_hist_alpha[j]()); }
  if(is_uhf) {
    sch(F_beta() = 0);
    for(int j = 0; j < idim; j++) { sch(F_beta() += X_final[j] * fock_hist_beta[j]()); }
  }
  // sch(F() += 0.5 * D()); //level shift
  sch.execute();
}

template std::tuple<double, double> exachem::scf::SCFIter::scf_iter_body<double>(
  ExecutionContext& ec, ChemEnv& chem_env, ScalapackInfo& scalapack_info, const int& iter,
  SCFVars& scf_vars, TAMMTensors& ttensors, EigenTensors& etensors
#if defined(USE_GAUXC)
  ,
  GauXC::XCIntegrator<Matrix>& gauxc_integrator
#endif
);

template std::tuple<std::vector<int>, std::vector<int>, std::vector<int>>
exachem::scf::SCFIter::compute_2bf_taskinfo<double>(
  ExecutionContext& ec, ChemEnv& chem_env, const SCFVars& scf_vars, const bool do_schwarz_screen,
  const std::vector<size_t>& shell2bf, const Matrix& SchwarzK, const size_t& max_nprim4,
  TAMMTensors& ttensors, EigenTensors& etensors, const bool cs1s2);

template void exachem::scf::SCFIter::compute_2bf<double>(
  ExecutionContext& ec, ChemEnv& chem_env, ScalapackInfo& scalapack_info, const SCFVars& scf_vars,
  const bool do_schwarz_screen, const std::vector<size_t>& shell2bf, const Matrix& SchwarzK,
  const size_t& max_nprim4, TAMMTensors& ttensors, EigenTensors& etensors, bool& is_3c_init,
  const bool do_density_fitting, double xHF);

template void exachem::scf::SCFIter::compute_2bf_ri<double>(
  ExecutionContext& ec, ChemEnv& chem_env, ScalapackInfo& scalapack_info, const SCFVars& scf_vars,
  const std::vector<size_t>& shell2bf, TAMMTensors& ttensors, EigenTensors& etensors,
  bool& is_3c_init, double xHF);

template void exachem::scf::SCFIter::init_ri<double>(ExecutionContext& ec, ChemEnv& chem_env,
                                                     ScalapackInfo& scalapack_info,
                                                     const SCFVars& scf_vars,
                                                     EigenTensors& etensors, TAMMTensors& ttensors);

template void exachem::scf::SCFIter::compute_2c_ints<double>(ExecutionContext& ec,
                                                             ChemEnv&          chem_env,
                                                             EigenTensors&     etensors,
                                                             const SCFVars&    scf_vars,
                                                             TAMMTensors&      ttensors);

template void exachem::scf::SCFIter::compute_3c_ints<double>(ExecutionContext&   ec,
                                                             ChemEnv&            chem_env,
                                                             const SCFVars&      scf_vars,
                                                             Tensor<TensorType>& xyZ);

template void exachem::scf::SCFIter::compute_2bf_ri_direct<double>(
  ExecutionContext& ec, ChemEnv& chem_env, const SCFVars& scf_vars,
  const std::vector<size_t>& shell2bf, TAMMTensors& ttensors, EigenTensors& etensors,
  const Matrix& SchwarzK);
