/*
 * ExaChem: Open Source Exascale Computational Chemistry Software.
 *
 * Copyright 2024 Pacific Northwest National Laboratory, Battelle Memorial Institute.
 *
 * See LICENSE.txt for details
 */

#include "input_options.hpp"

void SCFOptions::print() {
  std::cout << std::defaultfloat;
  std::cout << std::endl << "SCF Options" << std::endl;
  std::cout << "{" << std::endl;
  std::cout << " charge            = " << charge << std::endl;
  std::cout << " multiplicity      = " << multiplicity << std::endl;
  std::cout << " level shift       = " << lshift << std::endl;
  std::cout << " tol_int           = " << tol_int << std::endl;
  std::cout << " tol_sch           = " << tol_sch << std::endl;
  std::cout << " tol_lindep        = " << tol_lindep << std::endl;
  std::cout << " conve             = " << conve << std::endl;
  std::cout << " convd             = " << convd << std::endl;
  std::cout << " diis_hist         = " << diis_hist << std::endl;
  std::cout << " AO_tilesize       = " << AO_tilesize << std::endl;
  std::cout << " writem            = " << writem << std::endl;
  std::cout << " damp              = " << damp << std::endl;
  std::cout << " n_lindep          = " << n_lindep << std::endl;
  if(!moldenfile.empty()) {
    std::cout << " moldenfile        = " << moldenfile << std::endl;
    // std::cout << " n_lindep = " << n_lindep <<  std::endl;
  }

  std::cout << " scf_type          = " << scf_type << std::endl;

  txt_utils::print_bool(" direct_df        ", direct_df);

  if(!xc_type.empty() || snK) {
    std::cout << " DFT " << std::endl << " {" << std::endl;
    std::cout << "  xc_type           = [ ";
    for(auto xcfunc: xc_type) { std::cout << " \"" << xcfunc << "\","; }
    std::cout << "\b ]" << std::endl;

    std::cout << "  xc_grid_type      = " << xc_grid_type << std::endl;
    std::cout << "  xc_pruning_scheme = " << xc_pruning_scheme << std::endl;
    std::cout << "  xc_weight_scheme  = " << xc_weight_scheme << std::endl;
    std::cout << "  xc_exec_space     = " << xc_exec_space << std::endl;
    std::cout << "  xc_basis_tol      = " << xc_basis_tol << std::endl;
    std::cout << "  xc_batch_size     = " << xc_batch_size << std::endl;
    std::cout << "  xc_lb_kernel      = " << xc_lb_kernel << std::endl;
    std::cout << "  xc_mw_kernel      = " << xc_mw_kernel << std::endl;
    std::cout << "  xc_int_kernel     = " << xc_int_kernel << std::endl;
    std::cout << "  xc_red_kernel     = " << xc_red_kernel << std::endl;
    std::cout << "  xc_lwd_kernel     = " << xc_lwd_kernel << std::endl;
    if(xc_radang_size.first > 0 && xc_radang_size.second > 0) {
      std::cout << "  xc_radang_size    = " << xc_radang_size.first << ", " << xc_radang_size.second
                << std::endl;
    }
    else { std::cout << "  xc_rad_quad       = " << xc_rad_quad << std::endl; }

    txt_utils::print_bool("  snK              ", snK);
    std::cout << "  xc_snK_etol       = " << xc_snK_etol << std::endl;
    std::cout << "  xc_snK_ktol       = " << xc_snK_ktol << std::endl;
    std::cout << " }" << std::endl;
  }

  if(scalapack_np_row > 0 && scalapack_np_col > 0) {
    std::cout << " scalapack_np_row  = " << scalapack_np_row << std::endl;
    std::cout << " scalapack_np_col  = " << scalapack_np_col << std::endl;
    if(scalapack_nb > 1) std::cout << " scalapack_nb      = " << scalapack_nb << std::endl;
  }
  std::cout << " restart_size      = " << restart_size << std::endl;
  txt_utils::print_bool(" restart          ", restart);
  txt_utils::print_bool(" debug            ", debug);
  if(restart) txt_utils::print_bool(" noscf            ", noscf);
  // txt_utils::print_bool(" sad         ", sad);
  if(mulliken_analysis || mos_txt || mo_vectors_analysis.first) {
    std::cout << " PRINT {" << std::endl;
    if(mos_txt) std::cout << std::boolalpha << "  mos_txt             = " << mos_txt << std::endl;
    if(mulliken_analysis)
      std::cout << std::boolalpha << "  mulliken_analysis   = " << mulliken_analysis << std::endl;
    if(mo_vectors_analysis.first) {
      std::cout << "  mo_vectors_analysis = [" << std::boolalpha << mo_vectors_analysis.first;
      std::cout << "," << mo_vectors_analysis.second << "]" << std::endl;
    }
    std::cout << " }" << std::endl;
  }
  std::cout << "}" << std::endl << std::flush;
} // END of SCFOptions::print

void CommonOptions::print() {
  std::cout << std::defaultfloat;
  std::cout << std::endl << "Common Options" << std::endl;
  std::cout << "{" << std::endl;
  std::cout << " maxiter    = " << maxiter << std::endl;
  std::cout << " basis      = " << basis << " ";
  std::cout << gaussian_type;
  std::cout << std::endl;
  if(!dfbasis.empty()) std::cout << " dfbasis    = " << dfbasis << std::endl;
  if(!basisfile.empty()) std::cout << " basisfile  = " << basisfile << std::endl;
  std::cout << " geom_units = " << geom_units << std::endl;
  txt_utils::print_bool(" debug     ", debug);
  if(!file_prefix.empty()) std::cout << " file_prefix    = " << file_prefix << std::endl;
  std::cout << "}" << std::endl;
}
