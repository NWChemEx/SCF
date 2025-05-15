// /*
//  * Copyright 2022 NWChemEx-Project
//  *
//  * Licensed under the Apache License, Version 2.0 (the "License");
//  * you may not use this file except in compliance with the License.
//  * You may obtain a copy of the License at
//  *
//  * http://www.apache.org/licenses/LICENSE-2.0
//  *
//  * Unless required by applicable law or agreed to in writing, software
//  * distributed under the License is distributed on an "AS IS" BASIS,
//  * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  * See the License for the specific language governing permissions and
//  * limitations under the License.
//  */

// #include "gauxc.hpp"
// #include <gauxc/xc_integrator.hpp>
// #include <gauxc/xc_integrator/impl.hpp>
// #include <gauxc/xc_integrator/integrator_factory.hpp>

// namespace scf::xc::gauxc {

// /**
//  *  All default grids are product quadratures consisting of NR Mura-Knowles
//  (MK)
//  *  radial quadratures and NA Lebedev (L) spherical quadratures.
//  *
//  *  Fine(Grid)
//  *    - 75 (MK) x 302 (L) grid for each atom
//  *    - Target Accuracy 1e-6
//  *  UltraFine(Grid)
//  *    - 99 (MK) x 590 (L) grid for each atom
//  *    - Target Accuracy 1e-8
//  *  SuperFine(Grid)
//  *    - 175 (MK) x 974 (L) grid for atoms with Z <= 2
//  *    - 250 (MK) x 974 (L) grid for atoms with Z >= 3
//  *    - Target Accuracy 1e-10
//  */
// static utilities::CaseInsensitiveMap<xc_grid_type> atomic_grid_types = {
//   {"Fine", xc_grid_type::FineGrid},
//   {"UltraFine", xc_grid_type::UltraFineGrid},
//   {"SuperFine", xc_grid_type::SuperFineGrid},
//   {"FineGrid", xc_grid_type::FineGrid},
//   {"UltraFineGrid", xc_grid_type::UltraFineGrid},
//   {"SuperFineGrid", xc_grid_type::SuperFineGrid}};

// #define GAUXC_NWX_XC_PAIR(FUNC) {#FUNC, ExchCXX::Functional::FUNC}

// std::map<std::string, ExchCXX::Functional> nwx_to_exchcxx = {
//   GAUXC_NWX_XC_PAIR(SVWN3), GAUXC_NWX_XC_PAIR(SVWN5),
//   GAUXC_NWX_XC_PAIR(BLYP), GAUXC_NWX_XC_PAIR(B3LYP), GAUXC_NWX_XC_PAIR(PBE),
//   GAUXC_NWX_XC_PAIR(revPBE), GAUXC_NWX_XC_PAIR(PBE0)};

// #undef GAUXC_NWX_XC_PAIR

// using vxc_type = simde::aos_xc_e_aos;

// // XC Integration
// MODULE_CTOR(GauXC_XC) {
//     satisfies_property_type<exc_type>();

//     add_submodule<XCQuadratureBatches>("Quadrature Batches")
//       .set_description("Generate XC Quadrature Batches");

//     add_submodule<integral_shape_pt>("Tensor Shape")
//       .set_description("Determines the shape of the resulting tensor");

//     add_input<std::string>("XC Grid")
//       .set_description("Specification of the Atomic Grid")
//       .set_default("UltraFine");

//     /// TODO: Replace with information from the runtime
//     add_input<bool>("On GPU").set_default(false);
// }

// MODULE_RUN(GauXC_XC) {
//     auto [braket] = exc_type::unwrap_inputs(inputs);

//     auto bra = braket;
//     auto op  = braket;
//     auto ket = braket;

//     // Extract P into replicated Eigen matrix
//     const auto& P       = xc.at<2>().value();
//     const auto& P_basis = xc.at<2>().basis_set();
//     if(bra != P_basis)
//         throw std::runtime_error("GauXC: inconsistent basis P/XC");

//     const size_t nbf = bra.basis_set().n_aos();
//     Eigen::MatrixXd P_eigen(nbf, nbf);
//     {
//         // Scope vector
//         auto P_vector = tensorwrapper::tensor::to_vector(P);
//         if(P_vector.size() != (nbf * nbf)) {
//             throw std::runtime_error("Vectorized Density Has Wrong Size");
//         }
//         std::copy(P_vector.begin(), P_vector.end(), P_eigen.data());
//     }

//     // Generate Partitioned Molecular Quadrature
//     auto grid_string = inputs.at("XC Grid").value<std::string>();
//     auto grid_spec   = atomic_grid_types.at(grid_string);
//     auto& lb_mod     = submods.at("Quadrature Batches");
//     const auto& lb   = lb_mod.run_as<XCQuadratureBatches>(
//       bra.basis_set(), GauXC::PruningScheme::Unpruned,
//       GauXC::RadialQuad::MuraKnowles, grid_spec);

//     // Create XC functional
//     auto func_spec = nwx_to_exchcxx.at(xc.at<0>());
//     GauXC::functional_type func(ExchCXX::Backend::builtin, func_spec,
//                                 ExchCXX::Spin::Unpolarized);

//     // Create XCIntegrator instance
//     auto on_gpu = inputs.at("On GPU").value<bool>();
//     auto ex_space =
//       (on_gpu) ? GauXC::ExecutionSpace::Device : GauXC::ExecutionSpace::Host;
//     GauXC::XCIntegratorFactory<Eigen::MatrixXd> integrator_factory(
//       ex_space, "Replicated", "Default", "Default", "Default");
//     auto integrator = integrator_factory.get_instance(func, lb);

//     // Do the XC integration
//     std::vector<simde::type::ao_basis_set> bases = {bra.basis_set(),
//                                                     ket.basis_set()};
//     auto shape = submods.at("Tensor Shape").run_as<integral_shape_pt>(bases);
//     auto [EXC, VXC_eigen] = integrator.eval_exc_vxc(P_eigen);
//     auto VXC = tensorwrapper::tensor::eigen_to_tensor_wrapper(VXC_eigen,
//     shape);

//     // Wrap results
//     auto rv = results();
//     rv      = exc_type::wrap_results(rv, EXC);
//     return vxc_type::wrap_results(rv, VXC);
// }

// } // namespace scf::xc::gauxc
