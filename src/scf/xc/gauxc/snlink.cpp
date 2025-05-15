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
// #include <tensorwrapper/tensor/backends/tiledarray.hpp>

// namespace scf::xc::gauxc {

// using k_type = simde::aos_k_e_aos;

// // K Integration
// MODULE_CTOR(snLinK) {
//     satisfies_property_type<k_type>();

//     add_submodule<XCQuadratureBatches>("Quadrature Batches")
//       .set_description("Generate XC Quadrature Batches");

//     add_input<std::string>("Integration Grid")
//       .set_description("Specification of the Atomic Grid")
//       .set_default("UltraFine");

//     /// TODO: Replace with information from the runtime
//     add_input<bool>("On GPU").set_default(false);
// }

// MODULE_RUN(snLinK) {
//     const auto&& [braket] = k_type::unwrap_inputs(inputs);
//     const auto& bra       = braket.bra().ao_basis_set();
//     const auto& k_e       = braket.op();
//     const auto& ket       = braket.ket().ao_basis_set();

//     simde::type::tensor K;
//     // Early Return
//     if(k == simde::type::el_scf_k{}) {
//         auto rv = results();
//         return k_type::wrap_results(rv, K);
//     }

//     if(bra != ket) throw std::runtime_error("GauXC: bra must be equal to
//     ket");

//     // Extract P into replicated Eigen matrix
//     const auto& P       = k.at<1>().value();
//     const auto& P_basis = k.at<1>().basis_set();
//     if(bra != P_basis)
//         throw std::runtime_error("GauXC: inconsistent basis P/K");

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
//     auto grid_string = inputs.at("Integration Grid").value<std::string>();
//     auto grid_spec   = atomic_grid_types.at(grid_string);
//     auto& lb_mod     = submods.at("Quadrature Batches");
//     const auto& lb   = lb_mod.run_as<XCQuadratureBatches>(
//       bra.basis_set(), GauXC::PruningScheme::Unpruned,
//       GauXC::RadialQuad::MuraKnowles, grid_spec);

//     // Create XCIntegrator instance
//     auto on_gpu   = inputs.at("On GPU").value<bool>();
//     auto ex_space = (on_gpu) ?
//                       GauXC::ExecutionSpace::Device :
//                       GauXC::ExecutionSpace::Host; // TODO handle CUDA/HIP
//     GauXC::functional_type func;                   // Dummy functional
//     GauXC::XCIntegratorFactory<Eigen::MatrixXd> integrator_factory(
//       ex_space, "Replicated", "Default", "Default", "Default");
//     auto integrator = integrator_factory.get_instance(func, lb);

//     // Do the K integration
//     GauXC::IntegratorSettingsSNLinK sn_link_settings;
//     auto K_eigen = integrator.eval_exx(P_eigen, sn_link_settings);
//     K            = tensorwrapper::tensor::eigen_to_tensor_wrapper(K_eigen);

//     /// TODO: shape needs to be handle by eigen_to_tensor_wrapper
//     auto& K_ta = tensorwrapper::tensor::backends::unwrap_ta(K);
//     auto& P_ta = tensorwrapper::tensor::backends::unwrap_ta(P);
//     K_ta       = TA::retile(K_ta, P_ta.trange());

//     // Wrap results
//     auto rv = results();
//     return k_type::wrap_results(rv, K);
// }

// } // namespace scf::xc::gauxc
