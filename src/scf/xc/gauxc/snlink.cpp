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
#include "utilities.hpp"
#include <gauxc/xc_integrator.hpp>
#include <gauxc/xc_integrator/impl.hpp>
#include <gauxc/xc_integrator/integrator_factory.hpp>

namespace scf::xc::gauxc {

using k_type = simde::aos_k_e_aos;

// K Integration
MODULE_CTOR(snLinK) {
    satisfies_property_type<k_type>();

    add_submodule<XCQuadratureBatches>("Quadrature Batches")
      .set_description("Generate XC Quadrature Batches");

    /// TODO: Replace with information from the runtime
    add_input<bool>("On GPU").set_default(false);
}

MODULE_RUN(snLinK) {
    const auto&& [braket] = k_type::unwrap_inputs(inputs);
    const auto& bra       = braket.bra();
    const auto& P         = braket.op().rhs_particle();
    const auto& ket       = braket.ket();

    if(bra != ket || P.basis_set() != bra)
        throw std::runtime_error("GauXC: bra must be equal to ket");

    // Extract P into replicated Eigen matrix
    auto P_eigen = tw_to_eigen<double>(P.value());

    // Generate Partitioned Molecular Quadrature
    auto& lb_mod   = submods.at("Quadrature Batches");
    const auto& lb = lb_mod.run_as<XCQuadratureBatches>(bra.ao_basis_set());

    // Create XCIntegrator instance
    auto on_gpu   = inputs.at("On GPU").value<bool>();
    auto ex_space = (on_gpu) ?
                      GauXC::ExecutionSpace::Device :
                      GauXC::ExecutionSpace::Host; // TODO handle CUDA/HIP
    GauXC::functional_type func;                   // Dummy functional
    GauXC::XCIntegratorFactory<Eigen::MatrixXd> integrator_factory(
      ex_space, "Replicated", "Default", "Default", "Default");
    auto integrator = integrator_factory.get_instance(func, lb);

    // Do the K integration
    GauXC::IntegratorSettingsSNLinK sn_link_settings;
    auto K_eigen = integrator.eval_exx(P_eigen, sn_link_settings);
    auto K       = eigen_to_tw(K_eigen, get_runtime());

    // Wrap results
    auto rv = results();
    return k_type::wrap_results(rv, K);
}

} // namespace scf::xc::gauxc
