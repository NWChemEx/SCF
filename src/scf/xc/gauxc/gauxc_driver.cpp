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

using chemist::qm_operator::xc_functional;

#define GAUXC_NWX_XC_PAIR(FUNC) \
    if(in == xc_functional::FUNC) return ExchCXX::Functional::FUNC

ExchCXX::Functional nwx_to_exchcxx(xc_functional in) {
    GAUXC_NWX_XC_PAIR(SVWN3);
    else GAUXC_NWX_XC_PAIR(SVWN5);
    else GAUXC_NWX_XC_PAIR(BLYP);
    else GAUXC_NWX_XC_PAIR(B3LYP);
    else GAUXC_NWX_XC_PAIR(PBE);
    else GAUXC_NWX_XC_PAIR(revPBE);
    else GAUXC_NWX_XC_PAIR(PBE0);
    else throw std::out_of_range("No such functional");
}

#undef GAUXC_NWX_XC_PAIR

// XC Integration
MODULE_CTOR(GauXCDriver) {
    satisfies_property_type<XCDriver>();

    add_submodule<XCQuadratureBatches>("Quadrature Batches")
      .set_description("Generate XC Quadrature Batches");

    /// TODO: Replace with information from the runtime
    add_input<bool>("On GPU").set_default(false);
}

MODULE_RUN(GauXCDriver) {
    const auto& [nwx_func, aos, P] = XCDriver::unwrap_inputs(inputs);
    const size_t nbf               = aos.size();

    auto P_eigen = tw_to_eigen<double>(P);

    // Generate Partitioned Molecular Quadrature
    const auto& lb =
      submods.at("Quadrature Batches").run_as<XCQuadratureBatches>(aos);

    // Create XC functional
    auto func_spec = nwx_to_exchcxx(nwx_func);
    GauXC::functional_type func(ExchCXX::Backend::builtin, func_spec,
                                ExchCXX::Spin::Unpolarized);

    // Create XCIntegrator instance
    auto on_gpu = inputs.at("On GPU").value<bool>();
    auto ex_space =
      (on_gpu) ? GauXC::ExecutionSpace::Device : GauXC::ExecutionSpace::Host;
    GauXC::XCIntegratorFactory<Eigen::MatrixXd> integrator_factory(
      ex_space, "Replicated", "Default", "Default", "Default");
    auto integrator = integrator_factory.get_instance(func, lb);

    // Do the XC integration
    auto [EXC_float, VXC_eigen] = integrator.eval_exc_vxc(P_eigen);

    auto VXC = eigen_to_tw<double>(VXC_eigen, get_runtime());
    simde::type::tensor EXC(EXC_float);

    // Wrap results
    auto rv = results();
    return XCDriver::wrap_results(rv, EXC, VXC);
}

} // namespace scf::xc::gauxc
