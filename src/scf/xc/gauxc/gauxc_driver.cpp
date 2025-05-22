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
#include <gauxc/xc_integrator.hpp>
#include <gauxc/xc_integrator/impl.hpp>
#include <gauxc/xc_integrator/integrator_factory.hpp>

namespace scf::xc::gauxc {

/**
 *  All default grids are product quadratures consisting of NR Mura-Knowles
 (MK)
 *  radial quadratures and NA Lebedev (L) spherical quadratures.
 *
 *  Fine(Grid)
 *    - 75 (MK) x 302 (L) grid for each atom
 *    - Target Accuracy 1e-6
 *  UltraFine(Grid)
 *    - 99 (MK) x 590 (L) grid for each atom
 *    - Target Accuracy 1e-8
 *  SuperFine(Grid)
 *    - 175 (MK) x 974 (L) grid for atoms with Z <= 2
 *    - 250 (MK) x 974 (L) grid for atoms with Z >= 3
 *    - Target Accuracy 1e-10
 */
static utilities::CaseInsensitiveMap<xc_grid_type> atomic_grid_types = {
  {"Fine", xc_grid_type::FineGrid},
  {"UltraFine", xc_grid_type::UltraFineGrid},
  {"SuperFine", xc_grid_type::SuperFineGrid},
  {"FineGrid", xc_grid_type::FineGrid},
  {"UltraFineGrid", xc_grid_type::UltraFineGrid},
  {"SuperFineGrid", xc_grid_type::SuperFineGrid}};

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

    add_input<std::string>("XC Grid")
      .set_description("Specification of the Atomic Grid")
      .set_default("UltraFine");

    /// TODO: Replace with information from the runtime
    add_input<bool>("On GPU").set_default(false);
}

MODULE_RUN(GauXCDriver) {
    const auto& [nwx_func, aos, P] = XCDriver::unwrap_inputs(inputs);

    tensorwrapper::allocator::Eigen<double> alloc(get_runtime());
    const size_t nbf = aos.size();
    Eigen::MatrixXd P_eigen(nbf, nbf);
    {
        auto P_vector = alloc.rebind(P.buffer());
        auto begin    = P_vector.get_immutable_data();
        auto end      = begin + (nbf * nbf);
        std::copy(begin, end, P_eigen.data());
    }

    // Generate Partitioned Molecular Quadrature
    auto grid_string = inputs.at("XC Grid").value<std::string>();
    auto grid_spec   = atomic_grid_types.at(grid_string);
    auto& lb_mod     = submods.at("Quadrature Batches");
    const auto& lb   = lb_mod.run_as<XCQuadratureBatches>(
      aos, GauXC::PruningScheme::Unpruned, GauXC::RadialQuad::MuraKnowles,
      grid_spec);

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

    tensorwrapper::shape::Smooth shape{nbf, nbf};
    tensorwrapper::layout::Physical layout(shape);
    auto pbuffer = alloc.allocate(layout);
    for(std::size_t i = 0; i < nbf; ++i) {
        for(std::size_t j = 0; j < nbf; ++j)
            pbuffer->set_elem({i, j}, VXC_eigen(i, j));
    }
    simde::type::tensor VXC(shape, std::move(pbuffer));
    simde::type::tensor EXC(EXC_float);

    // Wrap results
    auto rv = results();
    return XCDriver::wrap_results(rv, EXC, VXC);
}

} // namespace scf::xc::gauxc
