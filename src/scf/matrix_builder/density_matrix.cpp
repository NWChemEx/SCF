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

#include "matrix_builder.hpp"

namespace scf::matrix_builder {
namespace {

const auto desc = R"(
)";

}

using pt = simde::aos_rho_e_aos<simde::type::cmos>;

MODULE_CTOR(DensityMatrix) {
    description(desc);
    satisfies_property_type<pt>();
    add_input<double>("cutoff")
      .set_description(
        "The cutoff for considering a state as part of the ensemble.")
      .set_default(1E-16);
}

MODULE_RUN(DensityMatrix) {
    const auto&& [braket] = pt::unwrap_inputs(inputs);
    const auto& bra_aos   = braket.bra();
    const auto& op        = braket.op();
    const auto& ket_aos   = braket.ket();
    const auto& cutoff    = inputs.at("cutoff").value<double>();

    const auto& mos     = op.orbitals();
    const auto& c       = mos.transform();
    const auto& weights = op.weights();
    auto n_aos          = c.logical_layout().shape().as_smooth().extent(0);

    // TODO: Compare the aos and make sure they're all the same

    // Step 1: Figure out which orbitals we need to grab
    using size_type = std::size_t;
    std::vector<size_type> participants;
    for(size_type i = 0; i < weights.size(); ++i)
        if(std::fabs(weights[i]) >= cutoff) participants.push_back(i);

    // For now we require that participants == [0, participants.size())
    for(size_type i = 0; i < participants.size(); ++i)
        if(participants[i] != i)
            throw std::runtime_error(
              "Pleas shuffle your orbitals so that the ensemble is a "
              "contiguous slice of the orbitals.");

    // Step 2: Grab the orbitals in the ensemble
    using allocator_type    = tensorwrapper::allocator::Eigen<double, 2>;
    using buffer_type       = typename allocator_type::eigen_buffer_type;
    using eigen_tensor_type = typename buffer_type::data_type;

    const auto& c_buffer = allocator_type::rebind(c.buffer());
    const auto& c_eigen  = c_buffer.value();

    using slice_t = Eigen::array<Eigen::Index, 2>;
    slice_t offsets{0, 0};
    slice_t extents{n_aos, Eigen::Index(participants.size())};
    eigen_tensor_type slice = c_eigen.slice(offsets, extents);

    // Step 3: CC_dagger
    using index_pair_t = Eigen::IndexPair<int>;
    Eigen::array<index_pair_t, 1> modes{index_pair_t(1, 1)};
    eigen_tensor_type p_eigen = slice.contract(slice, modes);
    tensorwrapper::shape::Smooth p_shape(n_aos, n_aos);
    tensorwrapper::layout::Physical l(p_shape);
    buffer_type p_buffer(p_eigen, l);
    simde::type::tensor p(p_shape, p_buffer);

    auto rv = results();
    return pt::wrap_results(rv, p);
}

} // namespace scf::matrix_builder