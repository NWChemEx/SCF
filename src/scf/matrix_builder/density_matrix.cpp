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
#include <unsupported/Eigen/CXX11/Tensor>

namespace scf::matrix_builder {
namespace {

const auto desc = R"(
)";

struct Kernel {
    Kernel(std::size_t n_aos, std::size_t n_occ) :
      m_n_aos(n_aos), m_n_occ(n_occ) {}
    std::size_t m_n_aos;
    std::size_t m_n_occ;

    template<typename FloatType>
    auto operator()(const std::span<FloatType>& c) {
        constexpr auto rmajor = Eigen::RowMajor;
        constexpr auto edynam = Eigen::Dynamic;
        using clean_type      = std::decay_t<FloatType>;
        using tensor_type = Eigen::Matrix<clean_type, edynam, edynam, rmajor>;
        using const_map_type = Eigen::Map<const tensor_type>;

        tensorwrapper::shape::Smooth p_shape{m_n_aos, m_n_aos};
        auto pp_buffer =
          tensorwrapper::buffer::make_contiguous<clean_type>(p_shape);

        // Step 2: Grab the orbitals in the ensemble
        const_map_type c_eigen(c.data(), m_n_aos, m_n_aos);
        auto slice = c_eigen.block(0, 0, m_n_aos, m_n_occ);

        // Step 3: CC_dagger
        tensor_type p_eigen = slice * slice.transpose();
        for(std::size_t i = 0; i < m_n_aos; ++i) {
            for(std::size_t j = 0; j < m_n_aos; ++j) {
                pp_buffer.set_elem({i, j}, p_eigen(i, j));
            }
        }
        return simde::type::tensor(p_shape, std::move(pp_buffer));
    }
};

} // namespace

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
    // TODO: Compare the aos and make sure they're all the same
    // const auto& bra_aos   = braket.bra();
    // const auto& ket_aos   = braket.ket();
    const auto& op     = braket.op();
    const auto& cutoff = inputs.at("cutoff").value<double>();

    const auto& mos     = op.orbitals();
    const auto& c       = mos.transform();
    const auto& weights = op.weights();
    auto n_aos          = c.logical_layout().shape().as_smooth().extent(0);

    // Step 1: Figure out which orbitals we need to grab
    using size_type = std::size_t;
    std::vector<size_type> participants;
    for(size_type i = 0; i < weights.size(); ++i)
        if(std::fabs(weights[i]) >= cutoff) participants.push_back(i);

    // For now we require that participants == [0, participants.size())
    for(size_type i = 0; i < participants.size(); ++i)
        if(participants[i] != i)
            throw std::runtime_error(
              "Please shuffle your orbitals so that the ensemble is a "
              "contiguous slice of the orbitals.");

    // TODO: The need to dispatch like this goes away when TW supports slicing
    using tensorwrapper::buffer::visit_contiguous_buffer;
    Kernel k(n_aos, participants.size());
    const auto& c_buffer = tensorwrapper::buffer::make_contiguous(c.buffer());
    auto p               = visit_contiguous_buffer(k, c_buffer);
    auto rv              = results();
    return pt::wrap_results(rv, p);
}

} // namespace scf::matrix_builder
