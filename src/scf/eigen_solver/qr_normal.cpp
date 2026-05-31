/*
 * Copyright 2026 NWChemEx-Project
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

#include "eigen_solver.hpp"
#include "qr_eigen_helpers.hpp"
#include <limits>
#include <simde/simde.hpp>

namespace scf::eigen_solver {

namespace {
template<typename FloatType>
struct Tolerance {
    static constexpr double value =
      10 * std::numeric_limits<FloatType>::epsilon();
};
#ifdef ENABLE_SIGMA
template<typename FloatType>
struct Tolerance<tensorwrapper::types::uncertain_type<FloatType>> {
    using float_type = tensorwrapper::types::uncertain_type<FloatType>;
    using value_t    = typename float_type::value_t;
    static constexpr double value =
      10 * std::numeric_limits<value_t>::epsilon();
};

template<typename FloatType>
struct Tolerance<tensorwrapper::types::interval_type<FloatType>> {
    using value_t =
      typename tensorwrapper::types::interval_type<FloatType>::value_t;
    static constexpr double value =
      10 * std::numeric_limits<value_t>::epsilon();
};
#endif

struct Kernel {
    std::size_t m_n_rows;
    std::size_t m_n_cols;

    using tensor_t = simde::type::tensor;
    using return_t = std::pair<tensor_t, tensor_t>;

    Kernel(std::size_t n_rows, std::size_t n_cols) :
      m_n_rows(n_rows), m_n_cols(n_cols) {}

    template<typename FloatType>
    return_t operator()(const std::span<FloatType>& A) {
        using clean_t = std::decay_t<FloatType>;
        const auto n  = m_n_rows;
        if(m_n_rows != m_n_cols) {
            throw std::runtime_error("QRNormal: matrix must be square");
        }
        const double tol = Tolerance<clean_t>::value;
        auto [evals, evecs] =
          detail::symmetric_qr_eigen<clean_t>(A, n, tol, 1000);

        using tensorwrapper::utilities::make_tensor;
        auto values  = make_tensor({n}, evals);
        auto vectors = make_tensor({n, n}, evecs);
        return std::make_pair(values, vectors);
    }
};
} // namespace

using pt = simde::EigenSolve;

const auto desc = R"(
 Eigen Solve via QR
 ---------------------------------

 Symmetric eigen solve by QR iteration using Householder QR
 decomposition at each step.
 )";

MODULE_CTOR(QRNormal) {
    description(desc);
    satisfies_property_type<pt>();
}

MODULE_RUN(QRNormal) {
    auto&& [A] = pt::unwrap_inputs(inputs);

    using tensorwrapper::buffer::make_contiguous;
    const auto& A_buffer = make_contiguous(A.buffer());
    const auto& A_shape  = A_buffer.shape();
    Kernel k(A_shape.extent(0), A_shape.extent(1));
    using tensorwrapper::buffer::visit_contiguous_buffer;
    auto [values, vectors] = visit_contiguous_buffer(k, A_buffer);

    auto rv = results();
    return pt::wrap_results(rv, values, vectors);
}

} // namespace scf::eigen_solver
