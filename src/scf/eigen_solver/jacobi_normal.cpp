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
#include "jacobi_eigen_helpers.hpp"
#include <limits>
#include <simde/simde.hpp>

namespace scf::eigen_solver {

namespace {

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
            throw std::runtime_error("JacobiNormal: matrix must be square");
        }
        double tol = []() {
            if constexpr(tensorwrapper::types::is_uq_type_v<clean_t>) {
                using value_t = typename clean_t::value_t;
                return value_t(10) * std::numeric_limits<value_t>::epsilon();
            } else {
                return 10 * std::numeric_limits<clean_t>::epsilon();
            }
        }();

        const auto max_sweeps =
          tensorwrapper::types::is_uq_type_v<clean_t> ? 2000 * n : 50 * n;
        auto [evals, evecs] =
          detail::symmetric_jacobi_eigen<clean_t>(A, n, tol, max_sweeps);

        using tensorwrapper::utilities::make_tensor;
        auto values  = make_tensor({n}, evals);
        auto vectors = make_tensor({n, n}, evecs);
        return std::make_pair(values, vectors);
    }
};
} // namespace

using pt = simde::EigenSolve;

const auto desc = R"(
 Eigen Solve via Jacobi
 ---------------------------------

 Symmetric eigen solve by cyclic Jacobi rotations using Eigen's
 JacobiRotation class.
 )";

MODULE_CTOR(JacobiNormal) {
    description(desc);
    satisfies_property_type<pt>();
}

MODULE_RUN(JacobiNormal) {
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
