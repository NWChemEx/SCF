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

#include "eigen_solver.hpp"
#include <Eigen/Eigen>
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

    template<typename FloatType0, typename FloatType1>
    return_t operator()(const std::span<FloatType0>& A,
                        const std::span<FloatType1>& B) {
        throw std::runtime_error(
          "EigenGeneralized Kernel: Mixed float types not supported");
    }

    template<typename FloatType>
    return_t operator()(const std::span<FloatType>& A,
                        const std::span<FloatType>& B) {
        using clean_t = std::decay_t<FloatType>;
        // Convert to Eigen buffers

        // Wrap the tensors in Eigen::Map objects to avoid copy
        const auto* pA = A.data();
        const auto* pB = B.data();
        auto rows      = m_n_rows;
        auto cols      = m_n_cols;

        constexpr auto rmajor = Eigen::RowMajor;
        constexpr auto edynam = Eigen::Dynamic;
        using clean_type      = std::decay_t<FloatType>;
        using matrix_type = Eigen::Matrix<clean_type, edynam, edynam, rmajor>;
        using map_type    = Eigen::Map<const matrix_type>;

        map_type A_map(pA, rows, cols);
        map_type B_map(pB, rows, cols);

        // Compute
        Eigen::GeneralizedSelfAdjointEigenSolver<matrix_type> ges(A_map, B_map);
        auto eigen_values  = ges.eigenvalues();
        auto eigen_vectors = ges.eigenvectors();

        // Wrap in TensorWrapper Tensor
        tensorwrapper::shape::Smooth vector_shape{rows};
        tensorwrapper::shape::Smooth matrix_shape{rows, cols};
        tensorwrapper::layout::Physical vector_layout(vector_shape);
        tensorwrapper::layout::Physical matrix_layout(matrix_shape);

        using tensorwrapper::buffer::make_contiguous;

        auto pvalues_buffer  = make_contiguous<clean_t>(vector_shape);
        auto pvectors_buffer = make_contiguous<clean_t>(matrix_shape);

        for(decltype(rows) i = 0; i < rows; ++i) {
            pvalues_buffer.set_elem({i}, eigen_values(i));
            for(decltype(cols) j = 0; j < cols; ++j) {
                pvectors_buffer.set_elem({i, j}, eigen_vectors(i, j));
            }
        }

        simde::type::tensor values(vector_shape, std::move(pvalues_buffer));
        simde::type::tensor vectors(matrix_shape, std::move(pvectors_buffer));
        return std::make_pair(values, vectors);
    }
};
} // namespace

using pt = simde::GeneralizedEigenSolve;

const auto desc = R"(
Generalized Eigen Solve via Eigen
---------------------------------

TODO: Write me!!!
)";

MODULE_CTOR(EigenGeneralized) {
    description(desc);
    satisfies_property_type<pt>();
}

MODULE_RUN(EigenGeneralized) {
    auto&& [A, B] = pt::unwrap_inputs(inputs);

    using tensorwrapper::buffer::make_contiguous;
    const auto& A_buffer = make_contiguous(A.buffer());
    const auto& B_buffer = make_contiguous(B.buffer());
    const auto& A_shape  = A_buffer.shape();
    Kernel k(A_shape.extent(0), A_shape.extent(1));
    using tensorwrapper::buffer::visit_contiguous_buffer;
    auto [values, vectors] = visit_contiguous_buffer(k, A_buffer, B_buffer);

    auto rv = results();
    return pt::wrap_results(rv, values, vectors);
}

} // namespace scf::eigen_solver
