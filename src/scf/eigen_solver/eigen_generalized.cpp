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

    // Convert to Eigen buffers
    tensorwrapper::allocator::Eigen<double> allocator(get_runtime());
    const auto& eigen_A = allocator.rebind(A.buffer());
    const auto& eigen_B = allocator.rebind(B.buffer());

    // Wrap the tensors in Eigen::Map objects to avoid copy
    const auto* pA      = eigen_A.data();
    const auto* pB      = eigen_B.data();
    const auto& shape_A = eigen_A.layout().shape().as_smooth();
    auto rows           = shape_A.extent(0);
    auto cols           = shape_A.extent(1);
    Eigen::Map<const Eigen::MatrixXd> A_map(pA, rows, cols);
    Eigen::Map<const Eigen::MatrixXd> B_map(pB, rows, cols);

    // Compute
    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> ges(A_map, B_map);
    auto eigen_values  = ges.eigenvalues();
    auto eigen_vectors = ges.eigenvectors();

    // Wrap in TensorWrapper Tensor
    tensorwrapper::shape::Smooth vector_shape{rows};
    tensorwrapper::shape::Smooth matrix_shape{rows, cols};
    tensorwrapper::layout::Physical vector_layout(vector_shape);
    tensorwrapper::layout::Physical matrix_layout(matrix_shape);

    auto pvalues_buffer  = allocator.allocate(vector_layout);
    auto pvectors_buffer = allocator.allocate(matrix_layout);

    for(auto i = 0; i < rows; ++i) {
        pvalues_buffer->at(i) = eigen_values(i);
        for(auto j = 0; j < cols; ++j) {
            pvectors_buffer->at(i, j) = eigen_vectors(i, j);
        }
    }

    simde::type::tensor values(vector_shape, std::move(pvalues_buffer));
    simde::type::tensor vectors(matrix_shape, std::move(pvectors_buffer));

    auto rv = results();
    return pt::wrap_results(rv, values, vectors);
}

} // namespace scf::eigen_solver