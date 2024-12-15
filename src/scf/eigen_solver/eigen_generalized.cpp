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
    using matrix_alloc_t = tensorwrapper::allocator::Eigen<double, 2>;
    using vector_alloc_t = tensorwrapper::allocator::Eigen<double, 1>;
    const auto& eigen_A  = matrix_alloc_t::rebind(A.buffer());
    const auto& eigen_B  = matrix_alloc_t::rebind(B.buffer());

    // Wrap the tensors in Eigen::Map objects to avoid copy
    const auto* pA      = eigen_A.value().data();
    const auto* pB      = eigen_B.value().data();
    const auto& shape_A = eigen_A.layout().shape().as_smooth();
    auto rows           = shape_A.extent(0);
    auto cols           = shape_A.extent(1);
    Eigen::Map<const Eigen::MatrixXd> A_map(pA, rows, cols);
    Eigen::Map<const Eigen::MatrixXd> B_map(pB, rows, cols);

    // Compute
    Eigen::GeneralizedEigenSolver<Eigen::MatrixXd> ges(A_map, B_map);
    auto eigen_values  = ges.eigenvalues();
    auto eigen_vectors = ges.eigenvectors();

    // Wrap in TensorWrapper Tensor
    tensorwrapper::shape::Smooth vector_shape{rows};
    tensorwrapper::shape::Smooth matrix_shape{rows, cols};
    tensorwrapper::layout::Physical vector_layout(vector_shape);
    tensorwrapper::layout::Physical matrix_layout(matrix_shape);

    using matrix_buffer = typename matrix_alloc_t::eigen_buffer_type;
    using vector_buffer = typename vector_alloc_t::eigen_buffer_type;

    typename vector_buffer::data_type vector_tensor(rows);
    typename matrix_buffer::data_type matrix_tensor(rows, cols);
    for(auto i = 0; i < rows; ++i) {
        vector_tensor(i) = eigen_values(i).real();
        for(auto j = 0; j < cols; ++j) {
            matrix_tensor(i, j) = eigen_vectors(i, j).real();
        }
    }

    vector_buffer values_buffer(vector_tensor, vector_layout);
    matrix_buffer vectors_buffer(matrix_tensor, matrix_layout);

    simde::type::tensor values(vector_shape, values_buffer);
    simde::type::tensor vectors(matrix_shape, vectors_buffer);

    auto rv = results();
    return pt::wrap_results(rv, values, vectors);
}

} // namespace scf::eigen_solver