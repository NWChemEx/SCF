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
#include "submodule_request.hpp"
#include "wtf/fp/float_view.hpp"
#include <sigma/sigma.hpp>
#include <simde/simde.hpp>
#include <tensorwrapper/tensorwrapper.hpp>

using pt = simde::GeneralizedEigenSolve;

namespace scf::eigen_solver {
namespace {
const auto desc = R"(
Generalized Eigen Solve via Ball Arithmetic
-------------------------------------------
https://www.texmacs.org/joris/ball/ball.html
TODO: Write me!!!
H = TMT - L
n = number of rows (or columns)
k = max[max(1/(L_jj -Lii)), max(1/L_ii)]
Omega_n = n by n matrix 0+/-1
eta = 6 sqrt(n) k ||H||

T_ball = T(1 + eta Omega_n)
L_ball = B(1, eta) * L

For the generalized problem, AC=BCL, we'll try:

H = TAC - TBCL

where T = C^T
)";

template<typename UQType>
auto wrap_subdiagonalization(const tensorwrapper::buffer::Contiguous& A_buf,
                             const tensorwrapper::buffer::Contiguous& B_buf,
                             pluginplay::SubmoduleRequest& eigen_solver_mod) {
    using namespace tensorwrapper::buffer;
    using namespace tensorwrapper::shape;

    auto A_data = get_raw_data<UQType>(A_buf);
    auto B_data = get_raw_data<UQType>(B_buf);

    auto A_shape = A_buf.shape().make_smooth();
    auto B_shape = B_buf.shape().make_smooth();
    assert(A_shape.rank() == 2);
    assert(B_shape.rank() == 2);
    auto n_rows = A_shape.extent(0);
    auto n_cols = A_shape.extent(1);
    assert(n_rows == B_shape.extent(0));
    assert(n_cols == B_shape.extent(1));
    Smooth shape({n_rows, n_cols});

    std::vector<double> A_median(shape.size());
    std::vector<double> B_median(shape.size());

    for(std::size_t i = 0; i < shape.size(); ++i) {
        A_median[i] = A_data[i].median();
        B_median[i] = B_data[i].median();
    }

    Contiguous A_median_buffer(std::move(A_median), shape);
    Contiguous B_median_buffer(std::move(B_median), shape);
    tensorwrapper::Tensor A_tensor(shape, std::move(A_median_buffer));
    tensorwrapper::Tensor B_tensor(shape, std::move(B_median_buffer));

    return eigen_solver_mod.run_as<pt>(A_tensor, B_tensor);
};

template<typename UQType>
auto convert_to_uq(const tensorwrapper::Tensor& values,
                   const tensorwrapper::Tensor& vectors) {
    using float_type = typename UQType::value_t;
    auto values_buf  = make_contiguous(values.buffer());
    auto vectors_buf = make_contiguous(vectors.buffer());
    auto n           = values_buf.shape().extent(0);
    tensorwrapper::shape::Smooth vector_shape{n};
    tensorwrapper::shape::Smooth matrix_shape{n, n};
    using tensorwrapper::buffer::make_contiguous;
    auto uq_values_buffer  = make_contiguous<UQType>(vector_shape);
    auto uq_vectors_buffer = make_contiguous<UQType>(matrix_shape);
    using wtf::fp::float_cast;
    for(std::size_t i = 0; i < n; ++i) {
        auto vi = float_cast<float_type>(values_buf.get_elem({i}));
        uq_values_buffer.set_elem({i}, UQType(vi));
        for(std::size_t j = 0; j < n; ++j) {
            auto vij = float_cast<float_type>(vectors_buf.get_elem({i, j}));
            uq_vectors_buffer.set_elem({i, j}, UQType(vij));
        }
    }
    tensorwrapper::Tensor uq_values(vector_shape, std::move(uq_values_buffer));
    tensorwrapper::Tensor uq_vectors(matrix_shape,
                                     std::move(uq_vectors_buffer));
    return std::make_pair(uq_values, uq_vectors);
}

template<typename UQType>
auto compute_residual(const tensorwrapper::Tensor& A,
                      const tensorwrapper::Tensor& B,
                      const tensorwrapper::Tensor& C,
                      const tensorwrapper::Tensor& L) {
    using label_type = tensorwrapper::Tensor::label_type;
    using namespace tensorwrapper::buffer;
    label_type ij("i,j");
    label_type ji("j,i");
    label_type jk("j,k");
    label_type ik("i,k");
    label_type j("j");

    tensorwrapper::Tensor TA, TAC;
    TA(ik)  = C(ji) * A(jk);
    TAC(ik) = TA(ij) * C(jk);

    tensorwrapper::Tensor TB, TBC;
    TB(ik)  = C(ji) * B(jk);
    TBC(ik) = TB(ij) * C(jk);

    // TODO: Replace when batch contraction works...
    auto TBC_buf     = make_contiguous(TBC.buffer());
    auto L_buf       = make_contiguous(L.buffer());
    auto TBCL_shape  = TBC_buf.shape().make_smooth();
    auto TBCL_buffer = make_contiguous<UQType>(TBCL_shape);

    auto n_rows = TBCL_shape.extent(0);
    auto n_cols = TBCL_shape.extent(1);
    using wtf::fp::float_cast;
    for(std::size_t i = 0; i < n_rows; ++i) {
        auto Li = float_cast<UQType>(L_buf.get_elem({i}));
        for(std::size_t j = 0; j < n_cols; ++j) {
            auto TBCij = float_cast<UQType>(TBC_buf.get_elem({i, j}));
            TBCL_buffer.set_elem({i, j}, TBCij * Li);
        }
    }
    tensorwrapper::Tensor TBCL(TBCL_shape, std::move(TBCL_buffer));

    tensorwrapper::Tensor H;
    H(ij) = TAC(ij) - TBCL(ij);
    return H;
}

template<typename UQType>
auto compute_eta(const tensorwrapper::Tensor& H) {
    auto H_buffer = make_contiguous(H.buffer());
    auto shape    = H_buffer.shape();
    auto n_rows   = shape.extent(0);

    using float_type = typename UQType::value_t;
    using wtf::fp::float_cast;
    float_type norm = 0.0;
    for(std::size_t i = 0; i < n_rows; ++i) {
        float_type col_norm = 0.0;
        for(std::size_t j = 0; j < n_rows; ++j) {
            auto Hji    = H_buffer.get_elem({j, i});
            auto Hji_uq = float_cast<UQType>(Hji);
            auto Hji_abs =
              std::max(std::abs(Hji_uq.lower()), std::abs(Hji_uq.upper()));
            col_norm += Hji_abs * Hji_abs;
        }
        col_norm = std::sqrt(col_norm / n_rows);
        norm     = std::max(norm, col_norm);
    }
    return 6.0 * std::sqrt(n_rows) * norm;
}

template<typename UQType>
auto t_ball(const tensorwrapper::Tensor& C,
            const typename UQType::value_t& eta) {
    using namespace tensorwrapper::buffer;
    UQType pm1(-1.0, 1.0);
    auto shape = make_contiguous(C.buffer()).shape().make_smooth();
    std::vector<UQType> eta_omega(shape.size(), eta * pm1);
    Contiguous eta_omega_buffer(std::move(eta_omega), shape);
    tensorwrapper::Tensor eta_omega_tensor(shape, std::move(eta_omega_buffer));
    tensorwrapper::Tensor T_ball, COmega;
    COmega("i,k") = C("i,j") * eta_omega_tensor("j,k");
    T_ball("i,j") = C("i,j") + COmega("i,j");
    return T_ball;
}

template<typename UQType>
auto l_ball(const tensorwrapper::Tensor& L,
            const typename UQType::value_t& eta) {
    using namespace tensorwrapper::buffer;
    UQType one_pm_eta(1.0 - eta, 1.0 + eta);
    auto shape = make_contiguous(L.buffer()).shape().make_smooth();
    std::vector<UQType> one_pm_eta_data(shape.size(), one_pm_eta);
    Contiguous one_pm_eta_buffer(std::move(one_pm_eta_data), shape);
    tensorwrapper::Tensor one_pm_eta_tensor(shape,
                                            std::move(one_pm_eta_buffer));
    tensorwrapper::Tensor L_ball;
    L_ball("i") = L("i") * one_pm_eta_tensor("i");
    return L_ball;
}

} // namespace

MODULE_CTOR(BallGeneralized) {
    description(desc);
    satisfies_property_type<pt>();

    add_submodule<pt>("Eigen Solve");
}

MODULE_RUN(BallGeneralized) {
    auto&& [A, B] = pt::unwrap_inputs(inputs);

    using uq_type = tensorwrapper::types::idouble;
    using namespace tensorwrapper::buffer;
    using namespace tensorwrapper::shape;

    auto A_buf = make_contiguous(A.buffer());
    auto B_buf = make_contiguous(B.buffer());
    auto shape = A_buf.shape().make_smooth();

    // N.b., wrap_subdiagonalization will verify shapes match
    auto n = shape.extent(0);

    Smooth vector_shape{n};
    Smooth matrix_shape{n, n};

    auto eigen_solver_mod = submods.at("Eigen Solve");
    auto [values, vectors] =
      wrap_subdiagonalization<uq_type>(A_buf, B_buf, eigen_solver_mod);

    // Here we need to convert the Eigen values and vectors to UQ type
    auto [uq_values, uq_vectors] = convert_to_uq<uq_type>(values, vectors);

    auto H = compute_residual<uq_type>(A, B, uq_vectors, uq_values);

    auto eta_val = compute_eta<uq_type>(H);

    auto C_ball = t_ball<uq_type>(uq_vectors, eta_val);
    auto L_ball = l_ball<uq_type>(uq_values, eta_val);

    auto rv = results();
    return pt::wrap_results(rv, L_ball, C_ball);
}
} // namespace scf::eigen_solver
