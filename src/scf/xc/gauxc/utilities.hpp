/*
 * Copyright 2025 NWChemEx-Project
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

#pragma once
#include <Eigen/Dense>
#include <tensorwrapper/tensorwrapper.hpp>

namespace scf::xc::gauxc {

template<typename FloatType>
auto tw_to_eigen(const tensorwrapper::Tensor& t) {
    auto n = t.logical_layout().shape().as_smooth().extent(0);
    auto m = t.logical_layout().shape().as_smooth().extent(1);
    Eigen::Matrix<FloatType, Eigen::Dynamic, Eigen::Dynamic> t_eigen(n, m);

    auto t_vector = tensorwrapper::buffer::make_contiguous(t.buffer());
    auto t_data   = t_vector.get_immutable_data();
    auto begin = wtf::buffer::contiguous_buffer_cast<const FloatType>(t_data);
    auto end   = begin.data() + (n * m);
    std::copy(begin.data(), end, t_eigen.data());
    return t_eigen;
}

template<typename FloatType>
auto eigen_to_tw(
  const Eigen::Matrix<FloatType, Eigen::Dynamic, Eigen::Dynamic>& t_eigen,
  const parallelzone::runtime::RuntimeView& rt) {
    auto n = static_cast<std::size_t>(t_eigen.rows());
    auto m = static_cast<std::size_t>(t_eigen.cols());

    tensorwrapper::shape::Smooth shape{n, m};
    auto pbuffer = tensorwrapper::buffer::make_contiguous<FloatType>(shape);
    for(std::size_t i = 0; i < n; ++i) {
        for(std::size_t j = 0; j < m; ++j)
            pbuffer.set_elem({i, j}, t_eigen(i, j));
    }
    return tensorwrapper::Tensor(std::move(shape), std::move(pbuffer));
}

} // namespace scf::xc::gauxc
