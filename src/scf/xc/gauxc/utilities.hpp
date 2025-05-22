#pragma once
#include <tensorwrapper/tensorwrapper.hpp>

namespace scf::xc::gauxc {

template<typename FloatType>
auto tw_to_eigen(const tensorwrapper::Tensor& t) {
    auto n = t.logical_layout().shape().as_smooth().extent(0);
    auto m = t.logical_layout().shape().as_smooth().extent(1);
    Eigen::Matrix<FloatType, Eigen::Dynamic, Eigen::Dynamic> t_eigen(n, m);

    const auto& rt = t.buffer().allocator().runtime();
    tensorwrapper::allocator::Eigen<FloatType> alloc(rt);

    auto t_vector = alloc.rebind(t.buffer());
    auto begin    = t_vector.get_immutable_data();
    auto end      = begin + (n * m);
    std::copy(begin, end, t_eigen.data());
    return t_eigen;
}

template<typename FloatType>
auto eigen_to_tw(
  const Eigen::Matrix<FloatType, Eigen::Dynamic, Eigen::Dynamic>& t_eigen,
  const parallelzone::runtime::RuntimeView& rt) {
    tensorwrapper::allocator::Eigen<FloatType> alloc(rt);
    auto n = static_cast<std::size_t>(t_eigen.rows());
    auto m = static_cast<std::size_t>(t_eigen.cols());

    tensorwrapper::shape::Smooth shape{n, m};
    tensorwrapper::layout::Physical layout(shape);
    auto pbuffer = alloc.allocate(layout);
    for(std::size_t i = 0; i < n; ++i) {
        for(std::size_t j = 0; j < m; ++j)
            pbuffer->set_elem({i, j}, t_eigen(i, j));
    }
    return tensorwrapper::Tensor(std::move(shape), std::move(pbuffer));
}

} // namespace scf::xc::gauxc