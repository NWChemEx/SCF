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

#include <pluginplay/submodule_request.hpp>
#include <simde/simde.hpp>
#include <stdexcept>
#include <tensorwrapper/buffer/contiguous.hpp>

namespace scf::eigen_solver {

namespace {

const auto desc = R"(
Eigen Solve Driver
------------------

Routes normal eigenvalue problems to the appropriate backend based on the
floating-point type stored in the input matrix. Plain ``float``/``double``
matrices use the deterministic Eigen backend; uncertain and interval types use
ball-arithmetic certification.
)";

constexpr const char* kNoneSubmodule      = "none";
constexpr const char* kUncertainSubmodule = "uncertain";
constexpr const char* kIntervalSubmodule  = "interval";

using pt = simde::EigenSolve;

struct Router {
    using return_t = std::pair<simde::type::tensor, simde::type::tensor>;

    const simde::type::tensor& m_A;
    pluginplay::type::submodule_map& m_submods;

    Router(const simde::type::tensor& A,
           pluginplay::type::submodule_map& submods) :
      m_A(A), m_submods(submods) {}

    template<typename FloatType>
    return_t operator()(const std::span<FloatType>&) {
        using clean_t = std::decay_t<FloatType>;
        if constexpr(tensorwrapper::types::is_interval_v<clean_t>) {
            return m_submods.at(kIntervalSubmodule).run_as<pt>(m_A);
        } else if constexpr(tensorwrapper::types::is_uncertain_v<clean_t>) {
            return m_submods.at(kUncertainSubmodule).run_as<pt>(m_A);
        } else {
            return m_submods.at(kNoneSubmodule).run_as<pt>(m_A);
        }
    }
};

} // namespace

MODULE_CTOR(EigenSolveDriver) {
    description(desc);
    satisfies_property_type<pt>();

    add_submodule<pt>(kNoneSubmodule);
    add_submodule<pt>(kUncertainSubmodule);
    add_submodule<pt>(kIntervalSubmodule);
}

MODULE_RUN(EigenSolveDriver) {
    auto&& [A] = pt::unwrap_inputs(inputs);

    using tensorwrapper::buffer::make_contiguous;
    using tensorwrapper::buffer::visit_contiguous_buffer;
    const auto& A_buffer = make_contiguous(A.buffer());
    Router router(A, submods);
    auto [values, vectors] = visit_contiguous_buffer(router, A_buffer);

    auto rv = results();
    return pt::wrap_results(rv, values, vectors);
}

} // namespace scf::eigen_solver
