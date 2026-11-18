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

#pragma once
#include <chemist/grid/grid_class.hpp>
#include <pluginplay/pluginplay.hpp>
#include <simde/evaluate_braket/evaluate_braket.hpp>
#include <simde/types.hpp>

namespace scf::property_types {

/** @brief Evaluates the quantity of another property type by quadrature.
 *
 *  Modules satisfying this property type evaluate the quantity computed by
 *  @p BasePT by quadrature on a grid they are given, rather than on a grid
 *  they build themselves. This lets the caller choose the grid, *e.g.*, to
 *  evaluate the same quadrature on a finer grid to estimate the error of the
 *  grid. The inputs are those of @p BasePT plus the grid, and the results are
 *  those of @p BasePT.
 *
 *  @tparam BasePT The property type whose quantity is evaluated by
 *                 quadrature, *e.g.*, simde::aos_xc_e_aos for the XC
 *                 potential.
 */
template<typename BasePT>
DECLARE_TEMPLATED_PROPERTY_TYPE(Quadrature, BasePT);

template<typename BasePT>
TEMPLATED_PROPERTY_TYPE_INPUTS(Quadrature, BasePT) {
    using grid_type = const chemist::Grid&;
    auto rv         = BasePT::inputs();
    auto rv0        = rv.template add_field<grid_type>("Grid");
    rv0["Grid"].set_description("The grid the quadrature is evaluated on");
    return rv0;
}

template<typename BasePT>
TEMPLATED_PROPERTY_TYPE_RESULTS(Quadrature, BasePT) {
    return BasePT::results();
}

/// Property type for evaluating the XC energy by quadrature
using xc_energy_quadrature =
  Quadrature<simde::eval_braket<simde::type::rscf_wf, simde::type::XC_e_type,
                                simde::type::rscf_wf>>;

/// Property type for evaluating the AO matrix of the XC potential by quadrature
using xc_potential_quadrature = Quadrature<simde::aos_xc_e_aos>;

} // namespace scf::property_types
