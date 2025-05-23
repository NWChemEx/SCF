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

#include "scf_modules.hpp"
#include <cmath>
#include <simde/simde.hpp>

namespace scf {

using pt = simde::charge_charge_interaction;

const auto desc = R"(
Charge-Charge Interaction via Coulomb's Law
-------------------------------------------

Let :math:`Q` be a set of :math:`N_Q` point charges such that the :math:`i`-th
point charge has charge :math:`q_i` and position :math:`\mathbf{r_i}`. The
electrostatic potential at the point :math:`\mathbf{r}` generated by :math:`Q`,
:math:`E_Q(\mathbf{r})` is given by:

.. math::

   E_Q(\mathbf{r}) = \sum_{i=1}^{N_Q} 
                     \frac{q_i}{\left|\mathbf{r} - \mathbf{r_i}\right|}

Let :math:`S` be a set of :math:`N_S` point charges that is disjoint with 
:math:`Q`, then the interaction of :math:`S` with :math:`Q`, :math:`V_{SQ}` is 
given by:

.. math::

   V_{SQ} =& \sum_{i=1}^{N_S} q_i E_Q(\mathbf{r_i})\\
          =& \sum_{i=1}^{N_S}\sum_{j=1}^{N_Q} 
             \frac{q_i q_j}{\left|\mathbf{r_i} - \mathbf{r_j}\right|}

This module will compute :math:`V_{SQ}`, including in the common scenario where
:math:`S` equals :math:`Q` (in which case terms for which the denominator is
zero are skipped).
)";

MODULE_CTOR(CoulombsLaw) {
    description(desc);
    satisfies_property_type<pt>();
}

MODULE_RUN(CoulombsLaw) {
    const auto& [qlhs, qrhs] = pt::unwrap_inputs(inputs);

    // This algorithm only works as long as qlhs and qrhs are disjoint or
    // identical
    double e = 0.0;
    for(const auto&& qi : qlhs) {
        for(const auto&& qj : qrhs) {
            if(qi == qj) break;
            auto dx = qi.x() - qj.x();
            auto dy = qi.y() - qj.y();
            auto dz = qi.z() - qj.z();

            e += qi.charge() * qj.charge() /
                 std::sqrt(dx * dx + dy * dy + dz * dz);
        }
    }

    auto rv = results();
    return pt::wrap_results(rv, e);
}

} // namespace scf