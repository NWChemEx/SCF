.. Copyright 2026 NWChemEx-Project
..
.. Licensed under the Apache License, Version 2.0 (the "License");
.. you may not use this file except in compliance with the License.
.. You may obtain a copy of the License at
..
.. http://www.apache.org/licenses/LICENSE-2.0
..
.. Unless required by applicable law or agreed to in writing, software
.. distributed under the License is distributed on an "AS IS" BASIS,
.. WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
.. See the License for the specific language governing permissions and
.. limitations under the License.

.. _xc_inputs_uq_design:

######################################
Uncertainty in the XC Inputs on a Grid
######################################

This page describes how SCF estimates the uncertainty in the quantities the
XC modules evaluate at the grid points: today the density, :math:`\rho_i`, and
the values of the atomic orbitals (AOs), :math:`\phi_{\mu i}`; for GGAs and
beyond also, *e.g.*, the gradient of the density. An error model estimates the
uncertainty of each value, and the value is returned with that uncertainty
attached, as a UQ type. How that uncertainty reaches the XC energy and
potential, and how the error of the grid itself is handled, is described in
:ref:`xc_uq_design`.

In the notation of :ref:`xc_uq_design`, this page covers the input errors,
:math:`\delta E^{Inputs}` and :math:`\delta V^{Inputs}`: the uncertainty of the
values *at* the grid points. The error of the grid, *i.e.*, the screening of
the quadrature and its incompleteness, is attached to the XC energy and
potential instead.

As in :ref:`xc_uq_design`, the UQ type is assumed to be an asymmetric interval
(see :ref:`xc_uq_asymmetric_intervals`): it stores a lower bound, an upper
bound, and a mean, which need not be the midpoint. For the density at a grid
point, the mean is :math:`\rho_i` and the bounds are
:math:`\rho_i - \delta\rho_i^-` and :math:`\rho_i + \delta\rho_i^+`.

***************************************
What are the Density and AOs on a Grid?
***************************************

Let :math:`\vec{r}_i` be the :math:`i`-th of :math:`N_g` grid points,
:math:`\phi_\mu` the :math:`\mu`-th of :math:`N_b` AOs, and
:math:`P_{\mu\nu}` the total (*i.e.*, summed over spin) AO density matrix. The
density at :math:`\vec{r}_i` is:

.. math::

   \rho_i \equiv \rho\left(\vec{r}_i\right) =
      \sum_{\mu\nu}^{N_b} \phi_{\mu i} P_{\mu\nu} \phi_{\nu i}
      = \boldsymbol{\phi}_i^T \mathbf{P} \boldsymbol{\phi}_i

where :math:`\phi_{\mu i} = \phi_\mu\left(\vec{r}_i\right)` is an element of
the AO collocation matrix, :math:`\boldsymbol{\Phi}`, and
:math:`\boldsymbol{\phi}_i` is its :math:`i`-th column.

In SCF the collocation matrix is computed by the ``AOs on a grid`` module,
which satisfies ``simde::AOCollocationMatrix``. It evaluates every AO at every
grid point. The density is computed by the ``Density2Grid`` module, which
satisfies ``simde::EDensityCollocationMatrix``. Its inputs are the grid and the
density (a ``simde::type::e_density`` and the AO basis set). The matrix the
``e_density`` holds is not :math:`\mathbf{P}`: it is the restricted density
matrix, :math:`\mathbf{D}`, formed from products of the occupied orbitals
without a spin factor. Assuming restricted orbitals, ``Density2Grid`` forms
:math:`\mathbf{P} = 2\mathbf{D}`, gets the collocation matrix from its
``"AOs on a grid"`` submodule, forms
:math:`\mathbf{X} = \mathbf{P}\boldsymbol{\Phi}`, and takes the column-wise dot
product of :math:`\boldsymbol{\Phi}` and :math:`\mathbf{X}`. Its result is a
tensor with one element per grid point. Error models that work with the
density matrix must likewise double the stored matrix, and its uncertainty,
before using it as :math:`\mathbf{P}`.

The XC modules use both: the density as the argument of the functional, and
the collocation matrix to assemble the AO matrix of the potential.

Why Inputs Lose Precision
=========================

The value of :math:`\rho_i` is uncertain for four reasons. The first two
introduce error when :math:`\rho_i` is evaluated. The third is uncertainty
inherited from :math:`\mathbf{P}`.

#. **Floating-point error.** :math:`\rho_i` is a sum of
   :math:`N_b^2` terms of both signs, since off-diagonal elements of
   :math:`\mathbf{P}` can be negative. Where the terms nearly cancel, the sum
   is much smaller than its largest terms and if care is not taken the small
   terms will be lost in the rounding error.
#. **Screening.** Production evaluators skip AOs whose value at a point is
   below a threshold, and skip blocks of :math:`\mathbf{P}` that are small.
   This means that even if the :math:`i`-th grid point is retained,
   :math:`\rho_i` may contain an error due to dropping contributions from
   the collocation matrix and/or :math:`\mathbf{P}`.
#. **Uncertainty already in the AO quantities.** The AOs and the AO density
   matrix may already carry uncertainty. This uncertainty propagates to the
   collocation matrix and the :math:`\rho_i` computed from it.
#. **Uncertain grid points.** Grid points may be uncertain, e.g., due to the
   floating point arithmetic used to compute their weights and locations.

*****************
Statement of Need
*****************

In :ref:`xc_uq_design` the uncertainty in the XC energy and potential due to
their inputs follows from evaluating them with UQ types. That requires the
inputs to carry their uncertainty, *i.e.*, an error model for each input.
None exists today. Sources of error were listed above.

**************
Considerations
**************

#. Same inputs as the modules they replace

   - The modules that attach the uncertainty must be drop-in replacements for
     ``Density2Grid`` and ``AOs on a grid``, so that the XC modules do not
     change. They therefore satisfy ``simde::EDensityCollocationMatrix`` and
     ``simde::AOCollocationMatrix``.

#. The evaluator is a choice

   - The algorithm that evaluates the quantity (``Density2Grid`` today, a
     screened or GPU evaluator later) is chosen independently of the algorithm
     that estimates its uncertainty.

#. The error model is a choice

   - Rigorous bounds, cheap bounds, and estimates all have uses. As in
     :ref:`xc_uq_design`, each source of uncertainty has its own error model,
     which is a submodule of the module that attaches the uncertainty, and
     the user chooses one by binding that submodule.

#. Uncertainty at the grid points only

   - The error models estimate the uncertainty of the values at the grid
     points, not their variation between grid points, nor the effect of
     dropping grid points. Both belong to the error of the grid (see
     :ref:`xc_uq_design`).

#. Sound physical constraints

   - The density is non-negative. Any uncertain density returned must keep its
     lower bound at or above zero, by intersecting with
     :math:`[0, \infty)`, which can only remove values that can not be the
     true one. With an asymmetric interval this leaves the mean unchanged.

#. Only sources of uncertainty know the uncertain type

   - As in :ref:`qcup_scf_design`, the uncertain type comes from a
     ``UQ Initializer`` submodule rather than an input of its own. As in
     :ref:`xc_uq_design`, the error models create the uncertainty, so they
     are the modules that call ``UQ Initializer``; the wrappers only receive
     results that are already uncertain.

#. Ordinary DFT is unaffected

   - Nothing changes unless the user binds the new modules.

****************************
Grid Input UQ Implementation
****************************

.. note::

   The initial implementation of QCUP-DFT will assume that the AOs and the AO
   density matrix have no uncertainty. Thus the modules in this section are
   not actually implemented yet, but the design is here to show how they are
   envisioned to work.

SCF will add two modules, ``UQ Density2Grid`` and ``UQ AOs on a grid``, that
follow the pattern of ``UQ XC Energy`` and ``UQ XC Potential`` in
:ref:`xc_uq_design`. Each satisfies the same property type as the module it
wraps, runs that module through an evaluator submodule, passes the result to
one error-model submodule per source of uncertainty, and returns the uncertain
result the error models produce. They are ordinary modules, not PluginPlay
drivers because they only manage their inputs and their submodules.

The error models satisfy the ``XCUncertain<U>`` property type of
:ref:`xc_uq_design`, with ``U`` the property type of the evaluator, *i.e.*,
``XCUncertain<simde::EDensityCollocationMatrix>`` and
``XCUncertain<simde::AOCollocationMatrix>``. Their inputs are the inputs of
``U`` followed by the result of the evaluator, so the evaluator is not run a
second time. Their results are the result of the evaluator with uncertainties
applied, *i.e.*, each error model calls ``UQ Initializer``, keeps the
evaluator's result as the mean, and attaches its deviations. The density error
models also pin the lower bound at zero.

**********************
Further Considerations
**********************

- **GGAs and beyond.** GGAs also need
  :math:`\nabla\rho_i = 2\sum_{\mu\nu}\nabla\phi_{\mu i}P_{\mu\nu}\phi_{\nu i}`
  on the grid, and the gradient of the collocation matrix. The same
  strategies apply, with the gradient of the collocation matrix in place of
  one factor of it. Meta-GGAs add the kinetic-energy density, which is
  non-negative like the density and needs its lower bound pinned the same way.
