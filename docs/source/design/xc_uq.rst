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

.. _xc_uq_design:

##########################################
Uncertainty in the XC Energy and Potential
##########################################

This page describes how SCF tracks the uncertainty in the exchange-correlation
(XC) energy, :math:`E^{XC}`, and the AO matrix of the XC potential,
:math:`V^{XC}`, when they are evaluated by quadrature. It is part of the
uncertainty-quantified SCF effort described in :ref:`qcup_scf_design`, and
extends it from Hartree-Fock to Kohn-Sham DFT.

The uncertainty has two kinds of sources. The quantities evaluated at the grid
points, *e.g.*, the density, are uncertain; how their uncertainty is estimated
is described in :ref:`xc_inputs_uq_design`. The quadrature itself is also
approximate, because the grid is incomplete. This page covers how both reach
:math:`E^{XC}` and :math:`V^{XC}`, and how the second is estimated.

.. _xc_uq_asymmetric_intervals:

********************************
Assumption: Asymmetric Intervals
********************************

Throughout this page and :ref:`xc_inputs_uq_design`, the UQ type is assumed to
be an **asymmetric interval**. It stores a lower bound, an upper bound, and a
mean (the best estimate of the value), and the mean need not be the midpoint of
the bounds. For a value :math:`x`, the bounds are :math:`x - \delta x^-` and
:math:`x + \delta x^+`, and the deviations
:math:`\delta x^-, \delta x^+ \ge 0` need not be equal.

The reason for this assumption is that the density can not be negative, so its
lower bound must be pinned at zero whenever an error model's estimate would go
below it. An asymmetric interval allows this without changing the mean. A
symmetric interval stores only the bounds: pinning
:math:`[\rho_i - r_i, \rho_i + r_i]` at zero gives :math:`[0, \rho_i + r_i]`,
whose midpoint, :math:`(\rho_i + r_i)/2`, is larger than :math:`\rho_i`, so the
best estimate of the density would be biased upward. The error of a quadrature
also tends to have a known direction, which an asymmetric interval can
represent. ``sigma::Interval`` stores only the bounds, so this type does not
exist yet.

*************************************
Errors in the XC Energy and Potential
*************************************

Let :math:`\vec{r}_i` and :math:`w_i` be the point and weight of the
:math:`i`-th of :math:`N_g` grid points, :math:`\phi_{\mu i}` the value of the
:math:`\mu`-th atomic orbital (AO) at :math:`\vec{r}_i`, and :math:`\rho_i` the
density at :math:`\vec{r}_i`. For an LDA with energy density :math:`f` and
:math:`v = \partial f / \partial \rho`, the XC energy and potential of the
functional, :math:`\tilde{E}^{XC}` and :math:`\tilde{V}^{XC}`, are:

.. math::

   \tilde{E}^{XC} \equiv&\int f\left[\rho(\vec{r})\right] d\vec{r}\\
   \tilde{V}^{XC}_{\mu\nu} \equiv&
            \int v\left[\rho(\vec{r})\right]
                 \phi_\mu(\vec{r}) \phi_\nu(\vec{r}) d\vec{r}.

These are the exact integrals of the chosen (approximate) functional, not the
exact XC energy and potential, which are unknown. In practice, the integrals
are evaluated by quadrature. Let :math:`E^{Quad}` and :math:`V^{Quad}` be the
values obtained by exactly applying a quadrature rule:

.. math::

   E^{Quad} =&\sum_i w_i f\left[\rho_i\right]\\
   V^{Quad}_{\mu\nu} =& \sum_i w_i \phi_{\mu i} v\left[\rho_i\right]
                         \phi_{\nu i}.

Of note these equations assume:

#. The inputs (*e.g.*, :math:`\rho_i`) are evaluated at the grid points and
   are assumed to be exact.
#. The sums are over all grid points in the quadrature rule, *i.e.*, no terms
   are screened out.

Subject to these assumptions, the errors in the quadratures are then:

.. math::

    \delta E^{Quad} =& E^{Quad} - \tilde{E}^{XC}\\
    \delta V^{Quad}_{\mu\nu} =& V^{Quad}_{\mu\nu} - \tilde{V}^{XC}_{\mu\nu}.

If the sums are screened, terms are dropped from them and we actually compute
:math:`E^{Scr}` and :math:`V^{Scr}`:

.. math::

   E^{Scr} =&
    \sum_{i\in\mathcal{K}} w_i f\left[\rho_i\right]\\
   V^{Scr}_{\mu\nu} =&
    \sum_{i\in\mathcal{K}_{\mu\nu}} w_i \phi_{\mu i} v\left[\rho_i\right]
                                       \phi_{\nu i},

where :math:`\mathcal{K}` is the set of grid points that are kept in the
energy, and :math:`\mathcal{K}_{\mu\nu}` is the set of grid points that are
kept for the AO pair :math:`\mu\nu` in the potential. The latter allows for
screening that drops individual AO pairs at a grid point (*e.g.*, because
:math:`\left|\phi_{\mu i}\phi_{\nu i}\right|` is small), not only whole grid
points. The screening errors relative to full quadrature are:

.. math::

    \delta E^{Scr} =& E^{Scr} - E^{Quad}\\
    \delta V^{Scr}_{\mu\nu} =& V^{Scr}_{\mu\nu} - V^{Quad}_{\mu\nu}.


Similarly, if the inputs are approximated we actually compute :math:`E^{Inputs}`
and :math:`V^{Inputs}` and the errors from using approximate inputs relative to
full quadrature are:

.. math::

    \delta E^{Inputs} =& E^{Inputs} - E^{Quad}\\
    \delta V^{Inputs}_{\mu\nu} =& V^{Inputs}_{\mu\nu} - V^{Quad}_{\mu\nu}.

Note that the total error will in general not simply be the sum of the three
errors because of cross terms. For example, if both screening and approximate
inputs are used, then screening out a term will also change the uncertainty from
the input uncertainty to that term. Other errors, such as floating-point
rounding, can be built up using the same framework, *i.e.*, as the difference
between the value computed with the error and the value computed without it.

Throughout this page, :math:`E^{XC}` and :math:`V^{XC}` denote the values that
are actually computed, *i.e.*, with whatever screening and approximate inputs
are used, and :math:`\tilde{E}^{XC}` and :math:`\tilde{V}^{XC}` always denote
the exact integrals, *i.e.*, the reference values.

*****************
Statement of Need
*****************

QCUP-SCF attaches an uncertainty to every approximately computed quantity by
storing it as a UQ type, and lets the uncertainty propagate through the rest of
the calculation by doing arithmetic with those types. For the ERIs, an error
model attaches the uncertainty, and the uncertainty in the Coulomb and exchange
matrices, the Fock matrix, and the energy follows from the arithmetic.

As just motivated, the uncertainty in the computed :math:`E^{XC}` and
:math:`V^{XC}` has three sources:

#. **Uncertainty in the inputs.** :math:`\rho_i` and :math:`\phi_{\mu i}`
   (and, for GGAs and beyond, further quantities such as
   :math:`\nabla\rho_i`) are uncertain. If they are UQ types, the XC modules
   evaluate the quadratures with UQ types, and the inputs' uncertainty
   propagates through the arithmetic without further consideration. Error
   models for attaching the uncertainty to the inputs are described in
   :ref:`xc_inputs_uq_design`. Without UQ types, we would need an error model
   to attach the uncertainty to :math:`E^{XC}` and
   :math:`V^{XC}` directly.
#. **Screening of the grid.** The quadrature sums are screened, so some terms
   are dropped. This is an error in the quadrature, but does not stem from the
   uncertainty in any of the inputs, so it needs an error model of its own,
   which attaches uncertainty to :math:`E^{XC}` and
   :math:`V^{XC}` directly.
#. **Incompleteness of the grid.** Even with exact inputs and the full
   quadrature, the quadrature differs from the integral. This error is attached
   to :math:`E^{XC}` and :math:`V^{XC}` by its own error model.

Error models for the first source are described in :ref:`xc_inputs_uq_design`.
The purpose of this page is to design the infrastructure and error models for
the second and third error sources. For simplicity we will refer to the latter
two uncertainty sources as **the error of the grid**.

**************
Considerations
**************

#. Each error source is counted once

   - Uncertainty in the inputs is attached to the inputs; uncertainty from
     the grid is attached to :math:`E^{XC}` and :math:`V^{XC}`.
     The error models must not double-count either source of uncertainty.

#. Same interface as the XC modules

   - The module for attaching uncertainty to :math:`E^{XC}` and
     :math:`V^{XC}` is a wrapper around an XC module. It has the same
     property types, so it can be leveraged wherever the XC module was
     leveraged.

#. As additional approximations are made, the module must be extendible so that
   it can attach uncertainty for those sources as well.

#. The error model is a choice

   - As with the ERIs, there are several ways to estimate the error from the
     grid. The module should facilitate choosing among them.

#. The uncertain type is chosen in one place

   - As in :ref:`qcup_scf_design`, only the modules that create uncertainty
     know the uncertain type, and they get it from a ``UQ Initializer``
     submodule. Here, those are the error models.

#. Ordinary DFT is unaffected

   - Uncertainty should be done in wrappers over ordinary DFT or rely on
     templated DFT kernels to propagate uncertainty.
   - Code meant for use with traditional DFT should not need to know about
     uncertainty propagation.

********************
XC UQ Implementation
********************

SCF adds two modules, ``UQ XC Energy`` and ``UQ XC Potential``, that are
respectively responsible for computing the XC energy and XC potential with
uncertainty. They follow the same pattern, so this section only describes the
potential module.

As part of this design we introduce a new property type ``Quadrature<T>``. This
PT will live in the SCF repo for now. ``T`` here will be the property type of
the quantity that is being evaluated by quadrature, *e.g.*, for the XC
potential it is
``simde::aos_xc_e_aos``. In addition to taking the inputs ``T`` needs,
``Quadrature<T>`` will also take the grid. SCF's XC energy and potential modules
will be updated accordingly so that they call modules satisfying
``Quadrature<T>``.

The ``UQ XC Potential`` module satisfies the same property type as an XC
potential module (i.e., ``simde::aos_xc_e_aos``). The initial implementation
will focus on attaching the uncertainty caused by using an incomplete grid. The
current control flow is as follows:

- Unwrap inputs and pass them to the wrapped XC potential module.
- Call a submodule to estimate the uncertainty from the incompleteness error
  and attach the uncertainty.
- Call submodules for other sources of uncertainty, if any (eventually
  screening). Module's design dictates whether the uncertainty is attached to
  the results of the wrapped XC potential module or to the uncertainty from the
  incompleteness error.
- Return the final UQ results.

``XCUncertain<Quadrature<T>>`` will be the property type that estimates the
uncertainty from the incompleteness of the grid. The template parameter ``T``
is the property type of the module being evaluated by quadrature, i.e., for the
XC potential it is ``simde::aos_xc_e_aos``. ``XCUncertain<U>`` will be a PT that
lives in SCF that takes the inputs of ``U`` and the results of calling a module
satisfying ``U``. The ``Uncertain<V>`` property type in SimDE should eventually
be updated to replace ``XCUncertain<U>``, but is out of scope for this design.
The result of calling a module
satisfying ``XCUncertain<Quadrature<T>>`` will be the input results with
uncertainties applied, i.e., a matrix of the same size as the XC potential
matrix where each element is a UQ type object that uses the input result as the
mean of the returned result. It is the module which satisfies
``XCUncertain<Quadrature<T>>`` that is responsible for calling
``UQ Initializer``.

The QCUP-SCF driver will be in charge of swapping out the error models. We
will add submodules to the QCUP-SCF driver for each source of uncertainty, e.g.,
a submodule satisfying ``XCUncertain<Quadrature<simde::aos_xc_e_aos>>`` will be
added to QCUP-SCF and targeted at the ``UQ XC Potential``'s submodule. Users
can thus control the error models from the driver level.

To summarize, the submodules of the ``UQ XC Potential`` module are:

.. list-table::
   :header-rows: 1
   :widths: 30 30 40

   * - Hook name
     - Default key
     - Purpose
   * - ``XC Potential``
     - ``SVWN3 Potential``
     - The module that evaluates the quadrature.
   * - ``Incompleteness Error Model``
     - ``Reference Grid Error``
     - The error model that estimates the uncertainty from the incompleteness
       of the grid.

The ``Reference Grid Error`` module refers to the "Reference Grid" error model
introduced in the next section.

***************************
Incompleteness Error Models
***************************

The strategies below estimate the incompleteness errors,
:math:`\delta E^{Quad}` and, element-wise, :math:`\delta V^{Quad}`. It is
assumed that no screening occurs in the grid. The value from the ordinary grid
is kept as the mean, and the estimate becomes the deviations.

Reference Grid
--------------

Evaluate the quadrature a second time on a finer grid, giving
:math:`E^{Quad}_{\mathrm{ref}}` (or :math:`V^{Quad}_{\mathrm{ref}}`), and use
the difference, :math:`d = E^{Quad}_{\mathrm{ref}} - E^{Quad}`, as an estimate
of :math:`-\delta E^{Quad}`. Since the direction of the error is known, the
deviation on the side of the reference value is :math:`|d|`, and the deviation
on the other side can be zero or a fraction of :math:`|d|`. This mirrors the
``Analytic Error`` ERI error model and its benchmark ERIs. It is an estimate,
not a bound, and costs one evaluation on the finer grid, including the density
and AOs on that grid.

The strategy is implemented as an error model, ``Reference Grid Error`` (see
`XC UQ Implementation`_) that satisfies the ``XCUncertain<Quadrature<T>>``
property type. It has a submodule, ``Reference XC Module``, which satisfies
``Quadrature<T>``. The module gets the original grid as an input and thus
should be able to determine a finer grid from it. This requires that the grid
has identifiers for its radial and angular quadratures. The module then passes
the inputs and the finer grid to the ``Reference XC Module`` submodule and
returns the difference. This requires that the ``Reference Grid Error`` module
also have a submodule for generating grids.

Grid Extrapolation
------------------

Evaluate on a sequence of grids of increasing size and estimate the remaining
error from how quickly the values converge. This is more expensive than a
single reference grid, but gives an estimate of the error of the finest grid
rather than the difference between two grids. It is a tool for choosing the
grid and for validating the reference-grid strategy. It could be implemented
as a series of successive calls to the ``Reference Grid Error`` module, but it
is not yet implemented.

Verified Cell Enclosure
-----------------------

Partition space into regions, :math:`\Omega_i`, of known volume, and let
:math:`\rho\left(\Omega_i\right)` be an enclosure of every value the density
takes in :math:`\Omega_i`. Then

.. math::

   \tilde{E}^{XC} \in \sum_i \left|\Omega_i\right|
                f\left[\rho\left(\Omega_i\right)\right]
            + \int_{\text{outside}} f\left[\rho\left(\vec{r}\right)\right]
              d\vec{r},

evaluated with UQ types, where the last term covers space outside all of the
regions and must be bounded separately, *e.g.*, from a bound on the density
far from the molecule. The enclosures :math:`\rho\left(\Omega_i\right)` are
computed by this error model itself, not taken from the inputs of
:ref:`xc_inputs_uq_design`, which only carry the uncertainty of the values at
the grid points; taking them from the inputs would count the variation of the
density across each region twice. This converges more slowly than traditional
quadrature and is thus expected to give a wider bound.

A Priori Estimates
------------------

The error of the radial and angular quadratures can be estimated from their
known properties, *e.g.*, the highest degree of spherical harmonic an angular
grid integrates exactly, together with how quickly the integrand's expansion
in spherical harmonics decays around each atom. This is cheap, but needs
analysis per quadrature, and is the least developed of the strategies.

Summary
-------

.. list-table::
   :header-rows: 1
   :widths: 25 20 25 30

   * - Strategy
     - Bound or estimate
     - Extra cost
     - Notes
   * - Reference grid
     - estimate
     - 1 evaluation on a finer grid
     - direction of the error is known
   * - Grid extrapolation
     - estimate
     - several evaluations
     - for validation and choosing grids
   * - Verified cell enclosure
     - bound
     - 1 evaluation with region enclosures
     - wide; LDA only as described
   * - A priori estimates
     - estimate
     - small
     - needs analysis per quadrature

**********************
Further Considerations
**********************

- **Rounding in the functional.** Evaluating :math:`f`, :math:`v`, and the sums
  over grid points also rounds. This is only accounted for if the arithmetic of
  the UQ type rounds outward. It is expected to be small compared with the
  other sources.
- **Thresholds in the functionals.** The XC modules need a policy for points
  whose uncertain density reaches below a density threshold, or zero. Pinning
  the lower bound of the density at zero keeps the density sound but does not
  bound :math:`\partial v/\partial\rho`. Dropping such points, as LibXC does,
  must itself be counted as an error, and is part of the error of the grid.
- **Uncertain grids.** ``GridFromIntegratorXX`` can make grids with uncertain
  points and weights. Uncertain points make the inputs uncertain and are
  covered in :ref:`xc_inputs_uq_design`. Uncertain weights enter the sums
  directly, so their effect follows from the UQ arithmetic in the XC modules.
- **GGAs and beyond.** The inputs grow (see :ref:`xc_inputs_uq_design`), but the
  split between the inputs and the error of the grid is unchanged.
