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

.. _qcup_scf_design:

##################
Designing QCUP-SCF
##################

This page describes how the uncertainty-quantified SCF (QCUP-SCF) is exposed
to users. QCUP-SCF is a PluginPlay driver, whose design is described in the
`PluginPlay documentation`_. This page covers the parts specific to QCUP-SCF.
These parts span SCF, SimDE, Integrals, and ChemCache.

.. _PluginPlay documentation:
   https://nwchemex.github.io/PluginPlay/developer/design/driver.html

*****************
What is QCUP-SCF?
*****************

QCUP-SCF is a Hartree-Fock calculation in which numerical approximations, most
notably the screening of electron repulsion integrals (ERIs), are tracked as
uncertainties. Instead of ``double`` values, the calculation uses an uncertain
type, *e.g.*, intervals or affine arithmetic. The result is an
energy (and orbitals) together with a bound or estimate of its error.

QCUP-SCF runs the same SCF loop as an ordinary SCF calculation. It differs in
three ways:

#. **Which uncertain type is used.** This is chosen once and applies to the
   whole calculation.
#. **How uncertainty enters.** Error models attach uncertainty to quantities
   that are computed approximately. Today only the ERIs have an error model,
   and there are several to choose from.
#. **How uncertainty leaves.** After convergence, uncertainty is attached to
   the orbitals (from the uncertainty in the Fock matrix and from incomplete
   convergence) and to the energy (from incomplete convergence).

*****************
Statement of Need
*****************

QCUP-SCF is configured today by hand-written helper functions in the
application that runs the calculations. They make about twenty
``change_input`` and ``change_submod`` calls on modules from four plugins.
The table below lists them. In addition, the uncertainty attached after
convergence is hard-coded in ``scf_loop.cpp``. It is applied whenever the Fock
matrix holds an uncertain type, and it can not be changed or turned off.

.. list-table:: How QCUP-SCF is configured today
   :header-rows: 1
   :widths: 30 35 35

   * - Setting
     - Module(s)
     - Mechanism
   * - Uncertain type
     - ``Evaluate 2-Index BraKet``, ``Evaluate 4-Index BraKet``, ``Overlap``,
       ``Kinetic``, ``Nuclear`` (Integrals)
     - ``"UQ Type"`` input on each module
   * - Uncertain type and Taylor-model order
     - ``UQ Initializer`` (Integrals), used by ``UQ Atom Symm Blocked Driver``
     - ``"UQ Type"`` and ``"Order"`` inputs
   * - Uncertain SAD guess
     - ``sto-3g atomic density matrix`` (ChemCache)
     - ``"With UQ?"`` input
   * - Eigensolver that supports uncertain types
     - ``Loop``, ``Diagonalization Fock update`` (SCF)
     - ``"Diagonalizer"`` submodule set to ``Generalized eigensolve``
   * - ERIs with an error model
     - ``Four center J builder``, ``Four center K builder`` (Integrals)
     - ``"Four-center ERI"`` submodule set to
       ``UQ Atom Symm Blocked Driver``
   * - ERI error model
     - ``UQ Atom Symm Blocked Driver``
     - ``"ERI Error"`` submodule set to ``Analytic Error`` or
       ``Primitive Error Model``
   * - Error estimate of the primitive model
     - ``Primitive Error Model``
     - ``"Error estimate"`` input
   * - Signed analytic error
     - ``UQ Atom Symm Blocked Driver``
     - ``"Signed Error"`` input
   * - Averaging of the error
     - ``UQ Atom Symm Blocked Driver``
     - ``"Mean Type"`` input
   * - Reference ERIs for the analytic error model
     - ``Analytic Error``
     - ``"ERI4s"`` submodule set to a copy of ``ERI4`` with a threshold of
       :math:`10^{-25}`
   * - Screening threshold
     - ``ERI4``
     - ``"Threshold"`` input
   * - Memoization
     - every module except ``ERI4``
     - ``turn_off_memoization()``
   * - Orbital and energy uncertainty after convergence
     - ``Loop``
     - hard-coded, applied when the Fock matrix is uncertain

This has the problems described in the PluginPlay driver design: there is no
single place for the user to make choices, the choices are not first-class,
the ``ModuleManager``'s modules are changed in place, and every application
must repeat the wiring.

***********************
QCUP-SCF Considerations
***********************

#. One place for the user

   - The user selects the uncertain type and every error model through one
     object, without knowing where in the call graph each choice is used.

#. Every error model is a choice

   - The ERI error model is already a submodule. The uncertainty attached after
     convergence must become choices too. For now each has one implementation
     plus "no uncertainty."

#. No second SCF loop

   - QCUP-SCF must reuse ``Loop`` rather than duplicate it, *e.g.*, as a
     separate ``qcup_scf_loop.cpp``.

#. Ordinary SCF is unaffected

   - Loading the QCUP-SCF wiring must not change ordinary SCF calculations, and
     both must be usable in the same session.

#. Each shared module does one thing

   - ``UQ Initializer`` works because it only answers "how do I make an
     uncertain value?" A general "settings" module that many modules read is
     avoided, because a new setting may not be read by the modules it is meant
     for, and two settings may conflict without anything reporting it.

#. Only sources of uncertainty need to know the uncertain type

   - Most modules learn the type from their inputs. For example,
     ``Eigen Solve`` dispatches on the element type of the matrix it is given.
     Only the modules that create uncertain values from ``double`` values need
     to be told which type to use.

***********************
Alternatives Considered
***********************

Several ways of propagating the user's choices were considered. They are
described, with examples, in the PluginPlay driver design. In summary:

- **Hierarchical sub-drivers**, where each module forwards the settings of the
  modules below it, add inputs to modules that do not use them. They also
  require copying submodules on every call.
- **A shared settings module** violates "each shared module does one thing."
- **Promoted inputs** and **a wiring function** change the ``ModuleManager``'s
  modules in place.
- **A module that rewires itself**, *e.g.*, a ``qcup_scf_loop.cpp`` that
  selects ``Loop``'s submodules, must copy its call graph by hand on every
  call.
- **A separate recipe class** would add a second concept next to modules, with
  its own API for options.

The chosen design is a PluginPlay driver, ``QCUP SCF``, that wraps
``SCF Driver``.

*****************
QCUP-SCF Strategy
*****************

The strategy has three parts:

#. the uncertainty attached after convergence becomes submodules of ``Loop``,
#. the modules that create uncertainty get the uncertain type from
   ``UQ Initializer``, and
#. SCF provides a ``QCUP SCF`` driver that wires everything together.

Uncertainty After Convergence
=============================

``scf_loop.cpp`` currently does the following after convergence, but only when
the Fock matrix holds an uncertain type:

- diagonalizes the converged Fock matrix again and attaches first-order
  eigenvector uncertainty from the uncertainty in the Fock matrix
  (``attach_eigenvector_uncertainty``);
- inflates each orbital coefficient by the magnitude of the final change in the
  density matrix (``inflate_uncertainty_from``); and
- inflates the total energy by the magnitude of the final change in the energy
  (``inflate_uncertainty``).

These become two new property types. They are declared in SCF, in
``scf/property_types.hpp``, until their APIs have settled, and may move to
SimDE later:

.. list-table::
   :header-rows: 1
   :widths: 25 40 35

   * - Property type
     - Inputs
     - Result
   * - ``MOUncertainty``
     - Fock matrix :math:`\mathbf{F}`, overlap :math:`\mathbf{S}`, orbital
       energies :math:`\boldsymbol{\epsilon}`, orbital coefficients
       :math:`\mathbf{C}`, final density change :math:`\Delta\mathbf{P}`
     - :math:`\mathbf{C}` with uncertainty attached
   * - ``ConvergenceUncertainty``
     - energy :math:`E`, final energy change :math:`\Delta E`
     - :math:`E` with uncertainty attached

SCF provides two modules for each:

- ``"No MO uncertainty"`` and ``"No convergence uncertainty"`` return their
  input unchanged. They are the defaults, so ordinary SCF behaves as it does
  today.
- ``"First-order MO uncertainty"`` wraps the current eigenvector correction and
  the :math:`\Delta\mathbf{P}` inflation. ``"Energy change uncertainty"`` wraps
  the current :math:`\Delta E` inflation. The helper functions in
  ``eigen_solver/eigenvector_uncertainty.hpp`` and
  ``eigen_solver/inflate_uncertainty.hpp`` are reused unchanged.

``Loop`` gains the submodules ``"MO uncertainty"`` and
``"Convergence uncertainty"`` and always calls them after convergence.
``check_if_uqtype`` is removed, so QCUP-SCF is a different wiring of ``Loop``
and nothing more.

The first-order correction needs the eigenpairs of the converged Fock matrix,
which ``Loop`` does not otherwise have: the orbitals of the last iteration are
the eigenvectors of the previous, DIIS-extrapolated, Fock matrix. ``Loop``
therefore diagonalizes the converged Fock matrix once more and passes its
eigenpairs to ``MOUncertainty``. For ordinary SCF this costs one extra
diagonalization, and the returned orbitals are the eigenvectors of the
converged Fock matrix. The energy is unchanged.

The eigenvector correction stays a one-shot correction applied after
convergence. It must not be fed back into the SCF iterations, where it would
compound. :math:`\mathbf{S}` is passed to ``MOUncertainty`` so that a future
model can include the overlap's metric term. The current model omits that term
because the overlap is exact.

Where the Uncertain Type Enters
===============================

The modules that create uncertain values are:

- the integral modules (``ERI4`` through the ``Evaluate 2-Index BraKet`` and
  ``Evaluate 4-Index BraKet`` modules, ``Overlap``, ``Kinetic``, and
  ``Nuclear``),
- the SAD density (``sto-3g atomic density matrix``), and
- the error models, through ``UQ Atom Symm Blocked Driver``.

The one-electron integrals and the SAD density only need the type because
WTF does not convert ``double`` values to uncertain values in mixed
arithmetic. The overlap, for example, is exact. Its uncertain values have a
radius of zero.

Each of these modules gets a ``"UQ Initializer"`` submodule, replacing the
``"UQ Type"`` and ``"With UQ?"`` inputs. All of them are bound to the same
``UQ Initializer`` module by default, so the ``QCUP SCF`` driver sets the type
in one place. A module can opt out by binding a copy of ``UQ Initializer``,
*e.g.*, one that produces ``double`` values.

- ``UQKind`` gains ``none``, for ``double`` values, and ``UQ Initializer``'s
  ``"UQ Type"`` defaults to ``"none"``, so ordinary calculations are
  unchanged. ``simde::visit_uq_kind`` turns the kind of a factory into the
  type a module fills its buffer with.
- The ERIs are exact; their uncertainty is added by an ERI error model. Since
  all integral modules are the same class, Integrals registers a copy,
  ``No UQ Initializer``, whose type is always ``"none"``, and binds it to
  ``ERI2``, ``ERI3``, ``ERI4``, and the copies of ``ERI4``. The driver's
  overlay of ``UQ Initializer`` therefore never reaches them. The
  one-electron integrals and the ``Evaluate N-Index BraKet`` modules use
  ``UQ Initializer``.
- ``UQ Atom Symm Blocked Driver`` throws if its factory makes ``double``
  values, since it would silently drop the errors.

ChemCache and SCF can not depend on Integrals, so the ``UQInitializer``
property type and the ``UQFactory`` class move from Integrals to SimDE. The
``UQ Initializer`` module stays in Integrals. ChemCache adds Integrals as a
preferred provider and defaults the ``"UQ Initializer"`` submodule of its
atomic density matrices to it. Without Integrals, those modules need a
``UQ Initializer`` to be supplied. For the same reason, the
property type of the ERI error models, ``Uncertainty<ERI4>``, also moves to
SimDE, since the ``QCUP SCF`` driver in SCF has a submodule of that type.

If WTF learns to convert ``double`` values to uncertain values, the
one-electron integrals and the SAD density no longer need
``"UQ Initializer"``. This is a possible follow-up and not required.

The QCUP SCF Driver
===================

``QCUP SCF`` runs exactly what ``SCF Driver`` runs, and only wires it
differently. It is therefore a driver that satisfies ``simde::AOEnergy``, has
a submodule, ``"SCF"``, bound to ``SCF Driver``, and whose ``run_`` runs
that submodule with its own inputs. The ``"SCF"`` submodule is also what puts
``SCF Driver``'s call graph in the driver's call graph.

The driver lives in SCF, even though most of the modules it configures come
from Integrals. This follows the rules in the driver design for drivers and
plugins:

- every module the driver chooses from another plugin is a submodule of the
  driver, whose default comes from SCF's preferred providers, which already
  include Integrals and ChemCache;
- the defaults set in the driver's constructor only name SCF's own modules;
  and
- targets in other plugins use qualified keys, *e.g.*,
  ``{"integrals", "UQ Initializer"}``.

When Integrals is loaded, ``QCUP SCF`` is ready to run with the defaults
below. When it is not, the submodules whose defaults come from Integrals stay
unset, and running the driver fails with the list of those submodules, as it
would for any module.

Settings
--------

The user's choices are ordinary inputs and submodules of ``QCUP SCF``. Each is
forwarded to the module that uses it. All targets in Integrals are qualified
with ``"integrals"``.

.. list-table::
   :header-rows: 1
   :widths: 25 12 23 40

   * - Setting
     - Kind
     - Default
     - Forwarded to
   * - ``Diagonalizer``
     - submodule
     - ``Generalized eigensolve``
     - ``Loop`` and ``Diagonalization Fock update``: ``"Diagonalizer"``
   * - ``UQ Type``
     - input
     - ``"interval"``
     - ``UQ Initializer``: ``"UQ Type"``
   * - ``ERIs``
     - submodule
     - ``UQ Atom Symm Blocked Driver`` (Integrals)
     - ``Four center J builder`` and ``Four center K builder``:
       ``"Four-center ERI"``
   * - ``ERI Error Model``
     - submodule
     - ``Analytic Error`` (Integrals)
     - ``UQ Atom Symm Blocked Driver``: ``"ERI Error"``
   * - ``Benchmark ERIs``
     - submodule
     - ``Analytic Benchmark ERI4`` (Integrals)
     - ``Analytic Error``: ``"ERI4s"``
   * - ``ERI Error Estimate``
     - input
     - ``"Fine"``
     - ``Primitive Error Model``: ``"Error estimate"``
   * - ``Mean Type``
     - input
     - ``"none"``
     - ``UQ Atom Symm Blocked Driver``: ``"Mean Type"``
   * - ``ERI Threshold``
     - input
     - ``1.0E-6``
     - ``ERI4``: ``"Threshold"``
   * - ``MO Uncertainty Model``
     - submodule
     - ``First-order MO uncertainty``
     - ``Loop``: ``"MO uncertainty"``
   * - ``Convergence Uncertainty Model``
     - submodule
     - ``Energy change uncertainty``
     - ``Loop``: ``"Convergence uncertainty"``

The ERI module, the error model, and the benchmark ERIs used to be fixed
wiring or helper-function calls. They are submodules of the driver because
they come from Integrals. A side effect is that users can swap them too. The
eigensolver is a submodule because PluginPlay drivers have no fixed wiring:
every value a driver passes down is an input or submodule of the driver, with
the value as its default.

The error models and the uncertainty models are submodules, so the user's
choice is checked against their property types.

``ERI Error Estimate`` only matters when ``ERI Error Model`` is
``Primitive Error Model``. With any other error model,
``Primitive Error Model`` is not in the call graph, so the setting is inactive.
It is skipped and reported when the driver is locked, as described in the
driver design. Likewise, ``Benchmark ERIs`` is inactive unless the error
model is ``Analytic Error``.

Fixed Wiring
------------

- Memoization is turned off for every module except ``ERI4``. This is
  required for now.
- The cross-plugin defaults that SCF's plugin already sets
  (``Born-Oppenheimer approximation``, ``AO integral driver``, ``Overlap``,
  ``sto-3g SAD density``) apply to the copies as well, because the copies are
  made from the ``ModuleManager``'s modules after ``wire_defaults``.

QCUP-SCF needs no conditional wiring today. Choices that would need it, *e.g.*,
density fitting, can be added later with ``when`` blocks.

Example
-------

A sketch of the driver:

.. code-block:: c++

   DECLARE_DRIVER(QCUPSCF);

   DRIVER_CTOR(QCUPSCF) {
       description("Uncertainty-quantified SCF");
       satisfies_property_type<simde::AOEnergy>();
       add_submodule<simde::AOEnergy>("SCF").set_default("SCF Driver");

       add_submodule<simde::GeneralizedEigenSolve>("Diagonalizer")
         .set_default("Generalized eigensolve")
         .set_required_target("Loop", "Diagonalizer")
         .set_target("Diagonalization Fock update", "Diagonalizer");

       add_input<std::string>("UQ Type")
         .set_default(std::string("interval"))
         .set_required_target({"integrals", "UQ Initializer"}, "UQ Type");

       add_submodule<simde::ERI4>("ERIs")
         .set_required_target({"integrals", "Four center J builder"},
                              "Four-center ERI")
         .set_required_target({"integrals", "Four center K builder"},
                              "Four-center ERI");

       add_submodule<simde::Uncertainty<simde::ERI4>>("ERI Error Model")
         .set_required_target({"integrals", "UQ Atom Symm Blocked Driver"},
                              "ERI Error");

       add_input<std::string>("ERI Error Estimate")
         .set_default(std::string("Fine"))
         .set_target({"integrals", "Primitive Error Model"},
                     "Error estimate");

       using scf::property_types::MOUncertainty;
       add_submodule<MOUncertainty>("MO Uncertainty Model")
         .set_default("First-order MO uncertainty")
         .set_required_target("Loop", "MO uncertainty");

       // ... and so on for the remaining settings in the table

       wire_memoization(false, {{"integrals", "ERI4"}});
   }

   DRIVER_RUN(QCUPSCF) {
       const auto& [aos, sys] = simde::AOEnergy::unwrap_inputs(inputs);
       const auto e = submods.at("SCF").run_as<simde::AOEnergy>(aos, sys);
       auto rv      = results();
       return simde::AOEnergy::wrap_results(rv, e);
   }

SCF's plugin adds the driver and sets the defaults of the submodules that come
from Integrals:

.. code-block:: c++

   PLUGIN_CTOR(SCF) {
       // ... as today ...
       add_module<QCUPSCF>("QCUP SCF");

       default_submod_from("integrals", "QCUP SCF", "ERIs",
                           "UQ Atom Symm Blocked Driver");
       default_submod_from("integrals", "QCUP SCF", "ERI Error Model",
                           "Analytic Error");
       default_submod_from("integrals", "QCUP SCF", "Benchmark ERIs",
                           "Analytic Benchmark ERI4");
   }

The settings that every QCUP-SCF depends on use ``set_required_target``, so
a call graph in which they would do nothing is an error rather than a
silently ordinary SCF.

Using it from Python:

.. code-block:: python

   import integrals
   import pluginplay
   import scf
   import simde

   mm = pluginplay.ModuleManager()
   mm.load_plugin(scf.SCF())
   mm.load_plugin(integrals.Integrals())
   # ... and the other plugins SCF needs
   mm.wire_defaults()

   mm.change_input("QCUP SCF", "UQ Type", "thresholded affine")
   mm.change_submod("QCUP SCF", "ERI Error Model", "Primitive Error Model")
   mm.change_input("QCUP SCF", "ERI Error Estimate", "Fine")

   e = mm.run_as(simde.AOEnergy(), "QCUP SCF", aos, sys)
   e_plain = mm.run_as(simde.AOEnergy(), "SCF Driver", aos, sys)  # unaffected

``nwchemex.load_modules`` loads the same plugins and calls ``wire_defaults``,
so users of NWChemEx get ``QCUP SCF`` without doing anything.

To compare with another uncertain type, the user copies the driver, as for
any other module that has been run:

.. code-block:: python

   mm.copy_module("QCUP SCF", "QCUP SCF affine")
   mm.change_input("QCUP SCF affine", "UQ Type", "affine")

After ``mm.resolve_drivers("QCUP SCF")``, or a run,
``mm.at("QCUP SCF").inactive_settings()`` lists the settings which had no
effect, *e.g.*, ``ERI Error Estimate`` when the analytic error model is
chosen.

The calculations in an interaction energy can each use their own driver,
*e.g.*, by binding ``"QCUP SCF"`` as the SCF submodule of an
interaction-energy driver.

Mapping From Today's Configuration
==================================

Every setting in the table in `Statement of Need`_ has a home:

- the uncertain type becomes the ``UQ Type`` input of the driver, and the
  Taylor-model order is dropped (see `Further Considerations`_);
- the ``"UQ Type"`` and ``"With UQ?"`` inputs on the integral and SAD modules
  become ``"UQ Initializer"`` submodules (see
  `Where the Uncertain Type Enters`_);
- memoization becomes fixed wiring, and the eigensolver becomes the
  ``Diagonalizer`` submodule of the driver;
- the ERI module and the reference ERIs become the ``ERIs`` and
  ``Benchmark ERIs`` submodules of the driver, since they come from
  Integrals;
- the error model, its inputs, and the screening threshold become settings
  of the driver; and
- the uncertainty attached after convergence becomes the
  ``MO Uncertainty Model`` and ``Convergence Uncertainty Model`` submodules of
  the driver.

The application that runs the QCUP-SCF calculations then only calls
``change_input`` and ``change_submod`` on ``QCUP SCF``.

Work by Repository
==================

In the order the work must be done:

#. **PluginPlay**: add ``DriverBase`` and the driver macros, including
   copying a module together with every module below it, inactive-setting
   reporting, and Python support. The design document goes in the same PR.
#. **SimDE**: move ``UQInitializer``, ``UQFactory``, and ``Uncertainty`` from
   Integrals.
#. **Integrals**: use the moved property types; replace the ``"UQ Type"``
   inputs with ``"UQ Initializer"`` submodules; register
   ``No UQ Initializer`` and ``Analytic Benchmark ERI4``.
#. **ChemCache**: replace ``"With UQ?"`` with a ``"UQ Initializer"``
   submodule, defaulted from Integrals.
#. **SCF**: add the ``MOUncertainty`` and ``ConvergenceUncertainty`` property
   types, their modules, and the new ``Loop`` submodules; remove
   ``check_if_uqtype``; add the ``QCUP SCF`` driver and its defaults to SCF's
   plugin.

**********************
Further Considerations
**********************

- **Error-model features that are not on master.** The application that runs
  QCUP-SCF today uses a ``"Signed Error"`` input and the ``"Schwarz"`` and
  ``"SchwarzGF"`` error estimates. These are not on Integrals' master branch,
  where ``Primitive Error Model`` supports ``"Tolerance"``, ``"Coarse"``, and
  ``"Fine"``. ``QCUP SCF`` therefore has no ``Signed Error`` setting yet. It
  should be added when they are merged, subject to the next point.
- **Settings that depend on a submodule choice.** ``Signed Error`` is only
  valid with ``Analytic Error``, whose output is a signed difference.
  ``Primitive Error Model``'s output is an unsigned bound, so setting
  ``Signed Error`` with it gives wrong uncertainties. Since
  ``UQ Atom Symm Blocked Driver`` is in the call graph either way, the setting
  is not inactive and nothing reports the mistake. ``when`` blocks test the
  driver's inputs, not which module a submodule is bound to, so they can not
  express this either. The options are to have ``UQ Atom Symm Blocked Driver``
  check the combination itself, to register the signed variant as its own
  error model (*e.g.*, ``Signed Analytic Error``) and drop the input, or to
  extend ``when`` to test submodule bindings.
- **More error models.** New error models, *e.g.*, for the one-electron
  integrals or for DFT quadrature, are added as modules and selected with new
  forwarded submodules, without changing the driver's other settings.
- **Overlap uncertainty.** ``MOUncertainty`` omits the metric term of the
  generalized eigenproblem because the overlap is exact. If an error model is
  ever added to the overlap, ``"First-order MO uncertainty"`` must be extended.
- **Taylor models.** ``UQKind`` has no Taylor-model kind, and
  ``UQ Initializer`` no ``"Order"`` input, because Taylor models are not ready
  to be used yet. Adding them back means a kind and a factory that stores the
  order, an ``"Order"`` input on ``UQ Initializer``, and a ``Taylor Order``
  setting on ``QCUP SCF`` forwarded to it.
