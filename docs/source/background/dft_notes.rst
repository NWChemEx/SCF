.. Copyright 2025 NWChemEx-Project
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

#########
DFT Notes
#########

.. |e_xc| replace:: :math:`E^{XC}\left[\rho\left(\vec{r}\right)\right]`
.. |v_xc| replace:: :math:`V^{XC}\left[\rho\left(\vec{r}\right)\right]`
.. |rho| replace:: :math:`\rho\left(\vec{r}\right)`
.. |rho_i| replace::  :math:`\rho_i`
.. |drho| replace:: :math:`\bigtriangledown\rho\left(\vec{r}\right)`

Algorithmically, DFT will differ from Hartee-Fock primarily in the fact that we
need to compute two additional quantities:

1. the exchange-correlation (XC) energy, |e_xc|, and
2. the XC potential, |v_xc|,

both of which are functionals of the electron density of the system, |rho|.
By definition, |v_xc| is:

.. math::

   \newcommand{\density}{\rho\left(\vec{r}\right)}
   \newcommand{\exc}[1]{E^{XC}\left[#1\right]}
   \newcommand{\vxc}[1]{V^{XC}\left[#1\right]}

   \vxc{\density} \equiv
       \frac{\partial \exc{\density}}{\partial \density}

meaning |e_xc| is given by:

.. math::

   \exc{\density} = \int \vxc{\density}\density d\vec{r}.

Conceptually, |e_xc| is the XC contribution to the electronic energy and
|v_xc| is the potential needed to make the non-interacting system have the same
density as the real system.

To make progress, we introduce an ansantze for |e_xc|, namely we
introduce a quantity :math:`\epsilon` called the XC energy density which is the
XC energy per electron. The form of :math:`\epsilon` determines the "flavor" of
DFT, e.g., the difference between say PBE and BLYP is the form of
:math:`\epsilon`. Whereas |v_xc| is a functional of |rho| alone,
:math:`\epsilon` will in general depend on not only |rho|, but additional
properties of the system.

In terms of :math:`\epsilon`, |e_xc| is now written:

.. math::

   \newcommand{\edensity}[1]{\epsilon\left[#1\right]}

   \exc{\density} = \int \edensity{\density,\cdots}\density d\vec{r}.

XC functionals are typically classified by the parameters that :math:`\epsilon`
depends on. For example, :math:`\epsilon` depends on

- only |rho| in the local density approximation (LDA),
- |rho| and |drho| (the gradient of |rho|) in the generalized gradient
  approximation (GGA), and
- |rho|, |drho|, the Laplacian of |rho|, and the kinetic energy density in meta
  GGAs.

For an LDA we have:

.. math::

   \vxc{\density} =&
      \frac{\partial \exc{\density}}{\partial \density}\\
      =& \edensity{\density} +
         \density\frac{\partial \edensity{\density}}{\partial \density}

*******************
Introduction of AOs
*******************

In Kohn-Sham DFT we express the density in terms of atomic orbitals.

.. math::
   \newcommand{\bf}[1]{\phi_{#1}\left(\vec{r}\right)}

   \density = \sum_{\mu\nu} \bf{\mu}P_{\mu\nu}\bf{\nu}

where :math:`P_{\mu\nu}` is the :math:`\mu\nu`-th element of the atomic
density matrix. Inserting this expression into the above equation for |e_xc|
results in:

.. math::

   \exc{\density} = \sum_{\mu\nu} P_{\mu\nu}
      \int \bf{\mu}\edensity{\density, \cdots}\bf{\nu} d\vec{r}

and in the LDA |v_xc| looks like:

.. math::

   \vxc{\density} = \edensity{\density} +
         \sum_{\mu\nu}P_{\mu\nu}\bf{\mu}\bf{\nu}
         \frac{\partial \edensity{\density}}{\partial \density}

Usually we do not want |v_xc| in real space, but rather in AO space.

Analytic solutions for the above integrals are not known and |e_xc| and |v_xc|
must be evaluated by quadrature.

**************************
Introduction of Quadrature
**************************

In solving an integral by quadrature, we make the following approximation:

.. math::

   \int f(\vec{r}) d\vec{r} \approx \sum_{i=1}^{N_g} w_i f(\vec{r_i})

That is we define a quadrature rule :math:`\mathcal{Q}` which is a set of
:math:`N_g` pairs of the form :math:`\lbrace w_i, \vec{r}_i\rbrace`. Here,
:math:`w_i` and :math:`\vec{r_i}` are respectively the quadrature weight and
real-space location of the :math:`i`-th grid point.

At this point it is helpful to define:

.. math::

   \newcommand{\densityg}[1]{\rho_{#1}}
   \newcommand{\bfg}[1]{\phi_{#1}}

   \densityg{i}\equiv&\rho\left(\vec{r_i}\right)\\
   \bfg{\mu i}\equiv&\phi_{\mu}\left(\vec{r_i}\right)

which respectively are the values of the density and the :math:`\mu`-th AO
evaluated at the :math:`i`-th grid point. |rho_i| is then given by:

.. math::



   \densityg{i} =& \sum_{\mu\nu} \bfg{\mu i}P_{\mu \nu}\bfg{\nu i}\\
                =& \sum_{\mu}\bfg{\mu i}X_{\mu i}

where in the second line we defined the common intermediate (the collocation
matrix):

.. math::

   X_{\mu i} = \sum_{\nu} P_{\mu\nu}\bfg{\nu i}

Using :math:`\mathcal{Q}`, |e_xc| becomes:

.. math::

   \exc{\density{}} =& \sum_{\mu\nu} P_{\mu\nu}
      \sum_i^{N_g} w_i\bfg{\mu i}\edensity{\densityg{i}, \cdots}\bfg{\nu i}\\
       =&  \sum_i^{N_g} w_i \edensity{\densityg{i}, \cdots}\densityg{i}
