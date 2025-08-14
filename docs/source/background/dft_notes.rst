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
.. |edensity| replace:: :math:`f\left[\rho\left(\vec{r}\right), \cdots\right]`
.. |eparticle| replace:: :math:`\epsilon\left[\rho\left(\vec{r}\right), \cdots\right]`
.. |rho| replace:: :math:`\rho\left(\vec{r}\right)`
.. |rho_i| replace::  :math:`\rho_i`
.. |drho| replace:: :math:`\left|\bigtriangledown\rho\left(\vec{r}\right)\right|^2`
.. |dr| replace:: :math:`d\vec{r}`

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

or rearranging for |e_xc|:

.. math::

   \exc{\density} = \int \vxc{\density}\density d\vec{r}.

Conceptually, |e_xc| is the XC contribution to the electronic energy and
|v_xc| is the potential needed to make the non-interacting system have the same
density as the real system. Unfortunately, we do not know the analytic form of
|e_xc| (or equivalently |v_xc|) and likely never will.

To make progress, we need to introduce an ansantze for either |e_xc| or |v_xc|.
To do this we define a quantity |edensity| called the XC energy density
which is the XC energy per infinitesimal volume |dr|. In terms of
|edensity|, |e_xc| is written as:

.. math::

   \newcommand{\edensity}[1]{f\left[#1\right]}
   \newcommand{\eparticle}[1]{\epsilon\left[#1\right]}

   \exc{\density} = \int \edensity{\density,\cdots} d\vec{r}.

Related to |edensity| is a quantity |eparticle| which is the XC energy per
unit particle. The exact relationship between the two is:

.. math::

   \edensity{\density,\cdots} = \eparticle{\density,\cdots}\density.

It should be noted that, unlike |v_xc| and |e_xc|, |edensity| and |eparticle|
will in general depend on additional functionals beyond the density.

.. warning::

   Colloquially speaking, "a DFT functional" (e.g., when someone says they
   used "the BLYP DFT functional") can refer to the analytic form of |edensity|
   or |eparticle|. Making matters worse, :math:`\epsilon` is commonly used to
   denote both quantities. When comparing equations it is critical to
   distinguish between these two quantities.

XC functionals are typically classified by the parameters that |edensity| or
|eparticle| depend on. More specifically dependence on:

- only |rho| defines the local density approximation (LDA),
- |rho| and |drho| (the square of the gradient of |rho|) defines the
  generalized gradient approximation (GGA), and
- |rho|, |drho|, the Laplacian of |rho|, and the kinetic energy density defines
  a  meta GGA.

For an LDA we have:

.. math::

   \vxc{\density} =& \frac{\partial \exc{\density}}{\partial \density}\\
                  =& \frac{\partial \edensity{\density}}{\partial \density}

*******************
Introduction of AOs
*******************

In Kohn-Sham DFT we solve the Kohn-Sham equations in an orbital basis that is
obtained as a linear combination of atomic orbitals (AOs). Assume that there
are :math:`N_b` AOs and let :math:`\phi_\mu\left(\vec{r}\right)` be the
:math:`\mu`-th AO. The equation for |e_xc| remains unchanged other than |rho|
is now:

.. math::
   \newcommand{\bf}[1]{\phi_{#1}\left(\vec{r}\right)}

   \density = \sum_{\mu}^{N_b}\sum_{\nu}^{N_b} \bf{\mu}P_{\mu\nu}\bf{\nu}

where :math:`P_{\mu\nu}` is the :math:`\mu\nu`-th element of the atomic
density matrix. In the AO basis set the :math:`\mu\nu`-th element of |v_xc| is:

.. math::

   V^{XC}_{\mu\nu} = \int \bf{\mu} \vxc{\density} \bf{\nu} d\vec{r}.

In the LDA this becomes:

.. math::

   V^{XC}_{\mu\nu} = \int \bf{\mu}
     \frac{\partial \edensity{\density}}{\partial \density} \bf{\nu} d\vec{r}.

For most DFT functionals, analytic solutions for the above integrals are not
known and |e_xc| and |v_xc| must be evaluated by quadrature.

**************************
Introduction of Quadrature
**************************

In solving an integral by quadrature, we make the following approximation:

.. math::

   \int g(\vec{r}) d\vec{r} \approx \sum_{i=1}^{N_g} w_i g(\vec{r_i}).

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
evaluated at the :math:`i`-th grid point. Similarly, we define:

.. math::

   \newcommand{\edensityg}[1]{f_{#1}}
   \newcommand{\dedensitygdrho}[1]{f_{#1}^{\left(\rho\right)}}

   \edensityg{i}\equiv&\edensity{\densityg{i}}\\
   \dedensitygdrho{i}\equiv&
     \left.
       \frac{\partial \edensity{\density{}}}
         {\partial \density{}}
     \right|_{\density{}=\densityg{i}}

which are the energy density, and the "derivative of the energy density with
respect to the density" evaluated at |rho_i|.

Using these quantities, |rho_i| is then given by:

.. math::

   \densityg{i} =& \sum_{\mu}^{N_b}
      \sum_{\nu}^{N_b} \bfg{\mu i}P_{\mu \nu}\bfg{\nu i}\\
                =& \sum_{\mu}^{N_b} \bfg{\mu i}X_{\mu i}

where in the second line we defined the common intermediate (the collocation
matrix):

.. math::

   X_{\mu i} = \sum_{\nu}^{N_b} P_{\mu\nu}\bfg{\nu i}

Using :math:`\mathcal{Q}`, |e_xc| becomes:

.. math::

   \exc{\density{}} = \sum_i^{N_g} w_i\edensityg{i}

and |v_xc| becomes:

.. math::

   V_{\mu\nu}^{XC} =&
     \sum_i^{N_g} w_i \bfg{\mu i} \dedensitygdrho{i} \bfg{\nu i}\\
                   =& \sum_i^{N_g} w_i \bfg{\mu i} Z_{\nu i}

where we defined the intermediate:

.. math::

   Z_{\mu i} =\dedensitygdrho{i} \bfg{\mu i}.

***********************
As a Sparse Map Problem
***********************

While the last sections have described DFT as a tensor problem it's usually not
solved as one.  DFT is not usually treated as a tensor problem because:

- Large tensors. Grids minimally use about 1000 grid points per atom (higher-
  quality grids tend to be order 10,000) and most AO basis sets have order 10
  basis functions per atom. Tensors like :math:`\phi_{\mu i}` then have
  minimally "10,000 times number of atoms squared" elements, meaning the tensor
  for 100 atoms already requires gigabytes of memory.
- Sparsity. Most DFT quantities are local. So if basis functions for a tensor
  element are spatially far a part, the element is usually close to zero.

To describe the sparsity we introduce sparse maps. Given two basis sets,
:math:`A` and :math:`B`, the sparse map :math:`L` maps each basis
function in :math:`A` to a subset of the basis functions in :math:`B`. Assume
we have some tensor with elements :math:`T_{ab}` where :math:`a` indexes basis
functions in :math:`A` and :math:`b` indexes basis functions in :math:`B`.
For a given value of :math:`a`, the non-zero elements of :math:`T_{ab}` are
those such that :math:`b` is in  :math:`L(a)`.

In DFT, we use atom-centered grids and AOs. It is therefore common to define
sparse maps :math:`L(A\rightarrow i)` and :math:`L(A\rightarrow \mu)` which
respectively map atom indices to grid points and atom indices to AOs. Using
these maps the equation for the density becomes:

.. math::

   \densityg{i_A} = \sum_{\mu_A} \bfg{\mu_A i_A}X_{\mu_A i_A}

where an index like :math:`i_A` is shorthand for restricting the value of
:math:`i` to those afforded by the sparse map :math:`L(A\rightarrow i)`.
Applying the same logic to the other DFT quantities:

.. math::

   X_{\mu_A i_A} =& \sum_{\nu_A} P_{\mu_A\nu_A}\bfg{\nu_A i_A}\\
   Z_{\mu_A i_A} =& \dedensitygdrho{i_A} \bfg{\mu_A i_A}\\
   V_{\mu_A\nu_A}^{XC} =& \sum_{i_A} w_{i_A}\bfg{\mu_A i_A} Z_{\nu_A i_A}.

Finally, the equation for |e_xc| becomes:

.. math::

   \exc{\density{}} = \sum_{A}\sum_{i_A} w_{i_A}\edensityg{i_A}

Of note, for a given grid we expect the number of grid points associated with
an atom to be roughly constant. Similarly, for a given AO basis set we expect
the number of AOs associated with an atom to also be roughly constant. This
means that cost to form all quantities will scale linearly with the number of
atoms.
