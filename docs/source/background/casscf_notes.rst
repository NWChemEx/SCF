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

############
CASSCF Notes
############

Disclaimer I'm not 100% convinced that CASSCF belongs in this repo, but I'm
going to start the notes here for now.

********
Notation
********

- :math:`p,q,r,s` denote MOs belonging to any orbital space.
- :math:`i,j,k,l` denote inactive doubly occupied MOs.
- :math:`a,b,c,d` denote virtual MOs.
- :math:`v,w,x,y` denote MOs inside the active space.

*************************************
Parameterization of the Wave Function
*************************************

.. |ket_CASSCF| replace:: :math:`\ket{\mathbf{C},\mathbf{\kappa}}`

The CASSCF wave function, |ket_CASSCF|, is defined as:

.. math::

   \newcommand{\ketCASSCF}[2]{\ket{#1,#2}}
   \newcommand{\kvector}{\mathbf{\kappa}}
   \newcommand{\Cvector}{\mathbf{C}}

   \ketCASSCF{\Cvector}{\kvector} = e^{-\hat{\kappa}} \sum_I C_I \ket{I}

.. |kappa| replace:: :math:`\kvector`
.. |C| replace:: :math:`\Cvector`
.. |ketI| replace:: :math:`\ket{I}`

where |C| is a vector of configuration interaction (CI) coefficients such that
:math:`C_I` is the contribution of the :math:`I`-th multi-electron state
(usually Slater determinants, but sometimes configuration state functions; we
use "determinants" for brevity, but the majority of the discussion applies
equally to CSFs), |ketI|, to the wave function and the :math:`e^{-\hat{\kappa}}`
operator defines a rotation of the molecular orbitals (MO) in terms of a set of
rotation parameters, |kappa|, via:

.. math::
   \newcommand{\Epq}[1]{E^-_{#1}}

    \hat{\kappa} = \sum_{p>q} \kappa_{pq} \Epq{pq}.

.. |Epq| replace:: :math:`\Epq{pq}`

Here |Epq| is the antisymmetric singlet excitation operator.

.. |ket0| replace:: :math:`\ket{0}`

In practice, CASSCF is an iterative process and it is useful to
rewrite the current wave function in terms of our last best guess for the wave
function, |ket0|. In terms of |ket0| the CASSCF wave function is given by:

.. math::
   \newcommand{\cvector}{\mathbf{c}}

   \ketCASSCF{\cvector}{\kvector} = e^{-\hat{\kappa}}
      \frac{\ket{0} + \hat{P}\ket{\cvector}}
           {\sqrt{1 + \braket{\cvector \mid \hat{P} \mid \cvector}}}

.. |c| replace:: :math:`\cvector`
.. |C0| replace:: :math:`\mathbf{C^{(0)}}`
.. |Phat| replace:: :math:`\hat{P}`
.. |zeta0| replace:: :math:`\mathbf{\zeta^{(0)}}`

here:

.. math::
   \newcommand{\refC}[1]{C_{#1}^{(0)}}

   \ket{0} =& \sum_{I} \refC{I} \ket{I}\\
   \ket{\cvector} =& \sum_{I} c_I \ket{I}\\
   \hat{P} =& 1 - \ket{0}\bra{0}

and |ket0| is assumed normalized. Some notes:

- The expression comes from writing |C|, the new CI coefficients, as the sum of
  |C0| and |c|, where |C0| is the current CI vector and |c| contains the
  variational parameters. At convergence, |c| will be zero.
- In general |c| is not orthogonal to |C0|. To ensure we only update our CI
  coefficients with new information, we introduce the projection operator |Phat|
  that removes any component of |c| already contained in |C0|, i.e.,
  :math:`\braket{0 \mid \hat{P}\mid\cvector} = 0`.
- The determinants, |ketI|, are assumed orthornormal.
- The determinants are defined in terms of the current MOs, not the
  MOs obtained by applying |kappa|.

It is convenient to express the CASSCF equations in tensor form. To that end,
we think of |kappa| as a vector, and not a matrix as the notation would suggest
(this is easily accomplished by thinking of :math:`pq` as a single index and not
a pair of indices). We then define several vectors:

.. math::

   \newcommand{\newparams}{\mathbf{\zeta}}
   \newcommand{\oldparams}{\mathbf{\zeta^{(0)}}}
   \newcommand{\stepparams}{\mathbf{\lambda}}
   \newcommand{\supervector}[2]{\left(\begin{array}{c}
                     #1\\
                     #2
                   \end{array}\right)}

    \newparams  =& \oldparams + \mathbf{P}\stepparams\\
    \oldparams  =& \supervector{\Cvector^{(0)}}{\mathbf{0}}\\
    \stepparams =& \supervector{\cvector}{\kvector}

.. |newparams| replace:: :math:`\newparams`
.. |oldparams| replace:: :math:`\oldparams`
.. |stepparams| replace:: :math:`\stepparams`
.. |P| replace:: :math:`\mathbf{P}`

Here |newparams| are the parameters for the current wave function, |oldparams|
are the parameters for the previous wave function, and |stepparams| is the
parameter step needed to form the |newparams|. Some notes:

- The zero in |oldparams| comes from the fact that the current MOs are used to
  define |ket0|, meaning that the rotation to the current MOs is zero.
- The previous statement means that we zero out the rotation parameters each
  iterations, which is equivalent to setting a new reference state.
- |P| is the matrix representation of the |Phat|. |P| again ensures that only
  the component of |stepparams| orthogonal to |oldparams| is used to update the
  wave function parameters.

Since |P| is defined in terms of |ket0|, which in turn is defined in terms of
|oldparams|, |P| has a simple form:

.. math::

    \mathbf{P} =& \mathbf{1} - \mathbf{O}\\
               =& \mathbf{1} - \oldparams\oldparams^T

where :math:`\mathbf{O} = \oldparams\oldparams^T` is the projection onto the
current wave function.

******************
Problem Definition
******************

The CASSCF Energy, :math:`E`, and wave function are obtained by
minimizing the expectation value of the Hamiltonian with respect to both
|C| and |kappa|, i.e.,

.. math::

   E = \min_{\cvector,\kvector}
       \frac{\braket{\cvector,\kvector\mid\hat{H}\mid\cvector,\kvector}}
            {\braket{\cvector,\kvector \mid \cvector,\kvector}}

The gradient of this expression is:

.. math::

   \newcommand{\ederiv}[1]{E^{(#1)}}
   \newcommand{\egrad}{\mathbf{\ederiv{1}}}
   \newcommand{\cgrad}{\mathbf{{^{c}}\ederiv{1}}}
   \newcommand{\ograd}{\mathbf{{^{o}}\ederiv{1}}}

   \egrad          = \supervector{\cgrad}{\ograd}

where:

.. math::

   \ederiv{1}_{I}  =& 2\braket{I\mid\hat{H}\mid 0} - 2\refC{I} E\\
   \ederiv{1}_{pq} =& \braket{0\mid\left[\Epq{pq},\hat{H}\right]\mid 0}.

In the above equation superscripts "c" and "o" respectively denote the gradient
of the energy with respect to the CI coefficients and the MO rotation
parameters. Whether the derivative is with respect to the CI coefficients or the
MO rotation parameters is obvious in element-wise equations and superscripts are
suppressed to avoid clutter.

In practice, CASSCF is numerically more difficult to converge than traditional
SCF and second-order optimization methods, i.e., those that rely on both the
gradient and the Hessian, are usually used. To that end, the CASSCF Hessian is
given by:

.. math::
   \newcommand{\supermatrix}[4]{\left(\begin{array}{cc}
                     #1 & #2\\
                     #3 & #4
                   \end{array}\right)}
   \newcommand{\ccHess}{\mathbf{{^{cc}}\ederiv{2}}}
   \newcommand{\coHess}{\mathbf{{^{co}}\ederiv{2}}}
   \newcommand{\ocHess}{\mathbf{{^{oc}}\ederiv{2}}}
   \newcommand{\ooHess}{\mathbf{{^{oo}}\ederiv{2}}}
   \newcommand{\ehess}{\mathbf{\ederiv{2}}}

   \ehess = \supermatrix{\ccHess}{\coHess}{\ocHess}{\ooHess}

.. |coHess| replace:: :math:`\coHess`
.. |ocHess| replace:: :math:`\ocHess`

Since the order of derivatives does not matter, |coHess| is simply the transpose
of |ocHess| and the Hessian has three unique blocks. The elements of the unique
blocks are:

.. math::

   \ederiv{2}_{IJ} =& 2\braket{I\mid\hat{H} - E \mid J} -
                      \refC{I}\ederiv{1}{J} - \refC{J}\ederiv{1}_{I}\\
   \ederiv{2}_{I,pq} =& 2\braket{I\mid\left[\Epq{pq},\hat{H}\right]\mid 0} -
                        2\refC{I}\ederiv{1}_{pq}\\
   \ederiv{2}_{pq,rs} =& \frac{1}{2}\left(1+\hat{P}_{pq,rs}\right)
                         \braket{0\mid\left[\Epq{pq},
                         \left[\Epq{rs},\hat{H}\right]\right]\mid 0}

.. |Ppqrs| replace:: :math:`\hat{P}_{pq,rs}`

where |Ppqrs| permutes the pair index :math:`pq` with the pair index :math:`rs`.

Setting the gradients equal to zero yields the CASSCF stationary condition:

.. math::

   0 =& \braket{0\mid\left[\Epq{pq},\hat{H}\right]\mid 0}\\
   0 =& 2\braket{I\mid\hat{H}\mid 0} - 2\refC{I} E,

also known as the generalized Brillouin theorem. The goal of CASSCF is to find
|newparams| such that the generalized Brillouin theorem is satisfied.

***************************************************************
Solving for :math:`\mathbf{\lambda}`: Newton Eigenvector Method
***************************************************************

To solve the CASSCF equations, we note that the optimum parameters are obtained
when :math:`\stepparams = \textbf{0}` (equivalently when |newparams| is the
same as |oldparams|). We thus Taylor series expand the energy
with respect to |stepparams| about the point :math:`\stepparams = \textbf{0}`
(note that this is not the same as expanding about |oldparams|). The Taylor
series expansion is:

.. math::

   Q(\stepparams) = E + \egrad^T \stepparams +
                   \frac{1}{2} \stepparams^T \ehess \stepparams.

Introducing the augmented Hessian:

.. math::

    \newcommand{\augHess}{\mathbf{G}}

    \augHess = \ehess + \oldparams \egrad^T + \egrad \oldparams^T

the Talor series takes the simpler form:

.. math::

   Q(\newparams) = E + \frac{1}{2} \newparams^T \augHess \newparams.

(n.b. I don't really see why this works, but I'm just going to accept it for
now). The optimization will be carried out subject to two constraints. The first
is that |newparams| is normalized with respect to |oldparams|, i.e.,

.. math::

   \newparams^T\mathbf{O}\newparams = 1.

The second constraint is that the the new parameters are "close" to the old
parameters because we don't trust the Taylor series too far from the expansion
point. Defining "close" via the constant :math:`h`, this constraint is expressed
as:

.. math::

   \newparams^T\mathbf{P}\newparams = h^2.

(note as a projection operator |P| satisfies :math:`\mathbf{P}^2 = \mathbf{P}`).
Here the constraint is only defined in terms of the component of |newparams|
that is not contained in |oldparams| since that is the only part that will move
our new parameters. To solve the constrained optimization problem, we introduce
two Lagrange multipliers, :math:`\mu` and :math:`\nu`, and define the
Lagrangian:

.. math::

   L(\newparams,\mu,\nu) = E + \frac{1}{2} \newparams^T \augHess \newparams -
      \frac{\mu}{2}\left(\newparams^T\mathbf{O}\newparams - 1\right) -
      \frac{\nu}{2}\left(\newparams^T\mathbf{P}\newparams - h^2\right).

Setting the derivative of the Lagrangian with respect to |newparams| equal to
zero we obtain:

.. math::

   \augHess \newparams = \mu \mathbf{O}\newparams + \nu \mathbf{P}\newparams.

This is a generalized eigenvalue equation. To make this more apparent we define
a matrix:

.. math::

   \newcommand{\Tmatrix}{\mathbf{T}\left(\alpha\right)}

    \Tmatrix = \mathbf{O} - \alpha^2\mathbf{P}.

.. |Tmatrix| replace:: :math:`\Tmatrix`

where :math:`\alpha^2=\frac{\mu}{\nu}`. With this definition the generalized
eigenvalue equation becomes:

.. math::

   \augHess \newparams = \nu \Tmatrix \newparams.

To recast this as a standard eigenvalue equation we factor |Tmatrix| as:

.. math::

   \Tmatrix = \Tmatrix^{1/2} \Tmatrix^{1/2}

and multiply both sides by :math:`\Tmatrix^{-1/2}` to obtain:

.. math::

   \Tmatrix^{-1/2} \augHess \Tmatrix^{-1/2}
      \left(\Tmatrix^{1/2}\newparams\right)
      = \nu \left(\Tmatrix^{1/2}\newparams\right).

Usingthe fact that:

.. math::

   \Tmatrix^{n} = \mathbf{O} - \alpha^{2n}\mathbf{P}

we define:

.. math::

    \newcommand{\yvector}{\mathbf{\xi}\left(\alpha\right)}

    \augHess\left(\alpha\right) =&
       \alpha^2\Tmatrix^{-1/2} \augHess \Tmatrix^{-1/2}\\
                                =& \ehess + \alpha\oldparams\egrad^T +
                                   \alpha\egrad \oldparams^T\\
    \yvector                    =& \Tmatrix^{1/2}\newparams.

.. |yvector| replace:: :math:`\yvector`

Resulting in the standard eigenvalue equation:

.. math::

   \augHess\left(\alpha\right) \yvector = \mu \yvector.

Solving the eigenvalue equation for the lowest eigenvalue, :math:`\mu`, yields
the update for the ground state wave function (in general the eigenvector for
the :math:`n`-th lowest eigenvalue yields the update for the :math:`n`-th
state). From the definition of |yvector| we have:

.. math::
   \newparams   =& \Tmatrix^{-1/2}\yvector\\
                =& \oldparams + \alpha^{-1}\mathbf{P}\yvector\\
                =& \oldparams + \mathbf{P}\stepparams\\
    \stepparams =& \alpha^{-1}\yvector

*****************
Working Equations
*****************

.. math::

   \newcommand{\iFock}[1]{F^I_{#1}}
   \newcommand{\aFock}[1]{F^A_{#1}}
   \newcommand{\eri}[2]{\left(#1\mid#2\right)}

   \iFock{pq} =& h_{pq} + \sum_{i}\left[2\eri{pq}{ii} - \eri{pi}{qi}\right]\\
   \aFock{pq} =& \sum_{uv} \gamma_{uv}
                  \left[\eri{pq}{uv} - \frac{1}{2}\eri{pu}{qv}\right]\\

The gradient is given in terms of the generalized Fock matrix:

.. math::

   \ederiv{1}_{pq} = 2\left(F_{pq} - F_{qp}\right)\\
   F_{ip} = 2\iFock{pi} + 2\aFock{pi}\\
   F_{up} = \sum_{v}\iFock{pv}\gamma_{uv} + Q_{up}\\
   F_{ap} = \mathbf{0} \\
   Q_{up} = \sum_{vwx} \Gamma_{uvwx}\eri{pv}{wx}

*********
Algorithm
*********

- Input: set of MO coefficients and a CI vector.
- Transform from AOs to current MOs
- For the energy of the :math:`n`-the state solve for the :math:`n`-th
  |yvector|.
- Compute |newparams| via :math:`\oldparams + \alpha^{-1}\mathbf{P}\yvector`.

Before leaving this section, it is useful to note that the CASSCF Hessian
contains terms that depend on the gradients. Such contributions will be zero at
convergence. It is thus useful to write the Hessian as the sum of a matrix that
will not be zero at convergence and the terms that depend on the gradient.

We define the matrix:

.. math::

   \newcommand{\Khess}[1]{K^{(2)}_{#1}}
   \newcommand{\Kmatrix}{\mathbf{\Khess{}}}
   \newcommand{\ccK}{{^{cc}}\Kmatrix}
   \newcommand{\coK}{{^{co}}\Kmatrix}
   \newcommand{\ocK}{{^{oc}}\Kmatrix}
   \newcommand{\ooK}{{^{oo}}\Kmatrix}

   \Kmatrix = \supermatrix{\ccK}{\coK}{\ocK}{\ooK}

.. |Kmatrix| replace:: :math:`\Kmatrix`

to be the part of the Hessian that does NOT depend on the gradient. Like the
Hessian, |Kmatrix| has three unique blocks because the off-diagonals are related
by a transpose. The elements of the |Kmatrix| are:

..  math::

    \Khess{IJ}    =& 2\braket{I\mid\hat{H} - E \mid J}\\
    \Khess{I,pq}  =& 2\braket{I\mid\left[\Epq{pq},\hat{H}\right]\mid 0}\\
    \Khess{pq,rs} =& \frac{1}{2}\left(1+\hat{P}_{pq,rs}\right)
                     \braket{0\mid\left[\Epq{pq},
                     \left[\Epq{rs},\hat{H}\right]\right]\mid 0}\\

In terms |Kmatrix|, the Hessian can be written as:

.. math::
   \newcommand{\ebar}{\mathbf{\bar{E}^{(1)}}}

   \ehess = \Kmatrix - \ebar\oldparams^T - \oldparams\ebar^T

where:

.. math::

   \ebar = \supervector{\cgrad}{2\ograd}.
