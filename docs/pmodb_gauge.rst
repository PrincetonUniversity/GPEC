.. _pmodb_gauge:

******************************************
Perturbed field strength (pmodb) and gauge
******************************************

This page states what the columns of ``gpec_pmodb_n#.out`` are, how they are
related, and which of them a field-strength consumer needs.  Verified identities
are proved with the FortSym native zero test by
``docs/derivations/pmodb_gauge/derive.py``; the verdicts are stored in
``docs/derivations/pmodb_gauge/fortsym_verdicts.json`` and the labels (G1,
G5.5, ...) refer to that file.

.. contents:: On this page
   :local:
   :depth: 1


Coordinates and displacement
============================

Assume smooth nested surfaces with B>0 and nonzero J, chi1 and q.
Fourier inversions below additionally exclude nonzero resonant sources.

GPEC works in the DCON chart :math:`(\psi,\theta,\zeta)` with angles in turns,
any straight-field-line Jacobian :math:`\mathcal{J}(\psi,\theta)` selected by
``jac_type``, :math:`\mathcal{J}B^\theta=\chi_1`,
:math:`\mathcal{J}B^\zeta=\chi_1 q` and Fourier factor
:math:`\exp[2\pi i(m\theta-n\zeta)]`.  The Euler-Lagrange solution supplies
:math:`\xi^\psi`, its radial derivative and the surface displacement
(Glasser 2016, Eqs. 9 and 18)

.. math::

   \xi_s = (\boldsymbol\xi\times\mathbf B)_\psi = \chi_1(q\xi^\theta-\xi^\zeta),

the in-surface displacement across :math:`\mathbf B` (``xms``;
``gpec_xclebsch`` writes :math:`\xi^\alpha=\xi_s/\chi_1`).  A field-aligned
displacement changes neither :math:`\xi^\psi` nor :math:`\xi_s` (G1.5).  The
parallel part is fixed in ``gpeq_contra`` by
:math:`\nabla\cdot\boldsymbol\xi=0` (G2), which gives ``xmt``
:math:`=\mathcal J\xi^\theta` and ``xmz`` :math:`=\mathcal J\xi^\zeta`.

For every ``jac_type`` (G1, G1F) the Jacobian-weighted perturbed field of
``gpeq_sol`` is :math:`\mathcal J\,\nabla\times(\boldsymbol\xi\times\mathbf B)`.


The pmodb columns
=================

The displacement identifications in this section assume unregularized
consistent components. The literal lagb expression is valid for both modes;
the regularized closure defect is given by R4-R5 below.

``eulb``
   :math:`\delta B_E=\hat{\mathbf b}\cdot\delta\mathbf B`, the Eulerian
   field strength (G3).

``lagb``
   :math:`\delta B_E+\xi^\psi\partial_\psi B+\xi_s/(\chi_1 q)\,\partial_\theta B`.
   This is the Lagrangian field strength of the displacement
   :math:`\boldsymbol\xi_g=\boldsymbol\xi_\perp+\lambda\mathbf B` with
   :math:`\xi_g^\zeta=0` in the *working* chart (``jac_type``) (G5.1-G5.4).
   It is neither the :math:`\boldsymbol\xi_\perp`-only value (negative
   control N1) nor Park et al. (2009) Eq. 7 with the divergence-free
   :math:`\boldsymbol\xi` that GPEC constructs (N2; G5.7 gives that value as
   ``lagb`` :math:`+\,\xi^\zeta/q\,\partial_\theta B`).  ``jac_out`` only
   relabels the harmonics; it does not change this gauge.

``Bdivxprp``
   :math:`B\,\nabla\cdot\boldsymbol\xi_\perp`, with
   :math:`\boldsymbol\xi_\perp=\boldsymbol\xi-(\boldsymbol\xi\cdot\mathbf B)\mathbf B/B^2`
   built from ``xmt``, ``xmz`` and the covariant ``xvt``, ``xvz`` (G4.1-G4.3).

``Bkxprp``
   :math:`\boldsymbol\xi_\perp\cdot\nabla B+\mu_0p'\xi^\psi/B` (G4.4).  This is
   :math:`B\boldsymbol\kappa\cdot\boldsymbol\xi_\perp` only through
   equilibrium force balance: exactly,
   :math:`B\boldsymbol\kappa\cdot\boldsymbol\xi_\perp=\boldsymbol\xi_\perp\cdot\nabla B
   +\boldsymbol\xi_\perp\cdot[(\nabla\times\mathbf B)\times\mathbf B]/B` (C2).
   The PENTRC/GPEC.jl matrices ``smat``, ``tmat`` (``idcon_build``) form the
   same quantity times :math:`\mathcal J B^2` from :math:`\xi^\psi` and
   :math:`\xi^\alpha` (G7).

``llagb`` (``gpec_diagnostics_pmodb``, ``pmodbmn_flag``)
   :math:`\delta B_L(\boldsymbol\xi_\perp)=\delta B_E+\boldsymbol\xi_\perp\cdot\nabla B`
   (G5.6).  ``lllagb`` is :math:`-(`\ ``Bdivxprp``\ :math:`+`\ ``Bkxprp``\ :math:`)`.


Relations between the columns
=============================

Gauge of ``lagb`` (G5.5, G5.8):

.. math::

   \mathrm{lagb}=\delta B_L(\boldsymbol\xi_\perp)-\frac{\xi_\perp^\zeta}{q}\partial_\theta B,
   \qquad
   \xi_\perp^\zeta=-\frac{B_\theta\xi_s/\chi_1+qB_\psi\xi^\psi}{B_\theta+qB_\zeta}.

Kinematic identity (C1; any :math:`\boldsymbol\xi_\perp`, no force balance)
and its pmodb form (C2 with G4; the optional direct contraction in the
GPEC chart with :math:`\det g=\mathcal J^2`, G6, returns UNKNOWN on ae635c80 and is recorded as an open gap, not as proved):

.. math::

   \delta B_L(\boldsymbol\xi_\perp)=-B(\nabla\cdot\boldsymbol\xi_\perp
   +\boldsymbol\kappa\cdot\boldsymbol\xi_\perp)
   =-(\mathrm{Bdivxprp}+\mathrm{Bkxprp})-\boldsymbol\xi_\perp\cdot\mathbf F/B,
   \qquad \mathbf F=(\nabla\times\mathbf B)\times\mathbf B-\mu_0\nabla p.

Hence

.. math::

   \mathrm{lagb}+\mathrm{Bdivxprp}+\mathrm{Bkxprp}
   =-\frac{\xi_\perp^\zeta}{q}\partial_\theta B-\frac{\boldsymbol\xi_\perp\cdot\mathbf F}{B}.

``lagb`` :math:`=-(`\ ``Bdivxprp``\ :math:`+`\ ``Bkxprp``\ :math:`)` is
therefore not an identity.  The first term is a gauge term of first order in
the perturbation. Its size depends on :math:`B_\theta/(qB_\zeta)`, the
Boozer :math:`B_\psi` and :math:`\partial_\theta B` and the displacement spectrum; growth toward the edge is a hypothesis
for a particular equilibrium, not a general bound.  The second term vanishes for an exact equilibrium.

Regularisation (``reg_flag``, Park et al. 2009 Eq. 10).  ``gpeq_sol`` uses
:math:`R\,\partial_\psi\xi^\psi`, :math:`R=s^2/(s^2+\sigma^2)`,
:math:`s=m-nq`, in the field and in ``xms``; ``gpeq_contra`` multiplies
``xwt``, ``xwz`` by :math:`R` again.  On one harmonic with
:math:`\theta`-independent :math:`\mathcal J` (R1-R3):

.. math::

   \mathcal J D_{effective}=(1-R)(\mathcal J R\,\partial_\psi\xi^\psi+\mathcal J'\xi^\psi),\qquad
   \mathcal J\nabla\cdot\delta\mathbf B=2\pi i\chi_1 s(1-R)\,\partial_\psi\xi^\psi.

Here the displayed R1 formula is the *effective* divergence diagnostic,
using :math:`R\partial_\psi\xi^\psi` for its radial derivative. The actual
divergence of the reconstructed displacement is instead (R1a)

.. math::

   \mathcal J\nabla\cdot\boldsymbol\xi_{rec}
   =\mathcal J(1-R^2)\partial_\psi\xi^\psi+(1-R)\mathcal J'\xi^\psi.

Also :math:`\chi_1(q\xi_{rec}^\theta-\xi_{rec}^\zeta)=R\xi_s` (R4).
The diagnostic llagb also differs from the projection of the reconstructed
contravariant displacement by (R5)

.. math::

   \mathrm{llagb}_{code}-(\delta B_E+\boldsymbol\xi_{rec,\perp}\cdot\nabla B)
   =(1-R)\frac{\xi_s}{\chi_1 q}\partial_\theta B.

These formulas are for the diagonal model only. Near rational surfaces the
regularized field, displacement and effective derivative therefore do not
form one exact ideal-induction pair. This is a limitation of the smoothing
policy, not evidence of resistive physics. The unconditional smoothing in
``gpout_pmodb`` even with ``reg_flag=f`` is a separate verified output bug.



What consumers need
===================

PENTRC, GPEC.jl ``KineticForces``
   The drive of Logan et al. (2013) Eqs. 4-5 needs :math:`\delta B` and
   :math:`\nabla\cdot\boldsymbol\xi` of the *same* displacement; the
   combination is then invariant under :math:`\boldsymbol\xi\to\boldsymbol\xi+\lambda\mathbf B`.
   The default PENTRC path (``peq_file`` = ``gpec_xclebsch``, ``set_peq``)
   uses :math:`\delta B/B=-(\nabla\cdot\boldsymbol\xi_\perp+\boldsymbol\kappa\cdot\boldsymbol\xi_\perp)`
   and :math:`\nabla\cdot\boldsymbol\xi_\perp` from the ``smat``-``zmat``
   matrices, i.e. a consistent :math:`\boldsymbol\xi_\perp` pair.  The
   ``pmodb_file`` path (``read_pmodb``) pairs ``lagb`` with
   :math:`-(\mathrm{lagb}+\mathrm{Bkxprp})/B`, which is not the divergence of
   :math:`\boldsymbol\xi_g`.  Under :math:`\boldsymbol\xi\to\boldsymbol\xi+\lambda\mathbf B`
   the consistent drive changes by a total derivative along the field line
   (K3), while this pairing changes it by
   :math:`-\mu B\lambda\,\partial_l B/v_\parallel` (K4).

Bounce-averaged codes in perturbed Boozer or Hamada coordinates (NEO-RT, NEO-2, Shaing-type)
   need :math:`\delta B` in the perturbed straight-field-line coordinates.
   That follows from :math:`\delta B_L(\boldsymbol\xi_\perp)` (``llagb``, or
   :math:`-(\mathrm{Bdivxprp}+\mathrm{Bkxprp})` away from rational surfaces)
   and :math:`\nabla\cdot\boldsymbol\xi_\perp` by a magnetic differential
   equation for the parallel gauge.  ``lagb`` must not be used as
   :math:`\delta B_L(\boldsymbol\xi_\perp)`.


Reproducing the proofs
======================

::

   PYTHONPATH=$FORTSYM/python python3 docs/derivations/pmodb_gauge/derive.py

FortSym (https://github.com/lazy-fortran/fortsym) is needed only to rerun the
proofs, not to build or run GPEC.


Target gauges in GPEC units (verified Q1-Q8)
============================================

Let :math:`D=\partial_\theta+q\partial_\zeta`, so
:math:`\mathbf B\cdot\nabla=(\chi_1/\mathcal J)D`.
For a matching equilibrium chart :math:`\mathcal J B^a=H(\psi)`,
:math:`a=2` is Boozer and :math:`a=0` is Hamada. A displacement with
:math:`b_X=\delta B_L(\boldsymbol\xi_X)` and
:math:`d_X=\nabla\cdot\boldsymbol\xi_X` has target-gauge source

.. math::

   T_a=d_X+a b_X/B,\qquad
   \boldsymbol\xi_{target}=\boldsymbol\xi_X+\lambda\mathbf B,\qquad
   T_a(\boldsymbol\xi_{target})=T_a(\boldsymbol\xi_X)
      +\frac{\mathbf B\cdot\nabla(\lambda B^a)}{B^a}.

For nonaxisymmetric modes the target flux-function perturbation is zero.
With :math:`w=\chi_1\lambda/\mathcal J`, this gives

.. math::

   Dw=-T_a,\qquad w_m=\frac{i T_{a,m}}{2\pi(m-nq)},\qquad
   b_{target}=b_X+w\partial_\theta B.

The source is formed in real space before Fourier transformation. A
nonzero source at :math:`m=nq` has no periodic solution. A resonant zero
source leaves an arbitrary homogeneous shift, set to zero here.
The factor :math:`2\pi` is from angles in turns and cancels the same
factor in :math:`\partial_\theta B`. Relabeling equilibrium harmonics
with ``jac_out`` alone does not solve this equation.

For a PEST chart :math:`\mathcal J/R^2` is a flux function; its
linear perturbation is :math:`\nabla\cdot\boldsymbol\xi-2\xi_R/R`
(Q8). PEST output is not automatically either target gauge. Transform the
equilibrium chart first, then solve the target equation. The Hamada scalar
on Boozer equilibrium angles uses :math:`T_2=2 b_H/B`, matching
rmp_torque c3614334, Eq. 4.14. Q5 proves composition of parallel shifts;
Q4 proves the angular-map determinant. The native regression compares the
finite spectral operator directly with that consumer's kernel.


Verified defects and numerical attribution
==========================================

* ``gpout_pmodb`` smooths the radial derivative even with ``reg_flag=f``;
  use ``xmp1_mn`` selected by ``gpeq_sol``. No change for ``reg_flag=t``.
* ``gpdiag_pmodb`` uses a stale ``jac`` in its first theta loop. The local
  Jacobian must be assigned before projecting out the parallel component.
* PENTRC ``read_pmodb`` overwrites the supplied divergence with
  ``-(lagb+Bkxprp)``. Use the supplied positive ``Bdivxprp`` and
  ``Bkxprp`` to form ``b_perp=-(Bdivxprp+Bkxprp)`` instead.
* Both CGL torque paths recover ``kappa.xi_perp`` as
  ``-(dbob+divx)/2`` although their input convention is
  ``dbob=-(divx+kapx)``. Remove the extra half. K1 and K5 prove the
  double-adiabatic pressure work for real and complex perturbations.

Finding (ii) from the TC24 reconstruction is a missing division by B in
the pressure *regressor*, not a GPEC normalization error (P1). For the
retained ITER n=3 output the original fit gives 0.1677; the corrected
pressure term plus the gauge term fits coefficients (1,1,1) with relative
residual about 1e-4. This is numerical corroboration only; ``docs/derivations/pmodb_gauge/check_output.py``
reproduces the fit from the pmodb, xclebsch, diagnostic and IDCON ASCII tables.

Finding (i) combines the gauge term, effective-derivative smoothing,
finite spectral/spline error, and the equilibrium force residual. The
identity above separates them; it does not prove a quantitative allocation
of the originally reported 0.5-19 percent to any one cause. The explicit
pressure term in curvature has coefficient one. Edge growth alone does
not diagnose nonideal plasma physics. Finding (iii) is consistent with the
explicit radial derivative in divergence and its smoothing near rational
surfaces. No change to the smoothing model itself is proposed here.


Source map and literature
==========================

Verified code contracts (line numbers at Fortran base 3622c25c):
``gpec/idcon.f:93`` (chi1=2pi psio), ``:616`` (projection geometry),
``:677`` (pressure in Pa), ``:1026`` (curvature and divergence matrices),
``gpec/gpeq.f:131`` (regularization selection), ``:208`` (surface-displacement
closure), ``gpec/gpout.f:3041`` (eulb/lagb), ``:3050`` (divergence and
pressure curvature), ``gpec/gpdiag.f:1012`` (diagnostic projection),
``:1061`` (llagb), ``pentrc/inputs.f90:376`` (ASCII columns),
``:396`` (old gauge mixing), ``:856`` (perpendicular matrix pair),
``pentrc/torque.F90:379`` (old CGL half).

Julia base 0304fd54c: ``src/KineticForces/KineticForcesStructs.jl:437``
uses the correct perpendicular dbob convention;
``src/KineticForces/Torque.jl:295`` repeats the CGL half.
``src/PerturbedEquilibrium/FieldReconstruction.jl``'s xi_s is Glasser's
surface displacement, not xi^zeta. Eq. 18 describes its elimination in the
Euler-Lagrange system; Eq. 9 defines it.

The physical definitions were checked against Park 2009, PoP 16 056115,
Eqs. 7 and 10; Park 2011, PoP 18 110702, Eqs. 4, 10 and 11;
Park and Logan 2017, PoP 24 032505, Eq. 13 and p. 6; Logan and Park 2013,
PoP 20 122507, Eqs. 4 and 5; and Boozer and Nuehrenberg 2006,
PoP 13 102501, Eqs. 42 and 43. The latter's Delta includes twice the
curvature; its Eq. 43 is Eulerian and is not the perpendicular Lagrangian
scalar. These are source contracts, separate from the native algebra proofs.

FortSym limits: #84, #85 and #86 are present at ae635c80. The direct G6
contraction remains UNKNOWN (33 symbols, rational metric jets); no probe
is promoted to proof. #91 (coefficient extraction after zero substitution)
and #92 (radical power reduction) are not exercised by the accepted
polynomial/rational gauge proofs. C1/C2 supply the invariant kinematic
identity independently of G6.
