.. _pmodb_gauge:

******************************************
Perturbed field strength (pmodb) and gauge
******************************************

This page states what the columns of ``gpec_pmodb_n#.out`` are, how they are
related, and which of them a field-strength consumer needs.  Every identity
quoted here is proved with the FortSym native zero test by
``docs/derivations/pmodb_gauge/derive.py``; the verdicts are stored in
``docs/derivations/pmodb_gauge/fortsym_verdicts.json`` and the labels (G1,
G5.5, ...) refer to that file.

.. contents:: On this page
   :local:
   :depth: 1


Coordinates and displacement
============================

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
and its pmodb form (C2 with G4; the same statement written directly in the
GPEC chart with :math:`\det g=\mathcal J^2`, G6, exceeds the FortSym
polynomial-term limit and is recorded as an open gap, not as proved):

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
the perturbation; it is largest where :math:`B_\theta/(qB_\zeta)`, the
Boozer :math:`B_\psi` and :math:`\partial_\theta B` are large, i.e. toward
the edge.  The second term vanishes for an exact equilibrium.

Regularisation (``reg_flag``, Park et al. 2009 Eq. 10).  ``gpeq_sol`` uses
:math:`R\,\partial_\psi\xi^\psi`, :math:`R=s^2/(s^2+\sigma^2)`,
:math:`s=m-nq`, in the field and in ``xms``; ``gpeq_contra`` multiplies
``xwt``, ``xwz`` by :math:`R` again.  On one harmonic with
:math:`\theta`-independent :math:`\mathcal J` (R1-R3):

.. math::

   \mathcal J\nabla\cdot\boldsymbol\xi=(1-R)(\mathcal J R\,\partial_\psi\xi^\psi+\mathcal J'\xi^\psi),\qquad
   \mathcal J\nabla\cdot\delta\mathbf B=2\pi i\chi_1 s(1-R)\,\partial_\psi\xi^\psi.

Inside :math:`|m-nq|\lesssim\sigma` the regularised displacement is neither
divergence free nor does it produce a divergence-free field, so
``Bdivxprp`` (which carries the explicit :math:`R\,\partial_\psi\xi^\psi`)
and ``lagb`` (which uses only ``xms``) are not related by the identities
above on the resonant harmonics.


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
