.. _pmodb_gauge:

******************************************
Perturbed field strength (pmodb) and gauge
******************************************

This page states what the columns of ``gpec_pmodb_n#.out`` are, how they are
related, and which of them a field-strength consumer needs.

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
displacement changes neither :math:`\xi^\psi` nor :math:`\xi_s`, so
:math:`(\xi^\psi,\xi_s)` determine :math:`\boldsymbol\xi_\perp` completely.
The parallel part used elsewhere in GPEC is fixed in ``gpeq_contra`` by
:math:`\nabla\cdot\boldsymbol\xi=0`, which gives ``xmt``
:math:`=\mathcal J\xi^\theta` and ``xmz`` :math:`=\mathcal J\xi^\zeta`.
The pmodb columns do not use ``xmt`` or ``xmz``.

For every ``jac_type`` the Jacobian-weighted perturbed field of
``gpeq_sol`` is :math:`\mathcal J\,\nabla\times(\boldsymbol\xi\times\mathbf B)`.


The pmodb columns
=================

``eulb``
   :math:`\delta B_E=\hat{\mathbf b}\cdot\delta\mathbf B`, the Eulerian
   field strength.

``lagb``
   :math:`\delta B_E+\xi^\psi\partial_\psi B+\xi_s/(\chi_1 q)\,\partial_\theta B`.
   This is the Lagrangian field strength of the displacement
   :math:`\boldsymbol\xi_g=\boldsymbol\xi_\perp+\lambda\mathbf B` with
   :math:`\xi_g^\zeta=0` in the *working* chart (``jac_type``).
   It is neither the :math:`\boldsymbol\xi_\perp`-only value nor Park et
   al. (2009) Eq. 7 with the divergence-free :math:`\boldsymbol\xi` that
   GPEC constructs (that value is ``lagb``
   :math:`+\,\xi^\zeta/q\,\partial_\theta B`).  ``jac_out`` only
   relabels the harmonics; it does not change this gauge.

``Bdivxprp``
   :math:`B\,\nabla\cdot\boldsymbol\xi_\perp`.  ``gpout_pmodb`` forms
   :math:`\mathcal J\boldsymbol\xi_\perp` directly from :math:`\xi^\psi` and
   :math:`\xi_s`,

   .. math::

      \mathcal J\xi_\perp^\theta=\chi_1^2\left(-f_2\,\xi^\psi+f_3\,\xi_s/\chi_1\right),
      \qquad
      \mathcal J\xi_\perp^\zeta=q\,\mathcal J\xi_\perp^\theta-\mathcal J\xi_s/\chi_1,

   with the equilibrium coefficients
   :math:`f_2=(g_{\psi\theta}+qg_{\psi\zeta})/(\mathcal JB^2)` and
   :math:`f_3=(g_{\theta\zeta}+qg_{\zeta\zeta})/(\mathcal JB^2)` (``eqfun``
   columns 2 and 3, covariant metric).  The same projection gives the
   ``xmat``-``zmat`` matrices of ``idcon_build``.

``Bkxprp``
   :math:`\boldsymbol\xi_\perp\cdot\nabla B+\mu_0p'\xi^\psi/B`.  This is
   :math:`B\boldsymbol\kappa\cdot\boldsymbol\xi_\perp` only through
   equilibrium force balance: exactly,
   :math:`B\boldsymbol\kappa\cdot\boldsymbol\xi_\perp=\boldsymbol\xi_\perp\cdot\nabla B
   +\boldsymbol\xi_\perp\cdot[(\nabla\times\mathbf B)\times\mathbf B]/B`.
   The PENTRC matrices ``smat``, ``tmat`` (``idcon_build``) form the
   same quantity times :math:`\mathcal J B^2` from :math:`\xi^\psi` and
   :math:`\xi^\alpha`.

``llagb`` (``gpec_diagnostics_pmodb``, ``pmodbmn_flag``)
   :math:`\delta B_L(\boldsymbol\xi_\perp)=\delta B_E+\boldsymbol\xi_\perp\cdot\nabla B`
   ``lllagb`` is :math:`-(`\ ``Bdivxprp``\ :math:`+`\ ``Bkxprp``\ :math:`)`,
   and ``Bkxprp`` :math:`=` ``llagb`` :math:`-` ``eulb`` :math:`+\mu_0p'\xi^\psi/B`.


Relations between the columns
=============================

Gauge of ``lagb``:

.. math::

   \mathrm{lagb}=\delta B_L(\boldsymbol\xi_\perp)-\frac{\xi_\perp^\zeta}{q}\partial_\theta B,
   \qquad
   \xi_\perp^\zeta=-\frac{B_\theta\xi_s/\chi_1+qB_\psi\xi^\psi}{B_\theta+qB_\zeta}.

Kinematic identity (any :math:`\boldsymbol\xi_\perp`, no force balance)
and its pmodb form:

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
the perturbation; its size depends on :math:`B_\theta/(qB_\zeta)`,
:math:`B_\psi`, :math:`\partial_\theta B` and the displacement spectrum.
The second term vanishes for an exact equilibrium.  In the output, finite
Fourier truncation and spline error add to both.


Regularisation
==============

With ``reg_flag`` (Park et al. 2009 Eq. 10), ``gpeq_sol`` replaces the
radial derivative by :math:`R\,\partial_\psi\xi^\psi`,
:math:`R=s^2/(s^2+\sigma^2)`, :math:`s=m-nq`, in the field and in ``xms``.
The regularised field is then not divergence free:

.. math::

   \mathcal J\nabla\cdot\delta\mathbf B=2\pi i\chi_1 s(1-R)\,\partial_\psi\xi^\psi.

``Bdivxprp`` uses the same derivative ``xmp1`` as the field, and
:math:`\boldsymbol\xi_\perp` is built from the same ``xms``, so the
regularisation cancels between the columns and the relation above holds
unchanged with the regularised quantities.
Near rational surfaces the regularised field and displacement do not form
an exact ideal-induction pair.  This is a limitation of the smoothing
policy, not evidence of resistive physics.  With ``reg_flag=f`` the
columns are unregularised and singular at rational surfaces.


What consumers need
===================

PENTRC
   The drive of Logan et al. (2013) Eqs. 4-5 needs :math:`\delta B` and
   :math:`\nabla\cdot\boldsymbol\xi` of the *same* displacement; the
   combination is then invariant under
   :math:`\boldsymbol\xi\to\boldsymbol\xi+\lambda\mathbf B` up to a total
   derivative along the field line.  Both PENTRC input paths use the
   :math:`\boldsymbol\xi_\perp` pair
   :math:`\delta B/B=-(\nabla\cdot\boldsymbol\xi_\perp+\boldsymbol\kappa\cdot\boldsymbol\xi_\perp)`
   and :math:`\nabla\cdot\boldsymbol\xi_\perp`: the default path
   (``peq_file`` = ``gpec_xclebsch``, ``set_peq``) from the ``smat``-``zmat``
   matrices, and the ``pmodb_file`` path (``read_pmodb``) from
   :math:`-(\mathrm{Bdivxprp}+\mathrm{Bkxprp})` and ``Bdivxprp``.  The
   ``lagb`` column is not read.  Regenerate ``pmodb_file`` inputs with GPEC
   1.6 or later.

Bounce-averaged codes in perturbed Boozer or Hamada coordinates (NEO-RT, NEO-2, Shaing-type)
   need :math:`\delta B` in the perturbed straight-field-line coordinates.
   That follows from :math:`\delta B_L(\boldsymbol\xi_\perp)` (``llagb``, or
   :math:`-(\mathrm{Bdivxprp}+\mathrm{Bkxprp})` away from rational surfaces)
   and :math:`\nabla\cdot\boldsymbol\xi_\perp` by a magnetic differential
   equation for the parallel gauge (next section).  ``lagb`` must not be
   used as :math:`\delta B_L(\boldsymbol\xi_\perp)`.


Target gauges in GPEC units
===========================

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
linear perturbation is :math:`\nabla\cdot\boldsymbol\xi-2\xi_R/R`.
PEST output is not automatically either target gauge. Transform the
equilibrium chart first, then solve the target equation.


History
=======

Before GPEC 1.6.0, ``Bdivxprp`` and ``Bkxprp`` subtracted a
Fourier-truncated parallel displacement, ``gpout_pmodb`` regularised the
divergence even with ``reg_flag=f``, ``gpdiag_pmodb`` used a stale Jacobian,
PENTRC's ``read_pmodb`` paired ``lagb`` with an inconsistent divergence, and
both PENTRC CGL paths halved the curvature term.  These were fixed in
1.6.0 (#340-#343).


Source map and literature
=========================

Code: ``gpeq_sol`` and ``gpeq_contra`` (``gpec/gpeq.f``, regularisation and
the divergence-free closure), ``gpout_pmodb`` (``gpec/gpout.f``, the pmodb
columns), ``gpdiag_pmodb`` (``gpec/gpdiag.f``, ``llagb`` and ``lllagb``),
``idcon_build`` (``gpec/idcon.f``, ``smat``-``zmat``), ``read_pmodb`` and
``set_peq`` (``pentrc/inputs.f90``), and ``direct_run`` (``equil/direct.f``,
``eqfun``).

The physical definitions were checked against Park 2009, PoP 16 056115,
Eqs. 7 and 10; Park 2011, PoP 18 110702, Eqs. 4, 10 and 11;
Park and Logan 2017, PoP 24 032505, Eq. 13 and p. 6; Logan and Park 2013,
PoP 20 122507, Eqs. 4 and 5; and Boozer and Nuehrenberg 2006,
PoP 13 102501, Eqs. 42 and 43. The latter's Delta includes twice the
curvature; its Eq. 43 is Eulerian and is not the perpendicular Lagrangian
scalar.
