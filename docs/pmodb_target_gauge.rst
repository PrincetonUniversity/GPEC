Optional perturbed-coordinate field strength
=============================================

With ``pmodb_flag=t`` and ``reg_flag=f``, NetCDF ``b_lag_perp`` supplies
``eulb+xi_perp.grad B``. It is the perpendicular-displacement value called
``llagb`` in ``gpdiag_pmodb``. It is omitted for regularized solutions, whose reconstructed surface
displacement need not equal xms/chi1. Existing ASCII columns keep their format.
This scalar alone is not the perturbed-Boozer input to NEO-RT or NEO-2.

Set ``pmodb_gauge_flag=t`` in ``GPEC_OUTPUT`` for an additional NetCDF
field ``b_lag_boozer`` or ``b_lag_hamada``. This requires ``reg_flag=f``
and a matching ``jac_type='boozer'`` or ``'hamada'``. The gauge is selected
by the working chart; ``jac_out`` selects only the harmonic representation.
The field has units Tesla and the same complex harmonic layout as ``b_lag``.
The output's ``gauge`` attribute records the physical gauge separately from
its ``jacobian`` representation. Use Boozer for NEO-type scalar-only
consumers and Hamada for Shaing-type formulas with that gauge contract.

Verified derivation
-------------------

GPEC angles are in turns and modes are exp(2 pi i(m theta-n zeta)). Write
``D=partial_theta+q partial_zeta``, so ``B.grad=(chi1/J)D``. On the matching
base ``J B^a=H(psi)``, with a=2 for Boozer or a=0 for Hamada. Under the
parallel shift ``xi_target=xi_X+lambda B``, ideal induction is unchanged,
``b_target=b_X+lambda B.grad B`` and ``d_target=d_X+B.grad lambda``.
The target condition for a nonaxisymmetric mode is ``d_target+a b_target/B=0``.
The product rule gives

.. math::

   T=d_X+a b_X/B,\qquad
   T_{target}=T+\frac{\mathbf B\cdot\nabla(\lambda B^a)}{B^a},\qquad
   w=\chi_1\lambda/J,\qquad
   Dw=-T,\qquad w_m=\frac{iT_m}{2\pi(m-nq)},\qquad
   b_{target}=b_X+w\partial_\theta B.

The implementation forms T in real space, truncates to the existing GPEC
poloidal spectrum, inverts each nonresonant mode, and forms the scalar in
real space before the existing output-coordinate transform. Identities are
exact; numerical output retains GPEC's finite spectral/spline resolution.
A nonzero resonant source is rejected because no periodic solution exists.
At a resonant zero source the free homogeneous shift is set to zero.
Numerically, abs(m-nq)<1e-10 is classified as resonant; a source coefficient
below 1e-12 times the largest source coefficient is treated as zero there.
The near-resonant result can be large: this feature applies no smoothing.

The input pair is ``b_X=b_lag_perp`` and the ordinary perpendicular
coordinate divergence, using the derivative selected by gpeq_sol. Regularized fields/derivatives are deliberately
excluded because they need not form one ideal-induction displacement.

Verification and scope
-----------------------

``gen/pmodb_gauge/prove.py`` proves Q1-Q8 and K5/P1 using FortSym's native
zero test, with no probe fallback (verified at ae635c80). The source, complex inverse
coefficient and scalar correction are separately proved and emitted by
``gen/pmodb_gauge/generate.py`` to ``gpec/pmodb_gauge_kernel.f90``.
Regenerate and check the resulting diff; ordinary builds require no CAS.

``regression/test/pmodb_target_gauge.f90`` exercises the actual operator:
parallel displacement gives zero in both gauges, a zero target source
leaves the scalar unchanged, and a nonzero resonant source fails. It also
writes a Hamada scalar transformed on Boozer equilibrium angles.
``regression/test/pmodb_target_compare.py`` compares it with rmp_torque's
c3614334/#55 kernel (external source supplied through PYTHONPATH).

For the Hamada scalar already relabeled onto Boozer equilibrium angles,
T=2 b_H/B; both paths use the same inverse-field-line operator. Equilibrium
relabeling must precede this conversion. PEST output is not supported as a
target; its J/R^2 constraint also involves radial geometric displacement.
This feature does not infer a new equilibrium chart from a scalar spectrum.

Source contracts: ``gpec/gpeq.f:131`` selects the derivative;
``gpec/gpdiag.f:1061`` defines the existing perpendicular scalar;
``gpec/gpout.f`` computes it before output transformation. The scalar-only
consumer formulas are those derived in rmp_torque c3614334,
``docs/boozer_gauge_ntv_derivation.md`` Eqs. 4.3-4.8 and 4.13-4.15.
Park 2011 PoP 18 110702 Eqs. 10-11 explain why a general kinetic action
with an explicit Jacobian term has a different input contract.
