.. _sign_conventions:

****************************
Sign Conventions Reference
****************************

This page is a comprehensive reference for the sign conventions used
throughout GPEC. Understanding these conventions is essential for
interpreting outputs, interfacing with other codes, and correctly
setting up kinetic calculations.

.. contents:: On this page
   :local:
   :depth: 2


Coordinate System
=================

GPEC uses right-handed magnetic coordinates :math:`(\psi, \theta, \zeta)` with
Fourier decomposition :math:`\exp(im\theta - in\phi)`.

Poloidal Flux :math:`\psi`
--------------------------

- Normalized from 0 (magnetic axis) to 1 (plasma boundary).
- :math:`\psi_0 = \psi_{\mathrm{bry}} - \psi_{\mathrm{axis}}` is forced positive
  on read (``read_eq_efit`` in ``equil/read_eq.f``). If the EQDSK has
  :math:`\psi_{\mathrm{bry}} < \psi_{\mathrm{axis}}`, both :math:`\psi_0` and
  the 2D flux array are sign-flipped.
- Occasionally :math:`\rho = \sqrt{\psi}` is used as a radius-like variable.

Poloidal Angle :math:`\theta`
-----------------------------

- "Upward outboard" convention, always. :math:`\theta` increases by 1 (not
  :math:`2\pi`) going once around the poloidal direction.
- This is fixed regardless of helicity or working coordinate choice.

Toroidal Coordinate :math:`\zeta` and :math:`\phi`
---------------------------------------------------

- The ignorable toroidal coordinate is :math:`\zeta = \phi/(2\pi) + \nu(\psi,\theta)`,
  where :math:`\nu` is a single-valued function that depends on the working
  coordinate system. PEST coordinates have :math:`\nu = 0`.
- :math:`\phi` is effectively **CCW** (counter-clockwise from above) for
  left-handed (LH) configurations, but **CW** (clockwise) for right-handed
  (RH) configurations.

Working Coordinate Options
--------------------------

Controlled by ``jac_type`` in ``equil.in``. The Jacobian is
:math:`J \propto B_p^{p_{bp}} \, B^{p_b} \, R^{-p_r}`:

.. list-table::
   :header-rows: 1
   :widths: 20 20 20 20

   * - Name
     - ``power_bp``
     - ``power_b``
     - ``power_r``
   * - Hamada (default)
     - 0
     - 0
     - 0
   * - PEST
     - 0
     - 0
     - 2
   * - Boozer
     - 0
     - 2
     - 0
   * - Equal-arc
     - 1
     - 0
     - 0

See also the coordinate discussion in :doc:`outputs`.


Helicity and Handedness
=======================

Helicity is computed in the ``gpec_main`` program (``gpec/gpec.f``):

.. code-block:: fortran

   ipd = 1.0
   btd = 1.0
   IF(ip_direction=="negative") ipd = -1.0
   IF(bt_direction=="negative") btd = -1.0
   helicity = ipd * btd

- ``ip_direction`` and ``bt_direction`` are set in ``coil.in``.
  "positive" means CCW viewed from above; "negative" means CW.
- **helicity = +1**: right-handed (RH) --- :math:`B_t` and :math:`I_p` in the
  same direction.
- **helicity = -1**: left-handed (LH) --- :math:`B_t` and :math:`I_p` opposed.
- The helicity value is stored in the ``gpec_control_output`` netcdf file.


:math:`F = R B_\phi`
=====================

The poloidal-current function :math:`F = R B_\phi` from the Grad-Shafranov
equation is **forced positive** via ``ABS()`` (``read_eq_efit`` in
``equil/read_eq.f``):

.. code-block:: fortran

   sq_in%fs(:,1) = ABS(sq_in%fs(:,1))

The code always works with :math:`|F|`. Any information about the sign of
:math:`B_t` must come from the ``bt_direction`` setting in ``coil.in``.


Safety Factor :math:`q`
========================

- Defined as :math:`q = (\mathbf{B} \cdot \nabla\zeta) / (\mathbf{B} \cdot \nabla\theta)`.
- **Not forced positive.** Read directly from the EFIT g-file without sign
  manipulation (unlike :math:`F`).
- In standard tokamak operation, EFIT provides :math:`q > 0`.
- Can be adjusted via ``newq0`` in ``equil.in``, which modifies :math:`F`
  to match while preserving the Grad-Shafranov solution
  (``direct_run`` in ``equil/direct.f``).


Mode Numbers :math:`m` and :math:`n`
=====================================

Toroidal Mode Number :math:`n`
------------------------------

- Set as a **positive integer** ``nn`` in ``dcon.in``. GPEC runs one toroidal
  harmonic at a time; results for multiple :math:`n` can be superposed.

Poloidal Mode Range
-------------------

The poloidal mode spectrum spans ``mlow`` to ``mhigh``, computed in the
``dcon`` program (``dcon/dcon.F``):

.. code-block:: fortran

   mlow  = MIN(nn*qmin, zero) - 4 - delta_mlow
   mhigh = nn*qmax + delta_mhigh

``delta_mlow`` and ``delta_mhigh`` (set in ``dcon.in``) widen the range
beyond the resonant modes.

Why Positive :math:`m` Is Always Resonant
-----------------------------------------

Resonant surfaces are found by ``sing_find`` in ``dcon/sing.f``.
It performs a binary search for flux surfaces where :math:`m = n \cdot q`.
Since :math:`n > 0` (by convention) and :math:`q > 0` (standard tokamak):

.. math::

   m_{\mathrm{res}} = \mathrm{NINT}(n \cdot q) > 0

**Resonant modes always have positive** :math:`m`. Negative-:math:`m` modes
are always non-resonant. This is by design: one can always plot the
:math:`m = 2` displacement profile and see resonant behavior at :math:`q = 2`,
while :math:`m = -2` is always non-resonant.


Spectrum Output Sign Conventions
================================

Real-Space Output
-----------------

For the real-space representation decomposed in :math:`\exp(-in\phi)` with
CCW :math:`\phi`, GPEC takes the complex conjugate for RH configurations.
This is implemented throughout ``gpec/gpout.f`` as:

.. code-block:: fortran

   -helicity * AIMAG(quantity)

For the full Fourier representation :math:`\exp(im\theta - in\phi)`, the
conjugate operation also flips up and down (not just the toroidal direction).

Interfacing with SURFMN
-----------------------

SURFMN expands in :math:`\exp(-im\theta - in\phi)` and always uses CCW
:math:`\phi`. To convert:

.. code-block:: python

   m_surfmn = helicity * m_gpec
   b_surfmn = real(b_m) - 1j * helicity * imag(b_m)

For LH configurations: only the sign of :math:`m` is flipped. For RH: :math:`m`
is unchanged but the complex conjugate is taken.

Interfacing with VACUUM
-----------------------

The VACUUM code uses CCW :math:`\phi` and downward outboard :math:`\theta`.
GPEC uses the complex conjugate of RH configurations when interfacing with
VACUUM.


Rotation Velocity Conventions (PENTRC)
======================================

:math:`\omega_E` (E x B Rotation)
----------------------------------

- Column 6 of the PENTRC kinetic profile file: :math:`\omega_E` in rad/s.
- Read by the ``read_kin`` subroutine in ``pentrc/inputs.f90``.
- **Sign convention**: positive :math:`\omega_E` means rotation in the
  direction of the toroidal coordinate :math:`\zeta`.
- Since :math:`\phi` direction depends on helicity, positive :math:`\omega_E`
  is effectively **co-current for RH plasmas** and **counter-current for
  LH plasmas**.

Diamagnetic Frequencies
-----------------------

Computed in ``read_kin`` (``pentrc/inputs.f90``) and ``tpsi``
(``pentrc/torque.F90``):

.. math::

   \omega_{*n} = -\frac{2\pi \, T_i}{e \, Z_i \, \chi_1 \, n_i} \frac{dn_i}{d\psi_n}

.. math::

   \omega_{*T} = -\frac{2\pi}{e \, Z_i \, \chi_1} \frac{dT_i}{d\psi_n}

where :math:`\chi_1 = 2\pi \psi_0` with :math:`\psi_0` the boundary poloidal flux.
The negative signs mean that a positive (outward-increasing) density or
temperature gradient yields a negative diamagnetic frequency.

Total Toroidal Rotation
-----------------------

.. code-block:: fortran

   wphi = welec + wdian + wdiat    ! tpsi in pentrc/torque.F90

The total toroidal rotation frequency is the sum of the E x B, density
diamagnetic, and temperature diamagnetic contributions.

Rotation Scaling Parameters
---------------------------

``pentrc.in`` provides two knobs:

- ``wefac``: direct multiplier on the :math:`\omega_E` profile.
- ``wpfac``: scales the total rotation :math:`\omega_\phi = \omega_E + \omega_{*n} + \omega_{*T}`
  by indirectly adjusting :math:`\omega_E`.

Energy Integral Resonance
-------------------------

In the PENTRC energy integral (``xintgrnd`` in ``pentrc/energy.f90``), the
resonance denominator involves:

.. math::

   n \omega_E + \ell_{\mathrm{eff}} \omega_b \sqrt{x} + n \omega_D x

where :math:`\omega_b` is the bounce frequency divided by :math:`x`,
:math:`\omega_D` is the magnetic precession frequency,
:math:`\ell_{\mathrm{eff}} = \ell - \sigma n q` is the effective bounce harmonic,
and :math:`x = E/T` is the normalized energy. The sign of :math:`\omega_E`
determines the direction of resonance in velocity space.


COCOS Compatibility
===================

GPEC does **not** use or reference the COCOS (Coordinate Convention Standard)
system. The conventions described in this document are GPEC-native and predate
COCOS. Users interfacing with COCOS-aware codes must manually translate between
conventions.


Quick Reference
===============

.. list-table::
   :header-rows: 1
   :widths: 30 50 20

   * - Quantity
     - Convention
     - Forced?
   * - :math:`\psi` (poloidal flux)
     - Normalized 0 (axis) to 1 (edge)
     - Yes, positive
   * - :math:`\theta` (poloidal angle)
     - Upward outboard
     - Fixed
   * - :math:`\phi` (toroidal angle)
     - CCW for LH, CW for RH
     - By helicity
   * - :math:`F = R B_\phi`
     - Always positive
     - Yes, ABS()
   * - :math:`q` (safety factor)
     - From EFIT, typically positive
     - No
   * - :math:`n` (toroidal mode)
     - Always positive
     - By convention
   * - Resonant :math:`m`
     - Always positive (since :math:`n > 0`, :math:`q > 0`)
     - By construction
   * - helicity
     - +1 RH, -1 LH
     - Computed
   * - :math:`\omega_E`
     - Positive = direction of :math:`\zeta`
     - No


Source Code References
======================

- **Poloidal flux sign**: ``read_eq_efit`` in ``equil/read_eq.f``
- **F = R*Bt forced positive**: ``read_eq_efit`` in ``equil/read_eq.f``
- **q definition**: ``dcon/README``
- **Helicity computation**: ``gpec_main`` program in ``gpec/gpec.f``
- **ip/bt direction**: ``input/coil.in``
- **Poloidal mode range**: ``dcon`` program in ``dcon/dcon.F``
- **Resonant surface finder**: ``sing_find`` in ``dcon/sing.f``
- **Output sign flips**: ``gpec/gpout.f`` (many locations, search ``helicity``)
- **omega_E input**: ``read_kin`` in ``pentrc/inputs.f90``
- **Diamagnetic frequencies**: ``tpsi`` in ``pentrc/torque.F90``
- **Energy integral**: ``xintgrnd`` in ``pentrc/energy.f90``
- **SURFMN interface**: ``docs/outputs.rst``
