PENT Namelist Inputs
********************

PENT_INPUT
==========

**kinetic_file** = "/p/gpec/users/nlogan/data/a10/kinetics/kinetic_a10_MARSbenchmark.dat"
  1 header line followed by column data for psi_n, ni(m^-3), ne(m^-3), ti(eV), te(eV), omega_E(rad/s).

**o1fun_file** = "ipec_order1_fun_n1.bin"
  IPEC output generated using ntv_flag or ascii B(psi,theta).

**o1mn_file** = "ipec_order1_n1.bin"
  IPEC output generated using ntv_flag, or ascii B(m,n).

**idconfile** = "euler.bin"


**o1native** = .TRUE.
  o1fun_file and o1mn_file are IPEC output bin files.

**imass** = 2
  Mass of main species in atomic mass units (0 for electrons).

**icharge** = 1
  Charge of main species in elementry charge units.

**zimp** = 6
  Charge of impurity species in elementry charge units.

**zmass** = 12
  Mass of impurity species in atomic mass units.

**nl** = 6
  Number of bounce harmonics calculated (-nl<=l<=nl).

**collision** = "harmonic"
  Collision operator. Choose from zero,small,krook, or harmonic.



PENT_CONTROL
============

**wefac** = 1
  Multiplier applied to omega_E profile.

**wdfac** = 1
  Multiplier applied to oemga_D.

**wpfac** = 1


**ptfac** = 0.001
  Fractional cutoff in Lambda space near trapped-passing boundry.

**kolim_flag** = .FALSE.
  Infinite rotation limit.

**neorot_flag** = .FALSE.
  Approximate numerator of resonant operator using neoclassical rotation.

**divxprp_flag** = .TRUE.
  Include arc-length change in perturbed Action.

**xlsode_flag** = .TRUE.
  Use LSODE package to integrate in energy space.

**lmdalsode_flag** = .TRUE.
  Use LSODE package to integrate in pitch angle.

**lsode_rtol** = 1e-09
  Relative tolarance in LSODE.

**lsode_atol** = 1e-12
  absolute tolarance in LSODE.

**xmax** = 72
  Upper limit of energy integration.

**ximag** = 0
  Step contour of integration off of real axis to avoid poles when collision is zero.

**nx** = 128
  Number of energy spline points if xlsode_flag false.

**nt** = 128
  Number of theta points used in bounce averaging.

**nlmda** = 128
  Number of Lambda points used in splines.

**maskpsi** = 2
  Step size taken by PENT on dcon psi grid (determined by idconfile).



PENT_OUTPUT
===========

**passing_flag** = .FALSE.
  Calculate torque and energy contributed by passing particles.

**gen_flag** = .TRUE.
  Perform full calculation in general geometry.

**rlar_flag** = .TRUE.
  Perform Reduced LAR approximate calculation.

**clar_flag** = .FALSE.
  Perform Circular LAR approximate calculation.

**kinetic_flag** = .FALSE.
  Output kinetic profiles.

**bounce_flag** = .FALSE.
  Output sample bounce averaging profiles.

**pitch_flag** = .FALSE.
  Output sample pitch angle profiles.

**energy_flag** = .FALSE.
  Output sample energy profiles.

**psiout(1)** = 0.1
  Include this surface in sampled outputs.

**psiout(10)** = 1.0
  Include this surface in sampled outputs.



PENT_ADMIN
==========

**fmnl_flag** = .FALSE.
  Do not change unless you know what you are doing.

**diag_flag** = .FALSE.
  Do not change unless you know what you are doing.



