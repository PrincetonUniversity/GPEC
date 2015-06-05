IPEC Namelist Inputs
********************

IPEC_INPUT
==========

**jac_in** = "pest"
  Coordinate system of 3D field input.  Choose from hamada, pest, boozer, equal_arc or other

**jsurf_in** = 1
  Set to 1 for surface weighted spectrum  Common for surfmn inputs

**tmag_in** = 0
  0 for ordinary toroidal angle 1 for magnetic coordinate angle

**data_flag** = .TRUE.
  Set true to read a 3D field spectrum from the file specified by the infile variable and apply it to the plasma boundary.

**data_type** = "surfmn"
  Type of vacuum 3D error field input file.  Choose from 'vac3d' or 'surfmn'.

**infile** = "/p/gpec/users/jpark/data/3d_field/d3d/surfmn.out.cossin.144182.03750_Intrins"
  Path to vacuum 3D field input file of type specified by data_type.  This file must specify the 3D vacuum field on the plasma boundary as determined in DCON.

**nmin** = 0
  The minimum toroidal mode number included in the infile data (typically 0).

**nmax** = 64
  The maximum toroidal mode number included in the infile data.

**mmin** = -64
  The minimum poloidal mode number included in the infile data (typically -64).

**mmax** = 64
  The maximum poloidal mode number included in the infile data (typically 64).

**harmonic_flag** = .FALSE.
  If true Tesla amplitudes must be assigned to variables cosmn(int) and sinmn(int), where int is the applied 3D field m  harmonic on the plasma surface. The toroidal mode number of the applied field is set by the DCON mode being analyzed.

**cosmn(3)** = 0.0001


**displacement_flag** = .FALSE.
  Take all boundary inputs as normal displacements rather than magnetic fields

**fixed_boundary_flag** = .FALSE.
  Sets the plasma boundary as defined in the equilibrium input file fixed in all IPEC calculations.

**coil_flag** = .TRUE.




IPEC_CONTROL
============

**resp_index** = 0


**sing_spot** = 0.0005


**reg_flag** = .TRUE.
  Regularize unintegrable pmodb in psi

**reg_spot** = 0.05
  dB ~ (m-nq)^2/[(m-nq)^2+reg_spot^2)] if reg_flag

**chebyshev_flag** = .TRUE.
  Set true for chebyshev polynomials of the 3D perturbations in (\psi,artheta).

**nche** = 20
  Number of coefficients used in the chebyshev polynomials (typicaly ~30).



IPEC_OUTPUT
===========

**jac_out** = "hamada"
  Coordinates of IPEC output files.  Choose from hamada, pest, boozer, equal_arc or other

**jsurf_out** = 0
  Set to 1 for a surface weighted spectrum b_{mn}=\oint\left(\delta\mathbf{B}\cdot\hat{n}ight)\left(	heta,\phiight)e^{-i\left(m	heta-n\phiight)}J\cdot|\nabla\psi|d	heta d\phi (an invariant flux on rational surfaces). Set to 0 for unweighted spectrum b_{mn}=\oint\left(\delta\mathbf{B}\cdot\hat{n}ight)\left(	heta,\phiight)e^{-i\left(m	heta-n\phiight)}d	heta d\phi. 

**tmag_out** = 1
  0 for ordinary toroidal angle. 1 for magnetic coordinate angle

**resp_flag** = .FALSE.
  Output ipec_response_n#.out file.

**singcoup_flag** = .FALSE.


**singfld_flag** = .FALSE.
  Output the q, \Psi_{N}, resonant \mathbf{B}, singular current, island width and Chirikov parameter at rational surfaces inside the plasma, including the plasma response.

**vsingfld_flag** = .FALSE.
  Output the \psi_{N}, resonant \mathbf{B}, singular current, island width and Chirikov parameter at rational surfaces inside the plasma, using vacuum calculations.  Requires coil_flag=True.

**pmodb_flag** = .FALSE.
  Plasma |\delta B| on flux surfaces.

**xbnormal_flag** = .FALSE.
  Output the normal components of the displacement, covarient perturbed magnetic field, and contravarient perturbed field (area weighted -> can compare between coordinates)  with plasma response included.

**vbnormal_flag** = .FALSE.
  Output the normal components of the displacement, covarient perturbed magnetic field, and contravarient perturbed field (area weighted -> can compare between coordinates) from vacuum calculations.  Requires coil_flag=True.

**fun_flag** = .FALSE.
  Outputs in functional format (on a coordinate grid) by back transforming the harmonic information (WARNING: this will result in larger data output).

**flux_flag** = .FALSE.
  New flux surfaces \psi(r+\delta r,z+\delta z,\phi) for visualization of perturbation.

**eqbrzphi_flag** = .FALSE.
  Output the original equilibrium field on a (r,z,\phi) grid.

**brzphi_flag** = .FALSE.
  Output the (r,z,\phi) components of the perturbed magnetic fields everywhere on the (r,z,\phi) grid.  - IPEC 1 only outputs perturbed magnetic fields valid inside the plasma on the (r,z,\phi) grid.  - IPEC 2.0 and later includes the plasma response and vacuum fields when the coil_flag is the vacuum input.

**xrzphi_flag** = .FALSE.
  Output the plasma displacement on a (r,z,\phi) grid.

**vbrzphi_flag** = .FALSE.
  Output the field due to the surface current defining the IPEC final solution boundary condition at the plasma surface (original external field boundary condition and plasma response) on an r,z grid. This does not represent a true field.

**vvbrzphi_flag** = .FALSE.
  Outputs the field due to the surface current defining the external 3D field boundary condition on a (r,z,\phi) grid. This field does not represent a true field.

**bin_flag** = .TRUE.
  Output binary files for use with the xdraw command when appropriate.

**bin_2d_flag** = .TRUE.
  Output binary files containing 2D contours for use with the xdraw command  when appropriate.

**vsbrzphi_flag** = .FALSE.


**ss_flag(7)** = .FALSE.


**ss_flag(8)** = .FALSE.


**xbrzphifun_flag** = .FALSE.


**arzphifun_flag** = .FALSE.


**ntv_flag** = .TRUE.
  Output first order quantities to bin files for interfacing with PENT code



IPEC_DIAGNOSE
=============

**div_flag** = .FALSE.


**radvar_flag** = .TRUE.




