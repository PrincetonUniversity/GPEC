"""
Script

Dictionary of tooltips explaining GPEC inputs.
"""

equil = {
'eq_type' : 
"Most common is 'efit'",
'jac_type' : 
"Coordinate system used in computations.\n  Choose from hamada, pest, boozer, equal_arc or other",
'psi_low'  : 
"Lower limit of normalized poloidal flux",
'psi_high'  : 
"Upper limit of normalized poloidal flux.\n  May be made variable by sas_flag in dcon namelist",
'newq0' : 
'Force q(0) as high as 1.2 to avoid internal kink instability'
}


dcon = {
"sas_flag" : 
"Move last flux surface inside of last rational surface to avoid peeling instability.",
"dmlim" : 
"Fractional distance between last rational surfaces if sas_flag true",
"nn" : 
"Single toroidal mode under consideration",
"delta_mlow" : 
"Range of poloidal mode numbers below minimum nq considered",
"delta_mhigh" : 
"Range of poloidal mode numbers above maximum nq considered",
"singfac_min" : 
"Fractional distance from singular surfaces approached by ODE solver",
"ucrit" : 
"Limit on fractional discrepency between large and small solutions (renormalizes matrix)",
"crossover" : ""
}

ipec = {
"jac_in" : 
"Coordinate system of 3D field input.\n  Choose from hamada, pest, boozer, equal_arc or other",
"power_rin" : 
"Manually defined magnetic coordinates power \alpha such that \mathbf{B}\cdot\\nabla\theta\propto R^{\alpha}r^{\beta}/B^{\gamma}B_{p}^{\delta}.",
"power_rcin" : 
"Manually defined magnetic coordinates power \beta such that \mathbf{B}\cdot\\nabla\theta\propto R^{\alpha}r^{\beta}/B^{\gamma}B_{p}^{\delta}.",
"power_bin" : 
"Manually defined magnetic coordinates power \gamma such that \mathbf{B}\cdot\\nabla\theta\propto R^{\alpha}r^{\beta}/B^{\gamma}B_{p}^{\delta}.",
"power_bpin" : 
"Manually defined magnetic coordinates power \delta such that \mathbf{B}\cdot\\nabla\theta\propto R^{\alpha}r^{\beta}/B^{\gamma}B_{p}^{\delta}.",
"jsurf_in" : 
"Set to 1 for surface weighted spectrum\n  Common for surfmn inputs",
"tmag_in" : 
"0 for ordinary toroidal angle \n1 for magnetic coordinate angle",
"data_flag" :
"Set true to read a 3D field spectrum from the file specified by the infile variable and apply it to the plasma boundary.",
"data_type" :
"Type of vacuum 3D error field input file.\n  Choose from 'vac3d' or 'surfmn'.",
"infile" :
"Path to vacuum 3D field input file of type specified by data_type.\n  This file must specify the 3D vacuum field on the plasma boundary as determined in DCON.",
"nmin" :
"The minimum toroidal mode number included in the infile data (typically 0).",
"nmax" :
"The maximum toroidal mode number included in the infile data.",
"mmin" :
"The minimum poloidal mode number included in the infile data (typically -64).",
"mmax" :  
"The maximum poloidal mode number included in the infile data (typically 64).",
"harmonic_flag" :
"If true Tesla amplitudes must be assigned to variables cosmn(int) and sinmn(int), where int is the applied 3D field m  harmonic on the plasma surface. The toroidal mode number of the applied field is set by the DCON mode being analyzed." ,
"fixed_boundary_flag" : 
"Sets the plasma boundary as defined in the equilibrium input file fixed in all IPEC calculations.",
"displacement_flag" : 
"Take all boundary inputs as normal displacements rather than magnetic fields",

"reg_flag" : 
"Regularize unintegrable pmodb in psi",
"reg_spot" : 
"dB ~ (m-nq)^2/[(m-nq)^2+reg_spot^2)] if reg_flag",
"chebyshev_flag" :
"Set true for chebyshev polynomials of the 3D perturbations in (\psi,\vartheta).",
"nche" :
"Number of coefficients used in the chebyshev polynomials (typicaly ~30).",

"jac_out" : 
"Coordinates of IPEC output files.\n  Choose from hamada, pest, boozer, equal_arc or other",

"jsurf_out" :
"Set to 1 for a surface weighted spectrum b_{mn}=\oint\left(\delta\mathbf{B}\cdot\hat{n}\right)\left(\theta,\phi\right)e^{-i\left(m\theta-n\phi\right)}J\cdot|\\nabla\psi|d\theta d\phi (an invariant flux on rational surfaces). Set to 0 for unweighted spectrum b_{mn}=\oint\left(\delta\mathbf{B}\cdot\hat{n}\right)\left(\theta,\phi\right)e^{-i\left(m\theta-n\phi\right)}d\theta d\phi. ",

"tmag_out" : 
"0 for ordinary toroidal angle. 1 for magnetic coordinate angle",

"nr" :
'Sets the r-resolution of the data output on a (r,z,\phi) grid.',

"nz" :
"Sets the z-resolution of the data output on a (r,z,\phi) grid.",

"fun_flag" :
"Outputs in functional format (on a coordinate grid) by back transforming the harmonic information (WARNING: this will result in larger data output).",

"flux_flag" :
"New flux surfaces \psi(r+\delta r,z+\delta z,\phi) for visualization of perturbation.",

"bin_flag" :
"Output binary files for use with the xdraw command when appropriate.",

"bin_2d_flag" :
"Output binary files containing 2D contours for use with the xdraw command  when appropriate.",

"resp_flag" :
"Output ipec_response_n#.out file.",

"singfld_flag" :
"Output the q, \Psi_{N}, resonant \mathbf{B}, singular current, island width and Chirikov parameter at rational surfaces inside the plasma, including the plasma response.",

"vsingfld_flag" :
"Output the \psi_{N}, resonant \mathbf{B}, singular current, island width and Chirikov parameter at rational surfaces inside the plasma, using vacuum calculations.\n  Requires coil_flag=True.",

"pmodb_flag" :
"Plasma |\delta B| on flux surfaces.",

"xbnormal_flag" :
"Output the normal components of the displacement, covarient perturbed magnetic field, and contravarient perturbed field (area weighted -> can compare between coordinates)  with plasma response included.",

"vbnormal_flag" :
"Output the normal components of the displacement, covarient perturbed magnetic field, and contravarient perturbed field (area weighted -> can compare between coordinates) from vacuum calculations.\n  Requires coil_flag=True.",

"eqbrzphi_flag" :
"Output the original equilibrium field on a (r,z,\phi) grid.",

"brzphi_flag" :
"Output the (r,z,\phi) components of the perturbed magnetic fields everywhere on the (r,z,\phi) grid.\n  - IPEC 1 only outputs perturbed magnetic fields valid inside the plasma on the (r,z,\phi) grid.\n  - IPEC 2.0 and later includes the plasma response and vacuum fields when the coil_flag is the vacuum input.",

"xrzphi_flag" :
"Output the plasma displacement on a (r,z,\phi) grid.",

"vbrzphi_flag" :
"Output the field due to the surface current defining the IPEC final solution boundary condition at the plasma surface (original external field boundary condition and plasma response) on an r,z grid. This does not represent a true field.",

"vpbrzphi_flag" :
"OBSOLETE.\n  - For IPEC 1. this output the plasma response field in the vacuum region on a (r,z,\phi) grid.",

"vvbrzphi_flag" :
"Outputs the field due to the surface current defining the external 3D field boundary condition on a (r,z,\phi) grid. This field does not represent a true field.",

"ntv_flag" : 
"Output first order quantities to bin files for interfacing with PENT code"
}


coil = {
"ceq_type" :
"This should be identical to the eq_type specified in equil.in",
"machine" :
"As of IPEC 3.00 supported machines include 'nstx', 'd3d', 'kstar', and 'iter'.",
"ip_direction" :
"Set as 'positive' (default) or 'negative' for CCW or CW from a top down view respectively.",
'bt_direction' :
'Set as "positive" (default) or "negative" for CCW or CW from a top down view respectively.',
'coil_num' :
'Total number of coil sets to be activated.',
'coil_name(1)' :
'Array values should be specified for each of the coil arrays to be activated. For example "coil_name(1) = c". The supported coil sets (and number of coils in each) as of IPEC 3.00 include,\n  - for nstx: "rwmef" (6), "tfef" (12), "pf5ef" (2), "ppu" (12), "ppl" (12), "psu" (12), "psl" (12), "vpu" (12), "vpl" (12), "vsu" (12), "vsl" (12), "hhfw" (24)\n  - for d3d: "iu" (6), "il" (6), "c" (6), "tbm_solenoid" (1), "tbm_racetrack" (1)\n  - for kstar: "fecu" (4), "fecm" (4), "fecl" (4)\n  -  for iter: "efcu" (6), "efcm" (6), "efcl" (6), "bl2u" (9), "bl2m" (9), "bl2l" (9), "avvu" (9), "avvm" (9), "avvl" (9)',
"coil_cur(1,1)" :
"Array of Ampere coil current for the nc-th coil in the set corresponding to coil_name(nc) (default  = 0).",
}

pent = {
"kinetic_file" : "1 header line followed by column data \nfor psi_n, ni(m^-3), ne(m^-3), ti(eV), te(eV), omega_E(rad/s).",
"o1fun_file" : "IPEC output generated using ntv_flag or ascii B(psi,theta).",
"o1mn_file" : "IPEC output generated using ntv_flag, or ascii B(m,n).",
"o1native" : "o1fun_file and o1mn_file are IPEC output bin files.",
"imass" : "Mass of main species in atomic mass units (0 for electrons).",
"icharge" : "Charge of main species in elementry charge units.",
"zimp" : "Charge of impurity species in elementry charge units.",
"zmass" : "Mass of impurity species in atomic mass units.",
"nl" : "Number of bounce harmonics calculated (-nl<=l<=nl).",
"collision" : "Collision operator. \nChoose from zero,small,krook, or harmonic.",

"wefac" : "Multiplier applied to omega_E profile.",
"wdfac" : "Multiplier applied to oemga_D.",
"ptfac" : "Fractional cutoff in Lambda space near trapped-passing boundry.",
"kolim_flag" : "Infinite rotation limit.",
"neorot_flag" : "Approximate numerator of resonant operator using neoclassical rotation.",
"divxprp_flag" : "Include arc-length change in perturbed Action.",
"xlsode_flag" : "Use LSODE package to integrate in energy space.",
"lmdalsode_flag" : "Use LSODE package to integrate in pitch angle.",
"lsode_rtol" : "Relative tolarance in LSODE.",
"lsode_atol" : "absolute tolarance in LSODE.",
"xmax" : "Upper limit of energy integration.",
"ximag" : "Step contour of integration off of real axis to avoid poles when collision is zero.",
"nx" : "Number of energy spline points if xlsode_flag false.",
"nt" : "Number of theta points used in bounce averaging.",
"nlmda" : "Number of Lambda points used in splines.",
"maskpsi" : "Step size taken by PENT on dcon psi grid (determined by idconfile).",

"passing_flag" : "Calculate torque and energy contributed by passing particles.",
"gen_flag" : "Perform full calculation in general geometry.",
"rlar_flag" : "Perform Reduced LAR approximate calculation.",
"clar_flag" : "Perform Circular LAR approximate calculation.",
"cgl_flag"  : "Fluid Chew-Goldberger-Low solution (omega_E -> inf).",
"kinetic_flag" : "Output kinetic profiles.",
"bounce_flag" : "Output sample bounce averaging profiles.",
"pitch_flag" : "Output sample pitch angle profiles.",
"energy_flag" : "Output sample energy profiles.",
"psiout(1)" : "Include this surface in sampled outputs.",
"psiout(2)" : "Include this surface in sampled outputs.",
"psiout(3)" : "Include this surface in sampled outputs.",
"psiout(4)" : "Include this surface in sampled outputs.",
"psiout(5)" : "Include this surface in sampled outputs.",
"psiout(6)" : "Include this surface in sampled outputs.",
"psiout(7)" : "Include this surface in sampled outputs.",
"psiout(8)" : "Include this surface in sampled outputs.",
"psiout(9)" : "Include this surface in sampled outputs.",
"psiout(10)" : "Include this surface in sampled outputs.",
"psiout" : "Include these 10 surfaces in sampled outputs.",

"fmnl_flag" : "Do not change unless you know what you are doing.",
"diag_flag" : "Do not change unless you know what you are doing."
}


alltips = {}

for k,v in equil.iteritems():
	alltips[k] = v
for k,v in dcon.iteritems():
	alltips[k] = v
for k,v in ipec.iteritems():
	alltips[k] = v
for k,v in coil.iteritems():
	alltips[k] = v
for k,v in pent.iteritems():
	alltips[k] = v
