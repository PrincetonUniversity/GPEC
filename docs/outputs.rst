Outputs and Post-Processing
***************************

Output Files
============

The IPEC outputs are entirely specified by flags (bool types t or f) set in the IPEC_output section of IPEC.in. All outputs are ASCII files, and can be grouped into two major categories.

ASCII File Outputs
------------------

A number of the flag options in IPEC.in instruct IPEC to output ASCII file data. Some of these outputs are always available. Some, however, require a input error field instead of a hard coded harmonic_flag call. Both groups are listed in detail here.

Outputs Always Available

These IPEC outputs can always be obtained from a equil.bin file output from DCON. 

IPEC_response_n#.out

    **Flag** resp_flag

    **Info** Energy for DCON eigenmodes and stability indices. Eigenvalues and eigenvectors for vacuum and plasma inductance (virtual casing currents to fields), plasma permeability (external fields to total fields), and plasma reluctance.

IPEC_singcoup_matrix_n#.out

    **Flag** singcoup_flag

    **Info** The coupling matrix to resonant fields, coupling matrix to singular currents, and to island half-widths is given for each rational surface within the plasma (q=2, 3, etc) for each surface the real and imaginary coupling constants are given for each poloidal mode number on the control surface.

IPEC_singcoup_svd_n#.out

    **Flag** singcoup_flag

    **Info** The SVD singular values (s) and eigen vectors for each coupling matrix in IPEC_singcoup_matrix_n#.out. Large s corresponds to large amplification, with the largest (most important mode) listed at the top. The results should be dotted with the unweighted normal field spectrum to give physical meaning.

Outputs Available When Error Field is Provided

These outputs are only available when an external error field file is provided as an input to IPEC. This means IPEC.in must have the data_flag turned on and a infile specified.

IPEC_control_n#.out

    **Flag** 

    **Info** The Plasma response for an external perturbation on the control surface. This includes the vacuum energy, surface energy, plasma energy, real and imaginary vacuum input \mathbf{B}_{in} and total field on plasma boundary\mathbf{B}_{out}as a function of poloidal mode number.

IPEC_singfld_n#.out

    **Flag** singfld_flag

    **Info** The \Psi_{N}, total resonant \mathbf{B} (real and imaginary), singular current (real and imaginary), island half width in units of normalized flux and Chirikov parameter at rational surface inside the plasma. 

    **Flag** singcoup_flag

    **Info** Additional section showing the overlap field and overlap percentage for each eigenmode in the singcoup_svd output.

IPEC_pmod_n#.out

    **Flag** pmodb_flag

    **Info** Eulerian and Lagrangian \left|\mathbf{B}\right|(real and imaginary) for each poloidal mode number at each value of \Psi_{N} output. This output is necessary for NTV post processing.

IPEC_xbnormal_n#.out

    **Flag** xbnormal_flag

    **Info** The normal components of the displacement, magnetic field without the plasma response, and magnetic field with the plasma response included for each poloidal mode number at each value of \Psi_{N} output.??

IPEC_*rzphi_n#.out

A number of output files have a similar structure. Here the * in the file name is replaced by the appropriate leading letters of the corresponding flag. For example the xrzphi_flag for n=1 creates a IPEC_xrzphi_n1.out file. some common properties of these files are:

• real and imaginary components: Output files contain two dimensional data on an \left(r,z\right) grid for a single toroidal harmonic. To translate into three dimensions, perform the transformationB\left(r,z,\phi\right)=B_{real}\left(r,z\right)\cos\left(n\phi\right)+B_{imag}\left(r,z\right)\sin\left(n\phi\right)
 

• l parameter: 1 designates points inside plasma, 0 points in vacuum, -1 points near/on surface (singularity)

    **Flag** eqbrzphi_flag

    **Info** The original equilibrium field on a \left(r,z,\phi\right) grid.

    **Flag** brzphi_flag

    **Info** The \left(r,z,\phi\right)components of the perturbed magnetic fields inside the plasma on the \left(r,z,\phi\right) grid.

    **Flag** xrzphi_flag

    **Info** The displacement on the \left(r,z,\phi\right) grid.

    **Flag** vbrzphi_flag

    **Info** The false perturbed magnetic field in the vacuum region on the \left(r,z,\phi\right) grid calculated using the IPEC boundary surface current composed of both the vacuum component and the plasma response.

    **Flag** vpbrzphi_flag

    **Info** The true perturbed magnetic field in the vacuum region on the \left(r,z,\phi\right) grid due to the plasma response alone calculated from the plasma response surface condition.

    **Flag** vvbrzphi_flag

    **Info** The false perturbed magnetic field in the vacuum region on the \left(r,z,\phi\right) grid calculated using the IPEC boundary surface current from the external fields alone.

    **Flag** ssbrzphi_flag

    **Info** The false perturbed magnetic field in the vacuum region on the \left(r,z,\phi\right) grid calculated using the IPEC boundary surface current from the external fields alone.

Binary File Outputs
-------------------

These files are designed for quick and easy visualization of results using the xdraw command. For more details on using xdraw see the devoted section on this page.

xbnormal.bin

    **Flag** bin_flag

    **Info** The normal displacement and magnetic field as functions of \Psi_{N} for xdraw.

xbnormal_2d.bin

    **Flag** bin_2d_flag

    **Info** Contour profiles of the normal displacement and magnetic field in (R,z) for xdraw.

pflux_re(im)_2d.bin

    **Flag** bin_2d_flag

    **Info** Contour profiles of the real (imaginary) perturbed flux in (R,z) for xdraw.

bnormal_spectrum.bin

    **Flag** bin_flag

    **Info** Surfmn type contours of the normal perturbed magnetic fields as a function of poloidal harmonic number and \Psi_{N}.

xdraw
=====

The binary IPEC outputs can be viewed using the commandxdraw filenamewhere filename is one of the .bin files created by IPEC (“.bin” excluded). This is a quick way to view results immediately as they are produced. The xdraw tool provides a highly interactive environment that takes keystroke inputs to change plot display options, navigate plots, display single or multiple responses at once, do limited post processing (get a gradient, or ratio), and save figures. For a full list of the command options, enter the xdraw environment and press “k”.

IDL NTV Post Processing
=======================

The calculation of neoclassical toroidal viscosity (NTV) requires knowledge of many plasma properties. These the 3D magnetic field within the plasma \delta B
  given by IPEC, as well as external variables such as pressure profiles, velocity profiles, and more that are not included in IPEC itself. To calculate an estimate of the NTV torque it is necessary to specify additional profile information with either measured data or analytical formulae. Either way, it should be emphasized that the calculation is a post-processing of IPEC data that is not self consistent with the IPEC solution. This is a reasonable estimate as long as the torque is not great enough to significantly effect the equilibrium.

A review of the NTV calculation and definitions of the involved variables can be found in Ref. [Park, 2009]. For the detailed mathematics, although important, are not the subject of this manual. The important highlights for the user to keep in mind are the following: 

• Necessary profiles include the toroidal rotationV_{\phi}, pressure P, and magnetic field B.

– Force balance and the neoclassical value of the poloidal velocity are used to calculate the radial electric field.qn\left(E_{r}-v_{\phi}B_{\theta}+v_{\theta}B_{\phi}\right)	=	\frac{\partial P}{\partial r}  v_{\theta}	=	\frac{c_{p}}{eB_{\theta}}\frac{\partial T_{i}}{\partial r}\approx\frac{1.17}{eB_{\theta}}\frac{\partial T_{i}}{\partial r}
 

– Specifically, the experimental profiles used are ion temperature T_{i}, electron temperature T_{e}, electron density n_{e}, impurity density n_{iz}, the effective charge Z_{eff}, and toroidal rotation V_{\phi}. Optionally the E\times B rotation and deuterium rotation can be specified as well.

• The integral over \kappa^{2}=\left[E-\mu B_{0}\left(1-\epsilon\right)\right]/2\mu B_{0}\epsilon is approximated away. ?

• There are multiple distinct regimes of collisional and methods of dealing with them. In the post processing these are output as variables with the prefixes

– CP \rightarrow Collisional Plateau regimes transport estimated by Shaing [Zhu2006, and references therein].

– IN \rightarrow Generalized combination of the 1/\nu and \nu regimes. This is the most trusted result at the time of this manual's writing.

Run Instructions
----------------

The first step in an NTV calculation for any given shot is to run the corresponding IPEC calculation. The output files necessary are the control, singfld, and pmodb, files. These should use the functional outputs. 

Having obtained the required output files, the user should enter the IDL environment and run one of the two NTV routines. These routines are 'total_energy_torque.pro' and 'simple_energy_torque.pro'. The first requires full experimental data and should be initialized in the package script found at '/u/jpark/analysis/lib/init_routines_jkp.idl', while the second defines analytic profiles for the necessary variables in order to calculate results from IPEC outputs and can be used on its own. 

An example of the run process within the IDL environment is shown here::

  IDL> @/u/jpark/analysis/lib/init_routines_jkp.idl
  IDL> gfile='g145117.03600_d3d_kinetic' ;efit equilibrium
  IDL> kdir = 'samplepath' ;directory of GA kinetic file 
  IDL> path1 = 'samplepath' ;directory of equilibrium file
  IDL> path2 = 'samplepath' ;directory of IPEC outputs 
  IDL> numstr = '145117.03600'	;shot.time
  IDL> t = total_energy_torque(nn=3,/adata,/atime, kdir=kdir,gfile=gfile,/diamag,/magpre,/counter,path1=path1,path2=path2,nl=8,numstr=numstr)

The result can be saved using::

  IDL> save,file='t145117.03600_d3d_kinetic_ievenc1_counter.tsav',t

or the user's choice of save file path. The result holds a huge amount of information. The structure can be found using the command::

  IDL> help, /str, tor 
  IDL> help, /str, t.TOT.L.NTVetc.

This example uses a number of specific call words, and there are many more that can be examined by looking into the nuts and bolts of the script 'total_energy_torque.pro'. So as not to leave the reader hanging, we will quickly review the calls used above. The keywords '/adata' and '/atime' tell the routine that it will be using data from Andrea Garofalo, the structure of which has been hard-coded into the routine and would need to be edited for any new data format. The keyword '/counter' changes the sign of some of Andrea Garofalo's data to be consistent with other conventions. The keywords '/diamag' and '/magpre' tell the routine to use the diamagnetic rotation and magnetic precession respectively, and should be called for almost all non-debugging calculation purposes. Finally, the variable 'nl' specifies the number of bounce harmonics to include in the calculation. There are many other arguments that could be passed to the total_energy_torque such as specifying a separate machine, or individual profiles. However, the user will be required to familiarize themselves with the internal structures of the function and the data before extending themselves to use its full functionality.
