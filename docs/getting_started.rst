.. _getting_started:


***************
Getting started
***************

.. _installing:

.. include:: ../install/README

.. _Run_Process:

Running using OMFIT
====================

`OMFIT <http://gafusion.github.io/OMFIT-source/>`_ is the recommended way to run GPEC for new users. Public releases of OMFIT are available on most major fusion clusters, or it can be downloaded to your personal computer.

Using OMFIT guides users through a clean workflow, providing many powerful tools to facilitate the necessary preparation of inputs. It provides the interface, for example, with EFIT and profile fitting tools on many machines. It also enables users to submit GPEC jobs remotely, quickly visualize outputs, perform powerful python post-processing, and share projects with others.

A complete description of the OMFIT workflow manager GPEC module can be found in `this google doc tutorial <https://docs.google.com/document/d/1qSUjJZYmET8-X08rRGBT6L_9D9fzEN_zCV75OX_ZbZk/edit?usp=sharing>`_.

Manual Run Process
==================

This section outlines the basic run process with reference to some of the most common input-output manipulation in GPEC runs. This is not meant to be a complete coverage of the many options and capabilities available to the user during this process. For a complete list of the input and output control variables please see the Namelist Variables section.

Public copies of the GPEC releases are available at /p/gpec on portal. The user should usually choose the latest release, lets call it <release>. To run GPEC, a user should first navigate to a directory in which they wish to run and write outputs. Once in the run-directory the user needs (at a minimum) the dcon, equil, vac, and gpec namelist files, which can be copied from /p/gpec/<release>/input for example. The user can then follow these simple steps to complete a GPEC run:

1. Edit the dcon.in file

  - DCON will only run for a single toroidal mode number. Set the toroidal mode number for the entire run here using the crucial parameter nn. GPEC outputs for multiple toroidal mode numbers can be added together in the post-processing.

2. Edit the equil.in file

  - Specify the specific plasma shot of interest by setting the initial equilibrium file (ex: EFIT output g-file).

  - Set the range of normalized flux to be included in the calculation.

3. $/p/gpec/<release>/bin/dcon

  - This creates the euler.bin file, which provides a orthonormal basis of force balance states for GPEC.

  - Once this is created, you can various effects on this equilibrium by changing the 3D field input in gpec.in without having to run DCON again.

4. Edit the gpec.in file

  - Specify the 3D error field input file to be used in GPEC to perturb the 2D equilibrium.

    - This must match the psihigh variable used in equil.in.

  - Set the flags to specify the desired outputs.

5. /p/gpec/<release>/bin/gpec

  - That’s it! Now go look at your output files.

Note that GPEC takes only two inputs. The file containing the 2D plasma equilibrium is specified in equil.in prior to the DCON run and the file defining the 3D error field on the plasma boundary surface is specified in GPEC.in prior to the the GPEC calculations.

.. _coordinates:

A Word On Coordinates
=====================

Keeping all the different coordinate systems used in a GPEC run can be intimidating to the first time user. Luckily, there is no reason to fear! GPEC easily switches between coordinate systems, and does not require consistency throughout the run process. That is, the inputs and outputs can have different coordinates. The coordinate system used for the DCON run specify the working coordinates for GPEC. The jac_in specified in gpec.in specifies the coordinate system used by the 3D error field input file, and jac_out specifies the coordinate system of the GPEC outputs. GPEC handles the transitions between these three coordinate systems all internally.

.. note:: The only real responsibility of the user is to make sure the GPEC input Jacobian is consistent with the 3D error field input file. However, the user should be aware that each coordinate system has different strengths. Hamada coordinates have the best convergence properties in general, but are poor for the inboard side of the plasma. PEST coordinates are user friendly but poor for outboard resolution. Boozer coordinates tend to find a middle ground between these, with mild convergence and fair resolution throughout the plasma.

The above requirements are all that is needed to understand the internal workings of GPEC's coordinate transforms. There is, however, a final question of the various “normal” toroidal angles as defined for a specific machine. GPEC does not take into account the sign of the toroidal field. The result is, that for machines such as DIII-D with “normal” toroidal field defined as negative (BCENTER<0 in efit equilibria) the phase of GPEC outputs will not match the machine-defined toroidal angle. To convert three dimensional GPEC outputs in such cases to machine coordinates it may be necessary to take the complex conjugate of a given amplitude and multiply equilibrium output toroidal fields by negative one.

.. _inputs:

Inputs
======

At this point it is obvious that GPEC calculations hinge on first acquiring two types of data. GPEC interfaces with the input data by simple reads of formatted text files. The format and file path's are specified in the .in files used in the run process, and are called on to initialize the ideal equilibrium solver.

The first input, equilibrium profile data, is specified in the equil.in file and can come in many different forms. By far the most common option is “efit”, in which case a two dimensional EFIT output g-file is used to specify the original axi-symmetric equilibrium. This equilibrium data and file format are widely used across multiple machines, and alternative formats are supported by GPEC when necessary. It is worth noting that significant differences between kinetic EFIT equilibria and the more simple non-kinetic EFITs has been noted. Obtaining the best 2D equilibrium data possible is encouraged for detailed calculations, however the casual user should not have need to concern themselves too heavily with this input source short of looking for the kinetic EFIT when possible.

The second input, the 3D “error field” data, is somewhat more specific in its requirements. The input data can be again read from a file, and accepts a small number of select formats. These files must specify the external field on the surface of the GPEC plasma boundary defined by the psihigh variable in equil.in. This means that the data must be pre-processed using knowledge of both the non-axisymmetric external field (applied and/or intrinsic), the same equilibrium used in the GPEC corresponding calculation, and the exact magnetic surface used in that calculation. This requirement of coordination by the user is an easy pitfall for the un-wary. If utilizing this method, be sure that all work is consistent! To avoid the pitfalls associated with the above method a in-house Biot-Savart module that contains various coil geometries has been added as of GPEC 2.00. Using the coil_flag in gpec.in will result in the reading of yet another input file “coil.in” and the fields from the coils/currents specified therein will be calculated on the appropriate surface. Finally, the harmonic_flag allows users to apply a set amplitude harmonic on the boundary. All the 3D field input may be used (or not) in any combination, and the total boundary condition will be a sum of those used.

.. _results:

Results
=======

The results of an GPEC run depend on the output flags selected in gpec.in prior to the run. Without getting into the specifics of each output, there are a few important properties any GPEC run that should always be kept in minds. Each GPEC run always returns its results in the form of consistently named .out files in the run-directory. This means that it will overwrite any previous run's results if they are not moved first. Always move any results worth keeping to a separate directory and keep track of which runs used which inputs. Outputs are expensive, so try to limit runs to only output the essential physics of interest each time.

An important property of GPEC results for post-processing is their linearity. GPEC calculates results for a single toroidal harmonic. When appropriate, results from multiple runs can simply be added (to obtain a total perturbed field for example). In the same spirit, the amplitude of a plasma response calculated by GPEC is linear with the error-field. Thus, plasma response fields can be multiplier by a constant to simulate changes in a driven 3D error field.
