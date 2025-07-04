
.. include:: ../install/README

.. _Run_Process:

Running using OMFIT
====================

`OMFIT <http://gafusion.github.io/OMFIT-source/>`_ is a framework for integrated modeling. It contains a GPEC module for preparing, running, and analyzing GPEC runs that is highly recommended for pure users (non-developers). OMFIT is the recommended way to run GPEC for new users. Public releases of OMFIT are available on most major fusion clusters, or it can be downloaded to your personal computer.

Using OMFIT guides users through a clean workflow, providing many powerful tools to facilitate the necessary preparation of inputs. It provides the interface, for example, with EFIT and profile fitting tools on many machines. It also enables users to submit GPEC jobs remotely, quickly visualize outputs, perform powerful python post-processing, and share projects with others.

A complete description of the OMFIT GPEC module can be found on the `module documentation page <http://gafusion.github.io/OMFIT-source/modules/mod_GPEC.html#mod-gpec>`_ and a thorough discussion of the workflow, inputs, outputs, and post-processing can be found in `the corresponding google doc tutorial <https://docs.google.com/document/d/1qSUjJZYmET8-X08rRGBT6L_9D9fzEN_zCV75OX_ZbZk/edit?usp=sharing>`_.



Manual Run Process
==================

This section outlines the basic run process with reference to some of the most common input-output manipulation in GPEC runs. This is not meant to be a complete coverage of the many options and capabilities available to the user during this process. For a complete list of the input and output control variables please see the :ref:`Inputs <inputs>` section.

.. note:: The use of $GPECHOME and executable commands such as `dcon`, `gpec` or `pentrc` assumes the appropriate environement variables have been set. The easiest way to do this is to use the module commands associated with one of the public installations.

To run GPEC, a user should first navigate to a directory in which they wish to run and write outputs. Once in the run-directory the user needs (at a minimum) the dcon, equil, vac, and gpec namelist files, which can be copied from $GPECHOME/input for example. The user can then follow these simple steps to complete a GPEC run:

1. Edit the dcon.in file

   - DCON will only run for a single toroidal mode number. Set the toroidal mode number for the entire run here using the crucial parameter nn. GPEC outputs for multiple toroidal mode numbers can be added together in the post-processing.

2. Edit the equil.in file

   - Specify the specific plasma shot of interest by setting the initial equilibrium file (ex: EFIT output g-file).

   - Set the range of normalized flux to be included in the calculation.

3. $GPECHOME/bin/dcon

   - This creates the euler.bin file, which provides a orthonormal basis of force balance states for GPEC.

   - Once this is created, you can various effects on this equilibrium by changing the 3D field input in gpec.in without having to run DCON again.

4. Edit the gpec.in file

   - Specify the 3D error field input file to be used in GPEC to perturb the 2D equilibrium.

     - This must match the psihigh variable used in equil.in.

   - Set the flags to specify the desired outputs.

5. $GPECHOME/bin/gpec

   - That’s it! Now go look at your output files.

Note that GPEC takes 2-3 inputs depending on whether kinetic effects are included. All calculations require a file containing the 2D plasma equilibrium, which is specified in equil.in prior to the DCON run. The most common equilibrium source is “efit”, in which case a two dimensional EFIT output g-file is used to specify the original axi-symmetric equilibrium. This equilibrium data and file format are widely used across multiple machines, but many alternative formats are also supported by GPEC. It is worth noting that significant differences between kinetic EFIT equilibria and the more simple non-kinetic EFITs has been noted. Obtaining the best 2D equilibrium data possible is encouraged for detailed calculations.

If kinetic effect are being included (via kin_flag in DCON or eny use of PENTRC), the kinetic profiles are also required. The `OMFIT tutorial <https://docs.google.com/document/d/1qSUjJZYmET8-X08rRGBT6L_9D9fzEN_zCV75OX_ZbZk/edit?usp=sharing>`_ details the necessary profiles and ascii file format required. Whenever possible, it is encouraged to use profiles consistent with the 2D equilibrium (especially if it is a kinetic EFIT).

The final input is the 3D “error field” data. The input data can be again read from a file, and accepts a small number of select formats. These files must specify the external field on the surface of the GPEC plasma boundary defined by the psihigh variable in equil.in. This means that the data must be pre-processed using knowledge of both the non-axisymmetric external field (applied and/or intrinsic), the same equilibrium used in the GPEC corresponding calculation, and the exact magnetic surface used in that calculation. This requirement of coordination by the user is an *easy pitfall for the un-wary*. If utilizing this method, be sure that all work is consistent! To avoid the pitfalls associated with the above method a in-house Biot-Savart module that contains various coil geometries has been added as of GPEC 2.00. Using the coil_flag in gpec.in will result in the reading of yet another input file “coil.in” and the fields from the coils/currents specified therein will be calculated on the appropriate surface. Finally, the harmonic_flag allows users to apply a set amplitude harmonic on the boundary. All the 3D field input may be used (or not) in any combination, and the total boundary condition will be a sum of those used.

Manual Examples
---------------

The GPEC repository contains a number of examples intended to provide both consistent regression testing and a useful reseource for users looking to manually run straightforward examples. These examples are located in $GPECHOME/docs/examples, and contain all the necessary inputs and namelists to run the appropriate suite of codes. Ideal and kinetic examples use dcon, gpec, and pentrc. Resistive examples contain the necessary inputs to run dcon, rdcon, rmatch and gpec.

Example console commands and output are given below.

- :ref:`DIII-D Ideal Example <DIIID_ideal_example>`
