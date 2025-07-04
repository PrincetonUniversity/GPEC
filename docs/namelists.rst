
Fully annotated namelist input files can be found in the $GPECHOME/input directory and any of the examples in $GPECHOME/docs/examples directories of any `public installation <installation>`.

The major input files requiring user interaction are as follows.

The EQUIL Namelist
===================

This namelist specifies the equilibrium inputs and how equilibrium quantities are splined for use in the various executables.

.. include:: ../input/equil.in
   :literal:

The DCON Namelist
==================

This namelist controls the behavior of the DCON code and determines its outputs.

.. include:: ../input/dcon.in
   :literal:

.. note:: Additional DCON documentation can be found in the pre-GPEC documentation :ref:`here <dcon>`.

The GPEC Namelist
==================

This namelist sets the inputs, controls the behavior, and determines the outputs of the GPEC code.

.. include:: ../input/gpec.in
   :literal:

The PENTRC Namelist
====================

This namelist sets the inputs, controls the behavior, and determines the outputs of the GPEC code.

.. include:: ../input/pentrc.in
   :literal:

