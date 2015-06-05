COIL Namelist Inputs
********************

COIL_CONTROL
============

**ceq_type** = "efit"
  This should be identical to the eq_type specified in equil.in

**machine** = "d3d"
  As of IPEC 3.00 supported machines include 'nstx', 'd3d', 'kstar', and 'iter'.

**ip_direction** = "positive"
  Set as 'positive' (default) or 'negative' for CCW or CW from a top down view respectively.

**bt_direction** = "negative"
  Set as "positive" (default) or "negative" for CCW or CW from a top down view respectively.

**coil_num** = 3
  Total number of coil sets to be activated.

**cmpsi** = 64


**cmtheta** = 480


**cmzeta** = 40


**coil_name(1)** = "iu"
  Array values should be specified for each of the coil arrays to be activated. For example "coil_name(1) = c". The supported coil sets (and number of coils in each) as of IPEC 3.00 include,  - for nstx: "rwmef" (6), "tfef" (12), "pf5ef" (2), "ppu" (12), "ppl" (12), "psu" (12), "psl" (12), "vpu" (12), "vpl" (12), "vsu" (12), "vsl" (12), "hhfw" (24)  - for d3d: "iu" (6), "il" (6), "c" (6), "tbm_solenoid" (1), "tbm_racetrack" (1)  - for kstar: "fecu" (4), "fecm" (4), "fecl" (4)  -  for iter: "efcu" (6), "efcm" (6), "efcl" (6), "bl2u" (9), "bl2m" (9), "bl2l" (9), "avvu" (9), "avvm" (9), "avvl" (9)

**coil_cur(1,1)** = 2178.0
  Array of Ampere coil current for the nc-th coil in the set corresponding to coil_name(nc) (default  = 0).

**coil_cur(1,2)** = 1577.0


**coil_cur(1,3)** = -601.0


**coil_cur(1,4)** = -2178.0


**coil_cur(1,5)** = -1577.0


**coil_cur(1,6)** = 601.0


**coil_name(2)** = "il"


**coil_cur(2,2)** = -2178.0


**coil_cur(2,3)** = -1577.0


**coil_cur(2,4)** = 601.0


**coil_cur(2,5)** = 2178.0


**coil_cur(2,6)** = 1577.0


**coil_cur(2,1)** = -601.0


**coil_name(3)** = "c"


**coil_cur(3,1)** = -1551.0


**coil_cur(3,2)** = 268.0


**coil_cur(3,3)** = 1744.0


**coil_cur(3,4)** = 1551.0


**coil_cur(3,5)** = -268.0


**coil_cur(3,6)** = -1744.0




COIL_OUTPUT
===========

**ipec_interface** = .TRUE.




