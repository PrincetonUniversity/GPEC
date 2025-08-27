# Rules for building sub-modules

# EQUIL
FCEQUIL = $(FC) $(FFLAGS) $(RECURSFLAG)
../equil/%.o: ../equil/%.f
	$(FCEQUIL) -c $<
../equil/%.o: ../equil/%.F
	$(FCEQUIL) -c $<
../equil/%.o: ../equil/%.f90
	$(FCEQUIL) -c $<
../equil/%.o: ../equil/%.F90
	$(FCEQUIL) -c $<

# LSODE
FCLSODE = $(FC) $(FFLAGS) $(RECURSFLAG) $(LEGACYFLAG) $(OMPFLAG) -I../equil
../lsode/%.o: ../lsode/%.f
	$(FCLSODE) -c $<
../lsode/%.o: ../lsode/%.F
	$(FCLSODE) -c $<

# ZLANGE
FCZLANGE = $(FC) $(FFLAGS) $(RECURSFLAG) $(LEGACYFLAG) -I../equil
../zlange/%.o: ../zlange/%.f
	$(FCZLANGE) -c $<
../zlange/%.o: ../zlange/%.F
	$(FCZLANGE) -c $<

# ZVODE
FCZVODE = $(FC) $(FFLAGS) $(RECURSFLAG) $(LEGACYFLAG) -I../equil
../zvode/%.o: ../zvode/%.f
	$(FCZVODE) -c $<
../zvode/%.o: ../zvode/%.F
	$(FCZVODE) -c $<

# ORBIT
FCORBIT = $(FC) $(FFLAGS) $(OMPFLAG) $(RECURSFLAG) -I$(MATHINC) -I../equil
../orbit/%.o: ../orbit/%.f
	$(FCORBIT) -c $<
../orbit/%.o: ../orbit/%.F
	$(FCORBIT) -c $<

# VACUUM
FCVACUUM = $(FC) $(FFLAGS) $(ZEROFLAG) $(RECURSFLAG) $(LEGACYFLAG) -I../equil
../vacuum/%.o: ../vacuum/%.f
	$(FCVACUUM) -c $<
../vacuum/%.o: ../vacuum/%.F
	$(FCVACUUM) -c $<

# PENTRC
FCPENTRC = $(FC) $(FFLAGS) $(OMPFLAG) $(RECURSFLAG) -I$(MATHINC) -I$(NETCDFINC) -I../equil -I../harvest
../pentrc/%.o: ../pentrc/%.f
	$(FCPENTRC) $(LEGACYFLAG) -c $<
../pentrc/%.o: ../pentrc/%.F
	$(FCPENTRC) $(LEGACYFLAG) -c $<
../pentrc/%.o: ../pentrc/%.f90
	$(FCPENTRC) -c $<
../pentrc/%.o: ../pentrc/%.F90
	$(FCPENTRC) -c $<

# DCON
FCDCON = $(FC) $(FFLAGS) $(OMPFLAG) $(RECURSFLAG) $(LEGACYFLAG) -I$(MATHINC) -I$(NETCDFINC) -I../equil -I../harvest -I../pentrc -I../vacuum
../dcon/%.o: ../dcon/%.f
	$(FCDCON) -c $<
../dcon/%.o: ../dcon/%.F
	$(FCDCON) -c $<

# RDCON
FCRDCON = $(FC) $(FFLAGS) $(OMPFLAG) $(LEGACYFLAG) -I$(MATHINC) -I$(NETCDFINC) -I../equil -I../vacuum
../rdcon/%.o: ../rdcon/%.f
	$(FCRDCON) -c $<
../rdcon/%.o: ../rdcon/%.F
	$(FCRDCON) -c $<

# STRIDE
FCSTRIDE = $(FC) $(FFLAGS) $(OMPFLAG) $(RECURSFLAG) $(LEGACYFLAG) -I$(MATHINC) -I$(NETCDFINC) -I../equil -I../vacuum
../stride/%.o: ../stride/%.f
	$(FCSTRIDE) -c $<
../stride/%.o: ../stride/%.F
	$(FCSTRIDE) -c $<

# MATCH
FCMATCH = $(FC) $(FFLAGS) $(OMPFLAG) $(RECURSFLAG) -I../equil
../match/%.o: ../match/%.f
	$(FCMATCH) -c $<
../match/%.o: ../match/%.F
	$(FCMATCH) -c $<

# RMATCH
FCRMATCH = $(FC) $(FFLAGS) $(RECURSFLAG) $(LEGACYFLAG) $(OMPFLAG) -I$(MATHINC) -I../equil
../rmatch/%.o: ../rmatch/%.f
	$(FCRMATCH) -c $<
../rmatch/%.o: ../rmatch/%.F
	$(FCRMATCH) -c $<

# SUM
FCSUM = $(FC) $(FFLAGS) -I../equil
../sum/%.o: ../sum/%.f
	$(FCSUM) -c $<
../sum/%.o: ../sum/%.F
	$(FCSUM) -c $<

# MULTI
FCMULTI = $(FC) $(FFLAGS) -I../equil
../multi/%.o: ../multi/%.f
	$(FCMULTI) -c $<
../multi/%.o: ../multi/%.F
	$(FCMULTI) -c $<

# COIL
FCCOIL = $(FC) $(FFLAGS) $(OMPFLAG) $(RECURSFLAG) -I$(MATHINC) -I$(NETCDFINC) -I../equil
../coil/%.o: ../coil/%.f
	$(FCCOIL) -c $<
../coil/%.o: ../coil/%.F
	$(FCCOIL) -c $<

# SLAYER
FCSLAYER = $(FC) $(FFLAGS) $(OMPFLAG) -I$(MATHINC) -I$(NETCDFINC) -I../equil
../slayer/%.o: ../slayer/%.f
	$(FCSLAYER) -c $<
../slayer/%.o: ../slayer/%.F
	$(FCSLAYER) -c $<

# GPEC
FCGPEC = $(FC) $(FFLAGS) $(OMPFLAG) $(RECURSFLAG) $(LEGACYFLAG) -I$(MATHINC) -I$(NETCDFINC) -I../equil -I../vacuum -I../harvest -I../pentrc -I../coil -I../slayer
../gpec/%.o: ../gpec/%.f
	$(FCGPEC) -c $<
../gpec/%.o: ../gpec/%.F
	$(FCGPEC) -c $<

# XDRAW
CCXDRAW = $(CC) $(X11_INC_DIR) $(CFLAGS) -DUNIX -O2
../xdraw/%.o: ../xdraw/%.c
	$(CCXDRAW) -c $<
../xdraw/%.o: ../xdraw/%.C
	$(CCXDRAW) -c $<