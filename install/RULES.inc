# Rules for building sub-modules

# EQUIL
FCEQUIL = $(FC) $(FFLAGS) $(RECURSFLAG)
../equil/%.o: ../equil/%.f
	$(FCEQUIL) -J../equil -c $< -o $@
../equil/%.o: ../equil/%.F
	$(FCEQUIL) -J../equil -c $< -o $@
../equil/%.o: ../equil/%.f90
	$(FCEQUIL) -J../equil -c $< -o $@
../equil/%.o: ../equil/%.F90
	$(FCEQUIL) -J../equil -c $< -o $@

# LSODE
FCLSODE = $(FC) $(FFLAGS) $(RECURSFLAG) $(LEGACYFLAG) $(OMPFLAG) -I../equil
../lsode/%.o: ../lsode/%.f
	$(FCLSODE) -J../lsode -c $< -o $@
../lsode/%.o: ../lsode/%.F
	$(FCLSODE) -J../lsode -c $< -o $@

# ZLANGE
FCZLANGE = $(FC) $(FFLAGS) $(RECURSFLAG) $(LEGACYFLAG) -I../equil
../zlange/%.o: ../zlange/%.f
	$(FCZLANGE) -J../zlange -c $< -o $@
../zlange/%.o: ../zlange/%.F
	$(FCZLANGE) -J../zlange -c $< -o $@

# ZVODE
FCZVODE = $(FC) $(FFLAGS) $(RECURSFLAG) $(LEGACYFLAG) -I../equil
../zvode/%.o: ../zvode/%.f
	$(FCZVODE) -J../zvode -c $< -o $@
../zvode/%.o: ../zvode/%.F
	$(FCZVODE) -J../zvode -c $< -o $@

# ORBIT
FCORBIT = $(FC) $(FFLAGS) $(OMPFLAG) $(RECURSFLAG) -I$(MATHINC) -I../equil
../orbit/%.o: ../orbit/%.f
	$(FCORBIT) -J../orbit -c $< -o $@
../orbit/%.o: ../orbit/%.F
	$(FCORBIT) -J../orbit -c $< -o $@

# VACUUM
FCVACUUM = $(FC) $(FFLAGS) $(ZEROFLAG) $(RECURSFLAG) $(LEGACYFLAG) -I../equil
../vacuum/%.o: ../vacuum/%.f
	$(FCVACUUM) -J../vacuum -c $< -o $@
../vacuum/%.o: ../vacuum/%.F
	$(FCVACUUM) -J../vacuum -c $< -o $@

# PENTRC
FCPENTRC = $(FC) $(FFLAGS) $(OMPFLAG) $(RECURSFLAG) -I$(MATHINC) -I$(NETCDFINC) -I../equil -I../harvest
../pentrc/%.o: ../pentrc/%.f
	$(FCPENTRC) $(LEGACYFLAG) -J../pentrc -c $< -o $@
../pentrc/%.o: ../pentrc/%.F
	$(FCPENTRC) $(LEGACYFLAG) -J../pentrc -c $< -o $@
../pentrc/%.o: ../pentrc/%.f90
	$(FCPENTRC) -J../pentrc -c $< -o $@
../pentrc/%.o: ../pentrc/%.F90
	$(FCPENTRC) -J../pentrc -c $< -o $@

# DCON
FCDCON = $(FC) $(FFLAGS) $(OMPFLAG) $(RECURSFLAG) $(LEGACYFLAG) -I$(MATHINC) -I$(NETCDFINC) -I../equil -I../harvest -I../pentrc -I../vacuum
../dcon/%.o: ../dcon/%.f
	$(FCDCON) -J../dcon -c $< -o $@
../dcon/%.o: ../dcon/%.F
	$(FCDCON) -J../dcon -c $< -o $@

# RDCON
FCRDCON = $(FC) $(FFLAGS) $(OMPFLAG) $(LEGACYFLAG) -I$(MATHINC) -I$(NETCDFINC) -I../equil -I../vacuum
../rdcon/%.o: ../rdcon/%.f
	$(FCRDCON) -J../rdcon -c $< -o $@
../rdcon/%.o: ../rdcon/%.F
	$(FCRDCON) -J../rdcon -c $< -o $@

# STRIDE
FCSTRIDE = $(FC) $(FFLAGS) $(OMPFLAG) $(RECURSFLAG) $(LEGACYFLAG) -I$(MATHINC) -I$(NETCDFINC) -I../equil -I../vacuum
../stride/%.o: ../stride/%.f
	$(FCSTRIDE) -J../stride -c $< -o $@
../stride/%.o: ../stride/%.F
	$(FCSTRIDE) -J../stride -c $< -o $@

# MATCH
FCMATCH = $(FC) $(FFLAGS) $(OMPFLAG) $(RECURSFLAG) -I../equil
../match/%.o: ../match/%.f
	$(FCMATCH) -J../match -c $< -o $@
../match/%.o: ../match/%.F
	$(FCMATCH) -J../match -c $< -o $@

# RMATCH
FCRMATCH = $(FC) $(FFLAGS) $(RECURSFLAG) $(LEGACYFLAG) $(OMPFLAG) -I$(MATHINC) -I../equil
../rmatch/%.o: ../rmatch/%.f
	$(FCRMATCH) -J../rmatch -c $< -o $@
../rmatch/%.o: ../rmatch/%.F
	$(FCRMATCH) -J../rmatch -c $< -o $@

# SUM
FCSUM = $(FC) $(FFLAGS) -I../equil
../sum/%.o: ../sum/%.f
	$(FCSUM) -J../sum -c $< -o $@
../sum/%.o: ../sum/%.F
	$(FCSUM) -J../sum -c $< -o $@

# MULTI
FCMULTI = $(FC) $(FFLAGS) -I../equil
../multi/%.o: ../multi/%.f
	$(FCMULTI) -J../multi -c $< -o $@
../multi/%.o: ../multi/%.F
	$(FCMULTI) -J../multi -c $< -o $@

# COIL
FCCOIL = $(FC) $(FFLAGS) $(OMPFLAG) $(RECURSFLAG) -I$(MATHINC) -I$(NETCDFINC) -I../equil
../coil/%.o: ../coil/%.f
	$(FCCOIL) -J../coil -c $< -o $@
../coil/%.o: ../coil/%.F
	$(FCCOIL) -J../coil -c $< -o $@

# SLAYER
FCSLAYER = $(FC) $(FFLAGS) $(OMPFLAG) -I$(MATHINC) -I$(NETCDFINC) -I../equil
../slayer/%.o: ../slayer/%.f
	$(FCSLAYER) -J../slayer -c $< -o $@
../slayer/%.o: ../slayer/%.F
	$(FCSLAYER) -J../slayer -c $< -o $@

# GPEC
FCGPEC = $(FC) $(FFLAGS) $(OMPFLAG) $(RECURSFLAG) $(LEGACYFLAG) -I$(MATHINC) -I$(NETCDFINC) -I../equil -I../vacuum -I../harvest -I../pentrc -I../coil -I../slayer
../gpec/%.o: ../gpec/%.f
	$(FCGPEC) -J../gpec -c $< -o $@
../gpec/%.o: ../gpec/%.F
	$(FCGPEC) -J../gpec -c $< -o $@

# XDRAW
CCXDRAW = $(CC) $(X11_INC_DIR) $(CFLAGS) -DUNIX -O2
../xdraw/%.o: ../xdraw/%.c
	$(CCXDRAW) -c $< -o $@
../xdraw/%.o: ../xdraw/%.C
	$(CCXDRAW) -c $< -o $@