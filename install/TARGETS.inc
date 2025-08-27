# Lower level targets

LIBDIR = ../lib

# Harvest library
harvest:
	cd ../; git submodule init; git submodule update
	cd ../harvest; $(MAKE) -j1 FC='$(FC)' CC='$(CC)' GACODE_ROOT= -f Makefile libharvest.a
	mkdir -p $(LIBDIR)
	cp -f ../harvest/libharvest.a $(LIBDIR)

# Equilibrium library
equil: harvest $(EQUIL_OBJS)
	mkdir -p $(LIBDIR)
	ar -r $(LIBDIR)/libequil.a $^

# LSODE integrator
lsode: $(LSODE_OBJS)
	mkdir -p $(LIBDIR)
	ar -r $(LIBDIR)/liblsode.a $^

# ZLANGE integrator
zlange: $(ZLANGE_OBJS)
	mkdir -p $(LIBDIR)
	ar -r $(LIBDIR)/libzlange.a $^

# ZVODE integrator
zvode: $(ZVODE_OBJS)
	mkdir -p $(LIBDIR)
	ar -r $(LIBDIR)/libzvode.a $^


# Orbit
ORBIT_LIBS = \
	-llsode \
	-lequil \
	$(MATHLIBS) \
	$(FCLIBS)
ORBIT_LIBDIRS = \
	-L$(LIBDIR) \
	-L$(MATHDIR) \
	-L$(FCDIR)
orbit: $(ORBIT_OBJS)
	mkdir -p ../bin
	$(FCORBIT) -o ../bin/orbit $^ $(ORBIT_LIBDIRS) $(ORBIT_LIBS) $(LDFLAGS)

# Vacuum code
vacuum: $(VACUUM_OBJS)
	mkdir -p $(LIBDIR)
	ar -r $(LIBDIR)/libvac.a $^

# PENTRC code
ifeq ($(wildcard ../pentrc/version.inc),)
pentrcversion: force
	@echo ">>> Creating version file"
	@echo "      CHARACTER(len=*), PARAMETER :: version ='"`git describe --tags`"'" | tee ../pentrc/version.inc
else ifeq ($(shell grep -F "'"`git describe --tags`"'" ../pentrc/version.inc),)
pentrcversion: force
	@echo ">>> Updating version file"
	@echo "      CHARACTER(len=*), PARAMETER :: version ='"`git describe --tags`"'" | tee ../pentrc/version.inc
else
pentrcversion:
	@echo ">>> Git version unchanged!"
endif
PENTRC_LIBS = \
	-llsode \
	-lequil \
	-lharvest \
	$(MATHLIBS) \
	$(NETCDF_LIBS) \
	$(FCLIBS) \
	$(NETCDF_EXTRA_LIBS)
PENTRC_LIBDIRS = \
	-L$(LIBDIR) \
	-L$(MATHDIR) \
	-L$(NETCDFDIR) \
	-L$(FCDIR)
pentrc: $(PENTRC_OBJS)
	mkdir -p ../bin
	$(FCPENTRC) -o ../bin/pentrc $^ $(PENTRC_LIBDIRS) $(PENTRC_LIBS) $(LDFLAGS)

# DCON code
DCON_LIBS = \
	-lpentrc \
	-llsode \
	-lequil \
	-lvac \
	-lharvest \
	$(MATHLIBS) \
	$(NETCDF_LIBS) \
	$(FCLIBS) \
	$(NETCDF_EXTRA_LIBS)
DCON_LIBDIRS = \
	-L$(LIBDIR) \
	-L$(MATHDIR) \
	-L$(NETCDFDIR) \
	-L$(FCDIR)
dcon: $(DCON_OBJS)
	mkdir -p ../bin
	$(FCDCON) -o ../bin/dcon $^ $(DCON_LIBDIRS) $(DCON_LIBS) $(LDFLAGS)

# RDCON code
ifeq ($(wildcard ../rdcon/version.inc),)
rdconversion: force
	@echo ">>> Creating version file"
	@echo "      CHARACTER(len=*), PARAMETER :: version ='"`git describe --tags`"'" | tee ../rdcon/version.inc
else ifeq ($(shell grep -F "'"`git describe --tags`"'" ../rdcon/version.inc),)
rdconversion: force
	@echo ">>> Updating version file"
	@echo "      CHARACTER(len=*), PARAMETER :: version ='"`git describe --tags`"'" | tee ../rdcon/version.inc
else
rdconversion:
	@echo ">>> Git version unchanged!"
endif
RDCON_LIBS = \
	-lpthread \
	-llsode \
	-lequil \
	-lvac \
	$(MATHLIBS) \
	$(NETCDF_LIBS) \
	$(FCLIBS) \
	$(NETCDF_EXTRA_LIBS)
RDCON_LIBDIRS = \
	-L$(LIBDIR) \
	-L$(MATHDIR) \
	-L$(NETCDFDIR) \
	-L$(FCDIR)
rdcon: $(RDCON_OBJS)
	mkdir -p ../bin
	$(FCRDCON) -o ../bin/rdcon $^ $(RDCON_LIBDIRS) $(RDCON_LIBS) $(LDFLAGS)

# STRIDE code
STRIDE_LIBS = \
	-lzvode \
	-lzlange \
	-llsode \
	-lequil \
	-lvac \
	$(MATHLIBS) \
	$(NETCDF_LIBS) \
	$(NETCDF_EXTRA_LIBS)
STRIDE_LIBDIRS = \
	-L$(LIBDIR) \
	-L$(MATHDIR) \
	-L$(NETCDFDIR)
stride: $(STRIDE_OBJS)
	mkdir -p ../bin
	$(FCSTRIDE) -o ../bin/stride $^ $(STRIDE_LIBDIRS) $(STRIDE_LIBS) $(LDFLAGS)

# match
MATCH_LIBS = \
	-lequil \
	-llsode \
	$(MATHLIBS) \
	$(FCLIBS)
MATCH_LIBDIRS = \
	-L$(LIBDIR) \
	-L$(MATHDIR) \
	-L$(FCDIR)
match: $(MATCH_OBJS)
	mkdir -p ../bin
	$(FCMATCH) -o ../bin/match $^ $(MATCH_LIBDIRS) $(MATCH_LIBS) $(LDFLAGS)

# rmatch
RMATCH_LIBS = \
	-lequil \
	-llsode \
	$(MATHLIBS) \
	$(NETCDF_LIBS) \
	$(NETCDF_EXTRA_LIBS) \
	$(FCLIBS)
RMATCH_LIBDIRS = \
	-L$(LIBDIR) \
	-L$(MATHDIR) \
	-L$(NETCDFDIR) \
	-L$(FCDIR)
rmatch: $(RMATCH_OBJS)
	mkdir -p ../bin
	$(FCRMATCH) -o ../bin/rmatch $^ $(RMATCH_LIBDIRS) $(RMATCH_LIBS) $(LDFLAGS)

# SLAYER
SLAYER_LIBS = \
	-llsode \
	-lequil \
	$(MATHLIBS) \
	$(NETCDF_LIBS) \
	$(NETCDF_EXTRA_LIBS)
SLAYER_LIBDIRS = \
	-L$(LIBDIR) \
	-L$(MATHDIR) \
	-L$(NETCDFDIR)
slayer:
	mkdir -p ../bin
	$(FCSLAYER) -o ../bin/slayer $^ $(SLAYER_LIBDIRS) $(SLAYER_LIBS) $(LDFLAGS)

# multi
MULTI_LIBS = \
	-lequil
MULTI_LIBDIRS = \
	-L$(LIBDIR)
multi: $(MULTI_OBJS)
	mkdir -p ../bin
	$(FCMULTI) -o ../bin/multi $^ $(MULTI_LIBDIRS) $(MULTI_LIBS) $(LDFLAGS)

# sum
SUM_LIBS = \
	-lequil \
	$(FCLIBS)
SUM_LIBDIRS = \
	-L$(LIBDIR) \
	-L$(FCDIR)
sum: $(SUM_OBJS)
	mkdir -p ../bin
	$(FCSUM) -o ../bin/sum $^ $(SUM_LIBDIRS) $(SUM_LIBS) $(LDFLAGS)

# Coil library
ifeq ($(wildcard ../coil/version.inc),)
coilversion: force
	@echo ">>> Creating version file"
	@echo "      CHARACTER(len=*), PARAMETER :: version ='"`git describe --tags`"'" | tee ../coil/version.inc
else ifeq ($(shell grep -F "'"`git describe --tags`"'" ../coil/version.inc),)
coilversion: force
	@echo ">>> Updating version file"
	@echo "      CHARACTER(len=*), PARAMETER :: version ='"`git describe --tags`"'" | tee ../coil/version.inc
else
coilversion:
	@echo ">>> Git version unchanged!"
endif
coil: $(COIL_OBJS)
	mkdir -p $(LIBDIR)
	ar -r $(LIBDIR)/libcoil.a $^

# GPEC
# decide if we need to update the version
# note gpec keeps track of its own version, which can differ from the synced dcon, pentrc, and stride ones
# this enables faster compiling during development
ifeq ($(wildcard ../gpec/version.inc),)
gpecversion: force
	@echo ">>> Creating version file"
	@echo "      CHARACTER(len=*), PARAMETER :: version ='"`git describe --tags`"'" | tee ../gpec/version.inc
else ifeq ($(shell grep -F "'"`git describe --tags`"'" ../gpec/version.inc),)
gpecversion: force
	@echo ">>> Updating version file"
	@echo "      CHARACTER(len=*), PARAMETER :: version ='"`git describe --tags`"'" | tee ../gpec/version.inc
else
gpecversion:
	@echo ">>> Git version unchanged!"
endif
GPEC_LIBS = \
	-lpentrc \
	-llsode \
	-lequil \
	-lvac \
	-lharvest \
	-lslayer \
	-lcoil \
	$(MATHLIBS) \
	$(NETCDF_LIBS) \
	$(NETCDF_EXTRA_LIBS)
GPEC_LIBDIRS = \
	-L$(LIBDIR) \
	-L$(MATHDIR) \
	-L$(NETCDFDIR)
gpec: $(GPEC_OBJS)
	mkdir -p ../bin
	$(FCGPEC) -o ../bin/gpec $^ $(GPEC_LIBDIRS) $(GPEC_LIBS) $(LDFLAGS)

# Misc
force: