# Lower level targets

LIBDIR = ../deps/lib
BINDIR = ../bin

# Harvest library
harvest: v
	cd ../; git submodule update --init ../deps/src/harvest
	cd ../deps/src/harvest; $(MAKE) -j1 FC='$(FC)' CC='$(CC)' GACODE_ROOT= -f Makefile libharvest.a
	mkdir -p $(LIBDIR)
	cp -f ../deps/src/harvest/libharvest.a $(LIBDIR)

# Equilibrium library
equil: v harvest $(EQUIL_OBJS)
	mkdir -p $(LIBDIR)
	ar -r $(LIBDIR)/libequil.a $(EQUIL_OBJS)

# LSODE integrator
lsode: v $(LSODE_OBJS)
	mkdir -p $(LIBDIR)
	ar -r $(LIBDIR)/liblsode.a $(LSODE_OBJS)

# ZLANGE integrator
zlange: v $(ZLANGE_OBJS)
	mkdir -p $(LIBDIR)
	ar -r $(LIBDIR)/libzlange.a $(ZLANGE_OBJS)

# ZVODE integrator
zvode: v $(ZVODE_OBJS)
	mkdir -p $(LIBDIR)
	ar -r $(LIBDIR)/libzvode.a $(ZVODE_OBJS)


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
orbit: v $(ORBIT_OBJS)
	mkdir -p $(BINDIR)
	$(FCORBIT) -o $(BINDIR)/orbit $(ORBIT_OBJS) $(ORBIT_LIBDIRS) $(ORBIT_LIBS) $(LDFLAGS)

# Vacuum code
vacuum: v $(VACUUM_OBJS) equil lsode
	mkdir -p $(LIBDIR)
	ar -r $(LIBDIR)/libvac.a $(VACUUM_OBJS)

# PENTRC code
ifeq ($(wildcard ../pentrc/version.inc),)
pentrcversion: v force
	@echo ">>> Creating version file"
	@echo "      CHARACTER(len=*), PARAMETER :: version ='"`git describe --tags`"'" | tee ../pentrc/version.inc
else ifeq ($(shell grep -F "'"`git describe --tags`"'" ../pentrc/version.inc),)
pentrcversion: v force
	@echo ">>> Updating version file"
	@echo "      CHARACTER(len=*), PARAMETER :: version ='"`git describe --tags`"'" | tee ../pentrc/version.inc
else
pentrcversion: v
	@echo ">>> Git version unchanged!"
endif
PENTRC_LIBS = \
	-llsode \
	-lequil \
	-lharvest \
	$(MATHLIBS) \
	$(NETCDFLIBS) \
	$(FCLIBS) \
	$(NETCDF_EXTRA_LIBS)
PENTRC_LIBDIRS = \
	-L$(LIBDIR) \
	-L$(MATHDIR) \
	-L$(NETCDFDIR) \
	-L$(FCDIR)
pentrc: v $(PENTRC_OBJS) lsode equil harvest
	mkdir -p $(BINDIR)
	$(FCPENTRC) -o $(BINDIR)/pentrc $(PENTRC_OBJS) $(PENTRC_LIBDIRS) $(PENTRC_LIBS) $(LDFLAGS)
	mkdir -p $(LIBDIR)
	ar -r $(LIBDIR)/libpentrc.a $(PENTRC_OBJS)

# DCON code
DCON_LIBS = \
	-lpentrc \
	-llsode \
	-lequil \
	-lvac \
	-lharvest \
	$(MATHLIBS) \
	$(NETCDFLIBS) \
	$(FCLIBS) \
	$(NETCDF_EXTRA_LIBS)
DCON_LIBDIRS = \
	-L$(LIBDIR) \
	-L$(MATHDIR) \
	-L$(NETCDFDIR) \
	-L$(FCDIR)
dcon: v $(DCON_OBJS) equil lsode vacuum pentrc harvest
	mkdir -p $(BINDIR)
	$(FCDCON) -o $(BINDIR)/dcon $(DCON_OBJS) $(DCON_LIBDIRS) $(DCON_LIBS) $(LDFLAGS)

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
rdconversion: force
	@echo ">>> Git version unchanged!"
endif
RDCON_LIBS = \
	-lpthread \
	-llsode \
	-lequil \
	-lvac \
	$(MATHLIBS) \
	$(NETCDFLIBS) \
	$(FCLIBS) \
	$(NETCDF_EXTRA_LIBS)
RDCON_LIBDIRS = \
	-L$(LIBDIR) \
	-L$(MATHDIR) \
	-L$(NETCDFDIR) \
	-L$(FCDIR)
rdcon: v $(RDCON_OBJS) equil lsode vacuum
	mkdir -p $(BINDIR)
	$(FCRDCON) -o $(BINDIR)/rdcon $(RDCON_OBJS) $(RDCON_LIBDIRS) $(RDCON_LIBS) $(LDFLAGS)

# STRIDE code
ifeq ($(wildcard ../stride/version.inc),)
strideversion: force
	@echo ">>> Creating version file"
	@echo "      CHARACTER(len=*), PARAMETER :: version ='"`git describe --tags`"'" | tee ../stride/version.inc
else ifeq ($(shell grep -F "'"`git describe --tags`"'" ../stride/version.inc),)
strideversion: force
	@echo ">>> Updating version file"
	@echo "      CHARACTER(len=*), PARAMETER :: version ='"`git describe --tags`"'" | tee ../stride/version.inc
else
strideversion: force
	@echo ">>> Git version unchanged!"
endif
STRIDE_LIBS = \
	-lzvode \
	-lzlange \
	-llsode \
	-lequil \
	-lvac \
	$(MATHLIBS) \
	$(NETCDFLIBS) \
	$(NETCDF_EXTRA_LIBS)
STRIDE_LIBDIRS = \
	-L$(LIBDIR) \
	-L$(MATHDIR) \
	-L$(NETCDFDIR)
stride: v $(STRIDE_OBJS) zlange zvode lsode equil vacuum
	mkdir -p $(BINDIR)
	$(FCSTRIDE) -o $(BINDIR)/stride $(STRIDE_OBJS) $(STRIDE_LIBDIRS) $(STRIDE_LIBS) $(LDFLAGS)

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
match: v $(MATCH_OBJS) equil lsode
	mkdir -p $(BINDIR)
	$(FCMATCH) -o $(BINDIR)/match $(MATCH_OBJS) $(MATCH_LIBDIRS) $(MATCH_LIBS) $(LDFLAGS)

# rmatch
RMATCH_LIBS = \
	-lequil \
	-llsode \
	$(MATHLIBS) \
	$(NETCDFLIBS) \
	$(NETCDF_EXTRA_LIBS) \
	$(FCLIBS)
RMATCH_LIBDIRS = \
	-L$(LIBDIR) \
	-L$(MATHDIR) \
	-L$(NETCDFDIR) \
	-L$(FCDIR)
rmatch: v $(RMATCH_OBJS) equil lsode
	mkdir -p $(BINDIR)
	$(FCRMATCH) -o $(BINDIR)/rmatch $(RMATCH_OBJS) $(RMATCH_LIBDIRS) $(RMATCH_LIBS) $(LDFLAGS)

# SLAYER
SLAYER_LIBS = \
	-llsode \
	-lequil \
	$(MATHLIBS) \
	$(NETCDFLIBS) \
	$(NETCDF_EXTRA_LIBS)
SLAYER_LIBDIRS = \
	-L$(LIBDIR) \
	-L$(MATHDIR) \
	-L$(NETCDFDIR)
slayer: v $(SLAYER_OBJS) lsode equil
	mkdir -p $(BINDIR)
	$(FCSLAYER) -o $(BINDIR)/slayer $(SLAYER_OBJS) $(SLAYER_LIBDIRS) $(SLAYER_LIBS) $(LDFLAGS)
	ar -r $(LIBDIR)/libslayer.a $(SLAYER_OBJS)

# multi
MULTI_LIBS = \
	-lequil
MULTI_LIBDIRS = \
	-L$(LIBDIR)
multi: v $(MULTI_OBJS) equil
	mkdir -p $(BINDIR)
	$(FCMULTI) -o $(BINDIR)/multi $(MULTI_OBJS) $(MULTI_LIBDIRS) $(MULTI_LIBS) $(LDFLAGS)

# sum
SUM_LIBS = \
	-lequil \
	$(FCLIBS)
SUM_LIBDIRS = \
	-L$(LIBDIR) \
	-L$(FCDIR)
sum: v $(SUM_OBJS) equil
	mkdir -p $(BINDIR)
	$(FCSUM) -o $(BINDIR)/sum $(SUM_OBJS) $(SUM_LIBDIRS) $(SUM_LIBS) $(LDFLAGS)

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
coilversion: force
	@echo ">>> Git version unchanged!"
endif
coil: v $(COIL_OBJS)
	mkdir -p $(LIBDIR)
	ar -r $(LIBDIR)/libcoil.a $(COIL_OBJS)

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
gpecversion: force
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
	$(NETCDFLIBS) \
	$(NETCDF_EXTRA_LIBS)
GPEC_LIBDIRS = \
	-L$(LIBDIR) \
	-L$(MATHDIR) \
	-L$(NETCDFDIR)
gpec: v $(GPEC_OBJS) pentrc lsode equil vacuum harvest slayer coil
	mkdir -p $(BINDIR)
	$(FCGPEC) -o $(BINDIR)/gpec $(GPEC_OBJS) $(GPEC_LIBDIRS) $(GPEC_LIBS) $(LDFLAGS)

# Xdraw
XDRAW_LIBS = \
	-lm \
	-lX11
XDRAW_LIBDIRS = \
	-L/usr/X11R6/lib64 \
	-L/opt/local/include \
	-L/opt/homebrew/include \
	$(X11_LIB_DIR)
xdraw: v xdraw_main pspack
xdraw_main: v $(XDRAW_OBJS1)
	mkdir -p $(BINDIR)
	$(CCXDRAW) -o $(BINDIR)/xdraw $(XDRAW_OBJS1) $(XDRAW_LIBDIRS) $(XDRAW_LIBS) $(LDFLAGS)
pspack: v $(XDRAW_OBJS2)
	mkdir -p $(BINDIR)
	$(CCXDRAW) -o $(BINDIR)/pspack $(XDRAW_OBJS2) $(XDRAW_LIBDIRS) $(XDRAW_LIBS) $(LDFLAGS)

# Misc
force: v