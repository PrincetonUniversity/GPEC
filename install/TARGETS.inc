#------------------------------------------------------------------------------
# TARGETS.inc
#
# Purpose:
#   Defines how to build libraries and executables from object files for each module.
#
# What it does:
#   - Specifies build commands for each library/executable (e.g. gpec, dcon, pentrc)
#   - Handles versioning for certain modules via git tags
#   - Sets up output directories (`../lib`, `../bin`)
#   - Lists linker flags and dependencies for each target
#
# How to use:
#   - To build a module, run `make <module>` (e.g. `make gpec`)
#   - To add a new target, copy a similar block and update for your module
#   - Most users do not need to modify this file directly
#------------------------------------------------------------------------------

# Setting installation directories.
# You can make with the following command to override the default installation
# directories. This is useful as it will automatically install dependencies in
# a custom location for later re-use.
# DEPSINSTALLDIR = /custom/path/to/libraries make
DEPSINSTALLDIR ?= ../deps
LIBDIR = ../deps/lib
BINDIR = ../bin

# Nice formatting for success messages
SUCCESS_COLOR = \033[1;92m
SUCCESS_RESET = \033[0m
define print_success
	@printf '%bSUCCESS:%b %s\n' '$(SUCCESS_COLOR)' '$(SUCCESS_RESET)' '$(1)'
endef

#------------------------------------------------------------------------------
# Submodule Dependencies
#------------------------------------------------------------------------------

# Harvest library
harvest: v
	cd ../; git submodule update --init deps/src/harvest
	cd ../deps/src/harvest; $(MAKE) -j1 FC='$(FC)' CC='$(CC)' GACODE_ROOT= -f Makefile libharvest.a
	mkdir -p $(LIBDIR)
	cp -f ../deps/src/harvest/libharvest.a $(LIBDIR)
	$(call print_success,harvest -> $(LIBDIR)/libharvest.a)

netcdf: v netcdf-c netcdf-fortran hdf5
	$(call print_success,All netcdf components compiled successfully)

netcdf-c: v hdf5
	cd ../; git submodule update --init deps/src/netcdf-c;
	if [ ! -f ../deps/src/netcdf-c/Makefile ]; then \
		cd ../deps/src/netcdf-c; \
			CPPFLAGS="-I$(realpath $(DEPSINSTALLDIR))/include" \
			LDFLAGS="-L$(realpath $(DEPSINSTALLDIR))/lib" \
			./configure --prefix=$(realpath $(DEPSINSTALLDIR)); \
	fi
	cd ../deps/src/netcdf-c; \
	$(MAKE); \
	$(MAKE) install
	$(call print_success,netcdf-c installed to $(DEPSINSTALLDIR))

netcdf-fortran: v netcdf-c
	cd ../; git submodule update --init deps/src/netcdf-fortran;
	if [ ! -f ../deps/src/netcdf-fortran/Makefile ]; then \
		cd ../deps/src/netcdf-fortran; \
			CPPFLAGS="-I$(realpath $(DEPSINSTALLDIR))/include" \
			LDFLAGS="-L$(realpath $(DEPSINSTALLDIR))/lib" \
			./configure --prefix=$(realpath $(DEPSINSTALLDIR)); \
	fi
	cd ../deps/src/netcdf-fortran; \
	$(MAKE); \
	$(MAKE) install
	$(call print_success,netcdf-fortran installed to $(DEPSINSTALLDIR))

hdf5: v
	cd ../; git submodule update --init deps/src/hdf5;
	if [ ! -f ../deps/src/hdf5/Makefile ]; then \
		cd ../deps/src/hdf5; \
		./configure --prefix=$(realpath $(DEPSINSTALLDIR)); \
	fi
	cd ../deps/src/hdf5; \
	$(MAKE); \
	$(MAKE) install;
	$(call print_success,hdf5 installed to $(DEPSINSTALLDIR))

OpenBLAS: v
	cd ../; git submodule update --init deps/src/OpenBLAS;
	cd ../deps/src/OpenBLAS; \
		PREFIX=$(realpath $(DEPSINSTALLDIR)) $(MAKE); \
		PREFIX=$(realpath $(DEPSINSTALLDIR)) $(MAKE) install;
	$(call print_success,OpenBLAS installed to $(DEPSINSTALLDIR))

# All dependencies
alldeps: v harvest netcdf OpenBLAS
	$(call print_success,All dependencies compiled successfully)

# Needed dependencies based on DEFAULTS.inc
neededdeps: v $(NEEDED_DEPS)
	$(call print_success,Needed dependencies compiled successfully)


#------------------------------------------------------------------------------
# Build Targets
#------------------------------------------------------------------------------


# Equilibrium library
equil: neededdeps harvest $(EQUIL_OBJS)
	mkdir -p $(LIBDIR)
	ar -r $(LIBDIR)/libequil.a $(EQUIL_OBJS)
	$(call print_success,libequil.a -> $(LIBDIR))

# LSODE integrator
lsode: neededdeps $(LSODE_OBJS)
	mkdir -p $(LIBDIR)
	ar -r $(LIBDIR)/liblsode.a $(LSODE_OBJS)
	$(call print_success,liblsode.a -> $(LIBDIR))

# ZLANGE integrator
zlange: neededdeps $(ZLANGE_OBJS)
	mkdir -p $(LIBDIR)
	ar -r $(LIBDIR)/libzlange.a $(ZLANGE_OBJS)
	$(call print_success,libzlange.a -> $(LIBDIR))

# ZVODE integrator
zvode: neededdeps $(ZVODE_OBJS)
	mkdir -p $(LIBDIR)
	ar -r $(LIBDIR)/libzvode.a $(ZVODE_OBJS)
	$(call print_success,libzvode.a -> $(LIBDIR))


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
orbit: neededdeps $(ORBIT_OBJS)
	rm -f ../orbit/orbit
	mkdir -p $(BINDIR)
	$(FCORBIT) -o ../orbit/orbit $(ORBIT_OBJS) $(ORBIT_LIBDIRS) $(ORBIT_LIBS) $(LDFLAGS)
	cp ../orbit/orbit $(BINDIR)
	$(call print_success,orbit -> $(BINDIR)/orbit)

# Vacuum code
vacuum: v neededdeps $(VACUUM_OBJS) equil lsode
	mkdir -p $(LIBDIR)
	ar -r $(LIBDIR)/libvac.a $(VACUUM_OBJS)
	$(call print_success,libvac.a -> $(LIBDIR))

# PENTRC code
../pentrc/version.inc:
	@ver="`git describe --tags`"; \
	if [ ! -f $@ ] || ! grep -F "'$$ver'" $@ >/dev/null 2>&1; then \
		echo ">>> Updating version file"; \
		echo "      CHARACTER(len=*), PARAMETER :: version ='$$ver'" | tee $@; \
	else \
		echo ">>> Git version unchanged!"; \
	fi
	$(call print_success,version file updated: $@)
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
pentrc: v neededdeps $(PENTRC_OBJS) lsode equil harvest
	mkdir -p $(BINDIR)
	rm -f ../pentrc/pentrc
	$(FCPENTRC) -o ../pentrc/pentrc $(PENTRC_OBJS) $(PENTRC_LIBDIRS) $(PENTRC_LIBS) $(LDFLAGS)
	cp ../pentrc/pentrc $(BINDIR)
	mkdir -p $(LIBDIR)
	ar -r $(LIBDIR)/libpentrc.a $(PENTRC_OBJS)
	$(call print_success,pentrc -> $(BINDIR)/pentrc and libpentrc.a -> $(LIBDIR))

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
dcon: v neededdeps $(DCON_OBJS) equil lsode vacuum pentrc harvest
	mkdir -p $(BINDIR)
	rm -f ../dcon/dcon
	$(FCDCON) -o ../dcon/dcon $(DCON_OBJS) $(DCON_LIBDIRS) $(DCON_LIBS) $(LDFLAGS)
	cp ../dcon/dcon $(BINDIR)
	$(call print_success,dcon -> $(BINDIR)/dcon)

# RDCON code
../rdcon/version.inc:
	@ver="`git describe --tags`"; \
	if [ ! -f $@ ] || ! grep -F "'$$ver'" $@ >/dev/null 2>&1; then \
		echo ">>> Updating version file"; \
		echo "      CHARACTER(len=*), PARAMETER :: version ='$$ver'" | tee $@; \
	else \
		echo ">>> Git version unchanged!"; \
	fi
	$(call print_success,version file updated: $@)
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
rdcon: v neededdeps $(RDCON_OBJS) equil lsode vacuum
	mkdir -p $(BINDIR)
	rm -f ../rdcon/rdcon
	$(FCRDCON) -o ../rdcon/rdcon $(RDCON_OBJS) $(RDCON_LIBDIRS) $(RDCON_LIBS) $(LDFLAGS)
	cp ../rdcon/rdcon $(BINDIR)
	$(call print_success,rdcon -> $(BINDIR)/rdcon)

# STRIDE code
../stride/version.inc:
	@ver="`git describe --tags`"; \
	if [ ! -f $@ ] || ! grep -F "'$$ver'" $@ >/dev/null 2>&1; then \
		echo ">>> Updating version file"; \
		echo "      CHARACTER(len=*), PARAMETER :: version ='$$ver'" | tee $@; \
	else \
		echo ">>> Git version unchanged!"; \
	fi
	$(call print_success,version file updated: $@)
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
stride: v neededdeps $(STRIDE_OBJS) zlange zvode lsode equil vacuum
	mkdir -p $(BINDIR)
	rm -f ../stride/stride
	$(FCSTRIDE) -o ../stride/stride $(STRIDE_OBJS) $(STRIDE_LIBDIRS) $(STRIDE_LIBS) $(LDFLAGS)
	cp ../stride/stride $(BINDIR)
	$(call print_success,stride -> $(BINDIR)/stride)

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
match: v neededdeps $(MATCH_OBJS) equil lsode
	mkdir -p $(BINDIR)
	rm -f ../match/match
	$(FCMATCH) -o ../match/match $(MATCH_OBJS) $(MATCH_LIBDIRS) $(MATCH_LIBS) $(LDFLAGS)
	cp ../match/match $(BINDIR)
	$(call print_success,match -> $(BINDIR)/match)

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
rmatch: v neededdeps $(RMATCH_OBJS) equil lsode
	mkdir -p $(BINDIR)
	rm -f ../rmatch/rmatch
	$(FCRMATCH) -o ../rmatch/rmatch $(RMATCH_OBJS) $(RMATCH_LIBDIRS) $(RMATCH_LIBS) $(LDFLAGS)
	cp ../rmatch/rmatch $(BINDIR)
	$(call print_success,rmatch -> $(BINDIR)/rmatch)

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
slayer: v neededdeps $(SLAYER_OBJS) lsode equil
	mkdir -p $(BINDIR)
	rm -f ../slayer/slayer
	$(FCSLAYER) -o ../slayer/slayer $(SLAYER_OBJS) $(SLAYER_LIBDIRS) $(SLAYER_LIBS) $(LDFLAGS)
	ar -r $(LIBDIR)/libslayer.a $(SLAYER_OBJS)
	cp ../slayer/slayer $(BINDIR)
	$(call print_success,slayer -> $(BINDIR)/slayer and libslayer.a -> $(LIBDIR))

# multi
MULTI_LIBS = \
	-lequil
MULTI_LIBDIRS = \
	-L$(LIBDIR)
multi: v neededdeps $(MULTI_OBJS) equil
	mkdir -p $(BINDIR)
	rm -f ../multi/multi
	$(FCMULTI) -o ../multi/multi $(MULTI_OBJS) $(MULTI_LIBDIRS) $(MULTI_LIBS) $(LDFLAGS)
	cp ../multi/multi $(BINDIR)
	$(call print_success,multi -> $(BINDIR)/multi)

# sum
SUM_LIBS = \
	-lequil \
	$(FCLIBS)
SUM_LIBDIRS = \
	-L$(LIBDIR) \
	-L$(FCDIR)
sum: v neededdeps $(SUM_OBJS) equil
	mkdir -p $(BINDIR)
	rm -f ../sum/sum
	$(FCSUM) -o ../sum/sum $(SUM_OBJS) $(SUM_LIBDIRS) $(SUM_LIBS) $(LDFLAGS)
	cp ../sum/sum $(BINDIR)
	$(call print_success,sum -> $(BINDIR)/sum)

# Coil library
../coil/version.inc:
	@ver="`git describe --tags`"; \
	if [ ! -f $@ ] || ! grep -F "'$$ver'" $@ >/dev/null 2>&1; then \
		echo ">>> Updating version file"; \
		echo "      CHARACTER(len=*), PARAMETER :: version ='$$ver'" | tee $@; \
	else \
		echo ">>> Git version unchanged!"; \
	fi
	$(call print_success,version file updated: $@)
coil: v neededdeps $(COIL_OBJS)
	mkdir -p $(LIBDIR)
	ar -r $(LIBDIR)/libcoil.a $(COIL_OBJS)
	$(call print_success,libcoil.a -> $(LIBDIR))

# GPEC
# decide if we need to update the version
# note gpec keeps track of its own version, which can differ from the synced dcon, pentrc, and stride ones
# this enables faster compiling during development
../gpec/version.inc:
	@ver="`git describe --tags`"; \
	if [ ! -f $@ ] || ! grep -F "'$$ver'" $@ >/dev/null 2>&1; then \
		echo ">>> Updating version file"; \
		echo "      CHARACTER(len=*), PARAMETER :: version ='$$ver'" | tee $@; \
	else \
		echo ">>> Git version unchanged!"; \
	fi
	$(call print_success,version file updated: $@)
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
gpec: v neededdeps $(GPEC_OBJS) pentrc lsode equil vacuum harvest slayer coil
	mkdir -p $(BINDIR)
	rm -f ../gpec/gpec
	$(FCGPEC) -o ../gpec/gpec $(GPEC_OBJS) $(GPEC_LIBDIRS) $(GPEC_LIBS) $(LDFLAGS)
	cp ../gpec/gpec $(BINDIR)
	$(call print_success,gpec -> $(BINDIR)/gpec)

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
	$(call print_success,Xdraw tools built successfully)
xdraw_main: v $(XDRAW_OBJS1)
	mkdir -p $(BINDIR)
	rm -f ../xdraw/xdraw
	$(CCXDRAW) -o ../xdraw/xdraw $(XDRAW_OBJS1) $(XDRAW_LIBDIRS) $(XDRAW_LIBS) $(LDFLAGS)
	cp ../xdraw/xdraw $(BINDIR)
	$(call print_success,xdraw -> $(BINDIR)/xdraw)
pspack: v $(XDRAW_OBJS2)
	mkdir -p $(BINDIR)
	rm -f ../xdraw/pspack
	$(CCXDRAW) -o ../xdraw/pspack $(XDRAW_OBJS2) $(XDRAW_LIBDIRS) $(XDRAW_LIBS) $(LDFLAGS)
	cp ../xdraw/pspack $(BINDIR)
	$(call print_success,pspack -> $(BINDIR)/pspack)
# Misc
force: v

mkbin: pentrc dcon rdcon stride match rmatch slayer multi sum gpec
	rm -f $(BINDIR)/*
	mkdir -p $(BINDIR)
	cp ../pentrc/pentrc $(BINDIR)
	cp ../dcon/dcon $(BINDIR)
	cp ../rdcon/rdcon $(BINDIR)
	cp ../stride/stride $(BINDIR)
	cp ../match/match $(BINDIR)
	cp ../rmatch/rmatch $(BINDIR)
	cp ../slayer/slayer $(BINDIR)
	cp ../multi/multi $(BINDIR)
	cp ../sum/sum $(BINDIR)
	cp ../gpec/gpec $(BINDIR)
	$(call print_success,All binaries copied to $(BINDIR))
