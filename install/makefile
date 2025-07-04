
include DEFAULTS.inc

# extend these variables to sub-makes (?= sets a default if user doesn't specify options)
FFLAGS ?= -g
CFLAGS ?= -g

export FFLAGS
export CFLAGS
export OMPFLAG
export LDFLAGS

# compile package

.IGNORE:

all: v harvest equil lsode zlange zvode orbit vacuum pentrc dcon match rdcon rmatch \
	multi sum xdraw slayer coil gpec stride mkbin

# fortran only
allf: v harvest equil lsode zlange zvode orbit vacuum pentrc dcon match rdcon rmatch \
	multi sum slayer coil gpec stride mkbin


# executables:

harvest:
	cd ../; git submodule init; git submodule update
	cd ../harvest; make FC='$(FC)' CC='$(CC)' GACODE_ROOT= -f Makefile libharvest.a
	mkdir -p ../lib
	cp -f ../harvest/libharvest.a ../lib/

xdraw:
	cd ../xdraw;  make

lsode:
	cd ../lsode;  make liblsode

zlange:
	cd ../zlange;  make libzlange

zvode:
	cd ../zvode;  make libzvode

equil:
	cd ../equil;  make libequil

vacuum:
	cd ../vacuum; make libvac

orbit:
	cd ../orbit;  make orbit

dcon:
	cd ../dcon;   make dcon

match:
	cd ../match;  make match

rdcon:
	cd ../rdcon;   make rdcon

rmatch:
	cd ../rmatch;  make rmatch

multi:
	cd ../multi;  make multi

sum:
	cd ../sum;    make sum

slayer:
	cd ../slayer; make slayer

coil:
	cd ../coil;   make libcoil

gpec:
	cd ../gpec;   make gpec

pentrc:
	cd ../pentrc; make pentrc

stride:
	cd ../stride; make stride

mkbin:
	mkdir -p ../bin
	cp -f ../orbit/orbit ../bin
	cp -f ../dcon/dcon ../bin
	cp -f ../match/match ../bin
	cp -f ../rdcon/rdcon ../bin
	cp -f ../rmatch/rmatch ../bin
	cp -f ../xdraw/xdraw ../bin
	cp -f ../multi/multi ../bin
	cp -f ../sum/sum ../bin
	cp -f ../slayer/slayer ../bin
	cp -f ../gpec/gpec ../bin
	cp -f ../pentrc/pentrc ../bin
	cp -f ../stride/stride ../bin

# print verbose info

v:
	@echo "------------------------------------------------------------"
	@echo "Defaults tailored for $(PLATFORM)"
	@echo "Compiling with FC = $(FC)"
	@echo "   using math libraries $(MATHLIBS)"
	@echo "   from $(MATHDIR)"
	@echo "   including $(MATHINC)"
	@echo "   and netcdf libraries $(NETCDFLIBS) $(NETCDF_C_LIBS)"
	@echo "   from $(NETCDFDIR) $(NETCDF_C_DIR)"
	@echo "   including $(NETCDFINC)"
	@echo "   and compiler libraries $(FCLIBS)"
	@echo "   from $(FCDIR)"
	@echo "   including $(FCINC)"
	@echo "   and using OMPFLAG $(OMPFLAG)"
	@echo
	@echo "Set environment FC,LAPACKHOME,NETCDFHOME to change these"
ifeq ("$(PLATFORM)","portal")
	@echo "   example: module load gcc4 lapack netcdf"
endif
	@echo "------------------------------------------------------------"

# cleaning

clean:
	cd ../harvest; make clean
	cd ../zlange;  make clean
	cd ../zvode;  make clean
	cd ../xdraw;  make clean
	cd ../lsode;  make clean
	cd ../equil;  make clean
	cd ../orbit;  make clean
	cd ../vacuum; make clean
	cd ../dcon;   make clean
	cd ../match;  make clean
	cd ../rdcon;  make clean
	cd ../rmatch; make clean
	cd ../multi;  make clean
	cd ../sum;    make clean
	cd ../slayer; make clean
	cd ../coil;   make clean
	cd ../gpec;   make clean
	cd ../pentrc; make clean
	cd ../stride; make clean

clear:
	cd ../docs/examples/run_ideal_example; rm -f *.out *.bin *.nc
	cd ../docs/examples/run_kinetic_example; rm -f *.out *.bin *.nc
	cd ../docs/examples/a10_ideal_example; rm -f *.out *.bin *.nc
	cd ../docs/examples/a10_kinetic_example; rm -f *.out *.bin *.nc
	cd ../docs/examples/a5_tearing_example; rm -f *.out *.bin *.nc
	cd ../docs/examples/DIIID_ideal_example; rm -f *.out *.bin *.nc
	cd ../docs/examples/DIIID_resistive_example; rm -f *.out *.bin *.nc *.dat
	cd ../docs/examples/DIIID_kinetic_example; rm -f *.out *.bin *.nc
	cd ../docs/examples/solovev_ideal_example; rm -f *.out *.bin *.nc
	cd ../docs/examples/solovev_kinetic_example; rm -f *.out *.bin *.nc
	cd ../docs/examples/solovev_resistive_example; rm -f *.out *.bin *.nc *.dat

realclean: clean
	clear
	rm -rf ../lib ../bin

print-% : ; $(info $* is a $(flavor $*) variable set to [$($*)]) @true
