
include DEFAULTS.inc

# extend these variables to sub-makes (?= sets a default if user doesn't specify options)
FFLAGS ?= -g
CFLAGS ?= -g

export FFLAGS
export CFLAGS
export OMPFLAG
export LDFLAGS

# compile package

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
	cd ../xdraw;  make -j 8

lsode:
	cd ../lsode;  make -j 8 liblsode

zlange:
	cd ../zlange;  make -j 8 libzlange

zvode:
	cd ../zvode;  make -j 8 libzvode

equil:
	cd ../equil;  make -j 8 libequil

vacuum:
	cd ../vacuum; make -j 8 libvac

orbit:
	cd ../orbit;  make -j 8 orbit

dcon:
	cd ../dcon;   make -j 8 dcon

match:
	cd ../match;  make -j 8 match

rdcon:
	cd ../rdcon;   make -j 8 rdcon

rmatch:
	cd ../rmatch;  make -j 8 rmatch

multi:
	cd ../multi;  make -j 8 multi

sum:
	cd ../sum;    make -j 8 sum

slayer:
	cd ../slayer; make -j 8 slayer

coil:
	cd ../coil;   make -j 8 libcoil

gpec:
	cd ../gpec;   make -j 8 gpec

pentrc:
	cd ../pentrc; make -j 8 pentrc

stride:
	cd ../stride; make -j 8 stride

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
