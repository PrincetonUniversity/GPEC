HOST := $(shell hostname)

ifneq (,$(findstring flux,$(HOST)))
#    include /p/transpusers/mgorelen/transp_main/environs/ppplflux/gnu/13.2.0/environ
    include /p/transpusers/mgorelen/transp_main/environs/ppplflux/gnu/13.2.0/makefile.mk
else
	include DEFAULTS.inc
endif


# extend these variables to sub-$(MAKE)s (?= sets a default if user doesn't specify options)
FFLAGS ?= -g
CFLAGS ?= -g

export FFLAGS
export CFLAGS
export OMPFLAG
export LDFLAGS
ifneq (,$(findstring flux,$(HOST)))
	FFLAGS := $(filter-out -c,$(FFLAGS))
endif
FFLAGS := $(filter-out -c,$(FFLAGS))
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
	cd ../harvest; $(MAKE) FC='$(FC)' CC='$(CC)' GACODE_ROOT= -f Makefile libharvest.a
	mkdir -p ../lib
	cp -f ../harvest/libharvest.a ../lib/

xdraw:
	cd ../xdraw;  $(MAKE)

lsode: equil
	cd ../lsode;  $(MAKE) liblsode

zlange: equil
	cd ../zlange;  $(MAKE) libzlange

zvode:
	cd ../zvode;  $(MAKE) libzvode

equil:
	cd ../equil;  $(MAKE) libequil

vacuum: equil
	cd ../vacuum; $(MAKE) libvac

orbit: equil lsode
	cd ../orbit;  $(MAKE) orbit

dcon: equil pentrc vacuum
	cd ../dcon;   $(MAKE) dcon

match: equil
	cd ../match;  $(MAKE) match

rdcon: equil
	cd ../rdcon;   $(MAKE) rdcon

rmatch: equil
	cd ../rmatch;  $(MAKE) rmatch

multi: equil
	cd ../multi;  $(MAKE) multi

sum: equil
	cd ../sum;    $(MAKE) sum

slayer: equil
	cd ../slayer; $(MAKE) slayer

coil: equil
	cd ../coil;   $(MAKE) libcoil

gpec: equil coil vacuum pentrc
	cd ../gpec;   $(MAKE) gpec

pentrc: equil
	cd ../pentrc; $(MAKE) pentrc

stride: equil vacuum
	cd ../stride; $(MAKE) stride

mkbin: orbit dcon match rdcon rmatch xdraw multi sum slayer gpec pentrc stride
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
	cd ../harvest; $(MAKE) clean
	cd ../zlange;  $(MAKE) clean
	cd ../zvode;  $(MAKE) clean
	cd ../xdraw;  $(MAKE) clean
	cd ../lsode;  $(MAKE) clean
	cd ../equil;  $(MAKE) clean
	cd ../orbit;  $(MAKE) clean
	cd ../vacuum; $(MAKE) clean
	cd ../dcon;   $(MAKE) clean
	cd ../match;  $(MAKE) clean
	cd ../rdcon;  $(MAKE) clean
	cd ../rmatch; $(MAKE) clean
	cd ../multi;  $(MAKE) clean
	cd ../sum;    $(MAKE) clean
	cd ../slayer; $(MAKE) clean
	cd ../coil;   $(MAKE) clean
	cd ../gpec;   $(MAKE) clean
	cd ../pentrc; $(MAKE) clean
	cd ../stride; $(MAKE) clean

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
