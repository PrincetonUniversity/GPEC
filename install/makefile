
include DEFAULTS.inc

# extend these variables to sub-makes
FFLAGS ?=
CFLAGS ?=
export FFLAGS
export CFLAGS

# compile package

.IGNORE:

all: v harvest lsode equil orbit vacuum pentrc dcon match multi sum \
	xdraw coil ipec mkbin

allf: v harvest lsode equil orbit vacuum pentrc dcon match multi sum \
	coil ipec mkbin


# executables:

harvest:
	cd ../; git submodule init; git submodule update
	cd ../harvest; make FC=$(FC) CC=$(CC) GACODE_ROOT= -f Makefile libharvest.a
	mkdir -p ../lib
	cp -f ../harvest/libharvest.a ../lib/

xdraw:
	cd ../xdraw;  make

lsode:
	cd ../lsode;  make liblsode

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

multi:
	cd ../multi;  make multi

sum:
	cd ../sum;    make sum

coil:
	cd ../coil;   make libcoil

ipec:  
	cd ../ipec;   make ipec

pentrc:  
	cd ../pentrc; make pentrc

mkbin:
	mkdir -p ../bin
	cp -f ../orbit/orbit ../bin
	cp -f ../dcon/dcon ../bin
	cp -f ../xdraw/xdraw ../bin
	cp -f ../match/match ../bin
	cp -f ../multi/multi ../bin
	cp -f ../sum/sum ../bin
	cp -f ../ipec/ipec ../bin
	cp -f ../pentrc/pentrc ../bin
	cp -f ../input/mput ../bin

# print verbose info

v:
	@echo "------------------------------------------------------------"
	@echo "Defaults tailored for $(PLATFORM)"
	@echo "Compiling with FC = $(FC)"
	@echo "   using math libraries $(MATHLIBS)"
	@echo "   from $(MATHDIR)"
	@echo "   including $(MATHINC)"
	@echo "   and netcdf libraries $(NETCDFLIBS)"
	@echo "   from $(NETCDFDIR)"
	@echo "   including $(NETCDFINC)"
	@echo "   and compiler libraries $(FCLIBS)"
	@echo "   from $(FCDIR)"
	@echo "   including $(FCINC)"
	@echo
	@echo "Set environement FC,LAPACKHOME,NETCDFHOME to change these"
ifeq ("$(PLATFORM)","portal")
	@echo "   example: module load gcc4 lapack netcdf"
endif
	@echo "------------------------------------------------------------"

# cleaning

clean:
	cd ../harvest; make clean
	cd ../xdraw;  make clean
	cd ../lsode;  make clean
	cd ../equil;  make clean
	cd ../orbit;  make clean
	cd ../vacuum; make clean
	cd ../dcon;   make clean
	cd ../match;  make clean
	cd ../multi;  make clean
	cd ../sum;    make clean
	cd ../coil;   make clean
	cd ../ipec;   make clean
	cd ../pentrc; make clean

realclean: clean
	rm -rf ../lib ../bin
