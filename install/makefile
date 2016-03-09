
include DEFAULTS.inc

# compile package

.IGNORE:

all: v harvest lsode equil orbit vacuum pentrc dcon match multi sum \
	xdraw coil ipec mkbin

# executables:

harvest:
	cd ../; git submodule init; git submodule update
	cd ../harvest; export GACODE_ROOT=;make CC=$(CC) -f Makefile libharvest.a
	mkdir -p ../lib
	cp -f ../harvest/libharvest.a ../lib/

xdraw:
	cd ../xdraw;  make -f makefile cflag=$(cflag)

lsode:
	cd ../lsode;  make -f makefile cflag=$(cflag)

equil:
	cd ../equil;  make -f makefile cflag=$(cflag)

orbit:
	cd ../orbit;  make -f makefile cflag=$(cflag)

vacuum:
	cd ../vacuum; make -f makefile cflag=$(cflag)

dcon:
	cd ../dcon;   make -f makefile cflag=$(cflag)

match:
	cd ../match;  make -f makefile cflag=$(cflag)

multi:
	cd ../multi;  make -f makefile cflag=$(cflag)

sum:
	cd ../sum;    make -f makefile cflag=$(cflag)

coil:
	cd ../coil;   make -f makefile cflag=$(cflag)

ipec:  
	cd ../ipec;   make -f makefile cflag=$(cflag)

pentrc:  
	cd ../pentrc; make -f makefile cflag=$(cflag)

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
	@echo "   using libraries $(MATHLIBS)"
	@echo "   from $(MATHHOME)"
	@echo "   and libraries $(NETCDFLIBS)"
	@echo "   from $(NETCDFHOME)"
	@echo
	@echo "Set environement FC,LAPACKHOME,NETCDFHOME to change these"
ifeq ("$(PLATFORM)","portal")
	@echo "   example: module load gcc4 lapack netcdf"
endif
	@echo "------------------------------------------------------------"

# cleaning

clean:
	cd ../xdraw;  make -f makefile clean
	cd ../lsode;  make -f makefile clean
	cd ../equil;  make -f makefile clean
	cd ../orbit;  make -f makefile clean
	cd ../vacuum; make -f makefile clean
	cd ../dcon;   make -f makefile clean
	cd ../match;  make -f makefile clean
	cd ../multi;  make -f makefile clean
	cd ../sum;    make -f makefile clean
	cd ../coil;   make -f makefile clean
	cd ../ipec;   make -f makefile clean
	cd ../pentrc; make -f makefile clean

realclean: clean
	rm -rf ../lib ../bin
