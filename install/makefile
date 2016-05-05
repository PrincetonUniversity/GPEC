
include DEFAULTS.inc

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
	cd ../xdraw;  make -f makefile cflag=$(cflag) options=$(options)

lsode:
	cd ../lsode;  make -f makefile cflag=$(cflag) options=$(options)

equil:
	cd ../equil;  make -f makefile cflag=$(cflag) options=$(options)

orbit:
	cd ../orbit;  make -f makefile cflag=$(cflag) options=$(options)

vacuum:
	cd ../vacuum; make -f makefile cflag=$(cflag) options=$(options)

dcon:
	cd ../dcon;   make -f makefile cflag=$(cflag) options=$(options)

match:
	cd ../match;  make -f makefile cflag=$(cflag) options=$(options)

multi:
	cd ../multi;  make -f makefile cflag=$(cflag) options=$(options)

sum:
	cd ../sum;    make -f makefile cflag=$(cflag) options=$(options)

coil:
	cd ../coil;   make -f makefile cflag=$(cflag) options=$(options)

ipec:  
	cd ../ipec;   make -f makefile cflag=$(cflag) options=$(options)

pentrc:  
	cd ../pentrc; make -f makefile cflag=$(cflag) options=$(options)

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
	@echo "   from $(MATHDIR)"
	@echo "   and libraries $(NETCDFLIBS)"
	@echo "   from $(NETCDFDIR)"
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
