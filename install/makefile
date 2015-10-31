#------------------------------------------------------------------
# General GNU make for portal, venus, or OSX platforms
# - Conditionals used to determine platform 
# - This defines platform specific defaults for FC,LAPACKHOME,NETCDFHOME
#
# On portal.pppl.gov (or portalr6.pppl.gov for more up-to-date libraries)
#   Must load compiler and basic libraries (These modify LD_LIBRARY_PATH).
#   Choose one of the following
#	$ module load lff95 lapack-blas netcdf
#	$ module load gcc4 lapack-blas netcdf
#	$ module load intel acml netcdf
#	Defaults will try to use gcc4
#
# On venus.gat.com
# - All paths are hardcoded at this point (08/2015).
#   -- We are waiting on a upgrade that will enable "module load" standardization
#
# On mac OSX
#   Prerequisites
#   - Install gfortran in your path.
#	  -- Available from macports or sourceforge
#   - Install lapack, blas, and netcdf libraries in usr/lib.
#
#
# If your platform is not properly included in this make file
# contact nlogan@pppl.gov with your operating system, compiler, and
# library installation details.
#
#------------------------------------------------------------------

# set platform dependent defaults

ifneq ($(shell echo $(HOSTNAME) | grep sunfire),)
	PLATFORM = portal
ifndef FC
ifndef F90
	FC = /usr/pppl/gcc/4.6.3/bin/gfortran4
else
	FC = $(F90)
endif
endif
ifndef LAPACKHOME
ifndef ACML_HOME
	LAPACKHOME = /usr/pppl/gcc/4.6-pkgs/lapack-3.3.1
else
	LAPACKHOME = $(ACML_HOME)
endif
endif
ifndef NETCDFHOME
	NETCDFHOME = /usr/pppl/gcc/4.6-pkgs/netcdf-4.1.2
endif
endif

ifneq ($(shell echo $(HOSTNAME) | grep venus),)
	PLATFORM = venus
ifndef FC
	FC = /task/imd/local64/bin/mpif90
endif
ifndef LAPACKHOME
	LAPACKHOME = /opt/pgi/linux86-64/12.10/lib #pgf
endif
ifndef NETCDFHOME
	NETCDFHOME = /task/imd/local64/lib #mpif90
endif
endif

ifneq ($(shell echo $(MACHTYPE) | grep apple),)
	PLATFORM = osx
ifndef FC
	FC = gfortran
endif
ifndef LAPACKHOME
	LAPACKHOME = /usr/lib
endif
ifndef NETCDFHOME
	NETCDFHOME = /usr/lib
endif
endif

# compile package

.IGNORE:

all: v harvest lsode equil orbit vacuum dcon match multi sum \
	xdraw coil ipec pentrc mkbin

# executables:

harvest:
	cd ../; git submodule init; git submodule update
	cd ../harvest; export GACODE_ROOT=;make CC=$(CC) -f Makefile libharvest.a
	mkdir -p ../lib
	cp -f ../harvest/libharvest.a ../lib/

xdraw:
	cd ../xdraw; make -f makefile FC=$(FC) LAPACKHOME=$(LAPACKHOME)

lsode:
	cd ../lsode; make -f makefile FC=$(FC) LAPACKHOME=$(LAPACKHOME)

equil:
	cd ../equil; make -f makefile FC=$(FC) LAPACKHOME=$(LAPACKHOME)

orbit:
	cd ../orbit; make -f makefile FC=$(FC) LAPACKHOME=$(LAPACKHOME)

vacuum:
	cd ../vacuum; make -f makefile FC=$(FC) LAPACKHOME=$(LAPACKHOME)

dcon:
	cd ../dcon; make -f makefile FC=$(FC) LAPACKHOME=$(LAPACKHOME)

match:
	cd ../match; make -f makefile FC=$(FC) LAPACKHOME=$(LAPACKHOME)

multi:
	cd ../multi; make -f makefile FC=$(FC) LAPACKHOME=$(LAPACKHOME)

sum:
	cd ../sum; make -f makefile FC=$(FC) LAPACKHOME=$(LAPACKHOME)

coil:
	cd ../coil; make -f makefile FC=$(FC) LAPACKHOME=$(LAPACKHOME)

ipec:  
	cd ../ipec; make -f makefile FC=$(FC) LAPACKHOME=$(LAPACKHOME)

pentrc:  
	cd ../pentrc; make -f makefile FC=$(FC) LAPACKHOME=$(LAPACKHOME)

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
	@echo "   using lapack libraries from $(LAPACKHOME)"
	@echo "   using netcdf libraries from $(NETCDFHOME)"
	@echo
	@echo "Set environement FC,LAPACKHOME,NETCDFHOME to change these"
ifeq ("$(PLATFORM)","portal")
	@echo "   example: module load gcc4 lapack netcdf"
endif
	@echo "------------------------------------------------------------"

# cleaning

clean:
	cd ../xdraw; make -f makefile FC=$(FC) LAPACKHOME=$(LAPACKHOME) clean
	cd ../lsode; make -f makefile FC=$(FC) LAPACKHOME=$(LAPACKHOME) clean
	cd ../equil; make -f makefile FC=$(FC) LAPACKHOME=$(LAPACKHOME) clean
	cd ../orbit; make -f makefile FC=$(FC) LAPACKHOME=$(LAPACKHOME) clean
	cd ../vacuum; make -f makefile FC=$(FC) LAPACKHOME=$(LAPACKHOME) clean
	cd ../dcon; make -f makefile FC=$(FC) LAPACKHOME=$(LAPACKHOME) clean
	cd ../match; make -f makefile FC=$(FC) LAPACKHOME=$(LAPACKHOME) clean
	cd ../multi; make -f makefile FC=$(FC) LAPACKHOME=$(LAPACKHOME) clean
	cd ../sum; make -f makefile FC=$(FC) LAPACKHOME=$(LAPACKHOME) clean
	cd ../coil; make -f makefile FC=$(FC) LAPACKHOME=$(LAPACKHOME) clean
	cd ../ipec; make -f makefile FC=$(FC) LAPACKHOME=$(LAPACKHOME) clean
	cd ../pentrc; make -f makefile FC=$(FC) LAPACKHOME=$(LAPACKHOME) clean

realclean: clean
	rm -rf ../lib ../bin
