# If using gfortran, the proper syntax is -frecursive.
# If using ifort, the proper syntax is -recursive.

ifneq (,$(findstring ifort,$(FC)))
    RECURSFLAG=-recursive
	LEGACYFLAG=
else ifneq (,$(findstring ifx,$(FC)))
	RECURSFLAG=-recursive
	LEGACYFLAG=
else
    RECURSFLAG=-frecursive
	LEGACYFLAG=-std=legacy
endif



include ../install/DEFAULTS.inc

IFLAGS = -I../equil -I../harvest -I$(MATHINC) -I$(NETCDFINC)
F90 = $(FC) $(FFLAGS) $(OMPFLAG) $(IFLAGS) $(RECURSFLAG)
export FFLAGS

.f.o:
	$(F90) $(LEGACYFLAG) -c $*.f

%.o: %.f90
	$(F90) -c $<

%.o: %.F90
	$(F90) -c $<

LIBDIR = ../lib

LIBS = \
	-llsode \
	-lequil \
	-lharvest \

OBJS = \
    params.o \
	utilities.o \
	dcon_interface.o \
	lsode1.o \
	lsode2.o \
	special.o \
    inputs.o \
	energy.o \
    pitch.o \
    torque.o \
	diagnostics.o \
	pentrc_interface.o \
	pentrc.o

all: harvest lsode equil pentrc

harvest:
	cd ../; git submodule init; git submodule update
	cd ../harvest; make FC=$(FC) CC=$(CC) GACODE_ROOT= -f Makefile libharvest.a
	mkdir -p ../lib
	cp -f ../harvest/libharvest.a ../lib/

lsode:
	cd ../lsode; make

equil:
	cd ../equil; make

ifeq ($(wildcard version.inc),)
version.inc: force
	@echo ">>> Creating version file"
	@echo "      CHARACTER(len=*), PARAMETER :: version ='"`git describe --tags`"'" | tee version.inc
else ifeq ($(shell grep -F "'"`git describe --tags`"'" version.inc),)
version.inc: force
	@echo ">>> Updating version file"
	@echo "      CHARACTER(len=*), PARAMETER :: version ='"`git describe --tags`"'" | tee version.inc
else
version.inc:
	@echo ">>> Git version unchanged!"
endif

force:

pentrc: $(OBJS)
	$(F90) -o pentrc $(OBJS) -L$(LIBDIR) $(LIBS) -L$(MATHDIR) $(MATHLIBS) -L$(NETCDFDIR) $(NETCDFLIBS) $(NETCDF_EXTRA_LIBS)
	mkdir -p ../lib
	ar -r ../lib/libpentrc.a *.o
	mkdir -p ../bin
	cp -f pentrc ../bin

# dependencies

lsode1.o: ../lsode/lsode.o
lsode2.o: ../lsode/lsode.o
params.o: version.inc ../equil/local_mod.mod
utilities.o: params.o
special.o: params.o utilities.o
dcon_interface.o: params.o utilities.o ../equil/spline_mod.mod ../equil/fspline_mod.mod ../equil/cspline_mod.mod ../equil/bicube_mod.mod
inputs.o:  params.o utilities.o special.o dcon_interface.o
energy.o:  params.o utilities.o lsode2.o
pitch.o:   params.o utilities.o lsode1.o energy.o special.o
torque.o:  params.o utilities.o lsode1.o special.o dcon_interface.o inputs.o energy.o pitch.o
diagnostics.o: params.o special.o energy.o pitch.o
pentrc_interface.o: params.o utilities.o special.o inputs.o energy.o pitch.o torque.o diagnostics.o
pentrc.o:  params.o  pentrc_interface.o

clean:
	rm -f *.o *.mod *.out *.bin pentrc *.original
