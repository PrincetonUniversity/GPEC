
include ../install/DEFAULTS.inc

IFLAGS = -I../equil -I../harvest -I../pentrc -I../vacuum -I$(MATHINC) -I$(NETCDFINC)
F90 = $(FC) $(FFLAGS) $(OMPFLAG) $(IFLAGS) $(RECURSFLAG) $(LEGACYFLAG)
export FFLAGS # for sub-makes

.f.o:
	$(F90) -c $*.f

.F.o:
	$(F90) -c $*.F

LIBDIR = ../lib

LIBS = \
	-lpentrc \
	-llsode \
	-lequil \
	-lvac \
	-lharvest \

OBJS = \
	debug.o \
	dcon_mod.o \
	mercier.o \
	bal.o \
	fourfit.o \
	sing.o \
	ode_output.o \
	free.o \
	ode.o \
	dcon_netcdf.o \
	resist.o \
	dcon.o

# targets

all: harvest lsode equil vacuum pentrc dcon

harvest:
	cd ../; git submodule init; git submodule update
	cd ../harvest; make FC=$(FC) CC=$(CC) GACODE_ROOT= -f Makefile libharvest.a
	mkdir -p ../lib
	cp -f ../harvest/libharvest.a ../lib/

lsode:
	cd ../lsode; make

equil:
	cd ../equil; make

vacuum:
	cd ../vacuum; make

pentrc:
	cd ../pentrc; make pentrc

# note we explicitly make pentrc to force sync the compile-time version variable
dcon: pentrc $(OBJS)
	$(F90) -o dcon $(OBJS) -L$(LIBDIR) $(LIBS) -L$(MATHDIR) $(MATHLIBS) -L$(NETCDFDIR) $(NETCDFLIBS) $(NETCDF_EXTRA_LIBS) $(LDFLAGS)
	mkdir -p ../bin
	cp -f dcon ../bin

clean:
	rm -f *.o *.mod *.out *.bin sum.dat *~ *.diff dcon *.original

# dependencies

dcon_mod.o: ../equil/spline_mod.mod ../equil/global_mod.mod \
	../equil/equil_mod.mod ../equil/equil_out_mod.mod ../pentrc/dcon_interface.mod
mercier.o: dcon_mod.o
bal.o: dcon_mod.o
fourfit.o: ../equil/fspline_mod.mod ../equil/global_mod.mod \
	../pentrc/pentrc_interface.mod ../pentrc/utilities.mod ../pentrc/torque.mod ../pentrc/inputs.mod \
	../pentrc/energy_integration.mod ../pentrc/pitch_integration.mod ../pentrc/dcon_interface.mod \
	dcon_mod.o
sing.o: fourfit.o
debug.o: ../equil/local_mod.mod
ode_output.o: dcon_mod.o sing.o debug.o
ode.o: ode_output.o free.o
dcon_netcdf.o: dcon_mod.o
free.o: ../equil/global_mod.mod dcon_netcdf.o ode_output.o
resist.o: dcon_mod.o
dcon.o: ../equil/equil_mod.mod ../equil/equil_out_mod.mod \
	../pentrc/inputs.o ../pentrc/pentrc_interface.o \
	bal.o mercier.o ode.o free.o resist.o
