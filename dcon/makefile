
include ../install/DEFAULTS.inc

FFLAGS = -vec_report0 -openmp
IFLAGS = -I../equil -I../harvest -I../pentrc -I$(MATHINC) -I$(NETCDFINC)
F90 = $(FC) $(FFLAGS) $(IFLAGS)
export FFLAGS # for sub-makes

.f.o:
	$(F90) -c $*.f

LIBDIR = ../lib

LIBS = \
	-llsode \
	-lequil \
	-lvac \
	-lharvest \
	-lpentrc

OBJS = \
	debug.o \
	dcon_mod.o \
	mercier.o \
	bal.o \
	fourfit.o \
	sing.o \
	ode_output.o \
	ode.o \
	dcon_netcdf.o \
	free.o \
	resist.o \
	dcon.o

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
	cd ../pentrc; make

dcon: $(OBJS)
	$(F90) -o dcon $(OBJS) -L$(LIBDIR) $(LIBS) -L$(MATHDIR) $(MATHLIBS) -L$(NETCDFDIR) $(NETCDFLIBS) -L$(NETCDF_C_DIR) $(NETCDF_C_LIBS) -L$(FCDIR) $(FCLIBS)
	mkdir -p ../bin
	cp -f dcon ../bin

# dependencies

dcon_mod.o: ../equil/spline_mod.mod ../equil/global_mod.mod \
	../equil/equil_mod.mod ../equil/equil_out_mod.mod
mercier.o: dcon_mod.o
bal.o: dcon_mod.o
fourfit.o: ../equil/fspline_mod.mod dcon_mod.o
sing.o: fourfit.o
ode_output.o: dcon_mod.o sing.o debug.o
ode.o: ode_output.o
dcon_netcdf.o: dcon_mod.o
free.o: ode.o dcon_netcdf.o
resist.o: dcon_mod.o
dcon.o: bal.o mercier.o ode.o free.o resist.o ../pentrc/inputs.o ../pentrc/dcon_interface.o ../pentrc/pentrc_interface.o

clean:
	rm -f *.o *.mod *.out *.bin sum.dat *~ *.diff dcon
