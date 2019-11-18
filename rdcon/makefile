include ../install/DEFAULTS.inc

FFLAGS = -qopenmp
IFLAGS = -I../equil -I$(MATHINC)
F90 = $(FC) $(FFLAGS) $(IFLAGS)
export FFLAGS # for sub-makes

.f.o:
	$(F90) -c $*.f

LIBDIR = ../lib

LIBS = \
	-lpthread \
	-llsode \
	-lequil \
	-lvac

OBJS = \
	debug.o \
	dcon_mod.o \
	toolbox.o \
	jacobi.o  \
	gal.o \
	mercier.o \
	bal.o \
	fourfit.o \
	sing1.o \
	sing.o \
	ode_output.o \
	ode.o \
	free.o \
	resist.o \
        lib_interface.o \
	dcon.o

all: rdcon

equil:
	cd ../equil; make

rdcon: $(OBJS)
	$(F90) -o rdcon $(OBJS) -L$(LIBDIR) $(LIBS) -L$(MATHDIR) $(MATHLIBS) -L$(FCDIR) $(FCLIBS)
	mkdir -p ../bin
	cp -f rdcon ../bin

# dependencies

dcon_mod.o: ../equil/spline_mod.mod ../equil/global_mod.mod \
	../equil/equil_mod.mod ../equil/equil_out_mod.mod
toolbox.o: dcon_mod.o
jacobi.o: ../equil/local.o
mercier.o: dcon_mod.o
bal.o: dcon_mod.o
fourfit.o: ../equil/fspline_mod.mod dcon_mod.o toolbox.o
sing1.o: fourfit.o
sing.o: sing1.o 
ode_output.o: sing.o debug.o
ode.o: ode_output.o
free.o: ode.o
resist.o: dcon_mod.o toolbox.o
gal.o: jacobi.o sing.o free.o
lib_interface.o: ../equil/local.o ../equil/global.o ../equil/utils.o ../equil/inverse.o
dcon.o: bal.o mercier.o ode.o resist.o gal.o lib_interface.o

clean:
	rm -f *.o *.mod *.out *.bin sum.dat *~ *.diff rdcon

rclean: clean
	cd ../equil; make clean
