# This is the makefile for DCON on LinuxIntel computers.

FFLAGS = -I../equil -g -vec_report0
F90 = ifort $(FFLAGS)

.f.o:
	$(F90) -c $*.f


LIBDIR = \
        -L../lib/LinuxIntel \
        -L /usr/local/lib \
        -L$(HOME)/lib \

LAPACKDIR = -L /u/zwang/numericalLib/lapack_3.4.2

LIBS = \
	$(LIBDIR) \
        -lpthread \
	-llsode \
	-lequil \
	-lvac  \
        $(LAPACKDIR) \
        -llapack \
        -lblas

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
	dcon.o

all: equil dcon

equil:
	cd ../equil; make

dcon: $(OBJS)
	$(F90) -o dcon $(OBJS) $(LIBS)

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
dcon.o: bal.o mercier.o ode.o resist.o gal.o

clean:
	rm -f *.o *.mod *.out *.bin sum.dat *~ *.diff dcon

rclean: clean
	cd ../equil; make clean
