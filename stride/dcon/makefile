# This is the makefile for STRIDE.

FFLAGS = -I../equil -xHost -O3 -qopenmp -mkl -diag-disable 10397 -warn all #-qopt-zmm-usage=high #-ipo
#FFLAGS = -I../equil -xHost -qopenmp -mkl -diag-disable 10397 -warn all -g -check all #-traceback #-qoptrpt
F90 = ifort $(FFLAGS)

.f.o:
	$(F90) -c $*.f

LIBDIR = \
	-L../lib \
	-L$(HOME)/lib

SUPERLUDIR = -L ../SuperLU/SuperLU_5.2.1/SRC
SUPERLUBLASDIR = -L ../SuperLU/

LIBS = \
	$(LIBDIR) \
	-lzvode \
	-lzlange \
	-llsode \
	-lequil \
	-lvac \
	$(SUPERLUDIR) \
	-lsuperlu \
	$(SUPERLUBLASDIR) \
	-lblas

OBJS = \
	debug.o \
	dcon_mod.o \
	bal.o \
	mercier.o \
	fourfit.o \
	sing.o \
	zvode1.o \
	ode.o \
	free.o \
	riccati.o \
	dcon.o

all: equil stride

equil:
	cd ../equil; make

stride: $(OBJS)
	$(F90) -o stride $(OBJS) $(LIBS)

# dependencies

dcon_mod.o: ../equil/spline_mod.mod ../equil/global_mod.mod \
	../equil/equil_mod.mod ../equil/equil_out_mod.mod
mercier.o: dcon_mod.o
bal.o: dcon_mod.o
fourfit.o: ../equil/fspline_mod.mod dcon_mod.o
sing.o: fourfit.o
zvode1.o:
riccati.o: sing.o
ode.o: sing.o debug.o zvode1.o ../equil/sparse.o riccati.o
free.o: ode.o
dcon.o: bal.o free.o mercier.o

clean:
	rm -f *.o *.mod *.out *.bin *.dat *~ *.diff stride *.log *.optrpt *.dci *__genmod.f90 fort.[0123456789]*
