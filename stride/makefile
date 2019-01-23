# This is the makefile for STRIDE.
include ../install/DEFAULTS.inc

FFLAGS = -qopenmp
# flags unique to stride (not to be clobbered by top level install/makefile options)
EXTRA_FFLAGS = -xHost -mkl -diag-disable 10397 #-03 #-g #-check all #-warn all #-qopt-zmm-usage=high #-ipo #-traceback #-qoptrpt
IFLAGS = -I../equil -I$(MATHINC) -I$(NETCDFINC)
F90 = $(FC) $(FFLAGS) $(EXTRA_FFLAGS) $(IFLAGS)

.f.o:
	$(F90) -c $*.f

LIBDIR = ../lib

SUPERLUDIR = ../superlu/SuperLU_5.2.1/SRC
SUPERLUBLASDIR = ../superlu/

LIBS = \
	-L$(LIBDIR) \
	-llsode \
	-lequil \
	-lvac \
	-L$(SUPERLUDIR) \
	-lsuperlu \
	-L$(SUPERLUBLASDIR) \
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

# targets

all: lsode equil vacuum superlu stride

lsode:
	cd ../lsode; make

equil:
	cd ../equil; make

vacuum:
	cd ../vacuum; make

superlu:
	$(SLU_COMMAND)

stride: $(OBJS)
	$(F90) -o stride $(OBJS) $(LIBS) -L$(MATHDIR) $(MATHLIBS)

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
