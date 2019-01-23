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
	-lzvode \
	-lzlange \
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
	stride_netcdf.o \
	free.o \
	riccati.o \
	stride.o

# targets

all: lsode zvode zlange superlu equil vacuum version stride

lsode:
	cd ../lsode; make

zvode:
	cd ../zvode; make

zlange:
	cd ../zlange; make

superlu:
	$(SLU_COMMAND)

equil:
	cd ../equil; make

vacuum:
	cd ../vacuum; make

stride: $(OBJS)
	$(F90) -o stride $(OBJS) $(LIBS) -L$(NETCDFDIR) $(NETCDFLIBS) -L$(NETCDF_C_DIR) $(NETCDF_C_LIBS)

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

clean:
	rm -f *.o *.mod *.out *.bin *.dat *~ *.diff stride *.log *.optrpt *.dci *__genmod.f90 fort.[0123456789]*

# dependencies

dcon_mod.o: ../equil/spline_mod.mod ../equil/global_mod.mod \
	../equil/equil_mod.mod ../equil/equil_out_mod.mod version.inc
mercier.o: dcon_mod.o
bal.o: dcon_mod.o
fourfit.o: ../equil/fspline_mod.mod dcon_mod.o
sing.o: fourfit.o
zvode1.o: ../zvode/zvode.o ../zlange/zlange.o
riccati.o: sing.o
ode.o: sing.o debug.o zvode1.o ../equil/sparse.o riccati.o
stride_netcdf.o: dcon_mod.o
free.o: ode.o stride_netcdf.o
stride.o: bal.o free.o mercier.o
