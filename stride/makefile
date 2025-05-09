# This is the makefile for STRIDE.

# If using gfortran, the proper syntax is -frecursive.
# If using ifort, the proper syntax is -recursive.

ifneq (,$(findstring ifort,$(FC)))
    RECURSFLAG=-recursive
else ifneq (,$(findstring ifx,$(FC)))
	RECURSFLAG=-recursive
else
    RECURSFLAG=-frecursive
endif

include ../install/DEFAULTS.inc

# flags unique to stride (not to be clobbered by top level install/makefile options)
IFLAGS = -I../equil -I../vacuum -I$(MATHINC) -I$(NETCDFINC)
F90 = $(FC) $(FFLAGS) $(OMPFLAG) $(IFLAGS) $(RECURSFLAG)

.f.o:
	$(F90) -c $*.f

.F.o:
	$(F90) -c $*.F

LIBDIR = ../lib

LIBS = \
	-L$(LIBDIR) \
	-lzvode \
	-lzlange \
	-llsode \
	-lequil \
	-lvac \
	-L$(MATHDIR) \
	$(MATHLIBS)

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

all: lsode zvode zlange equil vacuum stride

lsode:
	cd ../lsode; make

zvode:
	cd ../zvode; make

zlange:
	cd ../zlange; make

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
	rm -f *.o *.mod *.out *.bin *.dat *~ *.diff stride *.log *.optrpt *.dci *__genmod.f90 *.original

# dependencies

dcon_mod.o: ../equil/spline_mod.mod ../equil/global_mod.mod \
	../equil/equil_mod.mod ../equil/equil_out_mod.mod version.inc
mercier.o: dcon_mod.o
bal.o: dcon_mod.o
fourfit.o: ../equil/fspline_mod.mod dcon_mod.o
sing.o: fourfit.o
zvode1.o: ../zvode/zvode.o ../zlange/zlange.o
riccati.o: sing.o
ode.o: sing.o debug.o zvode1.o riccati.o
stride_netcdf.o: dcon_mod.o
free.o: ode.o stride_netcdf.o
stride.o: bal.o free.o mercier.o
