include ../install/DEFAULTS.inc

IFLAGS = -I../equil -I../vacuum -I$(MATHINC) -I$(NETCDFINC)
F90 = $(FC) $(FFLAGS) $(OMPFLAG) $(IFLAGS)
export FFLAGS # for sub-makes

.f.o:
	$(F90) -c $*.f

LIBDIR = ../lib

LIBS = \
	-lpthread \
	-llsode \
	-lequil \
	-lvac \

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
	dcon_netcdf.o \
	free.o \
	resist.o \
        lib_interface.o \
	dcon.o

all: lsode equil vacuum rdcon

lsode:
	cd ../lsode; make

equil:
	cd ../equil; make

vacuum:
	cd ../vacuum; make

rdcon: $(OBJS)
	$(F90) -o rdcon $(OBJS) -L$(LIBDIR) $(LIBS) -L$(MATHDIR) $(MATHLIBS) -L$(FCDIR) $(FCLIBS) -L$(NETCDFDIR) $(NETCDFLIBS) $(NETCDF_EXTRA_LIBS) $(LDFLAGS)
	mkdir -p ../bin
	cp -f rdcon ../bin

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
	rm -f *.o *.mod *.out *.bin sum.dat *~ *.diff rdcon *.original

rclean: clean
	cd ../equil; make clean

# dependencies

dcon_mod.o: ../equil/spline_mod.mod ../equil/global_mod.mod \
	../equil/equil_mod.mod ../equil/equil_out_mod.mod version.inc
toolbox.o: dcon_mod.o
jacobi.o: ../equil/local_mod.mod
mercier.o: dcon_mod.o
bal.o: dcon_mod.o
fourfit.o: ../equil/fspline_mod.mod dcon_mod.o toolbox.o
sing1.o: fourfit.o
sing.o: sing1.o
ode_output.o: sing.o debug.o
ode.o: ode_output.o
dcon_netcdf.o: dcon_mod.o
free.o: ode.o dcon_netcdf.o
resist.o: dcon_mod.o toolbox.o
gal.o: jacobi.o sing.o free.o
lib_interface.o: ../equil/local_mod.mod ../equil/global.o ../equil/utils.o ../equil/inverse.o
dcon.o: bal.o mercier.o ode.o resist.o gal.o lib_interface.o