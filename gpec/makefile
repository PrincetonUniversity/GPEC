
include ../install/DEFAULTS.inc

IFLAGS = -I../equil -I../harvest -I../pentrc -I../coil -I../slayer -I$(MATHINC) -I$(NETCDFINC)
F90 = $(FC) $(FFLAGS) $(OMPFLAG) $(IFLAGS)
export FFLAGS # for sub-makes

.f.o:
	$(F90) -c $*.f

LIBDIR = ../lib

LIBS = \
	-lpentrc \
	-llsode \
	-lequil \
	-lvac \
	-lharvest \
	-lslayer \
	-lcoil

OBJS = \
	gpglobal.o \
	ismath.o \
	idcon.o \
	rdcon.o \
	gpeq.o \
	gpresp.o \
	gpdiag.o \
	gpvacuum.o \
	gpout.o \
	gpec.o

# targets

all: pentrc coil gpec

# takes care of harvest, equil, lsode, etc.
pentrc:
	cd ../pentrc; make all

coil:
	cd ../coil; make


# decide if we need to update the version
# note gpec keeps track of its own version, which can differ from the synced dcon & pentrc one
# this enables faster compiling during development
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

gpec: $(OBJS)
	$(F90) -o gpec $(OBJS) -L$(LIBDIR) $(LIBS) -L$(MATHDIR) $(MATHLIBS) -L$(NETCDFDIR) $(NETCDFLIBS) $(NETCDF_EXTRA_LIBS)
	mkdir -p ../lib
	ar -r ../lib/libgpec.a *.o
	mkdir -p ../bin
	cp -f gpec ../bin

clean:
	rm -f *.o *.mod *.out *.bin gpec

# dependencies

gpglobal.o: ../equil/fspline_mod.mod ../equil/bicube_mod.mod
ismath.o: gpglobal.o  ../equil/local_mod.mod
idcon.o: gpglobal.o ismath.o
rdcon.o: gpglobal.o ismath.o idcon.o
gpeq.o: idcon.o ismath.o
gpresp.o: gpeq.o
gpvacuum.o : gpeq.o
gpdiag.o: gpresp.o gpvacuum.o
gpout.o: gpdiag.o gpresp.o gpvacuum.o ../coil/field.o ../coil/coil.o \
    ../equil/equil_mod.mod ../equil/equil_out_mod.mod \
    ../pentrc/pentrc_interface.mod ../pentrc/inputs.mod \
    ../slayer/gslayer_mod.mod \
	version.inc 
gpec.o: gpdiag.o gpout.o rdcon.o
