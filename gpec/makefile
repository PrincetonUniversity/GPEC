
include ../install/DEFAULTS.inc

FFLAGS = -vec_report0
IFLAGS = -I../equil -I../coil -I../harvest -I$(MATHINC) -I$(NETCDFINC)
F90 = $(FC) $(FFLAGS) $(IFLAGS)
export FFLAGS # for sub-makes

.f.o:
	$(F90) -c $*.f

LIBDIR = ../lib

LIBS = \
	-llsode \
	-lequil \
	-lvac \
	-lcoil \
	-lharvest

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

all: dcon coil gpec

dcon:
	cd ../dcon; make all

coil:
	cd ../coil; make

gpec: $(OBJS)
	$(F90) -o gpec $(OBJS) -L$(LIBDIR) $(LIBS) -L$(MATHDIR) $(MATHLIBS) -L$(NETCDFDIR) $(NETCDFLIBS) -L$(NETCDF_C_DIR) $(NETCDF_C_LIBS) -L$(FCDIR) $(FCLIBS)
	mkdir -p ../lib
	ar -r ../lib/libgpec.a *.o
	mkdir -p ../bin
	cp -f gpec ../bin

# dependencies

ismath.o : gpglobal.o
idcon.o: ismath.o
rdcon.o: gpglobal.o idcon.o ismath.o
gpeq.o: idcon.o ismath.o
gpresp.o: gpeq.o
gpvacuum.o : gpeq.o
gpdiag.o: gpresp.o gpvacuum.o
gpout.o: gpresp.o gpvacuum.o ../coil/field.o ../coil/coil.o
gpec.o: gpdiag.o gpout.o rdcon.o

clean:
	rm -f *.o *.mod *.out *.bin gpec
