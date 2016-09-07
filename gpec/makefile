
include ../install/DEFAULTS.inc

FFLAGS =
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
	ipglobal.o \
	ismath.o \
	idcon.o \
	rdcon.o \
	ipeq.o \
	ipresp.o \
	ipdiag.o \
	ipvacuum.o \
	ipout.o \
	gpec.o

all: dcon coil gpec

dcon:
	cd ../dcon; make all

coil:
	cd ../coil; make

gpec: $(OBJS)
	$(F90) -o gpec $(OBJS) -L$(LIBDIR) $(LIBS) -L$(MATHDIR) $(MATHLIBS) -L$(NETCDFDIR) $(NETCDFLIBS) -L$(FCDIR) $(FCLIBS)
	mkdir -p ../lib
	ar -r ../lib/libgpec.a *.o
	mkdir -p ../bin
	cp -f gpec ../bin

# dependencies

ismath.o : ipglobal.o
idcon.o: ismath.o
rdcon.o: ipglobal.o idcon.o ismath.o
ipeq.o: idcon.o ismath.o
ipresp.o: ipeq.o
ipvacuum.o : ipeq.o
ipdiag.o: ipresp.o ipvacuum.o
ipout.o: ipresp.o ipvacuum.o ../coil/field.o ../coil/coil.o
gpec.o: ipdiag.o ipout.o rdcon.o

clean:
	rm -f *.o *.mod *.out *.bin gpec
