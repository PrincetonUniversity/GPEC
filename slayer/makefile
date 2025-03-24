
include ../install/DEFAULTS.inc

IFLAGS = -I../equil -I$(MATHINC) -I$(NETCDFINC)
F90 = $(FC) $(FFLAGS) $(OMPFLAG) $(IFLAGS)
export FFLAGS

.f.o:
	$(F90) -c $*.f

LIBDIR = ../lib

LIBS = \
	-llsode \
	-lequil \
	-lzvode \

OBJS = \
       sglobal.o \
       params.o \
       delta.o \
       gslayer.o \
       slayer.o

all: lsode zvode equil slayer

lsode: 
	cd ../lsode; make

zvode: 
	cd ../zvode; make

equil:
	cd ../equil; make

slayer: $(OBJS)
	$(F90) -o slayer $(OBJS) -L$(LIBDIR) $(LIBS) -L$(MATHDIR) $(MATHLIBS) -L$(NETCDFDIR) $(NETCDFLIBS) $(NETCDF_EXTRA_LIBS)
	mkdir -p ../lib/
	ar -r ../lib/libslayer.a *.o
	mkdir -p ../bin
	cp -f slayer ../bin

# dependencies

params.o : sglobal.o
delta.o : sglobal.o ../equil/cspline_mod.mod
gslayer.o : sglobal.o delta.o
slayer.o : sglobal.o params.o delta.o
sglobal.o : ../equil/cspline_mod.mod

clean:
	rm -f *.o *.mod *.out *.bin slayer
