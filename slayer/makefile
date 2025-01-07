
include ../install/DEFAULTS.inc

IFLAGS = -I$(MATHINC) -I$(NETCDFINC)
F90 = $(FC) $(FFLAGS) $(OMPFLAG) $(IFLAGS)
export FFLAGS

.f.o:
	$(F90) -c $*.f

LIBDIR = ../lib

LIBS = \
	-llsode \
	-lzvode \

OBJS = \
       sglobal.o \
       params.o \
       delta.o \
       gslayer.o \
       slayer.o

all: lsode zvode slayer

lsode: 
	cd ../lsode; make

zvode: 
	cd ../zvode; make

slayer: $(OBJS)
	$(F90) -o slayer $(OBJS) -L$(LIBDIR) $(LIBS) -L$(MATHDIR) $(MATHLIBS) -L$(NETCDFDIR) $(NETCDFLIBS) $(NETCDF_EXTRA_LIBS)
	mkdir -p ../lib/
	ar -r ../lib/libslayer.a *.o
	mkdir -p ../bin
	cp -f slayer ../bin

# dependencies

params.o : sglobal.o
delta.o : sglobal.o
gslayer.o : sglobal.o delta.o
slayer.o : sglobal.o params.o delta.o

clean:
	rm -f *.o *.mod *.out *.bin slayer
