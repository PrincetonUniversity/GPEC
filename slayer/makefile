
include ../install/DEFAULTS.inc

IFLAGS = -I../equil -I../harvest -I../pentrc -I$(MATHINC) -I$(NETCDFINC)
F90 = $(FC) $(FFLAGS) $(OMPFLAG) $(IFLAGS)
export FFLAGS

.f.o:
	$(F90) -c $*.f

LIBDIR = ../lib

LIBS = \
	-llsode \
	-lpentrc \
	-lequil \
	-lharvest

OBJS = \
       sglobal.o \
       params.o \
       layerinputs.o \
       delta.o \
       gslayer.o \
       slayer.o

all: equil harvest lsode pentrc slayer

harvest:
	cd ../; git submodule init; git submodule update
	cd ../harvest; make FC=$(FC) CC=$(CC) GACODE_ROOT= -f Makefile libharvest.a
	mkdir -p ../lib
	cp -f ../harvest/libharvest.a ../lib/

equil:
	cd ../equil; make

lsode:
	cd ../lsode; make

pentrc:
	cd ../pentrc; make pentrc

slayer: $(OBJS)
	$(F90) -o slayer $(OBJS) -L$(LIBDIR) $(LIBS) -L$(MATHDIR) $(MATHLIBS) -L$(NETCDFDIR) $(NETCDFLIBS) $(NETCDF_EXTRA_LIBS)
	mkdir -p ../lib/
	ar -r ../lib/libslayer.a *.o
	mkdir -p ../bin
	cp -f slayer ../bin

# dependencies

params.o : sglobal.o
delta.o : sglobal.o
gslayer.o : sglobal.o layerinputs.o delta.o
slayer.o : sglobal.o params.o delta.o ../pentrc/grid.mod ../equil/spline_mod.mod

clean:
	rm -f *.o *.mod *.out *.bin slayers

