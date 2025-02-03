
include ../install/DEFAULTS.inc

IFLAGS = -I../equil -I$(MATHINC)
F90 = $(FC) $(FFLAGS) $(OMPFLAG) $(IFLAGS)
export FFLAGS

.f.o:
	$(F90) -c $*.f

LIBDIR = ../lib

LIBS = \
	-llsode \
	-lequil

OBJS = \
	orbit_mod.o \
	restart.o \
	field.o \
	inner.o \
	outer.o \
	inner_gc.o \
	outer_gc.o \
	orbit.o

all: lsode equil orbit

lsode:
	cd ../lsode; make

equil:
	cd ../equil; make

orbit: $(OBJS)
	$(F90) -o orbit $(OBJS) -L$(LIBDIR) $(LIBS) -L$(MATHDIR) $(MATHLIBS) -L$(FCDIR) $(FCLIBS)
	mkdir -p ../bin
	cp -f orbit ../bin

# dependencies

orbit_mod.o: ../equil/equil_mod.mod ../equil/equil_out_mod.mod
restart.o: orbit_mod.o
field.o: orbit_mod.o
inner.o: field.o restart.o
outer.o: field.o restart.o
inner_gc.o: field.o restart.o
outer_gc.o: field.o restart.o
orbit.o: inner.o outer.o inner_gc.o outer_gc.o

clean:
	rm -f *.o *.mod *.out *.bin orbit *.original
