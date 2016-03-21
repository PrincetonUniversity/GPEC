
include ../install/DEFAULTS.inc

FFLAGS = -I../equil
F90 = $(FC) $(options) $(FFLAGS)

.f.o:
	$(F90) -c $*.f

LIBDIR = \
	-L$(MATHDIR) \
	-L../lib

LIBS = \
	$(LIBDIR) \
	-llsode \
	-lequil \
	$(MATHLIBS)
	
OBJS = \
	orbit_mod.o \
	restart.o \
	field.o \
	inner.o \
	outer.o \
	inner_gc.o \
	outer_gc.o \
	orbit.o

all: v orbit

orbit: $(OBJS)
	$(F90) -o orbit $(OBJS) $(LIBS)
	mkdir -p ../bin
	cp -f orbit ../bin

v:
    @echo "cflag  = $(cflag)"


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
	rm -f *.o *.mod *.out *.bin orbit
