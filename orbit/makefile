# FC, LIBDIR and FCLIBS are set in the platform-specific included file
include ../install/make.inc.${GPEC_PLATFORM}

FFLAGS = -I../equil -g
F90 = $(FC) $(FFLAGS)

.f.o:
	$(F90) -c $*.f

LIBS = \
	$(LIBDIR) \
	-llsode \
	-lequil \
	$(FCLIBS)

OBJS = \
	orbit_mod.o \
	restart.o \
	field.o \
	inner.o \
	outer.o \
	inner_gc.o \
	outer_gc.o \
	orbit.o

orbit: $(OBJS)
	$(F90) -o orbit $(OBJS) $(LIBS)

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
