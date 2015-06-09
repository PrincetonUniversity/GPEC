# FC is set in the platform-specific included file
include ../install/make.inc.${GPEC_PLATFORM}

FFLAGS = -I../equil -g

F90 = $(FC) $(FFLAGS)

.f.o:
	$(F90) -c $*.f

OBJS = \
	coil.o \
	field.o 

libcoil: $(OBJS)
	mkdir -p ../lib
	ar -r ../lib/libcoil.a *.o

# dependencies

field.o: coil.o

clean:
	rm -f *.o *.mod
