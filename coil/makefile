# FC is set in the platform-specific included file
include ../install/make.inc

# hardcode pacakge directory variable?
#sed -i '0,:#REPLACE-WITH-PATH:{s:#REPLACE-WITH-PATH:'`pwd`':}' coil.f

FFLAGS = -I../equil

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
