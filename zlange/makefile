# This is the makefile for ZLANGE.

include ../install/DEFAULTS.inc

F90 = $(FC) $(FFLAGS) -I../globals

.f.o:
	$(F90) -c $*.f

OBJS = \
	zlange.o

libzlange: $(OBJS)
	mkdir -p ../lib
	ar -r ../lib/libzlange.a *.o

clean:
	rm -f *.o *.a *~
