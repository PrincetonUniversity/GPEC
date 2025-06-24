# This is the makefile for ZVODE.

include ../install/DEFAULTS.inc

F90 = $(FC) $(FFLAGS) $(RECURSFLAG) $(LEGACYFLAG) -I../equil

.f.o:
	$(F90) -c $*.f

OBJS = \
	zvode.o

libzvode: $(OBJS)
	mkdir -p ../lib
	ar -r ../lib/libzvode.a *.o

clean:
	rm -f *.o *.a *~
