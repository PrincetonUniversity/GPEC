# This is the makefile for ZVODE.

# If using gfortran, the proper syntax is -frecursive.
# If using ifort, the proper syntax is -recursive.

ifneq (,$(findstring ifort,$(FC)))
    EXTRAFLAGS=-recursive
else ifneq (,$(findstring ifx,$(FC)))
	EXTRAFLAGS=-recursive
else
    EXTRAFLAGS=-frecursive -std=legacy
endif

include ../install/DEFAULTS.inc

F90 = $(FC) $(FFLAGS) $(EXTRAFLAGS) -I../equil

.f.o:
	$(F90) -c $*.f

OBJS = \
	zvode.o

libzvode: $(OBJS)
	mkdir -p ../lib
	ar -r ../lib/libzvode.a *.o

clean:
	rm -f *.o *.a *~
