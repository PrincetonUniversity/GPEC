# This is the makefile for ZLANGE.
# If using gfortran, the proper syntax is -frecursive.
# If using ifort, the proper syntax is -recursive.

ifneq (,$(findstring ifort,$(FC)))
    RECURSFLAG=-recursive
    LEGACYFLAG=
else ifneq (,$(findstring ifx,$(FC)))
    RECURSFLAG=-recursive
    LEGACYFLAG=
else
    RECURSFLAG=-frecursive
    LEGACYFLAG=-std=legacy
endif

include ../install/DEFAULTS.inc

F90 = $(FC) $(FFLAGS) $(RECURSFLAG) $(LEGACYFLAG) -I../equil

.f.o:
	$(F90) -c $*.f

OBJS = \
	zlange.o

libzlange: $(OBJS)
	mkdir -p ../lib
	ar -r ../lib/libzlange.a *.o

clean:
	rm -f *.o *.a *~
