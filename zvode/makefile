# This is the makefile for ZVODE.

# If using gfortran, the proper syntax is -frecursive.
# If using ifort, the proper syntax is -recursive.

ifneq (,$(findstring ifort,$(FC)))
    RECURSFLAG=-recursive
else ifneq (,$(findstring ifx,$(FC)))
	RECURSFLAG=-recursive
else
    RECURSFLAG=-frecursive
endif

include ../install/DEFAULTS.inc

F90 = $(FC) $(FFLAGS) $(RECURSFLAG) -I../globals

.f.o:
	$(F90) -c $*.f

OBJS = \
	zvode.o

libzvode: $(OBJS)
	mkdir -p ../lib
	ar -r ../lib/libzvode.a *.o

clean:
	rm -f *.o *.a *~
