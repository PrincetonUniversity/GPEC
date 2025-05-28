# This is the makefile for VACUUM.

# If using gfortran, the proper syntax is -finit-local-zero.
# If using ifort, the proper syntax is -zero.

ifneq (,$(findstring ifort,$(FC)))
    ZEROFLAG=-zero
else
    ZEROFLAG=-finit-local-zero
endif

F90 = $(FC) $(FFLAGS) $(ZEROFLAG) -I../equil

.f.o:
	$(F90) -c $*.f

OBJS = \
	bzio.o \
	aleg_new.o \
	vacuum_global.o \
	vacuum_io.o \
	vacuum_ma.o \
	vacuum_math.o \
	vacuum_penn.o \
	vacuum_sprk2.o \
	vacuum_ut.o \
	vacuum_vac.o

libvac: $(OBJS)
	mkdir -p ../lib
	ar -r ../lib/libvac.a *.o

vacuum_diag.o: vacuum_global.o
vacuum_io.o: vacuum_global.o
vacuum_ma.o: vacuum_global.o
vacuum_math.o: vacuum_global.o
vacuum_sprk2.o: vacuum_global.o
vacuum_ut.o: vacuum_global.o
vacuum_vac.o: vacuum_global.o

clean:
	rm -f *.o *.mod *.original
