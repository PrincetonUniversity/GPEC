# This is the makefile for VACUUM.

FFLAGS = -O3 -qopenmp -xHost #-qopt-zmm-usage=high #-ipo

F90 = $(FC) $(FFLAGS)

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
	rm -f *.o *.mod
