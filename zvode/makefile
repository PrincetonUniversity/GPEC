# This is the makefile for ZVODE.

#FFLAGS = -qopt-report -xHost -g -check all -qopenmp
FFLAGS = -xHost -O3 -qopenmp #-qopt-zmm-usage=high #-ipo
F90 = ifort $(FFLAGS)

.f.o:
	$(F90) -c $*.f

OBJS = \
	zvode.o

libzvode: $(OBJS)
	mkdir -p ../lib
	ar -r ../lib/libzvode.a *.o

clean:
	rm -f *.o *.a *~
