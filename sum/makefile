
include ../install/DEFAULTS.inc

F90 = $(FC) $(FFLAGS) -I../globals 

.f.o:
	$(F90) -c $*.f

OBJS = \
	io.o \
	sum.o

LIBDIR = ../lib

LIBS = \
	-lglobals


sum: $(OBJS)
	$(F90) -o sum $(OBJS) -L$(LIBDIR) $(LIBS) -L$(FCDIR) $(FCLIBS) 
	mkdir -p ../bin
	cp -f sum ../bin

# dependencies

local.o: io.o
utils.o: local.o
sum.o: ../globals/utils.o ../globals/local_mod.mod

clean:
	rm -f *.o *.mod *.out *.bin sum *~ *.original
