
include ../install/DEFAULTS.inc

FFLAGS = 
F90 = $(FC) $(options) $(FFLAGS)

.f.o:
	$(F90) -c $*.f

OBJS = \
	io.o \
	local.o \
	utils.o \
	sum.o

sum: $(OBJS)
	$(F90) -o sum $(OBJS)  -L$(FCDIR) $(FCLIBS)
	mkdir -p ../bin
	cp -f sum ../bin

# dependencies

local.o: io.o
utils.o: local.o
sum.o: utils.o

clean:
	rm -f *.o *.mod *.out *.bin sum *~
