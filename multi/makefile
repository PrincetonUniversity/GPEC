# FC is set in the platform-specific included file
include ../install/make.inc

#FFLAGS = -g
FFLAGS = 
F90 = $(FC) $(FFLAGS)

.f.o:
	$(F90) -c $*.f

OBJS = \
	io.o \
	local.o \
	multi.o

multi: $(OBJS)
	$(F90) -o multi $(OBJS)

# dependencies

local.o: io.o
multi.o: local.o

clean:
	rm -f *.o *.mod *.out *.bin multi *~
