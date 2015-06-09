# FC is set in the platform-specific included file
include ../install/make.inc.${GPEC_PLATFORM}

FFLAGS = -g
F90 = $(FC) $(FFLAGS)

.f.o:
	$(F90) -c $*.f

OBJS = \
	io.o \
	local.o \
	utils.o \
	sum.o

sum: $(OBJS)
	$(F90) -o sum $(OBJS)

# dependencies

local.o: io.o
utils.o: local.o
sum.o: utils.o

clean:
	rm -f *.o *.mod *.out *.bin sum *~
