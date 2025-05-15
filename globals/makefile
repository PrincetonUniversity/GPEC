# If using gfortran, the proper syntax is -frecursive.
# If using ifort, the proper syntax is -recursive.

ifneq (,$(findstring ifort,$(FC)))
    RECURSFLAG=-recursive
else ifneq (,$(findstring ifx,$(FC)))
	RECURSFLAG=-recursive
else
    RECURSFLAG=-frecursive
endif

F90 = $(FC) $(FFLAGS) $(RECURSFLAG)

.f.o:
	$(F90) -c $*.f

%.o: %.f90
	$(F90) -c $<

OBJS = \
	io.o \
	machinelocal.o \
	utils.o \
	spline.o \
	bicube.o \
	cspline.o \
	fft.o \
	fspline.o \
	grid.o \

# targets

all : libequil

libglobals: $(OBJS)
	mkdir -p ../lib/
	ar -r ../lib/libglobals.a *.o

clean:
	rm -f *.o *.mod *.original

# dependencies

machinelocal.o : io.o
