# FC, LIBDIR and FCLIBS are set in the platform-specific included file
include ../install/make.inc.${GPEC_PLATFORM}

FFLAGS = -I../equil -ffree-line-length-0

F90 = $(FC) $(FFLAGS)

.f.o:
	$(F90) -c $*.f

%.o: %.f90
	$(F90) -c $<

LIBS = \
	$(LIBDIR) \
	-llsode \
	-lequil \
	$(FCLIBS)

OBJS = \
	dcon_interface.o \
	utilities.o \
    params.o \
	lsode1.o \
	lsode2.o \
	special.o \
    grid.o \
    inputs.o \
	energy.o \
    pitch.o \
    torque.o \
	diagnostics.o \
	pentrc_interface.o \
	pentrc.o

pentrc: $(OBJS)
	$(F90) -o pentrc $(OBJS) $(LIBS)
	mkdir -p ../lib
	ar -r ../lib/libpentrc.a *.o

# dependencies

lsode1.o:
lsode2.o:
params.o:
grid.o: 
utilities.o: params.o
special.o: params.o utilities.o
dcon_interface.o: params.o utilities.o
inputs.o:  params.o utilities.o special.o dcon_interface.o
energy.o:  params.o utilities.o lsode2.o
pitch.o:   params.o utilities.o lsode1.o energy.o special.o
torque.o:  params.o utilities.o special.o dcon_interface.o inputs.o energy.o pitch.o
diagnostics.o: params.o special.o energy.o pitch.o
pentrc.o:  utilities.o special.o energy.o pitch.o torque.o inputs.o diagnostics.o 
pentrc_interface.o: params.o utilities.o inputs.o energy.o pitch.o torque.o

clean:
	rm -f *.o *.mod *.out *.bin pentrc
