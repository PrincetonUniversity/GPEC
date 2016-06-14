
include ../install/DEFAULTS.inc

FFLAGS = -I../equil -I../harvest -I$(MATHINC) -I$(NETCDFINC)

F90 = $(FC) $(options) $(FFLAGS)

.f.o:
	$(F90) -c $*.f

%.o: %.f90
	$(F90) -c $<

LIBDIR = ../lib

LIBS = \
	-llsode \
	-lequil \
	-lharvest

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
	$(F90) -o pentrc $(OBJS) -L$(FCDIR) $(FCLIBS) -L$(MATHDIR) $(MATHLIBS) -L$(NETCDFDIR) $(NETCDFLIBS) -L$(LIBDIR) $(LIBS)
	mkdir -p ../lib
	ar -r ../lib/libpentrc.a *.o
	mkdir -p ../bin
	cp -f pentrc ../bin

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
pentrc_interface.o: params.o utilities.o special.o inputs.o energy.o pitch.o torque.o diagnostics.o
pentrc.o:  params.o pentrc_interface.o

clean:
	rm -f *.o *.mod *.out *.bin pentrc
