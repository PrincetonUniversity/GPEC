
FFLAGS = -I ../../equil -g
F90 = ifort $(FFLAGS)

LIBDIR = ../../lib

LIBS = \
	$(LIBDIR) \
	-lequil

OBJS = \
	deltag_queue.o \

all: equil deltag

equil:
	cd ../../equil; make

deltag: $(OBJS)
	$(F90) -c *.f90 
	$(F90) -o deltag $(OBJS) -L$(LIBS)

# dependencies

deltag_queue.o: ../../equil/local.o

clean:
	rm -f *.o *.mod *.out  *~ match
