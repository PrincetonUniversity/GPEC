
FFLAGS = -I ../../equil -g
F90 = ifort $(FFLAGS)

LIBDIR = ../../lib

LAPACKDIR = -L /u/zwang/numericalLib/lapack_3.4.2

LIBS = \
	$(LIBDIR) \
	-lequil \
     $(LAPACKDIR) \
     -llapack \
     -lblas


OBJS = \
	deltag_queue.o \
	deltag_coefficient_manage.o \

all: equil deltag

equil:
	cd ../../equil; make

deltag: $(OBJS)
	$(F90) -c *.f90 
	$(F90) -o deltag $(OBJS) -L$(LIBS)

# dependencies

deltag_queue.o: ../../equil/local.o
deltag_coefficient_manage.o: ../../equil/fspline.o

clean:
	rm -f *.o *.mod *.out  *~ match
