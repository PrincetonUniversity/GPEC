
IFLAGS = -I../equil

F90 = $(FC) $(FFLAGS) $(IFLAGS)

.f.o:
	$(F90) -c $*.f

OBJS = \
	coil.o \
	field.o 

libcoil: $(OBJS)
	mkdir -p ../lib
	ar -r ../lib/libcoil.a *.o

# dependencies

field.o: coil.o

clean:
	rm -f *.o *.mod
