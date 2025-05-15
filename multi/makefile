
F90 = $(FC) $(FFLAGS) -I../equil

.f.o:
	$(F90) -c $*.f 

OBJS = \
	io.o \
	multi.o

LIBDIR = ../lib

LIBS = \
	-lequil

multi: $(OBJS)
	$(F90) -o multi $(OBJS) -L$(LIBDIR) $(LIBS) 
	mkdir -p ../bin
	cp -f multi ../bin

# dependencies

local.o: io.o
multi.o: ../equil/local_mod.mod

clean:
	rm -f *.o *.mod *.out *.bin multi *~ *.original
