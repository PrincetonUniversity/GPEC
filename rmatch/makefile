# This is the makefile for DELTAR on LinuxIntel computers.

FFLAGS = -I ../equil -m64 -g
F90 = ifort $(FFLAGS)

.f.o:
	$(F90) -c $*.f

LIBDIR = \
	-L../lib/LinuxIntel \
	-L /usr/local/lib \
        -L$(HOME)/lib \

LAPACKDIR = -L /u/zwang/numericalLib/lapack_3.4.2

LIBS = \
	$(LIBDIR) \
	-llsode \
        $(LAPACKDIR) \
	-llapack \
        -lequil \
	-lblas \

OBJS = \
	io.o \
	local.o \
	gamma.o \
	deltar.o \
	inner.o  \
	jacobi.o \
	deltac.o \
	msing.o \
	match.o

all: equil match

equil:
	cd ../equil; make

match: $(OBJS)
	$(F90) -o match $(OBJS) $(LIBS)

# dependencies

msing.o: local.o
gamma.o: local.o
deltar.o: gamma.o
deltac.o: jacobi.o
inner.o: deltac.o
match.o: deltar.o deltac.o inner.o msing.o

clean:
	rm -f *.o *.mod *.out  *~ match
