include ../install/DEFAULTS.inc

FFLAGS = -g
IFLAGS = -I../equil -I$(MATHINC)
F90 = $(FC) $(FFLAGS) $(IFLAGS)
export FFLAGS # for sub-makes

.f.o:
	$(F90) -c $*.f

LIBDIR = ../lib

LIBS = \
	-llsode \
	-lequil

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

all: rmatch

equil:
	cd ../equil; make

rmatch: $(OBJS)
	$(F90) -o rmatch $(OBJS) -L$(LIBDIR) $(LIBS) -L$(MATHDIR) $(MATHLIBS) -L$(FCDIR) $(FCLIBS) -Wl,-rpath=$(MATHDIR)
	mkdir -p ../bin
	cp -f rmatch ../bin

# dependencies

msing.o: local.o
gamma.o: local.o
deltar.o: gamma.o
deltac.o: jacobi.o
inner.o: deltac.o
match.o: deltar.o deltac.o inner.o msing.o

clean:
	rm -f *.o *.mod *.out  *~ match
