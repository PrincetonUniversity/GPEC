include ../install/DEFAULTS.inc

IFLAGS = -I../globals -I../equil -I$(MATHINC)
F90 = $(FC) $(FFLAGS) $(OMPFLAG) $(IFLAGS)
export FFLAGS # for sub-makes

.f.o:
	$(F90) -c $*.f

LIBDIR = ../lib

LIBS = \
	-llsode \
	-lequil \
	-lglobals \

OBJS = \
	gamma.o \
	deltar.o \
	inner.o  \
	jacobi.o \
        inps.o \
        inpso.o \
	deltac.o \
	msing.o \
	match.o

all: rmatch

equil:
	cd ../equil; make

rmatch: $(OBJS)
	$(F90) -o rmatch $(OBJS) -L$(LIBDIR) $(LIBS) -L$(MATHDIR) $(MATHLIBS) -L$(FCDIR) $(FCLIBS) $(LDFLAGS)
	mkdir -p ../bin
	cp -f rmatch ../bin

# dependencies

msing.o: ../globals/local_mod.mod
gamma.o: ../globals/local_mod.mod
deltar.o: gamma.o
inps.o:  ../globals/local_mod.mod ../globals/utils.o
inpso.o: inps.o
deltac.o: jacobi.o inpso.o
inner.o: deltac.o ../globals/local_mod.mod
match.o: deltar.o deltac.o inner.o msing.o

clean:
	rm -f *.o *.mod *.out  *~ match *.original
