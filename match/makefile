
include ../install/DEFAULTS.inc

IFLAGS = -I../equil
F90 = $(FC) $(FFLAGS) $(OMPFLAG) $(IFLAGS) $(RECURSFLAG)
export FFLAGS # for sub-makes


.f.o:
	$(F90) -c $*.f

LIBDIR = ../lib

LIBS = \
	-lequil \
	-llsode \

OBJS = \
	qsort.o \
	gamma.o \
	deltar.o \
	root.o \
	match_mod.o \
	ideal.o \
	debug.o \
	matrix.o \
	resist.o \
	match.o

all: lsode equil

lsode:
	cd ../lsode; make

equil:
	cd ../equil; make

match: $(OBJS)
	$(F90) -o match $(OBJS) -L$(LIBDIR) $(LIBS) -L$(MATHDIR) $(MATHLIBS) -L$(FCDIR) $(FCLIBS) $(LDFLAGS)
	mkdir -p ../bin
	cp -f match ../bin

# dependencies

debug.o: ../equil/local_mod.mod
gamma.o: ../equil/local_mod.mod
deltar.o: gamma.o
root.o: ../equil/local_mod.mod
match_mod.o: deltar.o root.o ../equil/utils_mod.mod
ideal.o: match_mod.o
matrix.o: match_mod.o debug.o
resist.o: matrix.o
match.o: ideal.o resist.o

clean:
	rm -f *.o *.mod *.out *.bin match *.original
