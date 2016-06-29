include ../install/DEFAULTS.inc

FFLAGS = -I../equil
F90 = $(FC) $(options) $(FFLAGS)

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

rdcon: $(OBJS)
    $(F90) -o dcon $(OBJS) -L$(LIBDIR) $(LIBS) -L$(MATHDIR) $(MATHLIBS) -L$(FCDIR) $(FCLIBS)
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
