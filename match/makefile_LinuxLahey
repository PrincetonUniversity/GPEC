
FFLAGS = -I../equil
F90 = $(FC) $(FFLAGS)

.f.o:
	$(F90) -c $*.f

LIBDIR = \
	-L${LAPACKHOME} \
	-L../lib

LIBS = \
	$(LIBDIR) \
	-lequil \
	-llsode \
	-llapack \
	-lblas

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

match: $(OBJS)
	$(F90) -o match $(OBJS) $(LIBS)

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
	rm -f *.o *.mod *.out *.bin match
