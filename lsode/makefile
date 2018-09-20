
include ../install/DEFAULTS.inc

FFLAGS =  -qopenmp
F90 = $(FC) $(FFLAGS)

.f.o: 
	$(F90) -c $*.f

OBJS = \
	daxpy.o \
	dcfode.o \
	dcopy.o \
	ddot.o \
	dewset.o \
	dgbfa.o \
	dgbsl.o \
	dgefa.o \
	dgesl.o \
	dintdy.o \
	dprepj.o \
	dscal.o \
	dsolsy.o \
	dsrcom.o \
	dstode.o \
	dumach.o \
	dvnorm.o \
	idamax.o \
	iumach.o \
	ixsav.o \
	lsode.o \
	xerrwd.o \
	xsetf.o \
	xsetun.o

liblsode: $(OBJS)
	mkdir -p ../lib
	ar -r ../lib/liblsode.a *.o

clean:
	rm -f *.o *.a
