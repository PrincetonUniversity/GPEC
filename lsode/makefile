# FC is set in the platform-specific included file
include ../install/make.inc.${GPEC_PLATFORM}

FFLAGS = -g
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
	mkdir -p ../lib/LinuxLahey64
	ar -r ../lib/LinuxLahey64/liblsode.a *.o

clean:
	rm -f *.o *.a
