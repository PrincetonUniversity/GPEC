# FC, LIBDIR and FCLIBS are set in the platform-specific included file
include ../install/make.inc.${GPEC_PLATFORM}

FFLAGS = -I../equil -I../coil

F90 = $(FC) $(FFLAGS)

.f.o:
	$(F90) -c $*.f

LIBS = \
	$(LIBDIR) \
	-llsode \
	-lequil \
	-lvac \
	-lcoil \
	$(FCLIBS)

OBJS = \
	ipglobal.o \
	ismath.o \
	idcon.o \
	rdcon.o \
	ipeq.o \
	ipresp.o \
	ipdiag.o \
	ipvacuum.o \
	ipout.o \
	ipec.o

ipec: $(OBJS)
	$(F90) -o ipec $(OBJS) $(LIBS)
	mkdir -p ../lib/LinuxLahey64
	ar -r ../lib/LinuxLahey64/libipec.a *.o

# dependencies

ismath.o : ipglobal.o
idcon.o: ismath.o
rdcon.o: ipglobal.o idcon.o ismath.o
ipeq.o: idcon.o ismath.o
ipresp.o: ipeq.o
ipvacuum.o : ipeq.o
ipdiag.o: ipresp.o ipvacuum.o
ipout.o: ipresp.o ipvacuum.o ../coil/field.o ../coil/coil.o
ipec.o: ipdiag.o ipout.o rdcon.o

clean:
	rm -f *.o *.mod *.out *.bin ipec
