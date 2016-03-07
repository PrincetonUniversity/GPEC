
FFLAGS = -I../equil -I../coil -I../harvest -I$(NETCDFHOME)/include

F90 = $(FC) $(FFLAGS)

.f.o:
	$(F90) -c $*.f

LIBDIR = \
	-L$(LAPACKHOME) \
	-L$(NETCDFHOME)/lib \
	-L../lib

LIBS = \
	$(LIBDIR) \
	-llsode \
	-lequil \
	-lvac \
	-lcoil \
	-lharvest \
	-llapack \
	-lblas \
	-lnetcdf
#	-lnetcdff

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
	@echo $(PLATFORM)
	$(F90) -o ipec $(OBJS) $(LIBS)
	mkdir -p ../lib
	ar -r ../lib/libipec.a *.o

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
