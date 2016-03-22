
include ../install/DEFAULTS.inc

FFLAGS = -I../equil -I../coil -I../harvest -I$(NETCDFHOME)/include

F90 = $(FC) $(options) $(FFLAGS)

.f.o:
	$(F90) -c $*.f

LIBDIR = \
	-L$(MATHDIR) \
	-L$(NETCDFDIR) \
	-L../lib

LIBS = \
	$(LIBDIR) \
	-llsode \
	-lequil \
	-lvac \
	-lcoil \
	-lharvest \
	$(MATHLIBS) \
	$(NETCDFLIBS)

OBJS = \
	gpec_global.o \
	math.o \
	dcon_interface.o \
	rdcon_interface.o \
	peq.o \
	resp.o \
	diagnostic.o \
	vacuum_interface.o \
	output_control.o \
	output_profile.o \
	output_cylindrical.o \
	ipec.o

ipec: $(OBJS)
	$(F90) -o ipec $(OBJS) $(LIBS)
	mkdir -p ../lib
	ar -r ../lib/libipec.a *.o
	mkdir -p ../bin
	cp -f ipec ../bin

# dependencies

ismath.o : ipglobal.o
dcon_interface.o: ismath.o
rdcon_interface.o: ipglobal.o dcon_interface.o ismath.o
peq.o: dcon_interface.o ismath.o
resp.o: peq.o
vacuum_interface.o : peq.o
diagnostic.o: resp.o vacuum_interface.o
output_control.o: resp.o vacuum_interface.o ../coil/field.o ../coil/coil.o
output_profile.o: resp.o vacuum_interface.o ../coil/field.o ../coil/coil.o
output_cylindrical.o: resp.o vacuum_interface.o ../coil/field.o ../coil/coil.o
ipec.o: diagnostic.o output_control.o rdcon_interface.o

clean:
	rm -f *.o *.mod *.out *.bin ipec
