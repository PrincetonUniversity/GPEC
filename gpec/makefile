
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
	-lpentrc \
	$(MATHLIBS) \
	$(NETCDFLIBS)

OBJS = \
	gpec_global.o \
	gpec_math.o \
	dcon_interface.o \
	rdcon_interface.o \
	peq.o \
	response.o \
	vacuum_interface.o \
	diagnostics.o \
	output_control.o \
	output_profile.o \
	output_cylindrical.o \
	gpec.o

gpec: $(OBJS)
	$(F90) -o gpec $(OBJS) $(LIBS)
	mkdir -p ../lib
	ar -r ../lib/libgpec.a *.o
	mkdir -p ../bin
	cp -f gpec ../bin

# dependencies

gpec_math.o : gpec_global.o
dcon_interface.o: gpec_global.o
rdcon_interface.o: gpec_global.o
peq.o: gpec_global.o gpec_math.o dcon_interface.o
response.o: gpec_global.o peq.o dcon_interface.o
vacuum_interface.o : gpec_global.o gpec_math.o peq.o
diagnostic.o: gpec_global.o gpec_math.o peq.o dcon_interface.o vacuum_interface.o
output_control.o: gpec_global.o gpec_math.o peq.o dcon_interface.o vacuum_interface.o ../coil/field.o
output_profile.o: gpec_global.o gpec_math.o peq.o dcon_interface.o ../coil/field.o ../coil/coil.o
output_cylindrical.o: gpec_global.o gpec_math.o peq.o dcon_interface.o vacuum_interface.o diagnostic.o ../coil/field.o
gpec.o: gpec_global.o response.o peq.o dcon_interface.o rdcon_interface.o diagnostics.o output_control.o output_cylindrical.o output_profile.o

clean:
	rm -f *.o *.mod *.out *.bin gpec
