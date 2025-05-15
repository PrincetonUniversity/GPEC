# If using gfortran, the proper syntax is -frecursive.
# If using ifort, the proper syntax is -recursive.

ifneq (,$(findstring ifort,$(FC)))
    RECURSFLAG=-recursive
else ifneq (,$(findstring ifx,$(FC)))
	RECURSFLAG=-recursive
else
    RECURSFLAG=-frecursive
endif

F90 = $(FC) $(FFLAGS) $(RECURSFLAG) -I../globals

.f.o:
	$(F90) -c $*.f

OBJS = \
	global.o \
	inverse.o \
	direct.o \
	read_eq.o \
	lar.o \
	gsec.o \
	sol.o \
	equil.o \
	equil_out.o

# targets

all : libequil

libequil: $(OBJS)
	mkdir -p ../lib/
	ar -r ../lib/libequil.a *.o

clean:
	rm -f *.o *.mod *.original

# dependencies

global.o: ../globals/bicube.o ../globals/cspline.o
inverse.o: global.o ../globals/utils.o
direct.o: global.o
read_eq.o: inverse.o direct.o
lar.o: inverse.o
gsec.o: inverse.o
sol.o: direct.o
equil.o: read_eq.o lar.o gsec.o sol.o
equil_out.o: global.o
