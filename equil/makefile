include ../install/DEFAULTS.inc

F90 = $(FC) $(FFLAGS) $(RECURSFLAG)

.f.o:
	$(F90) -c $*.f

%.o: %.f90
	$(F90) -c $<


OBJS = \
	io.o \
	local.o \
	utils.o \
	spline.o \
	grid.o \
	bicube.o \
	cspline.o \
	fft.o \
	fspline.o \
	global.o \
	inverse.o \
	direct.o \
	read_eq.o \
	lar.o \
	gsec.o \
	sol.o \
	equil.o \
	equil_out.o \
	jacobi.o

# targets

all : libequil

libequil: $(OBJS)
	mkdir -p ../lib/
	ar -r ../lib/libequil.a *.o

clean:
	rm -f *.o *.mod *.original

# dependencies

local.o: io.o
utils.o: local.o
spline.o: local.o
grid.o: local.o spline.o
bicube.o: spline.o
cspline.o: local.o
fft.o: local.o
fspline.o: local.o spline.o cspline.o utils.o fft.o
global.o: bicube.o cspline.o
inverse.o: global.o utils.o grid.o
direct.o: global.o
read_eq.o: inverse.o direct.o
lar.o: inverse.o
gsec.o: inverse.o
sol.o: direct.o
equil.o: read_eq.o lar.o gsec.o sol.o
equil_out.o: global.o
