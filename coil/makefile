# This is the makefile for COIL

# If using gfortran, the proper syntax is -frecursive.
# If using ifort, the proper syntax is -recursive.

ifneq (,$(findstring ifort,$(FC)))
    RECURSFLAG=-recursive
else ifneq (,$(findstring ifx,$(FC)))
	RECURSFLAG=-recursive
else
    RECURSFLAG=-frecursive
endif

include ../install/DEFAULTS.inc

IFLAGS = -I../equil -I$(MATHINC) -I$(NETCDFINC) -I$(CRLIBMINC)
F90 = $(FC) $(FFLAGS) $(OMPFLAG) $(IFLAGS) $(RECURSFLAG)

.f.o:
	$(F90) -c $*.f

.F.o:
	$(F90) -c $*.F

OBJS = \
	coil.o \
	field.o 

# targets

all: libcoil

# decide if we need to update the version
# note gpec keeps track of its own version, which can differ from the synced dcon & pentrc one
# this enables faster compiling during development
ifeq ($(wildcard version.inc),)
version.inc: force
	@echo ">>> Creating version file"
	@echo "      CHARACTER(len=*), PARAMETER :: version ='"`git describe --tags`"'" | tee version.inc
else ifeq ($(shell grep -F "'"`git describe --tags`"'" version.inc),)
version.inc: force
	@echo ">>> Updating version file"
	@echo "      CHARACTER(len=*), PARAMETER :: version ='"`git describe --tags`"'" | tee version.inc
else
version.inc:
	@echo ">>> Git version unchanged!"
endif

force:

libcoil: $(OBJS)
	mkdir -p ../lib
	ar -r ../lib/libcoil.a *.o

clean:
	rm -f *.o *.mod ../lib/libcoil.a *.original


# dependencies

coil.o: version.inc
field.o: coil.o
