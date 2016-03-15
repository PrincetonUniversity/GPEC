#------------------------------------------------------------------
# General GNU make for portal, venus, or OSX platforms
# - Conditionals used to determine platform
# - This defines platform specific defaults for FC,LAPACKHOME,NETCDFHOME
#
# On portal.pppl.gov (or portalr6.pppl.gov for more up-to-date libraries)
#   Must load compiler and basic libraries (These modify LD_LIBRARY_PATH).
#   Choose one of the following
#	$ module load lff95 lapack-blas netcdf
#	$ module load gcc4 lapack-blas netcdf
#	$ module load intel acml netcdf
#	Defaults will try to use gcc4
#
# On venus.gat.com
# - All paths are hardcoded at this point (08/2015).
#   -- We are waiting on a upgrade that will enable "module load" standardization
#
# On mac OSX
#   Prerequisites
#   - Install gfortran in your path.
#	  -- Available from macports or sourceforge
#   - Install lapack, blas, and netcdf libraries in usr/lib.
#
#
# If your platform is not properly included in this make file
# contact nlogan@pppl.gov with your operating system, compiler, and
# library installation details.
#
#------------------------------------------------------------------

# general defaults
#------------------

# use acml or lapack based on what is loaded
ifndef ACML_HOME
    MATHHOME = $(LAPACKHOME)
    MATHLIBS = -llapack -lblas
else
    MATHHOME = $(ACML_HOME)
    MATHLIBS = -lacml
endif

# check if netcdfhome has netcdff library
ifneq ($(shell ls $(NETCDFHOME) | grep netcdff),)
    NETCDFLIBS = -lnetcdf -lnetcdff
else
    NETCDFLIBS = -lnetcdf
endif


# platform dependent defaults
#-----------------------------

# force default if necessary compiler/libraries not specified
ifneq ($(shell echo $(HOSTNAME) | grep sunfire),)
    PLATFORM = portal
    ifeq "$(or $(cflag), $(and $(FC), $(MATHHOME), $(NETCDFHOME)))" ""
        cflag = portal_gfortran
    endif
else ifneq ($(shell echo $(HOSTNAME) | grep venus),)
    PLATFORM = venus
    ifeq "$(or $(cflag), $(and $(FC), $(MATHHOME), $(NETCDFHOME)))" ""
        cflag = venus_mpi
    endif
else ifneq ($(shell echo $(MACHTYPE) | grep apple),)
    PLATFORM = osx
    ifeq "$(or $(cflag), $(and $(FC), $(MATHHOME), $(NETCDFHOME)))" ""
        cflag = osx_gfortran
    endif
endif

# set compiler and libraries based on cflag
#-----------------------------------------------

ifdef cflag
	ifeq ($(cflag),zwang_intel)
		PLATFORM = portal
		FC = /usr/pppl/intel/2013.1.3/bin/ifort
		CC = /usr/pppl/intel/2013.1.3/bin/icc
		MATHHOME = /u/zwang/numericalLib/lapack_3.4.2
		MATHLIBS = -llapack -lblas
		NETCDFHOME = /usr/pppl/intel/11-pkgs/netcdf-3.6.3
		NETCDFLIBS = -lnetcdf -lnetcdff
	else ifeq ($(cflag),portal_intel)
		PLATFORM = portal
		FC = /usr/pppl/intel/11.0/081/bin/intel64/ifort
		MATHHOME = /usr/pppl/intel/11-pkgs/acml-4.3.0/ifort64/lib
		MATHLIBS = -lacml
		NETCDFHOME = /usr/pppl/intel/11-pkgs/netcdf-3.6.3
		NETCDFLIBS = -lnetcdf -lnetcdff
	else ifeq ($(cflag),portalr6_intel)
		PLATFORM = portal
		FC = /usr/pppl/intel/2015.u1/bin/ifort
		MATHHOME = /usr/pppl/intel/2015-pkgs/acml-5.3.1/ifort64/lib
		MATHLIBS = -lacml
		NETCDFHOME = /usr/pppl/intel/2015-pkgs/netcdf-4.3.2
		NETCDFLIBS = -lnetcdf -lnetcdff
	else ifeq ($(cflag),portal_gfortran)
		PLATFORM = portal
		FC = /usr/pppl/gcc/4.6.3/bin/gfortran4
		MATHHOME = /usr/pppl/gcc/4.6-pkgs/lapack-3.3.1
		MATHLIBS = -llapack -lblas
		NETCDFHOME = /usr/pppl/gcc/4.6-pkgs/netcdf-4.1.2
		NETCDFLIBS = -lnetcdf -lnetcdff
	else ifeq ($(cflag),portal_lfc)
		PLATFORM = portal
		FC = /usr/pppl/lf64x/lf6480/bin/lfc
		MATHHOME = /usr/pppl/lff95/8.x86_64-pkgs/lapack-3.1.1.x86_64
		MATHLIBS = -llapack -lblas
		NETCDFHOME = /usr/pppl/lff95/8.x86_64-pkgs/netcdf-3.6.2.x86_64
		NETCDFLIBS = -lnetcdf
	else ifeq ($(cflag),venus_mpi)
        FC = /task/imd/local64/bin/mpif90
        MATHHOME = /opt/pgi/linux86-64/12.10/lib #pgf
        MATHLIBS = -llapack -lblas
        NETCDFHOME = /task/imd/local64/lib #mpif90
        NETCDFLIBS = -lnetcdf
	else ifeq ($(cflag),osx_gfortran)
        FC = gfortran
        MATHHOME = /usr/lib
        MATHLIBS = -llapack -lblas
        NETCDFHOME = /usr/lib
        NETCDFLIBS = -lnetcdf
	else
#$(error No corresponding option of cflag $(cflag) to set compile environment.)
	endif
endif
