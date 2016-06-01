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
    MATHDIR = $(LAPACKHOME)
else
    MATHDIR = $(ACML_HOME)
endif


# platform dependent defaults
#-----------------------------

# force default if necessary compiler/libraries not specified
ifneq ($(shell echo $(HOSTNAME) | grep sunfire),)
    PLATFORM = portal
    ifeq "$(or $(cflag), $(and $(FC), $(MATHDIR), $(NETCDFHOME)))" ""
        cflag = portal_lfc
    endif
else ifneq ($(shell echo $(HOSTNAME) | grep venus),)
    PLATFORM = venus
    ifeq "$(or $(cflag), $(and $(FC), $(MATHDIR), $(NETCDFHOME)))" ""
        cflag = venus_mpi
    endif
else ifneq ($(shell echo $(MACHTYPE) | grep apple),)
    PLATFORM = osx
    ifeq "$(or $(cflag), $(and $(FC), $(MATHDIR), $(NETCDFHOME)))" ""
        cflag = osx_gfortran
    endif
endif


# set compiler and libraries based on cflag
#-----------------------------------------------

ifdef cflag
	ifeq ($(cflag),zwang_intel)
	    # module load intel/2013.1.3 lapack netcdf
	    # status: success!
		PLATFORM = portal
		FC = /usr/pppl/intel/2013.1.3/bin/ifort
		CC = /usr/pppl/intel/2013.1.3/bin/icc
		MATHDIR = /u/zwang/numericalLib/lapack_3.4.2
		NETCDFHOME = /usr/pppl/intel/11-pkgs/netcdf-3.6.3
	else ifeq ($(cflag),portal_intel) # currently fails
	    #  module load intel acml netcdf
		PLATFORM = portal
		FC = /usr/pppl/intel/11.0/081/bin/intel64/ifort
		CC = /usr/pppl/intel/11.0/081/bin/intel64/icc
		MATHDIR = /usr/pppl/intel/11-pkgs/acml-4.3.0/ifort64
		NETCDFHOME = /usr/pppl/intel/11-pkgs/netcdf-3.6.3
	else ifeq ($(cflag),portal_intel13)
	    #  module load intel/2013.1.3 acml netcdf
	    # status: fails to compile with "for_init.c:(.text+0x8ce): undefined reference to `_intel_fast_memmove'"
		PLATFORM = portal
		FC = /usr/pppl/intel/2013.1.3/bin/ifort
		CC = /usr/pppl/intel/2013.1.3/bin/icc
		MATHDIR = /usr/pppl/intel/2013-pkgs/acml-5.3.1/ifort64
		#MATHDIR = /usr/pppl/lapack/3.3.1
		NETCDFHOME = /usr/pppl/intel/2013-pkgs/netcdf-4.3.2
		#NETCDFHOME = /usr/pppl/intel/11-pkgs/netcdf-3.6.3
	else ifeq ($(cflag),portalr6_intel) # currently fails
	    # module load intel/2015.u1 acml netcdf
		PLATFORM = portal
		FC = /usr/pppl/intel/2015.u1/bin/ifort
		CC = /usr/pppl/intel/2015.u1/bin/icc
		MATHDIR = /usr/pppl/intel/2015-pkgs/acml-5.3.1/ifort64
		NETCDFHOME = /usr/pppl/intel/2015-pkgs/netcdf-4.3.2
	else ifeq ($(cflag),portalr6_gfortran)
	    # module load gcc4/4.6.3; module load lapack-blas/3.3.1 netcdf/4.1.2
	    # status: success! (executor must have gcc4 modules loaded to use libs)
		PLATFORM = portal
		FC = /usr/pppl/gcc/4.6.3/bin/gfortran4
		CC = /usr/pppl/gcc/4.6.3/bin/gcc4
		MATHDIR = /usr/pppl/gcc/4.6-pkgs/lapack-3.3.1
		NETCDFHOME = /usr/pppl/gcc/4.6-pkgs/netcdf-4.1.2
	else ifeq ($(cflag),portal_lfc)
	    # module load lff95/8.00a_64bit lapack-blas/3.1.1 netcdf/3.6.2
		# status: success! (executor does not need same environment)
		PLATFORM = portal
		FC = /usr/pppl/lf64x/lf6480/bin/lfc
		CC = /usr/pppl/gcc/4.6.3/bin/gcc4
		MATHDIR = /usr/pppl/lff95/8.x86_64-pkgs/lapack-3.1.1.x86_64
		NETCDFHOME = /usr/pppl/lff95/8.x86_64-pkgs/netcdf-3.6.2.x86_64
	else ifeq ($(cflag),venus_mpi)
        FC = /task/imd/local64/bin/mpif90
        MATHDIR = /opt/pgi/linux86-64/12.10/lib #pgf
        NETCDFHOME = /task/imd/local64/lib #mpif90
	else ifeq ($(cflag),osx_gfortran)
        FC = gfortran
        MATHDIR = /usr/lib
        NETCDFHOME = /opt/local/lib
	else
#$(error No corresponding option of cflag $(cflag) to set compile environment.)
	endif
endif


# specifics of which libraries are used
#---------------------------------------

# tack on lib dir to home dirs if necessary
ifneq ($(wildcard $(MATHDIR)/lib/.),)
    MATHDIR := $(MATHDIR)/lib
endif
NETCDFDIR = $(NETCDFHOME)
ifneq ($(wildcard $(NETCDFHOME)/lib/.),)
    NETCDFDIR = $(NETCDFHOME)/lib
else
    NETCDFDIR = $(NETCDFHOME)
endif

# check if math dir has lapack or acml library
ifneq ($(shell ls $(MATHDIR) | grep liblapack),)
    MATHLIBS = -llapack -lblas
else
    MATHLIBS = -lacml
endif

# check if netcdf dir has netcdff library
ifneq ($(shell ls $(NETCDFDIR) | grep netcdff),)
    NETCDFLIBS = -lnetcdf -lnetcdff
else
    NETCDFLIBS = -lnetcdf
endif
