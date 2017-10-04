#------------------------------------------------------------------
# General GNU make for portal, venus, or OSX platforms
# - Conditionals used to determine platform
# - This defines platform specific defaults for FC,LAPACKHOME,NETCDFHOME
#
# On portal.pppl.gov (or portalr6.pppl.gov for more up-to-date libraries)
#   Must load compiler and basic libraries (These modify LD_LIBRARY_PATH).
#   Choose one of the following
#    $ module load gcc4 lapack-blas netcdf
#    $ module load lff95 lapack-blas netcdf
#    $ module load pgf netcdf acml/5.3.1/pgi64
#    $ module load intel netcdf acml/5.3.1/ifort64
#   Defaults will try to use gcc4
#
# On venus.gat.com
# - All paths are hardcoded at this point (08/2015).
#   -- We are waiting on a upgrade that will enable "module load" standardization
#
# On mac OSX
#   Prerequisites
#   - Install gfortran in your path.
#      -- Available from macports or sourceforge
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
else
    MATHHOME = $(ACML_HOME)
endif


# platform dependent defaults
#-----------------------------

# force default if necessary compiler/libraries not specified
ifneq ($(shell echo $(HOSTNAME) | grep sunfire),)
    PLATFORM = portal
    ifeq "$(or $(GPEC_PLATFORM), $(and $(FC), $(MATHHOME), $(NETCDFHOME)))" ""
        GPEC_PLATFORM = portal_intel
    endif
else ifneq ($(shell echo $(HOSTNAME) | grep venus),)
    PLATFORM = venus
    ifeq "$(or $(GPEC_PLATFORM), $(and $(FC), $(MATHHOME), $(NETCDFHOME)))" ""
        GPEC_PLATFORM = venus_mpi
    endif
else ifneq ($(shell echo $(MACHTYPE) | grep apple),)
    PLATFORM = osx
    ifeq "$(or $(GPEC_PLATFORM), $(and $(FC), $(MATHHOME), $(NETCDFHOME)))" ""
        GPEC_PLATFORM = osx_gfortran
    endif
endif


# set compiler and libraries based on GPEC_PLATFORM
#-----------------------------------------------

ifdef GPEC_PLATFORM
    export GPEC_PLATFORM # extends it to recursive makes
	ifeq ($(GPEC_PLATFORM),portal_intel)
	    # $ module load intel acml netcdf
        # make options=-vec_report0 2>&1 | tee log.txt
	# make options=-traceback -p for debugging memory leak
        # status: Success 06/14/2016, solovev_ideal_example DCON 9s, GPEC 3m 1s.
		PLATFORM = portal
		FC = /usr/pppl/intel/11.0/081/bin/intel64/ifort -O2		
		CC = /usr/pppl/intel/11.0/081/bin/intel64/icc -O2
		MATHHOME = /usr/pppl/intel/11-pkgs/acml-4.3.0/ifort64
		NETCDFHOME = /usr/pppl/intel/11-pkgs/netcdf-3.6.3
    else ifeq ($(GPEC_PLATFORM),portal_gfortran)
	    # $ module load gcc4 acml netcdf
        # status: Success 06/14/2016. solovev_ideal_example DCON 15s, GPEC 14m 32s.
		PLATFORM = portal
		FC = /usr/pppl/gcc/4.6.3/bin/gfortran4
		CC = /usr/pppl/gcc/4.6.3/bin/gcc4
		MATHHOME = /usr/pppl/gcc/4.6-pkgs/lapack-3.3.1
		NETCDFHOME = /usr/pppl/gcc/4.6-pkgs/netcdf-4.1.2
    else ifeq ($(GPEC_PLATFORM),portal_pgf)
	    # $ module load pgf acml netcdf; export ACML_HOME=$ACML_HOME/pgi64
        # $ make GPEC_PLATFORMS="'-tp px-64'" options="'-tp px-64'" 2>&1 | tee log.txt
        # status: Success 06/17/2016. solovev_ideal_example DCON 29s, GPEC 12m 41s.
		PLATFORM = portal
		FC = /usr/pppl/pgi/11.10/linux86-64/11.10/bin/pgf90
		CC = /usr/pppl/pgi/11.10/linux86-64/11.10/bin/pgcc
		MATHHOME = /usr/pppl/pgi/11-pkgs/acml-5.1.0/pgi64
		NETCDFHOME = /usr/pppl/pgi/11-pkgs/netcdf-4.1.1
	else ifeq ($(GPEC_PLATFORM),portal_lfc)
	    # $ module load lff95/8.00a_64bit lapack-blas/3.1.1 netcdf/3.6.2
		# status: Success 06/17/2016. solovev_ideal_example DCON 15s, GPEC 15m 36s.
		PLATFORM = portal
		FC = /usr/pppl/lf64x/lf6480/bin/lfc
		CC = /usr/pppl/gcc/4.6.3/bin/gcc4
		MATHHOME = /usr/pppl/lff95/8.x86_64-pkgs/lapack-3.1.1.x86_64
		NETCDFHOME = /usr/pppl/lff95/8.x86_64-pkgs/netcdf-3.6.2.x86_64
    else ifeq ($(GPEC_PLATFORM),portalr6_intel)
	    # $ module load intel acml netcdf
        # status: Success 06/14/2016. solovev_ideal_example DCON 15s, GPEC segfault in xclebsch.
		PLATFORM = portal
		FC = /usr/pppl/intel/2015.u1/bin/ifort
		CC = /usr/pppl/intel/2015.u1/bin/icc
		MATHHOME = /usr/pppl/intel/2015-pkgs/acml-5.3.1/ifort64
		NETCDFHOME = /usr/pppl/intel/2015-pkgs/netcdf-4.3.2
	else ifeq ($(GPEC_PLATFORM),zwang_intel)
	    # $ module load intel/2013.1.3 lapack netcdf
	    # status: Success.
		PLATFORM = portal
		FC = /usr/pppl/intel/2013.1.3/bin/ifort
		CC = /usr/pppl/intel/2013.1.3/bin/icc
		MATHHOME = /u/zwang/numericalLib/lapack_3.4.2
		NETCDFHOME = /usr/pppl/intel/11-pkgs/netcdf-3.6.3
	else ifeq ($(GPEC_PLATFORM),venus_mpi)
        FC = /task/imd/local64/bin/mpif90
        MATHHOME = /opt/pgi/linux86-64/12.10/lib #pgf
        NETCDFHOME = /task/imd/local64/lib #mpif90
	else ifeq ($(GPEC_PLATFORM),osx_gfortran)
        FC = gfortran
        MATHHOME = /usr/lib
        NETCDFHOME = /opt/local/lib
	else
        $(error No corresponding option of GPEC_PLATFORM $(GPEC_PLATFORM) to set compile environment.)
	endif
endif


# specifics of which libraries are used
#---------------------------------------

# tack on lib dir to home dirs if necessary
ifneq ($(wildcard $(MATHHOME)/lib/.),)
    MATHDIR = $(MATHHOME)/lib
else
    MATHDIR = $(MATHHOME)
endif
ifneq ($(wildcard $(NETCDFHOME)/lib/.),)
    NETCDFDIR = $(NETCDFHOME)/lib
else
    NETCDFDIR = $(NETCDFHOME)
endif

# tack on include dir to home dirs if necessary
ifneq ($(wildcard $(MATHHOME)/include/.),)
    MATHINC = $(MATHHOME)/include
else
    MATHINC = $(MATHHOME)
endif
ifneq ($(wildcard $(NETCDFHOME)/include/.),)
    NETCDFINC = $(NETCDFHOME)/include
else
    NETCDFINC = $(NETCDFHOME)
endif

# include compiler specific libs for intel fast moves
ifneq ($(wildcard $(F90HOME)/lib/intel64/.),)
    FCDIR = $(F90HOME)/lib/intel64
    #FCLIBS = -lirc -limf # outdated
    FCLIBS = -lintlc -limf
    FCINC = $(F90HOME)/include/intel64
else
    FCDIR = .
    FCLIBS = 
    FCINC = .
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
