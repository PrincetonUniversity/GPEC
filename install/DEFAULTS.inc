#------------------------------------------------------------------------------
# This file is used to dynamically set some default compiler options
# based on what libraries are available on the current machine.
#
# Prior to version 1.2, GPEC could be compiled using ifort, pgf, gfortran,
# or lfc. All public releases on portal were compiled using lfc, which
# required only default user path libraries.
#
# With the release of 1.2, the options restricted to more modern compilers
# compatible wth openmp.
# Machine specific instructions are below.
#
# Note that GPEC uses git to fetch/update the harvest repository as well as
# to label the version at compilation time. Thus, you must have a working git
# command in order to compile.
#
# On portal.pppl.gov
#   Users can run GPEC on portal using the public module. Note, portal now "hides" such modules behind a mod_* gateway module
#     $ module load mod_gpec gpec
#   Developers may want to load the development branch module to get the latest compiler settings if they've changed since the last release
#     $ module load mod_gpec gpec/0.0
#   Alternatively, we can manually load the 2018 intel compiler with supporting libraries,
#     $ module load intel/2018.u1 netcdf-fortran/4.4.4 acml/5.3.1/ifort64 git
#   Note that before the netcdf module was properly set up, we had to manually load supporting modules
#     $ module load intel/2018.u1 szip hdf hdf5-parallel/1.10.2 curl netcdf-c/4.6.0 netcdf-fortran/4.4.4 acml/5.3.1/ifort64 git
#   Any of these should properly set up your LD_LIBRARY_PATH and PATH. Executing,
#     $ make v
#   in this install directory should produce the following verbose message,
#
#   ------------------------------------------------------------
#   Defaults tailored for
#   Compiling with FC = /usr/pppl/intel/2018.u1/bin/ifort
#      using math libraries -lacml
#      from /usr/pppl/acml/5.3.1/ifort64/lib
#      including /usr/pppl/acml/5.3.1/ifort64/include
#      and netcdf libraries -lnetcdff -lnetcdf
#      from /usr/pppl/intel/2018-pkgs/netcdf-fortran-4.4.4/lib /usr/pppl/intel/2018-pkgs/netcdf-c-4.6.0/lib
#      including /usr/pppl/intel/2018-pkgs/netcdf-fortran-4.4.4/include
#      and compiler libraries -lintlc -limf
#      from /usr/pppl/intel/2018.u1/lib/intel64
#      including /usr/pppl/intel/2018.u1/include/intel64
#
#
#   Set environment FC,LAPACKHOME,NETCDFHOME to change these
#   ------------------------------------------------------------
#
#   Now compiling using
#    $ make
#   Should set all the right flags and run through.
#
#   Developers can also load the default intel/11.0.081 and associated libraries,
#     $ module load intel netcdf acml/5.3.1/ifort64 git
#   This should properly set up your LD_LIBRARY_PATH and PATH. Executing,
#     $ make v
#   in this install directory should produce the following verbose message.
#
#   ------------------------------------------------------------
#   Defaults tailored for
#   Compiling with FC = /usr/pppl/intel/11.0/081/bin/intel64/ifort
#      using math libraries -lacml
#      from /usr/pppl/acml/5.3.1/ifort64/lib
#      including /usr/pppl/acml/5.3.1/ifort64/include
#      and netcdf libraries -lnetcdf -lnetcdff
#      from /usr/pppl/intel/11-pkgs/netcdf-3.6.3/lib
#      including /usr/pppl/intel/11-pkgs/netcdf-3.6.3/include
#      and compiler libraries -lintlc -limf
#      from /usr/pppl/intel/11.0/081/lib/intel64
#      including /usr/pppl/intel/11.0/081/include/intel64
#
#   Set environment FC,LAPACKHOME,NETCDFHOME to change these
#   ------------------------------------------------------------
#
#   In this case, one should compile using
#    $ make FFLAGS="-openmp -vec_report0" 2>&1 | tee make.log
#
# On ukstar
#   Developers should use the default intel compiler and load the netcdf_intel module
#     $ module load netcdf/netcdf_intel
#     $ export F90HOME=/opt/intel/compilers_and_libraries_2018.1.163/linux
#     $ export FC=/opt/intel/compilers_and_libraries_2018.1.163/linux/bin/intel64/ifort
#   This should properly set up your LD_LIBRARY_PATH and PATH. Executing,
#     $ make v
#   in this install directory should produce the following verbose message.
#
#   ------------------------------------------------------------
#   Defaults tailored for
#   Compiling with FC = /opt/intel/compilers_and_libraries_2018.1.163/linux/bin/intel64/ifort
#      using math libraries -lmkl_rt
#      from /opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64
#      including /opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/include/intel64/lp64
#      and netcdf libraries -lnetcdff -lnetcdf
#      from /usr/local/netcdf/lib /usr/local/netcdf/lib
#      including /usr/local/netcdf/include
#      and compiler libraries -lirc -limf
#      from /opt/intel/compilers_and_libraries_2018.1.163/linux/lib/intel64
#      including /opt/intel/compilers_and_libraries_2018.1.163/linux/compiler/include/intel64
#
#   Set environment FC,LAPACKHOME,NETCDFHOME to change these
#   ------------------------------------------------------------
#
#   Now compile using
#    $ make OMPFLAG="-qopenmp" 2>&1 | tee make.log
#
# On iris.gat.com
#   Developers should use the ifort compiler and associated netcdf module
#     $ module load ifort/ifort2017 netcdf/4.4.1-mpich3.2-ifort2017
#   Then set one undefined environment variable
#     $ export F90HOME=/fusion/usc/opt/intel2018/compilers_and_libraries_2018.2.199/linux
#   Alternately, you can simply load the gpec module!
#     $ module load gpec
#   This should properly set up your LD_LIBRARY_PATH, PATH, and environment variables.
#    Either way, executing,
#     $ make v
#   in this install directory should produce the following verbose message.
#
#   ------------------------------------------------------------
#   Defaults tailored for
#   Compiling with FC = /fusion/usc/opt/intel2018/compilers_and_libraries_2018.2.199/linux/bin/intel64/ifort
#      using math libraries -lmkl_rt
#      from /fusion/usc/opt/intel2017//compilers_and_libraries/linux/mkl/lib/intel64
#      including /fusion/usc/opt/intel2017//compilers_and_libraries/linux/mkl/include/intel64/lp64
#      and netcdf libraries -lnetcdff -lnetcdf
#      from /fusion/usc/opt/netcdf/netcdf-4.4.1_mpich-3.2_ifort2017/lib /fusion/usc/opt/netcdf/netcdf-4.4.1_mpich-3.2_ifort2017/lib
#      including /fusion/usc/opt/netcdf/netcdf-4.4.1_mpich-3.2_ifort2017/include
#      and compiler libraries -lirc -limf
#      from /fusion/usc/opt/intel2018/compilers_and_libraries_2018.2.199/linux//lib/intel64
#      including /fusion/usc/opt/intel2018/compilers_and_libraries_2018.2.199/linux//compiler/include/intel64
#
#   Set environment FC,LAPACKHOME,NETCDFHOME to change these
#   ------------------------------------------------------------
#
#   Now compile using
#    $ make FFLAGS="-qopenmp" 2>&1 | tee make.log
#
# On toki0*.aug.ipp.mpg.de
#   Developers should use the ifort compiler as well as the associated impi and netcdf modules
#     $ module load intel/17.0 impi/4.1.3 netcdf-mpi/4.3.2 mkl/2017 git
#   Then set undefined environment variables
#     $ export F90HOME=/afs/@cell/common/soft/intel/ics2017.1/17.0/linux
#     $ export FC=/afs/@cell/common/soft/intel/ics2017.1/17.0/linux/bin/intel64/bin/ifort
#     $ export NETCDFHOME=/afs/rzg/common/soft/netcdf/4.3.2/amd64_sles11/intel/14.0/impi-4.1.3
#   Alternately, you can simply load the gpec module!
#     $ module use /afs/.ipp/amd64_sles15/modules.addon.2019/TOK
#     $ module load gpec
#   This should properly set up your LD_LIBRARY_PATH, PATH, and environment variables.
#   Either way, executing,
#     $ make v
#   in this install directory should produce the following verbose message.
#
#   nlogan@toki02 105)
#   ------------------------------------------------------------
#   Defaults tailored for
#   Compiling with FC = /afs/@cell/common/soft/intel/ics2017.1/17.0/linux/bin/intel64/bin/ifort
#      using math libraries -lmkl_rt
#      from /afs/@cell/common/soft/intel/ics2017.1/17.0/linux/mkl/lib/intel64
#      including /afs/@cell/common/soft/intel/ics2017.1/17.0/linux/mkl/include/intel64/lp64
#      and netcdf libraries -lnetcdff -lnetcdf
#      from /afs/rzg/common/soft/netcdf/4.3.2/amd64_sles11/intel/14.0/impi-4.1.3/lib /afs/rzg/common/soft/netcdf/4.3.2/amd64_sles11/intel/14.0/impi-4.1.3/lib
#      including /afs/rzg/common/soft/netcdf/4.3.2/amd64_sles11/intel/14.0/impi-4.1.3/include
#      and compiler libraries -lirc -limf
#      from /afs/@cell/common/soft/intel/ics2017.1/17.0/linux/compiler/lib/intel64
#      including /afs/@cell/common/soft/intel/ics2017.1/17.0/linux/compiler/include/intel64
#
#   Set environment FC,LAPACKHOME,NETCDFHOME to change these
#   ------------------------------------------------------------
#
#   Now compile using
#     $ make FFLAGS='-openmp' 2>&1 | tee make.log
#
# On MIT's mferws01 linux machine
#   Developers should use the gfortran compiler as well as the associated netcdf. Unfortunately, there are no modules here
#   Set environment variables manually as:
#     $ export FC=gfortran
#     $ export LAPACKHOME=/usr/lib/x86_64-linux-gnu/
#     $ export NETCDFHOME=/usr/lib/x86_64-linux-gnu/
#     $ export NETCDFINC=/usr/include
#   Executing,
#     $ make v
#   in this install directory should produce the following verbose message.
#   ------------------------------------------------------------
#   Defaults tailored for
#   Compiling with FC = gfortran
#      using math libraries -llapack -lblas
#      from /usr/lib/x86_64-linux-gnu/
#      including /usr/lib/x86_64-linux-gnu/
#      and netcdf libraries -lnetcdff -lnetcdf
#      from /usr/lib/x86_64-linux-gnu/ /usr/lib/x86_64-linux-gnu/
#      including /usr/lib/x86_64-linux-gnu/
#      and compiler libraries
#      from .
#      including .
#      and using OMPFLAG -fopenmp
#
#   Set environment FC,LAPACKHOME,NETCDFHOME to change these
#   ------------------------------------------------------------
#
#   Now compile using
#     $ make FFLAGS='' 2>&1 | tee make.log
#
# On NERSC Permutter
#   GPEC is installed at /global/common/software/m3739/perlmutter/
#   Load the linux modules
#     $ module use /global/common/software/m3739/perlmutter/modules
#     $ module load gpec
#   Executing,
#     $ make v
#   in this install directory should produce the following verbose message.
#   ------------------------------------------------------------
#   Defaults tailored for
#   Compiling with FC = gfortran
#      using math libraries -llapack -lblas
#      from /usr/lib64
#      including /usr/lib64
#      and netcdf libraries -lnetcdff -lnetcdf
#      from /opt/cray/pe/netcdf/4.9.0.9/gnu/12.3/lib /opt/cray/pe/netcdf/4.9.0.9/gnu/12.3/lib
#      including /opt/cray/pe/netcdf/4.9.0.9/gnu/12.3/include
#      and compiler libraries
#      from .
#      including .
#      and using OMPFLAG -fopenmp
#
#   Set environment FC,LAPACKHOME,NETCDFHOME to change these
#   ------------------------------------------------------------
#
#   Now compile using
#     $ make FFLAGS='-fallow-argument-mismatch' 2>&1 | tee make.log
#   
#
# On intel mac OSX
#   If one has installed the free gfortran compiler and built netcdf and lapack,
#   a typical set of paths looks like the following.
#     $ make v
#   ------------------------------------------------------------
#   Defaults tailored for
#   Compiling with FC = /usr/local/bin/gfortran
#      using math libraries -llapack -lblas
#      from /usr/lib
#      including /usr/lib
#      and netcdf libraries -lnetcdff -lnetcdf
#      from /opt/local/lib /opt/local/lib
#      including /opt/local/include
#      and compiler libraries
#      from .
#      including .
#
#   Set environment FC,LAPACKHOME,NETCDFHOME to change these
#   --------------------------------------------------------------
#
#   Now compile using
#     $ make FFLAGS='-fallow-argument-mismatch' 2>&1 | tee make.log
#
#
#   On ARM mac OSX using Apple's optimized accelerate library
#   The accelerate framework contains the necessary math libraries from
#   BLAS and LAPACK, so you should compile netcdf-c and netcdf-fortran
#   locally (or install them via homebrew). 
#   If you install netcdf-c & -fortran using homebrew, you  
#   can do `echo $HOMEBREW_PREFIX` to get the appropriate NETCDFHOME
#   variable, as it has changed recently. Example with:
#   CC = /opt/homebrew/bin/gcc-14
#   FC = /opt/homebrew/bin/gfortran-14
#   NETCDFHOME = /opt/homebrew
#   LAPACKHOME = /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/
#   
#   $ make v
#   ------------------------------------------------------------
#   Defaults tailored for 
#   Compiling with FC = /opt/homebrew/bin/gfortran-14
#      using math libraries -llapack -lblas
#      from /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/
#      including /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/
#      and netcdf libraries -lnetcdff -lnetcdf
#      from /opt/homebrew/lib /opt/homebrew/lib
#      including /opt/homebrew/include
#      and compiler libraries 
#      from .
#      including .
#      and using OMPFLAG -fopenmp
#   
#   Set environment FC,LAPACKHOME,NETCDFHOME to change these
#   ------------------------------------------------------------
#   
#   Now compile using
#     $ make FFLAGS='-fallow-argument-mismatch' 2>&1 | tee make.log
#
# Some common modifications of the above settings might include,
#     $ make FFLAGS="-qopt-report -g -O0 -check all"  # for debugging
#     $ make FFLAGS="-O3"     # for optimized version
#     $ make FFLAGS="-xHost -O3 -qopt-zmm-usage=high -ipo"  # for old STRIDE options
# When using modern gfortran, some additional flags needed may include
#     $ make FFLAGS="-fallow-argument-mismatch"  # to pass over old, inconsistent portions of code
#
# On macs using gfortran, one may need to add links to X11 in xdraw/makefile for xdraw to compile properly:
#     CFLAGS += -I/opt/X11/include
#     LDFLAGS += -L/opt/X11/lib
#     LIBS += -lX11
#
# If your platform is not properly included in this make file
# contact nikolas.logan@columbia.edu with your operating system, compiler, and
# library installation details.
#
#   ---------------------------------------------------------------
# On SNU PLARE's Linux machine
#   Developers can use intel compilers with the following modules
#   module load git compiler-rt/latest mkl/latest 
#   module load netcdf-fortran/intel/4.5.2
#   And Set environment variables manually as:
#     $ export /APP/enhpc/compiler/intel/oneapi/intelpython/python3.7/lib:/APP/netcdf-fortran-4.5.2_gnu/lib64:/usr/local/lib:/APP/MATLAB/R2020a/bin/glnxa64:$LD_LIBRARY_PATH
#  Executing,
#     $ make v
#   in this install directory should produce the following verbose message.
#   --------------------------------------------------------------
#   Defaults tailored for
#   Compiling with FC = /APP/enhpc/compiler/intel/oneapi/compiler/2021.3.0/linux/bin/intel64/ifort
#      using math libraries -lmkl_rt
#      from /APP/enhpc/compiler/intel/oneapi/mkl/2021.3.0/lib/intel64
#      including /APP/enhpc/compiler/intel/oneapi/mkl/2021.3.0/include/intel64/lp64
#      and netcdf libraries -lnetcdff -lnetcdf
#      from /APP/netcdf-fortran-4.5.2_icc/lib /APP/netcdf-fortran-4.5.2_icc/lib
#      including /APP/netcdf-fortran-4.5.2_icc/include
#      and compiler libraries -lirc -limf
#      from /APP/enhpc/compiler/intel/oneapi/compiler/2021.3.0/linux/compiler/lib/intel64
#      including /APP/enhpc/compiler/intel/oneapi/compiler/2021.3.0/linux/compiler/include/intel64
#      and using OMPFLAG -qopenmp
#   Set environment FC,LAPACKHOME,NETCDFHOME to change these
#-----------------------------------------------------------------
#
#   Now compile using
#    $ make
#
#-----------------------------------------------------------------
# general defaults
#-----------------------------------------------------------------

# The new modules on portal at PPPL now use LAPACK_HOME instead of LAPACKHOME
ifndef LAPACKHOME
   ifdef LAPACK_HOME
         LAPACKHOME = $(LAPACK_HOME)
   endif
endif

# If using the Intel compiler, the MKL library is a better choice for lapack+blas
# The ACML library is not being maintained anymore but we can leave it there for now

# use acml or lapack based on what is loaded
ifdef MKLROOT
    MATHHOME = $(MKLROOT)
  	MATHLIBS = -lmkl_rt
else ifdef OPENBLASHOME
    MATHHOME = $(OPENBLASHOME)
    MATHLIBS = -lopenblas
else ifdef LAPACKHOME
    MATHHOME = $(LAPACKHOME)
    MATHLIBS = -llapack -lblas
else ifdef ACML_HOME
    MATHHOME = $(ACML_HOME)
    MATHLIBS = -lacml
else
    $(info "Math libraries require MKLROOT, LAPACKHOME, OPENBLASHOME, or ACML_HOME")
endif

# specifics of which libraries are used
#---------------------------------------

# tack on lib and include directories to math home if necessary
ifneq ($(wildcard $(MATHHOME)/lib/intel64/.),)
    MATHDIR = $(MATHHOME)/lib/intel64
else ifneq ($(wildcard $(MATHHOME)/lib/.),)
    MATHDIR = $(MATHHOME)/lib
else ifneq ($(wildcard $(MATHHOME)/lib64/.),)
    MATHDIR = $(MATHHOME)/lib64
else
    MATHDIR = $(MATHHOME)
endif
ifneq ($(wildcard $(MATHHOME)/include/intel64/lp64/.),)
    MATHINC = $(MATHHOME)/include/intel64/lp64
else ifneq ($(wildcard $(MATHHOME)/include/.),)
    MATHINC = $(MATHHOME)/include
else
    MATHINC = $(MATHHOME)
endif

# tack on lib and include directories to netcdf home if necessary
# and catch the fact that PPPL split up fortran and c netcdf modules
ifdef NETCDF_FORTRAN_HOME
    ifneq ($(wildcard $(NETCDF_FORTRAN_HOME)/lib/.),)
        NETCDFDIR = $(NETCDF_FORTRAN_HOME)/lib
    else
        NETCDFDIR = $(NETCDF_FORTRAN_HOME)
    endif
    ifndef NETCDFINC
        ifneq ($(wildcard $(NETCDF_FORTRAN_HOME)/include/.),)
            NETCDFINC = $(NETCDF_FORTRAN_HOME)/include
        else
            NETCDFINC = $(NETCDF_FORTRAN_HOME)
        endif
    endif
else
    ifndef NETCDFHOME
        ifndef NETCDF_DIR
            $(info "Must set NETCDF_FORTRAN_HOME, NETCDFHOME, or NETCDF_DIR")
        else
            NETCDFHOME = $(NETCDF_DIR)
        endif
    endif
    ifneq ($(wildcard $(NETCDFHOME)/lib/.),)
        NETCDFDIR = $(NETCDFHOME)/lib
    else
        NETCDFDIR = $(NETCDFHOME)
    endif
    ifndef NETCDFINC
        ifneq ($(wildcard $(NETCDFHOME)/include/.),)
            NETCDFINC = $(NETCDFHOME)/include
        else
            NETCDFINC = $(NETCDFHOME)
        endif
    endif
endif
# check if the netcdf and netcdf libraries are split up or not
ifdef NETCDF_C_HOME
    NETCDFLIBS = -lnetcdff
    NETCDF_C_LIBS = -lnetcdf
    NETCDF_C_DIR = $(NETCDF_C_HOME)/lib
    NETCDF_EXTRA_LIBS = -L$(NETCDF_C_DIR) $(NETCDF_C_LIBS)
else
    ifneq ($(shell ls $(NETCDFDIR) | grep netcdff),)
        NETCDFLIBS = -lnetcdff
        NETCDF_C_LIBS = -lnetcdf
        NETCDF_C_DIR = $(NETCDFDIR)
        NETCDF_EXTRA_LIBS = -L$(NETCDF_C_DIR) $(NETCDF_C_LIBS)
    else
        NETCDFLIBS = -lnetcdf
        NETCDF_C_LIBS = -lnetcdf
        NETCDF_C_DIR = $(NETCDFDIR)
        NETCDF_EXTRA_LIBS =
    endif
endif

# Intel specifics
#----------------
# include compiler specific libs for intel fast moves
ifneq ($(wildcard $(F90HOME)/lib/intel64/.),)
    FCDIR = $(F90HOME)/lib/intel64
#    FCLIBS = -lirc -limf # irc.a needed for static compile
    FCLIBS = -lintlc -limf # intel changed irc to intlc ~2013?
    FCINC = $(F90HOME)/include/intel64
    NETCDF_EXTRA_LIBS += -L$(FCDIR) $(FCLIBS)
else ifneq ($(wildcard $(F90HOME)/compiler/lib/intel64/.),)
    FCDIR = $(F90HOME)/compiler/lib/intel64
    FCLIBS = -lirc -limf # irc.a needed for static compile
    FCINC = $(F90HOME)/compiler/include/intel64
    NETCDF_EXTRA_LIBS += -L$(FCDIR) $(FCLIBS)
else
    FCDIR = .
    FCLIBS =
    FCINC = .
endif

# XDRAW options
#----------------
ifdef X11_HOME
    X11_INC_DIR = -I$(X11_HOME)/include
    X11_LIB_DIR = -L$(X11_HOME)/lib
endif

#Different flags for different compilers
ifeq ($(findstring ifort,$(FC)),ifort)
	OMPFLAG ?= -qopenmp
    ZEROFLAG ?= -zero
    RECURSFLAG ?= -recursive
    LEGACYFLAG ?=
else ifeq ($(findstring ifx,$(FC)),ifx)
	OMPFLAG ?= -qopenmp
    ZEROFLAG ?= -zero
    RECURSFLAG ?= -recursive
    LEGACYFLAG ?=
else ifeq ($(findstring gfortran,$(FC)),gfortran)
	OMPFLAG ?= -fopenmp
    ZEROFLAG ?= -finit-local-zero
    RECURSFLAG ?= -frecursive
    LEGACYFLAG ?= -std=legacy
else ifeq ($(findstring pgfortran,$(FC)),pgfortran)
	OMPFLAG ?= -mp
    ZEROFLAG ?= -zeroinit
    RECURSFLAG ?= -Mrecursive
    LEGACYFLAG ?= -std=legacy
else
    ZEROFLAG ?= -zero
	OMPFLAG ?= -openmp
    RECURSFLAG ?= -recursive
    LEGACYFLAG ?= -std=legacy
endif