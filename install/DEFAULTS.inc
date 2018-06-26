#------------------------------------------------------------------------------
# This file is used to dynamically set some default compiler options
# based on what libraries are available on the current machine.
#
# Prior to version 1.2, GPEC could be compiled using ifort, pgf, gfortran,
# or lfc. All public releases on portal wheree compiled using lfc, which
# required only default user path libraries.
#
# With the release of 1.2, the options have been greatly simplified by
# restricting all compilation to intel compilers compatible wth openmpi.
# Machine specific instructions are below.
#
# Note that GPEC uses git to fetch/update the harvest repository as well as
# to label the version at compilation time. Thus, you must have a working git
# command in order to compile.
#
# On portal.pppl.gov
#   Developers should load the default intel/11.0.081 and associated libraries,
#     $ module load intel openmpi netcdf acml/5.3.1/ifort64 git
#   This should properly set up your LD_LIBRARY_PATH and PATH. Executing,
#     $ make v
#   in this install directory should produce the following verbose message,
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
#   Now compile using
#    $ make FFLAGS="-openmp -vec_report0"
#
# On ukstar
#   Developers should use the default intel compiler and load the netcdf_intel module
#     $ module load netcdf/netcdf_intel
#     $ export F90HOME=/opt/intel/compilers_and_libraries_2018.1.163/linux
#     $ export FC=/opt/intel/compilers_and_libraries_2018.1.163/linux/bin/intel64/ifort

#   This should properly set up your LD_LIBRARY_PATH and PATH. Executing,
#     $ make v
#   in this install directory should produce the following verbose message,
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
#    $ make FFLAGS="-qopenmp"
#
# On iris.gat.com
#   With the move to licensed intel compilers, OSX is no longer actively supported.
#
# On mac OSX
#   With the move to licensed intel compilers, OSX is no longer actively supported.
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
ifdef MKLROOT
    MATHHOME = $(MKLROOT)
  	MATHLIBS = -lmkl_rt
else ifdef LAPACKHOME
    MATHHOME = $(LAPACKHOME)
    MATHLIBS = -llapack -lblas
else ifdef ACML_HOME
    MATHHOME = $(ACML_HOME)
    MATHLIBS = -lacml
else
    echo "Math libraries require MKLROOT, LAPACKHOME, or ACML_HOME"
endif

# specifics of which libraries are used
#---------------------------------------

# tack on lib and include directories to math home if necessary
ifneq ($(wildcard $(MATHHOME)/lib/intel64/.),)
    MATHDIR = $(MATHHOME)/lib/intel64
else ifneq ($(wildcard $(MATHHOME)/lib/.),)
    MATHDIR = $(MATHHOME)/lib
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
    ifneq ($(wildcard $(NETCDF_FORTRAN_HOME)/include/.),)
        NETCDFINC = $(NETCDF_FORTRAN_HOME)/include
    else
        NETCDFINC = $(NETCDF_FORTRAN_HOME)
    endif
else
    ifndef NETCDFHOME
        ifndef NETCDF_DIR
            echo "Must set NETCDF_FORTRAN_HOME, NETCDFHOME, or NETCDF_DIR"
        else
            NETCDFHOME = $(NETCDF_DIR)
        endif
    endif
    ifneq ($(wildcard $(NETCDFHOME)/lib/.),)
        NETCDFDIR = $(NETCDFHOME)/lib
    else
        NETCDFDIR = $(NETCDFHOME)
    endif
    ifneq ($(wildcard $(NETCDFHOME)/include/.),)
        NETCDFINC = $(NETCDFHOME)/include
    else
        NETCDFINC = $(NETCDFHOME)
    endif
endif
# check if the netcdf and netcdf libraries are split up or not
ifdef NETCDF_C_HOME
    NETCDFLIBS = -lnetcdff
    NETCDF_C_LIBS = -lnetcdf
    NETCDF_C_DIR = $(NETCDF_C_HOME)/lib
else
    ifneq ($(shell ls $(NETCDFDIR) | grep netcdff),)
        NETCDFLIBS = -lnetcdff
        NETCDF_C_LIBS = -lnetcdf
        NETCDF_C_DIR = $(NETCDFDIR)
    else
        NETCDFLIBS = -lnetcdf
        NETCDF_C_LIBS = $(NETCDFDIR)
        NETCDF_C_DIR = -lnetcdf
    endif
endif

# Intel specifics
#----------------
# include compiler specific libs for intel fast moves
ifneq ($(wildcard $(F90HOME)/lib/intel64/.),)
    FCDIR = $(F90HOME)/lib/intel64
    FCLIBS = -lirc -limf # irc.a needed for static compile
    #FCLIBS = -lintlc -limf # intlc was new for intel?
    FCINC = $(F90HOME)/include/intel64
else ifneq ($(wildcard $(F90HOME)/compiler/lib/intel64/.),)
    FCDIR = $(F90HOME)/lib/intel64
    FCLIBS = -lirc -limf # irc.a needed for static compile
    FCINC = $(F90HOME)/compiler/include/intel64
else
    FCDIR = .
    FCLIBS =
    FCINC = .
endif
