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
# On portal.pppl.gov
#   Developers should load the default intel/11.0.081 and associated libraries,
#     $ module load intel openmpi netcdf acml/5.3.1/ifort64
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
# On venus.gat.com
#   <Proper compiler options need to be updated for v1.2>
#
# On mac OSX
#   With the move to licensed intel compilers, OSX is no longer supported.
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
