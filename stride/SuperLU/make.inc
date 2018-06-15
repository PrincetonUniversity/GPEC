############################################################################
#
#  Program:         SuperLU
#
#  Module:          make.inc
#
#  Purpose:         Top-level Definitions
#
#  Creation date:   May 10, 2015
#
#  Modified:	    
#		    
#
############################################################################
#
# ToDo: Make this work for shared libraries and out-of-source builds.
#
#  The name of the libraries to be created/linked to
#
SuperLUroot	= /p/gpec/users/aglasser/gpec_stride/stride/SuperLU
SUPERLULIB   	= $(SuperLUroot)/SuperLU_5.2.1/SRC/libsuperlu.a

TMGLIB       	= libtmglib.a

BLASLIB	  	= $(SuperLUroot)/libblas.a
LIBS		= $(SUPERLULIB) $(SuperLUroot)/libblas.a

#
#  The archiver and the flag(s) to use when building archive (library)
#  If your system has no ranlib, set RANLIB = echo.
#
ARCH         = ar
ARCHFLAGS    = cr
RANLIB       = ranlib

CC           = gcc
CFLAGS 	     = -O3 -DNDEBUG -DPRNTlevel=0 -DAdd_  
NOOPTS       = -O0
FORTRAN	     = ifort

LOADER       = $(CC)
LOADOPTS     = -g

#
# The directory in which Matlab is installed
#
MATLAB	     = /Applications/MATLAB_R2015b.app
