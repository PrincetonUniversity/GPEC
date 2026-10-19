#------------------------------------------------------------------------------
# SOURCES.inc
#
# Purpose:
#   Automatically finds all source files (.f, .f90, .F, .F90, .c) for each GPEC module.
#
# What it does:
#   - Uses `wildcard` to locate Fortran and C source files in each module directory
#   - Defines lists of object files for each module (EQUIL_OBJS, LSODE_OBJS, etc.)
#   - Keeps build up-to-date as files are added/removed
#
# How to use:
#   - Usually you won't need to edit this file! It detects new files automatically.
#   - To add a new source file, simply place it in the appropriate module directory.
#   - To add a new module, copy the pattern shown and create NEWMODULE_SOURCES/OBJS.
#   - If you are adding a c source and c is not included in your module's sources
#     section below, then add it! See vacuum as an example. 
#------------------------------------------------------------------------------

# EQUIL sources
EQUIL_SOURCES = $(wildcard ../equil/*.f ../equil/*.f90 ../equil/*.F ../equil/*.F90)
EQUIL_OBJS = $(EQUIL_SOURCES:.f=.o)
EQUIL_OBJS := $(EQUIL_OBJS:.f90=.o)
EQUIL_OBJS := $(EQUIL_OBJS:.F=.o)
EQUIL_OBJS := $(EQUIL_OBJS:.F90=.o)

# LSODE sources
LSODE_SOURCES = $(wildcard ../lsode/*.f ../lsode/*.f90 ../lsode/*.F ../lsode/*.F90)
LSODE_OBJS = $(LSODE_SOURCES:.f=.o)
LSODE_OBJS := $(LSODE_OBJS:.f90=.o)
LSODE_OBJS := $(LSODE_OBJS:.F=.o)
LSODE_OBJS := $(LSODE_OBJS:.F90=.o)

#ZLANGE sources
ZLANGE_SOURCES = $(wildcard ../zlange/*.f ../zlange/*.f90 ../zlange/*.F ../zlange/*.F90)
ZLANGE_OBJS = $(ZLANGE_SOURCES:.f=.o)
ZLANGE_OBJS := $(ZLANGE_OBJS:.f90=.o)
ZLANGE_OBJS := $(ZLANGE_OBJS:.F=.o)
ZLANGE_OBJS := $(ZLANGE_OBJS:.F90=.o)

# ZVODE sources
ZVODE_SOURCES = $(wildcard ../zvode/*.f ../zvode/*.f90 ../zvode/*.F ../zvode/*.F90)
ZVODE_OBJS = $(ZVODE_SOURCES:.f=.o)
ZVODE_OBJS := $(ZVODE_OBJS:.f90=.o)
ZVODE_OBJS := $(ZVODE_OBJS:.F=.o)
ZVODE_OBJS := $(ZVODE_OBJS:.F90=.o)

# ORBIT sources
ORBIT_SOURCES = $(wildcard ../orbit/*.f ../orbit/*.f90 ../orbit/*.F ../orbit/*.F90)
ORBIT_OBJS = $(ORBIT_SOURCES:.f=.o)
ORBIT_OBJS := $(ORBIT_OBJS:.f90=.o)
ORBIT_OBJS := $(ORBIT_OBJS:.F=.o)
ORBIT_OBJS := $(ORBIT_OBJS:.F90=.o)

# VACUUM sources
VACUUM_SOURCES = $(wildcard ../vacuum/*.f ../vacuum/*.f90 ../vacuum/*.F ../vacuum/*.F90 ../vacuum/*.c)
VACUUM_OBJS = $(VACUUM_SOURCES:.f=.o)
VACUUM_OBJS := $(VACUUM_OBJS:.f90=.o)
VACUUM_OBJS := $(VACUUM_OBJS:.F=.o)
VACUUM_OBJS := $(VACUUM_OBJS:.F90=.o)
VACUUM_OBJS := $(VACUUM_OBJS:.c=.o)

# PENTRC sources
PENTRC_SOURCES = $(wildcard ../pentrc/*.f ../pentrc/*.f90 ../pentrc/*.F ../pentrc/*.F90)
PENTRC_OBJS = $(PENTRC_SOURCES:.f=.o)
PENTRC_OBJS := $(PENTRC_OBJS:.f90=.o)
PENTRC_OBJS := $(PENTRC_OBJS:.F=.o)
PENTRC_OBJS := $(PENTRC_OBJS:.F90=.o)

# DCON sources
DCON_SOURCES = $(wildcard ../dcon/*.f ../dcon/*.f90 ../dcon/*.F ../dcon/*.F90)
DCON_OBJS = $(DCON_SOURCES:.f=.o)
DCON_OBJS := $(DCON_OBJS:.f90=.o)
DCON_OBJS := $(DCON_OBJS:.F=.o)
DCON_OBJS := $(DCON_OBJS:.F90=.o)

# RDCON sources
RDCON_SOURCES = $(wildcard ../rdcon/*.f ../rdcon/*.f90 ../rdcon/*.F ../rdcon/*.F90)
RDCON_OBJS = $(RDCON_SOURCES:.f=.o)
RDCON_OBJS := $(RDCON_OBJS:.f90=.o)
RDCON_OBJS := $(RDCON_OBJS:.F=.o)
RDCON_OBJS := $(RDCON_OBJS:.F90=.o)

# MATCH sources
MATCH_SOURCES = $(wildcard ../match/*.f ../match/*.f90 ../match/*.F ../match/*.F90)
MATCH_OBJS = $(MATCH_SOURCES:.f=.o)
MATCH_OBJS := $(MATCH_OBJS:.f90=.o)
MATCH_OBJS := $(MATCH_OBJS:.F=.o)
MATCH_OBJS := $(MATCH_OBJS:.F90=.o)

# RMATCH sources
RMATCH_SOURCES = $(wildcard ../rmatch/*.f ../rmatch/*.f90 ../rmatch/*.F ../rmatch/*.F90)
RMATCH_OBJS = $(RMATCH_SOURCES:.f=.o)
RMATCH_OBJS := $(RMATCH_OBJS:.f90=.o)
RMATCH_OBJS := $(RMATCH_OBJS:.F=.o)
RMATCH_OBJS := $(RMATCH_OBJS:.F90=.o)

# STRIDE sources
STRIDE_SOURCES = $(wildcard ../stride/*.f ../stride/*.f90 ../stride/*.F ../stride/*.F90)
STRIDE_OBJS = $(STRIDE_SOURCES:.f=.o)
STRIDE_OBJS := $(STRIDE_OBJS:.f90=.o)
STRIDE_OBJS := $(STRIDE_OBJS:.F=.o)
STRIDE_OBJS := $(STRIDE_OBJS:.F90=.o)

# SUM sources
SUM_SOURCES = $(wildcard ../sum/*.f ../sum/*.f90 ../sum/*.F ../sum/*.F90)
SUM_OBJS = $(SUM_SOURCES:.f=.o)
SUM_OBJS := $(SUM_OBJS:.f90=.o)
SUM_OBJS := $(SUM_OBJS:.F=.o)
SUM_OBJS := $(SUM_OBJS:.F90=.o)

# MULTI sources
MULTI_SOURCES = $(wildcard ../multi/*.f ../multi/*.f90 ../multi/*.F ../multi/*.F90)
MULTI_OBJS = $(MULTI_SOURCES:.f=.o)
MULTI_OBJS := $(MULTI_OBJS:.f90=.o)
MULTI_OBJS := $(MULTI_OBJS:.F=.o)
MULTI_OBJS := $(MULTI_OBJS:.F90=.o)

# SLAYER sources
SLAYER_SOURCES = $(wildcard ../slayer/*.f ../slayer/*.f90 ../slayer/*.F ../slayer/*.F90)
SLAYER_OBJS = $(SLAYER_SOURCES:.f=.o)
SLAYER_OBJS := $(SLAYER_OBJS:.f90=.o)
SLAYER_OBJS := $(SLAYER_OBJS:.F=.o)
SLAYER_OBJS := $(SLAYER_OBJS:.F90=.o)

# COIL sources
COIL_SOURCES = $(wildcard ../coil/*.f ../coil/*.f90 ../coil/*.F ../coil/*.F90)
COIL_OBJS = $(COIL_SOURCES:.f=.o)
COIL_OBJS := $(COIL_OBJS:.f90=.o)
COIL_OBJS := $(COIL_OBJS:.F=.o)
COIL_OBJS := $(COIL_OBJS:.F90=.o)

# GPEC sources
GPEC_SOURCES = $(wildcard ../gpec/*.f ../gpec/*.f90 ../gpec/*.F ../gpec/*.F90)
GPEC_OBJS = $(GPEC_SOURCES:.f=.o)
GPEC_OBJS := $(GPEC_OBJS:.f90=.o)
GPEC_OBJS := $(GPEC_OBJS:.F=.o)
GPEC_OBJS := $(GPEC_OBJS:.F90=.o)

# XDRAW sources
XDRAW_OBJS1 = \
	../xdraw/binread.o \
	../xdraw/event.o \
	../xdraw/gradient.o \
	../xdraw/limits.o \
	../xdraw/ps.o \
	../xdraw/redraw.o \
	../xdraw/setcolor.o \
	../xdraw/spline.o \
	../xdraw/spline1.o \
	../xdraw/xcontour.o \
	../xdraw/xdraw.o \
	../xdraw/xedit.o \
	../xdraw/xinit.o \
	../xdraw/xtools.o

XDRAW_OBJS2 = \
	../xdraw/ps.o \
	../xdraw/pspack.o

# All object files
ALL_OBJS = $(EQUIL_OBJS) $(LSODE_OBJS) $(ZLANGE_OBJS) $(ZVODE_OBJS) $(ORBIT_OBJS) \
		   $(VACUUM_OBJS) $(PENTRC_OBJS) $(DCON_OBJS) $(RDCON_OBJS) \
		   $(MATCH_OBJS) $(RMATCH_OBJS) $(STRIDE_OBJS) $(SUM_OBJS) $(MULTI_OBJS) \
		   $(SLAYER_OBJS) $(COIL_OBJS) $(GPEC_OBJS) $(XDRAW_OBJS1) $(XDRAW_OBJS2)