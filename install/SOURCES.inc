 # Auto-discover sources

EQUIL_SOURCES = $(wildcard ../equil/*.f ../equil/*.f90 ../equil/*.F ../equil/*.F90)
EQUIL_OBJS = $(EQUIL_SOURCES:.f=.o)
EQUIL_OBJS := $(EQUIL_OBJS:.f90=.o)

LSODE_SOURCES = $(wildcard ../lsode/*.f ../lsode/*.f90 ../lsode/*.F ../lsode/*.F90)
LSODE_OBJS = $(LSODE_SOURCES:.f=.o)
LSODE_OBJS := $(LSODE_OBJS:.f90=.o)

ZLANGE_SOURCES = $(wildcard ../zlange/*.f ../zlange/*.f90 ../zlange/*.F ../zlange/*.F90)
ZLANGE_OBJS = $(ZLANGE_SOURCES:.f=.o)
ZLANGE_OBJS := $(ZLANGE_OBJS:.f90=.o)

ZVODE_SOURCES = $(wildcard ../zvode/*.f ../zvode/*.f90 ../zvode/*.F ../zvode/*.F90)
ZVODE_OBJS = $(ZVODE_SOURCES:.f=.o)
ZVODE_OBJS := $(ZVODE_OBJS:.f90=.o)

ORBIT_SOURCES = $(wildcard ../orbit/*.f ../orbit/*.f90 ../orbit/*.F ../orbit/*.F90)
ORBIT_OBJS = $(ORBIT_SOURCES:.f=.o)
ORBIT_OBJS := $(ORBIT_OBJS:.f90=.o)

VACUUM_SOURCES = $(wildcard ../vacuum/*.f ../vacuum/*.f90 ../vacuum/*.F ../vacuum/*.F90)
VACUUM_OBJS = $(VACUUM_SOURCES:.f=.o)
VACUUM_OBJS := $(VACUUM_OBJS:.f90=.o)

PENTRC_SOURCES = $(wildcard ../pentrc/*.f ../pentrc/*.f90 ../pentrc/*.F ../pentrc/*.F90)
PENTRC_OBJS = $(PENTRC_SOURCES:.f=.o)
PENTRC_OBJS := $(PENTRC_OBJS:.f90=.o)

DCON_SOURCES = $(wildcard ../dcon/*.f ../dcon/*.f90 ../dcon/*.F ../dcon/*.F90)
DCON_OBJS = $(DCON_SOURCES:.f=.o)
DCON_OBJS := $(DCON_OBJS:.f90=.o)

MATCH_SOURCES = $(wildcard ../match/*.f ../match/*.f90 ../match/*.F ../match/*.F90)
MATCH_OBJS = $(MATCH_SOURCES:.f=.o)
MATCH_OBJS := $(MATCH_OBJS:.f90=.o)

RDCON_SOURCES = $(wildcard ../dcon/*.f ../dcon/*.f90 ../dcon/*.F ../dcon/*.F90)
RDCON_OBJS = $(RDCON_SOURCES:.f=.o)
RDCON_OBJS := $(RDCON_OBJS:.f90=.o)

RMATCH_SOURCES = $(wildcard ../rmatch/*.f ../rmatch/*.f90 ../rmatch/*.F ../rmatch/*.F90)
RMATCH_OBJS = $(RMATCH_SOURCES:.f=.o)
RMATCH_OBJS := $(RMATCH_OBJS:.f90=.o)

MULTI_SOURCES = $(wildcard ../multi/*.f ../multi/*.f90 ../multi/*.F ../multi/*.F90)
MULTI_OBJS = $(MULTI_SOURCES:.f=.o)
MULTI_OBJS := $(MULTI_OBJS:.f90=.o)

SUM_SOURCES = $(wildcard ../sum/*.f ../sum/*.f90 ../sum/*.F ../sum/*.F90)
SUM_OBJS = $(SUM_SOURCES:.f=.o)
SUM_OBJS := $(SUM_OBJS:.f90=.o)

SLAYER_SOURCES = $(wildcard ../slayer/*.f ../slayer/*.f90 ../slayer/*.F ../slayer/*.F90)
SLAYER_OBJS = $(SLAYER_SOURCES:.f=.o)
SLAYER_OBJS := $(SLAYER_OBJS:.f90=.o)

COIL_SOURCES = $(wildcard ../coil/*.f ../coil/*.f90 ../coil/*.F ../coil/*.F90)
COIL_OBJS = $(COIL_SOURCES:.f=.o)
COIL_OBJS := $(COIL_OBJS:.f90=.o)

GPEC_SOURCES = $(wildcard ../gpec/*.f ../gpec/*.f90 ../gpec/*.F ../gpec/*.F90)
GPEC_OBJS = $(GPEC_SOURCES:.f=.o)
GPEC_OBJS := $(GPEC_OBJS:.f90=.o)

STRIDE_SOURCES = $(wildcard ../stride/*.f ../stride/*.f90 ../stride/*.F ../stride/*.F90)
STRIDE_OBJS = $(STRIDE_SOURCES:.f=.o)
STRIDE_OBJS := $(STRIDE_OBJS:.f90=.o)

XDRAW_SOURCES = $(wildcard ../xdraw/*.f ../xdraw/*.f90 ../xdraw/*.F ../xdraw/*.F90)
XDRAW_OBJS = $(XDRAW_SOURCES:.f=.o)
XDRAW_OBJS := $(XDRAW_OBJS:.f90=.o)