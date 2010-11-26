#-------------------------------------------------------------------------------

HOST         := $(shell hostname -s)
OS           := $(shell uname)

DEBUG         = N
RELAT         = Y
INTERPOLATION = C
TVD           = Y
BNDRY         = N
PREC          = 16
DIMS          = 3
TEST          = WAVE
LOSS          = N

# compiler and options
#
FC            = gfortran
FFLAGS        = -O2
LD            = $(FC)
LDFLAGS       = $(FFLAGS)
LIBS          =

# include the host compilation options
#
ifeq ($(wildcard make.$(HOST)),make.$(HOST))
$(info )
$(info Using host compilation options from "$(wildcard make.$(HOST))".)
$(info )
include make.$(HOST)
else
$(info )
$(info Using default compilation options.)
$(info )
endif

ifeq ($(DEBUG),Y)
FFLAGS += -g -DDEBUG
endif
ifeq ($(RELAT),Y)
FFLAGS += -DRELAT
endif
ifeq ($(INTERPOLATION),L)
FFLAGS += -DTRILIN
endif
ifeq ($(INTERPOLATION),C)
FFLAGS += -DTRICUB
ifeq ($(TVD),Y)
FFLAGS += -DTVD
endif
endif
ifeq ($(BNDRY),Y)
FFLAGS += -DBNDRY
endif
FFLAGS += -DPREC=$(PREC)
FFLAGS += -DDIMS=$(DIMS)
ifeq ($(TEST),WAVE)
FFLAGS += -DTEST -DWTEST
endif
ifeq ($(TEST),ISLAND)
FFLAGS += -DTEST -DITEST
endif
ifeq ($(LOSS),Y)
FFLAGS += -DLOSS
endif

ifeq ($(TEST),N)
# HDF5 include and library paths
#
FFLAGS += -I$(HDF5DIR)/include
LIBS   += -L$(HDF5DIR)/lib
LIBS   += -lhdf5_fortran -lhdf5

# FITS library path
#
LIBS   += -L$(FITSDIR)/lib
LIBS   += -lcfitsio -lz
endif

#-------------------------------------------------------------------------------

.SUFFIXES:
.SUFFIXES: .F90 .o

.F90.o:
	$(FC) -c $(FFLAGS) $*.F90

#-------------------------------------------------------------------------------

name = paccel

default: $(name).x

sources = fields.F90 fitsio.F90 hdf5io.F90 paccel.F90 params.F90 particles.F90
objects = fields.o   fitsio.o   hdf5io.o   paccel.o   params.o   particles.o
files   = $(sources) makefile params.in license.txt

$(name).x: $(objects)
	$(LD) $(LDFLAGS) $(objects) -o $(name).x $(LIBS)

dist: $(files)
	tar czvf $(name).tar.gz $(files)

clean:
	rm -rf *.x *.o *.mod *.il

clean-all:
	rm -rf *.x *.o *.mod *.il *.dat *.bin *~

#-------------------------------------------------------------------------------

fields.o        : makefile fields.F90 params.o fitsio.o hdf5io.o
fitsio.o        : makefile fitsio.F90 params.o
hdf5io.o        : makefile hdf5io.F90 params.o
paccel.o        : makefile paccel.F90 fields.o params.o particles.o
params.o        : makefile params.F90
particles.o     : makefile particles.F90 fields.o params.o

#-------------------------------------------------------------------------------
