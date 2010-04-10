#-------------------------------------------------------------------------------

HOST    := $(shell hostname)
OS      := $(shell uname)

DEBUG   = N
INTERP  = C
TVD     = Y
BNDRY   = N
PREC    = 8
DIMS    = 2

FC      = ifort
FFLAGS  = -O3 -i4 -mtune=core2 -xS -I$(HDF5DIR)/include
LD      = $(FC)
LDFLAGS = $(FFLAGS)
LIBS    =

ifeq ($(DEBUG),Y)
FFLAGS += -g -DDEBUG
endif
ifeq ($(INTERP),L)
FFLAGS += -DTRILIN
endif
ifeq ($(INTERP),C)
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

# HDF5 include and library paths
#
FFLAGS += -I$(HDF5DIR)/include
LIBS   += -L$(HDF5DIR)/lib
LIBS   += -lhdf5_fortran -lhdf5

# FITS library path
#
LIBS   += -L$(FITSDIR)/lib
LIBS   += -lcfitsio -lz

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
