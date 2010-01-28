#-------------------------------------------------------------------------------

HOST    := $(shell hostname)
OS      := $(shell uname)

DEBUG   = N
PREC    = 8

FC      = ifort
FFLAGS  = -O3 -i4 -mtune=core2 -xS -I$(HDF5DIR)/include
LD      = $(FC)
LDFLAGS = $(FFLAGS)
LIBS    =

ifeq ($(DEBUG),Y)
FFLAGS += -g
endif
FFLAGS += -DPREC=$(PREC)

# HDF5 include and library paths
#
FFLAGS += -I$(HDF5DIR)/include
LIBS   += -L$(HDF5DIR)/lib
LIBS   += -lhdf5_fortran -lhdf5 -lsz

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

sources = fitsio.F90 hdf5io.F90 interpolation.F90 paccel.F90 params.F90
objects = fitsio.o   hdf5io.o   interpolation.o   paccel.o   params.o
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

paccel.o        : makefile paccel.F90 params.o hdf5io.o fitsio.o interpolation.o
params.o        : makefile params.F90
fitsio.o        : makefile fitsio.F90 params.o
hdf5io.o        : makefile hdf5io.F90 params.o
interpolation.o : makefile interpolation.F90 fitsio.o hdf5io.o params.o

#-------------------------------------------------------------------------------
