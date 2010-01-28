#-------------------------------------------------------------------------------

HOST    := $(shell hostname)
OS      := $(shell uname)

FC      = ifort
#FFLAGS  = -g -O3 -i4 -mtune=core2 -xS -fno-omit-frame-pointer -ipo -finline-limit=1000 -I$(HDF5DIR)/include
FFLAGS  = -g -O3 -i4 -I$(HDF5DIR)/include
LD      = $(FC)
LDFLAGS = $(FFLAGS)
LIBS    = -L$(HDF5DIR)/lib -lhdf5_fortran -lhdf5 -L$(FITSDIR)/lib -lcfitsio -lz

#-------------------------------------------------------------------------------

.SUFFIXES:
.SUFFIXES: .F90 .o

.F90.o:
	$(FC) -c $(FFLAGS) $*.F90

#-------------------------------------------------------------------------------

name = paccel

default: $(name).x

sources = $(name).F90 params.F90 fitsio.F90 hdf5io.F90 interpolation.F90
objects = $(name).o   params.o   fitsio.o   hdf5io.o   interpolation.o
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
