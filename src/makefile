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
.SUFFIXES: .f90 .o

.f90.o:
	$(FC) -c $(FFLAGS) $*.f90

#-------------------------------------------------------------------------------

name = paccel

default: $(name).x

sources = $(name).f90 params.f90 fitsio.f90 hdf5io.f90 interpolation.f90
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

paccel.o        : makefile paccel.f90 params.o hdf5io.o fitsio.o interpolation.o
params.o        : makefile params.f90
fitsio.o        : makefile fitsio.f90 params.o
hdf5io.o        : makefile hdf5io.f90 params.o
interpolation.o : makefile interpolation.f90 fitsio.o hdf5io.o params.o

#-------------------------------------------------------------------------------
