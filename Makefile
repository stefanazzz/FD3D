#
# Makefile for PML3d Finite differences
#

#FC = f90       # f90 compiler
FC = gfortran   # f90 compiler
FC7 = f77      # f77 compiler

OBJ  =  triffy.o StressPMLpack.o VeloPMLpack.o InitDist.o PMLMed.o drives.o drivev.o VectI.o

#FFLAGS = -O3 -pfa -64 -mp
#FFLAGS = -ffpe-trap=invalid,zero,overflow -g
#FFLAGS = -g
FFLAGS = -O3 -floop-nest-optimize

###############################
.SUFFIXES: .f90 .dec

all: triffy_run.x

DEC =  triffy.dec parameters.f
##############################
# Objects
##############################
StressPMLpack.o: StressPMLpack.f90 $(DEC)
	$(FC) $(FFLAGS) -c $<

VeloPMLpack.o: VeloPMLpack.f90 $(DEC)
	$(FC) $(FFLAGS) -c $<

InitDist.o: InitDist.f90 $(DEC)
	$(FC) $(FFLAGS) -c $<

VectI.o: VectI.f90 $(DEC)
	$(FC) $(FFLAGS) -c $<

PMLMed.o: PMLMed.f90 $(DEC)
	$(FC) $(FFLAGS) -c $<

drives.o: drives.f90 $(DEC)
	$(FC) $(FFLAGS) -c $<

drivev.o: drivev.f90 $(DEC)
	$(FC) $(FFLAGS) -c $<

triffy.o: triffy.f90 $(DEC)
	$(FC) $(FFLAGS) -c $<

triffy_run.x: $(OBJ) 
	$(FC) $(FFLAGS) -o $@ $(OBJ) 
 
clean:
	rm -f *.o
