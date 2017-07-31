# simple makefile for schro_1d program
# Copyright (C) Fabien Brieuc - 2017
FC = gfortran
FFLAGS = -O3 -fconvert='big-endian' -Wall
#FFLAGS = -O3 -fconvert='big-endian' -Wall -Wextra -fbounds-check
FLIBS = -L/usr/local/lib -llapack
EXE = schro_1d.exe
OBJ = schro_1d.f90

$(EXE): $(OBJ)
	$(FC) schro_1d.f90 -o $(EXE) $(FFLAGS) $(FLIBS)

clean:
	rm -f *.mod *.exe
