# Makefile for use on chmlinxx computers.
# Adapted from make file written by C. Symmonds for the MCE program which was inturn
# made with the aid of the fmkmf utility (Copyright (c) The University of Edinburgh)
# subsequent streamlining alterations have been made 

F90=gfortran
FFLAGS= -Wall -Wextra -ffast-math -ftree-vectorize -march=native -Ofast -O3 -fopenmp 
FLINKS= -llapack -lblas
VPATH=../src

# Link step

all: ZOMBIE

ZOMBIE:randgen_mod.o mod_types.o vars.o infnan_constants.o inf.o infnan.o alarrays.o operators.o readpars.o outputs.o electrons.o zom.o ham.o imgtp.o clean.o gradient_descent.o MainZombie.o #neural_network.o #gram_schmidt.o
		 $(F90) $(FFLAGS) -o $@ $^ $(FLINKS)


# Compile steps

randgen_mod.o:randgen_mod.f90
	 $(F90) $(FFLAGS) -c $<

mod_types.o:mod_types.f90
	 $(F90) $(FFLAGS) -c $<

infnan_constants.o:infnan_constants.f90
	 $(F90) $(FFLAGS) -c $<

inf.o:inf.f90 infnan_constants.o 
	 $(F90) $(FFLAGS) -c $<
	 
infnan.o:infnan.f90 infnan_constants.o inf.o
	 $(F90) $(FFLAGS) -c $<
	 
vars.o:vars.f90 mod_types.o  
	 $(F90) $(FFLAGS) -c $<

alarrays.o:alarrays.f90 vars.o mod_types.o 
	 $(F90) $(FFLAGS) -c $<

operators.o:operators.f90 vars.o alarrays.o mod_types.o 
	 $(F90) $(FFLAGS) -c $<

readpars.o:readpars.f90 vars.o mod_types.o 
	 $(F90) $(FFLAGS) -c $<

outputs.o:outputs.f90 vars.o mod_types.o 
	 $(F90) $(FFLAGS) -c $<

electrons.o:electrons.f90 vars.o alarrays.o
	 $(F90) $(FFLAGS) -c $<

zom.o:zom.f90 vars.o infnan.o mod_types.o randgen_mod.o
	 $(F90) $(FFLAGS) -c $<
	 
ham.o:ham.f90 mod_types.o vars.o alarrays.o  
	 $(F90) $(FFLAGS) -c $<

imgtp.o:imgtp.f90 vars.o alarrays.o mod_types.o 
	 $(F90) $(FFLAGS) -c $<

clean.o:clean.f90 vars.o alarrays.o zom.o ham.o outputs.o readpars.o mod_types.o 
	 $(F90) $(FFLAGS) -c $<

# neural_network.o:neural_network.f90 vars.o alarrays.o mod_types.o 
#	 $(F90) $(FFLAGS) -c $<

gradient_descent.o:gradient_descent.f90 vars.o alarrays.o ham.o imgtp.o outputs.o infnan.o zom.o mod_types.o operators.o randgen_mod.o #neural_network.o
	 $(F90) $(FFLAGS) -c $<

# gram_schmidt.o:gram_schmidt.f90 vars.o alarrays.o ham.o gradient_descent.o grad_d.o imgtp.o outputs.o infnan.o zom.o grad_calc.o
# 	 $(F90) $(FFLAGS) -c $<

MainZombie.o:MainZombie.f90 randgen_mod.o mod_types.o vars.o infnan_constants.o inf.o infnan.o operators.o alarrays.o readpars.o outputs.o electrons.o zom.o ham.o imgtp.o clean.o gradient_descent.o #neural_network.o #gram_schmidt.o
	 $(F90) $(FFLAGS) -c $<

# This entry allows you to type " make clean " to get rid of
# all object and module files

clean:
	rm -rf ZOMBIE f_{files,modd}* *.o *.mod *.exe *~ ../run/*.exe ../run/*~ ../src/*~