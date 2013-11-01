ROUTINES = kriging.o cokriging.o 
MODULES = params.o grid.o matrix.o regression.o correlation.o sensitivity.o mle.o  sequential_sampling.o error.o utils.o 
ANALYTICAL = analytical_functions.o analytical_solver.o
LDFLAGS=-I/usr/lib/lapack95_modules/ -I/usr/lib/fortarrange_header/ -lfort_arrange -llapack95 -llapack -lblas 
BINDIR=../bin

OPT = -O3
#PROFILE=-pg #comment out if you don't want profiling

default : f.f90 $(MODULES) $(ROUTINES) $(ANALYTICAL)
	gfortran $(PROFILE) -o ../foo f.f90 $(MODULES) $(ROUTINES) $(ANALYTICAL) $(LDFLAGS) -J$(BINDIR)

sensitivity : s.f90 $(MODULES) $(ROUTINES) $(ANALYTICAL)
	gfortran $(PROFILE) -o ../bar s.f90 $(MODULES) $(ROUTINES) $(ANALYTICAL) $(LDFLAGS) -J$(BINDIR)

modular_sequential : adaptive_sequential_sampling.f90 $(MODULES) $(ROUTINES) 
	gfortran $(PROFILE) -o ../spam adaptive_sequential_sampling.f90 $(MODULES) $(ROUTINES)

%.o : %.f90
	gfortran $(PROFILE) $(OPT) -c $< $(LDFLAGS) -J$(BINDIR)

clear :
	rm -f *.o ../foo ../bar
	rm -f $(BINDIR)/*.mod

sensitivity-test : sensitivity.o correlation.o regression.o params.o matrix.o
	gfortran -c sensitivity.f90 params.o correlation.o regression.o matrix.o $(LDFLAGS)
	gfortran -o sens ./unit-tests/sensitivity-unittests.f90 sensitivity.o params.o correlation.o regression.o matrix.o $(LDFLAGS) && ./sens
	
