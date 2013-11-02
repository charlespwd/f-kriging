ROUTINES = kriging.o cokriging.o sequential_sampling.o 
MODULES = params.o grid.o matrix.o regression.o correlation.o \
			 sensitivity.o mle.o error.o utils.o optimization.o linesearch.o 
ANALYTICAL = analytical_functions.o analytical_solver.o
LDFLAGS=-I/usr/lib/lapack95_modules/ -I/usr/lib/fortarrange_header/ -lfort_arrange -llapack95 -llapack -lblas 
BINDIR=../bin

OPT = -O3 
#PROFILE=-pg #comment out if you don't want profiling

default : f.f90 $(MODULES) $(ROUTINES) $(ANALYTICAL)
	gfortran $(OPT) $(PROFILE) -o ../foo f.f90 $(MODULES) $(ROUTINES) $(ANALYTICAL) $(LDFLAGS) -J$(BINDIR)

sensitivity : s.f90 $(MODULES) $(ROUTINES) $(ANALYTICAL)
	gfortran $(PROFILE) $(OPT) -o ../bar s.f90 $(MODULES) $(ROUTINES) $(ANALYTICAL) $(LDFLAGS) -J$(BINDIR)

modular : modular-test.f90 $(MODULES) $(ROUTINES) $(ANALYTICAL) 
	gfortran $(PROFILE) $(OPT) -o ../spam modular-test.f90 $(MODULES) $(ROUTINES) $(ANALYTICAL) $(LDFLAGS) -J$(BINDIR)

linesearch.o : fortLinesearch/linesearch.f90
	gfortran $(PROFILE) $(OPT) -c fortLinesearch/linesearch.f90 $(LDFLAGS) -J$(BINDIR)

%.o : %.f90
	gfortran $(PROFILE) $(OPT) -c $< $(LDFLAGS) -J$(BINDIR)

clear :
	rm -f *.o ../foo ../bar
	rm -f $(BINDIR)/*.mod

sensitivity-test : sensitivity.o correlation.o regression.o params.o matrix.o
	gfortran -c sensitivity.f90 params.o correlation.o regression.o matrix.o $(LDFLAGS)
	gfortran -o sens ./unit-tests/sensitivity-unittests.f90 sensitivity.o params.o correlation.o regression.o matrix.o $(LDFLAGS) && ./sens
	
