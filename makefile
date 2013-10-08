ROUTINES=optimize_theta_mle.o kriging.o cokriging.o analytical_solver.f90
CONSTRUCT=construct_beta.o construct_delta.o construct_DR.o construct_kriging_RS.o
FUNCTIONS= get_trace.o get_sigma2.o init_theta.o get_mse.o
MODULE=params.o analytical_functions.o grid.o utils.o matrix.o regression.o correlation.o sensitivity.o
LDFLAGS=-I/usr/lib/lapack95_modules/ -llapack95 -llapack -lblas 

OPT = -O3
#PROFILE=-pg #comment out if you don't want profiling

default : f.f90 $(MODULE) $(FUNCTIONS) $(CONSTRUCT) $(ROUTINES)
	gfortran $(PROFILE) -o foo f.f90 $(MODULE) $(TEST) $(FUNCTIONS) $(CONSTRUCT) $(ROUTINES) $(LDFLAGS) 

%.o : %.f90
	gfortran $(PROFILE) $(OPT) -c $< $(LDFLAGS)

clear :
	rm *.o *.mod foo

sensitivity : sensitivity.o
	gfortran -c sensitivity.f90 params.o correlation.o regression.o matrix.o $(LDFLAGS)
	gfortran -o sens ./unit-tests/sensitivity-unittests.f90 sensitivity.o params.o correlation.o regression.o matrix.o $(LDFLAGS) && ./sens
	
