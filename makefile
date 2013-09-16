ROUTINES=optimize_theta_mle.o kriging.o cokriging.o
CONSTRUCT=construct_beta.o construct_delta.o construct_f.o construct_fmat.o construct_R.o construct_DR.o construct_kriging_RS.o
FUNCTIONS=get_rxy.o invertr.o eye.o rescale.o get_trace.o get_sigma2.o init_theta.o
MODULE=params.o
TEST=dump.o
LDFLAGS=-I/usr/lib/lapack95_modules/ -llapack95 -llapack -lblas 

#OPT = -O1
PROFILE=-pg #comment out if you don't want profiling

default : f.f90 $(MODULE) $(TEST) $(FUNCTIONS) $(CONSTRUCT) $(ROUTINES)
	gfortran $(PROFILE) -o foo f.f90 $(MODULE) $(TEST) $(FUNCTIONS) $(CONSTRUCT) $(ROUTINES) $(LDFLAGS) 

%.o : %.f90
	gfortran $(PROFILE) $(OPT) -c $< $(LDFLAGS)

clear :
	rm *.o *.mod foo

