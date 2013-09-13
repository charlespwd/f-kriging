CONSTRUCT=construct_f.o construct_fmat.o construct_R.o construct_DR.o
FUNCTIONS=rxy.o 
MODULE=params.o
TEST=dump.o
LDFLAGS=-I/usr/lib/lapack95_modules/ -llapack95 -llapack -lblas 

PROFILE=-pg #comment out if you don't want profiling

default : f.f90 $(CONSTRUCT) $(FUNCTIONS) $(MODULE) $(TEST)
	gfortran $(PROFILE) -o foo f.f90 $(TEST) $(FUNCTIONS) $(CONSTRUCT) $(LDFLAGS) 

%.o : %.f90
	gfortran $(PROFILE) -O -c $<

%.mod : %.f90
	gfortran $(PROFILE) -O -c $<
