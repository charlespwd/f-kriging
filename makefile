CONSTRUCT=construct_beta.o construct_f.o construct_fmat.o construct_R.o construct_DR.o
FUNCTIONS=rxy.o invertr.o eye.o rescale.o
MODULE=params.o
TEST=dump.o
LDFLAGS=-I/usr/lib/lapack95_modules/ -llapack95 -llapack -lblas 

#OPT = -O1
#PROFILE=-pg #comment out if you don't want profiling

default : f.f90 $(MODULE) $(TEST) $(FUNCTIONS) $(CONSTRUCT)
	gfortran $(PROFILE) -o foo f.f90 $(MODULE) $(TEST) $(FUNCTIONS) $(CONSTRUCT) $(LDFLAGS) 

%.o : %.f90
	gfortran $(PROFILE) $(OPT) -c $< $(LDFLAGS)

clear :
	rm *.o *.mod foo

