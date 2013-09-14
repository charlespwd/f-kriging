MODULE PARAMS
   DOUBLE PRECISION, PARAMETER :: TOL = 1D-9 ! theta converged when 1x10-9
   DOUBLE PRECISION, PARAMETER :: CONDTOL = 1D15 
   INTEGER, PARAMETER :: MAXCOUNT = 400 ! Max iteration of MLE
   INTEGER, PARAMETER :: NUGGET = 9 ! Power of the eps added to the diagonal
   INTEGER, PARAMETER :: ORDER = 0 ! order of the regression 
   INTEGER, PARAMETER :: Pc = 2    ! power of the correlation (2=gaussian)
END MODULE

