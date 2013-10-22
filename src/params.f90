MODULE PARAMS
   DOUBLE PRECISION, PARAMETER :: TOL = 1D-9 ! theta converged when 1x10-9
   DOUBLE PRECISION, PARAMETER :: CONDTOL = 1.0d0 * (10.0d0 ** 15)
   DOUBLE PRECISION, PARAMETER :: Raug = 0.01d0
   INTEGER, PARAMETER :: MAXCOUNT = 10 ! Max iteration of MLE
   INTEGER, PARAMETER :: NUGGET = 8 ! Power of the eps added to the diagonal
   INTEGER, PARAMETER :: ORDER = 2 ! order of the regression 
   INTEGER, PARAMETER :: Pc = 2    ! power of the correlation (2=gaussian)
   DOUBLE PRECISION, PARAMETER :: PI = 4.D0*DATAN(1.D0)
   INTEGER :: MSEMODE
END MODULE

