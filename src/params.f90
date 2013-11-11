module PARAMS
   double precision, parameter :: TOL = 1D-9 ! theta converged when 1x10-9
   double precision, parameter :: CONDTOL = 1.0d0 * (10.0d0 ** 15)
   double precision, parameter :: Raug = 0.01d0
   integer :: MAXCOUNT = 10 ! Max iteration of MLE
   integer, parameter :: NUGGET = 8 ! Power of the eps added to the diagonal
   integer :: modify_count = 0
   integer, parameter :: Pc = 2    ! power of the correlation (2=gaussian)
   double precision, parameter :: PI = 4.D0*DATAN(1.D0)
end module

