!!
! this function creates computes the regression mat
! Right now, it's [1, x_1^1, x_2^1, x_1^2, x_2^2 ..
!
! input:
!  x: location of evaluation (f(x))
!
! output:
!  f(x): [1, x_1^1, ...]
SUBROUTINE construct_f(f,x,Order,D,nsnap)
   IMPLICIT NONE
   INTEGER :: Order, D, nsnap
   INTEGER :: oo,dd,i
   DOUBLE PRECISION :: f(1 + Order*D)
   DOUBLE PRECISION :: x(D)
   i = 2;
   f(1) = 1;
   DO oo=1,Order
      DO dd=1,D
         f(i) = x(dd) ** oo
         i = i + 1
      END DO
   END DO
END SUBROUTINE
   
