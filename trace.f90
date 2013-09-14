! computes the trace of a square matrix
! Arguments
!  A: square matrix, mxm
!  m: size of matrix
DOUBLE PRECISION FUNCTION trace(A,m)
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: m
   DOUBLE PRECISION, INTENT(IN) :: A(m,m)
   INTEGER :: ii
   trace = 0
   DO ii=1,m
      trace = trace + A(ii,ii)
   END DO
END FUNCTION
