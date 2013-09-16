! computes the trace of a square matrix
! Arguments
!  A: square matrix, mxm
!  m: size of matrix
DOUBLE PRECISION FUNCTION get_trace(A,m)
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: m
   DOUBLE PRECISION, INTENT(IN) :: A(m,m)
   INTEGER :: ii
   get_trace = 0
   DO ii=1,m
      get_trace = get_trace + A(ii,ii)
   END DO
END FUNCTION
