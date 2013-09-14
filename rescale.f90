!!
! This routine rescales the vector x in a range 
! of distance 1. As suggested in kriging
!
! ARGUMENTS:
!  X (inout) : the vector to be scaled
!  D (in) : # of dimensions
!  Ns (in) : # of snapshots
!  Xmax (in) : array of maximums per column of X (i.e. per dimension)
!  Xmin (in) : array of min per column of X (i.e. per dimension)
!  Xout (out), optional : if present, x isn't overwritten, the scaled
SUBROUTINE rescale(X, D, Ns, XMAX, XMIN)
   INTEGER, INTENT(IN) :: Ns, D
   DOUBLE PRECISION, INTENT(IN) :: XMAX(D), XMIN(D)
   DOUBLE PRECISION, INTENT(INOUT) :: X(Ns,D)
   INTEGER :: dd

   DO dd=1,D
      X(:,dd) = X(:,dd) / (XMAX(dd) - XMIN(dd)) 
   END DO
END SUBROUTINE

