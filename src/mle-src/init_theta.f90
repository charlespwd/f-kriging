! This routine initializes theta if it isn't already and 
! initializes the bounds on theta during the MLE
!
! Arguments:
!  theta(inout) : in : the initvalue of theta, default should be < 0
!                 out: the initvalue of theta, if not < 0 on input, unchanged.
!  bounds(out) :  [D,2], the first column contains the minimum thresholds per
!                 dimension, the second column contains the maximums
!  X : the set of snapshot positions, [Ns,D]
!  D : the number of dimensions of the domain
!  Ns : the number of snapshots
SUBROUTINE init_theta(theta,bounds,X,D,Ns)
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: D, Ns
   DOUBLE PRECISION, INTENT(IN) :: X(Ns,D)
   DOUBLE PRECISION, INTENT(INOUT) :: theta(D)
   DOUBLE PRECISION, INTENT(OUT) :: bounds(D,2)
   DOUBLE PRECISION :: xmax(D), xmin(D)
   INTEGER :: I_MIN = 1, I_MAX = 2
   INTEGER :: ii
   xmax = maxval(X,1)
   xmin = minval(X,1) 
   DO ii=1,D
      IF (theta(ii) < 0) THEN
         theta(ii) = log(5.0d0) / (xmax(ii) - xmin(ii))
      END IF
      bounds(ii,I_MIN) = (1.0d0/1000.0d0) * theta(ii)
      bounds(ii,I_MAX) = 1000.0d0 * theta(ii)
   END DO
END SUBROUTINE
