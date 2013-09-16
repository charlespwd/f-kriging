!!
! This function computes the spatial correlation matrix R vectorially
! 
! inputs:
!  R: the correlation matrix, NxN
!	theta: the correlation factor;
! 	X: the vector of snapshot positions
!  D: # of dimensions
!  Ns: # of snapshots
!  Pc: Power of correlation 
!
! outputs: 
! 	R	 
SUBROUTINE construct_R(R,theta,X,D,Ns,Pc)
   IMPLICIT NONE
   DOUBLE PRECISION :: get_rxy
   INTEGER :: D, Ns, Pc
   INTEGER :: ii, jj, dd
   DOUBLE PRECISION :: R(Ns,Ns),theta(D)
   DOUBLE PRECISION :: X(Ns,D)
   DOUBLE PRECISION :: tmp
   DO ii=1,Ns
      DO jj=ii,Ns
         ! inlined for performance, this gets called so much that the overhead
         ! isn't worth the better readability. I left the call there
         ! instead
!         R(ii,jj) = get_rxy(theta,X(ii,:),X(jj,:),D,Pc)
         tmp = 0.0d0
         DO dd=1,D
            tmp = tmp &
               - theta(dd) * (abs(X(ii,dd) - X(jj,dd)) ** Pc)
         END DO
         R(ii,jj) = exp(tmp)
         R(jj,ii) = R(ii,jj)
      END DO
   END DO
END SUBROUTINE
