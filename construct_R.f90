!!
! This function computes the spatial correlation matrix R vectorially
! 
! inputs:
!  R: the correlation matrix, NxN
!	theta: the correlation factor;
! 	snap_pos: the vector of snapshot positions
!  D: # of dimensions
!  Ns: # of snapshots
!  Pc: Power of correlation 
!
! outputs: 
! 	R	 
SUBROUTINE construct_R(R,theta,snap_pos,D,Ns,Pc)
   IMPLICIT NONE
   DOUBLE PRECISION :: rxy
   INTEGER :: D, Ns, Pc
   INTEGER :: ii, jj
   DOUBLE PRECISION :: R(Ns,Ns),theta(D)
   DOUBLE PRECISION :: snap_pos(Ns,D)
   DO ii=1,Ns
      DO jj=ii,Ns
         R(ii,jj) = rxy(theta,snap_pos(ii,:),snap_pos(jj,:),D,Pc)
         R(jj,ii) = R(ii,jj)
      END DO
   END DO
END SUBROUTINE
