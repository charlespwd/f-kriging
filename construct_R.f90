!!
! This function computes the spatial correlation matrix R vectorially
! 
! inputs:
!  R: the correlation matrix, NxN
!	theta: the correlation factor;
! 	snap_pos: the vector of snapshot positions
!  D: # of dimensions
!  nsnap: # of snapshots
!  Pc: Power of correlation 
!
! outputs: 
! 	R	 
SUBROUTINE construct_R(R,theta,snap_pos,D,nsnap,Pc)
   IMPLICIT NONE
   DOUBLE PRECISION :: rxy
   INTEGER :: D, Nsnap, Pc
   INTEGER :: ii, jj
   DOUBLE PRECISION :: R(nsnap,nsnap),theta(D)
   DOUBLE PRECISION :: snap_pos(nsnap,D)
   DO ii=1,nsnap
      DO jj=ii,nsnap
         R(ii,jj) = rxy(theta,snap_pos(ii,:),snap_pos(jj,:),D,Pc)
         R(jj,ii) = R(ii,jj)
      END DO
   END DO
END SUBROUTINE
