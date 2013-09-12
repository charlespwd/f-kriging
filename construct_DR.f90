!!
! This function computes the derrivative of the spatial correlation matrix 
!
! inputs:
!	R: the correlation matrix;
! 	snap_pos: the snapshot positions;
!  D: dimension of domain
!  nsnap: # of snapshots
!
! outputs: 
! 	DR(x,y,d);	 
SUBROUTINE construct_DR(DR,R,snap_pos,D,nsnap)
   USE PARAMS
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: nsnap, D
   DOUBLE PRECISION, INTENT(OUT) :: DR(nsnap,nsnap,D)
   DOUBLE PRECISION, INTENT(IN) :: R(nsnap,nsnap)
   DOUBLE PRECISION, INTENT(IN) :: snap_pos(nsnap,D)
   INTEGER :: ii, jj, kk
   
   DO ii=1,nsnap
      DO jj=ii,nsnap
         DO kk=1,D
            DR(ii,jj,kk) = 0 - R(ii,jj) * ((abs(snap_pos(jj,kk) - snap_pos(ii,kk))) ** Pc)
            DR(jj,ii,kk) = DR(ii,jj,kk)
         END DO
      END DO
   END DO
END SUBROUTINE
