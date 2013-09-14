!!
! This function computes the derrivative of the spatial correlation matrix 
!
! inputs:
!	R: the correlation matrix;
! 	snap_pos: the snapshot positions;
!  D: dimension of domain
!  Ns: # of snapshots
!
! outputs: 
! 	DR(x,y,d);	 
SUBROUTINE construct_DR(DR,R,snap_pos,D,Ns)
   USE PARAMS
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: Ns, D
   DOUBLE PRECISION, INTENT(OUT) :: DR(Ns,Ns,D)
   DOUBLE PRECISION, INTENT(IN) :: R(Ns,Ns)
   DOUBLE PRECISION, INTENT(IN) :: snap_pos(Ns,D)
   INTEGER :: ii, jj, kk
   
   DO ii=1,Ns
      DO jj=ii,Ns
         DO kk=1,D
            DR(ii,jj,kk) = 0 - R(ii,jj) * ((abs(snap_pos(jj,kk) - snap_pos(ii,kk))) ** Pc)
            DR(jj,ii,kk) = DR(ii,jj,kk)
         END DO
      END DO
   END DO
END SUBROUTINE
