!!
! A routine to inverse the correlation matrix
! it is semi intelligent as it corrects it by calling nuggetcorrect
! in the event the condition number is too large. 
! 
! ARGUMENTS:
!  R (inout) : the correlation matrix, gets modified if nuggetcorrected.
!  I (in) : NsxNs, the identity matrix
!  Rinv (out) : the inverse of the correlation matrix, ns x ns
!  Ns : number of snapshots, i.e. dimension of matrix;
SUBROUTINE invertR(R,I,Rinv,Ns)
   USE LA_PRECISION,ONLY:WP=>DP
   USE F95_LAPACK,ONLY:LA_GESVX
   IMPLICIT NONE 
   INTEGER, INTENT(IN) :: Ns
   DOUBLE PRECISION, INTENT(IN) :: I(Ns,Ns)
   DOUBLE PRECISION, INTENT(INOUT) :: R(Ns,Ns)
   DOUBLE PRECISION, INTENT(OUT) :: Rinv(Ns,Ns)
   DOUBLE PRECISION :: wR(Ns,Ns), wI(Ns,Ns)
   DOUBLE PRECISION :: RCOND
   INTEGER :: INFO

   wR = R ! we do not want R and I to be modified by LAPACK
   wI = I

   CALL LA_GESVX(wR,wI,Rinv,RCOND=RCOND,INFO=INFO)
   CALL NUGGETCORRECT(R,I,Rinv,Ns,RCOND ** (-1))
END SUBROUTINE

!! 
! nugget correct, this function is CALLed by inverse R
! When a matrix is close to singular, a small epsilon is added to the 
! diagonal of the correlation matrix to permit it to solve
SUBROUTINE nuggetcorrect(R,I,Rinv,Ns,cond)
   USE PARAMS,ONLY:NUGGET,CONDTOL
   USE LA_PRECISION,ONLY:WP=>DP
   USE F95_LAPACK,ONLY:LA_GESVX
   IMPLICIT NONE
   INTEGER,INTENT(IN) :: Ns
   DOUBLE PRECISION,INTENT(inout) :: Rinv(Ns,Ns), R(Ns,Ns) 
   DOUBLE PRECISION,INTENT(in) :: I(Ns,Ns) 
   DOUBLE PRECISION :: wI(Ns,Ns),wR(Ns,Ns)
   INTEGER :: ii,info
   DOUBLE PRECISION :: eps, cond,rcon

! could possibly optimize this by giving the factored L matrix since
! you need to calculate LA_GESVX to get the condition number. 
!   CALL LA_GESVX(R,I,Rinv,RCOND=rcon)
!   cond = rcon ** (-1) 
   IF (cond > condtol) THEN 
      WRITE(*,'(a25,E11.5)') 'applying nugget to cond= ',cond
      eps = (10 + Ns) * 10D0 ** (-NUGGET) 
      DO ii=1,Ns
        R(ii,ii) = R(ii,ii) + eps
      END DO

      ! we do not want R or I to be modified by the LAPACK subroutines
      wR = R
      wI = I

      CALL LA_GESVX(wR,wI,Rinv,RCOND=rcon,INFO=info) 
      WRITE(*,'(a25,E11.5)') 'cond after nugget= ', rcon**(-1)
   END IF
END SUBROUTINE

