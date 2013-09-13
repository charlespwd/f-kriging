!! 
! nugget correct, this function is CALLed by inverse R
! When a matrix is close to singular, a small epsilon is added to the 
! diagonal of the correlation matrix to permit it to solve
SUBROUTINE nuggetcorrect(R,Rinv,I,NS,cond)
   USE PARAMS,ONLY:NUGGET,CONDTOL
   USE LA_PRECISION,ONLY:WP=>DP
   USE F95_LAPACK,ONLY:LA_GESVX
   IMPLICIT NONE
   INTEGER,INTENT(IN) :: NS
   DOUBLE PRECISION,INTENT(inout) :: Rinv(NS,NS), R(NS,NS) 
   DOUBLE PRECISION,INTENT(in) :: I(NS,NS)  
   DOUBLE PRECISION :: wI(NS,NS),wR(NS,NS)
   INTEGER :: ii,info
   DOUBLE PRECISION :: eps, cond,rcon

! could possibly optimize this by giving the factored L matrix since
! you need to calculate LA_GESVX to get the condition number. 
!   CALL LA_GESVX(R,I,Rinv,RCOND=rcon)
!   cond = rcon ** (-1) 
   IF (cond > condtol) THEN 
      WRITE(*,'(a25,E11.5)') 'applying nugget to cond= ',cond
      eps = (10 + NS) * 10D0 ** (-NUGGET) 
      DO ii=1,NS
        R(ii,ii) = R(ii,ii) + eps
      END DO

      ! we do not want R or I to be modified by the LAPACK subroutines
      wR = R
      wI = I

      CALL LA_GESVX(wR,wI,Rinv,RCOND=rcon,INFO=info) 
      WRITE(*,'(a25,E11.5)') 'cond after nugget= ', rcon**(-1)
   END IF
END SUBROUTINE
