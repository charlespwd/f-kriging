SUBROUTINE nuggetcorrect(R,NSnap)
   USE PARAMS
   IMPLICIT NONE
   INTEGER,INTENT(IN) :: NSnap
   DOUBLE PRECISION,INTENT(inout) :: R(NSnap,NSnap)
   INTEGER :: ii
   DOUBLE PRECISION :: eps, cond
   
   ! estimate the condition number;
   
   IF (cond > condtol) 
      eps = (10 + nsnap) * 10D1 ** (-NUGGET) 
      DO ii=1:NSnap
         R(ii,ii) = R(ii,ii) + eps
      END DO
   END IF
END SUBROUTINE
