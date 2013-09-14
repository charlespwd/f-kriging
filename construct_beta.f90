! beta = (F'*Rinv*F) \ (F'*Rinv*Y);
SUBROUTINE construct_beta(beta,F,Rinv,Y,ORDER,D,Ns)
   USE LA_PRECISION, ONLY:WP=>DP
   USE F95_LAPACK, ONLY:LA_GESV,LA_GESVX
   implicit none 

   ! arguments
   INTEGER, INTENT(IN) :: D, Ns, Order
   DOUBLE PRECISION, INTENT(IN) :: F(Ns,1+ORDER*D),Rinv(Ns,Ns),Y(Ns,1)
   DOUBLE PRECISION, INTENT(OUT) :: beta(1+Order*D,1)
   
   ! work variables
   INTEGER :: Fdim
   DOUBLE PRECISION :: FTRinv(1+Order*D,Ns) 
   DOUBLE PRECISION :: LHS(1+Order*D,1+Order*D), RHS(1+Order*D,1)
   
   Fdim = 1+Order*D
   
   ! from BLAS
   ! F' * Rinv, [Ns * Fdim] ** T x [Ns x Ns] = [Fdim x Ns]
   CALL DGEMM('T','N',Fdim,Ns,Ns,1.0D0,F,Ns,Rinv,Ns,0.0D0,FTRinv,Fdim)

   ! LHS = F'*Rinv*F, [Fdim x Ns] * [Ns x Fdim] = [Fdim x Fdim]
   CALL DGEMM('N','N',Fdim,Fdim,Ns,1.0D0,FTRinv,Fdim,F,Ns,0.0D0,LHS,Fdim)

   ! RHS = F'*Rinv * Y, [Fdim x Ns]*[Ns x 1] = [Fdim x 1]
   CALL DGEMM('N','N',Fdim,1,Ns,1.0D0,FTRinv,Fdim,Y,Ns,0.0D0,RHS,Fdim)

   ! LHS \ RHS
   CALL LA_GESV(LHS,RHS)
   beta = RHS
   
   !beta = RHS
END SUBROUTINE
