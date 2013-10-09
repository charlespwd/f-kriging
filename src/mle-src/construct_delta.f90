!	delta(1,d) = -0.5*(trace_RDR(d) - (1.0 / sigma2) *...
!		((Y-F*beta)'*Rinv_DR(:,:,d)*Rinv*(Y-F*beta)));
SUBROUTINE construct_delta(delta,Rinv,Rinv_DR,trace_RDR,YmFb,sigma2,Order,D,Ns)
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: Order, D, Ns
	DOUBLE PRECISION, INTENT(IN) :: Rinv(Ns,Ns), Rinv_DR(Ns,Ns,d)  
   DOUBLE PRECISION, INTENT(IN) :: trace_RDR(D,1), YmFb(Ns,1)
   DOUBLE PRECISION, INTENT(IN) :: sigma2
   DOUBLE PRECISION, INTENT(OUT) :: delta(D,1)
   DOUBLE PRECISION ::  tmp(1,ns),tmp2(1,ns), res(1,1)
   INTEGER :: fdim, ii

   fdim = 1+Order*D
   
   DO ii=1,D
      ! tmp = YmFb**T * Rinv_DR(:,:,d) ; * Rinv * YmFb
      ! [ns,1]**T x [Ns,ns] = [1,ns]
      call DGEMM('T','n',1,ns,ns,1.0d0,YmFb,ns,Rinv_DR(:,:,ii),ns,0.0d0,tmp,1)

		! tmp2 = tmp * Rinv
		! [1,ns] x [ns,ns] = [1,ns] 
      call DGEMM('n','n',1,ns,ns,1.0d0,tmp,1,Rinv,ns,0.0d0,tmp2,1)
		
		! res = -1/sigma2 * tmp2 * YmFb + trace_RDR(d)
		! [1,1] = [1,ns] x [ns,1] 
		res(1,1) = trace_RDR(ii,1) 
		call DGEMM('n','n',1,1,ns,((-1.0d0)/sigma2),tmp2,1,YmFb,ns,1.0d0,res,1)
		delta(ii,1) = (-0.5d0) * res(1,1)
   END DO
END SUBROUTINE
