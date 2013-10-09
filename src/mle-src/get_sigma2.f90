!sigma2 = (Y-F*beta)' * Rinv * (Y - F*beta) / nsnap;
DOUBLE PRECISION FUNCTION get_sigma2(YmFb,Rinv,Order,D,Ns)
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: Order, D, Ns
   DOUBLE PRECISION, INTENT(IN) :: YmFb(Ns,1), Rinv(Ns,Ns)
   DOUBLE PRECISION ::  tmp(1,ns), res(1,1)
   INTEGER :: fdim

   fdim = 1+Order*D
   
   ! tmp = (Y-F*beta)**T * Rinv
   ! [1,ns] =[ns,1]**T x [ns,ns] 
   call DGEMM('T','n',1,ns,ns,1.0d0,YmFb,ns,Rinv,ns,0.0d0,tmp,1)

   ! res = 1/ns * tmp * YmFb
   ! [1,1] = [1,ns]x[ns,1] 
   call DGEMM('n','n',1,1,ns,(1.0d0/Ns),tmp,1,YmFb,ns,0.0d0,res,1)

   get_sigma2 = res(1,1)
END FUNCTION


   
