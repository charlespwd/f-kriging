SUBROUTINE construct_kriging_RS(YNEW,XNEW,MSE,Y,X,F,R,theta,D,Ns,NsNew)
   USE PARAMS, ONLY:Order,Pc
   IMPLICIT NONE

   ! ARGUMENTS
   INTEGER, INTENT(IN) :: D,Ns,NsNew
   DOUBLE PRECISION, INTENT(IN) :: X(Ns,D), XNEW(NsNew,D)
   DOUBLE PRECISION, INTENT(IN) :: Y(Ns,1)
   DOUBLE PRECISION, INTENT(IN) :: F(NS,1+Order*D)
   DOUBLE PRECISION, INTENT(IN) :: R(Ns,Ns)
   DOUBLE PRECISION, INTENT(IN) :: theta(D)
   DOUBLE PRECISION, INTENT(OUT) :: YNEW(NsNew,1)
   DOUBLE PRECISION, INTENT(OUT) :: MSE

   ! Functions
   DOUBLE PRECISION :: get_sigma2, get_rxy, get_mse

   ! Work variables
   DOUBLE PRECISION :: Rinv(Ns,Ns), I(Ns,Ns)
   DOUBLE PRECISION :: beta(1+Order*D,1) 
   DOUBLE PRECISION :: YmFb(Ns,1)
   DOUBLE PRECISION :: RinvYmFb(Ns,1)
   DOUBLE PRECISION :: sigma2
   DOUBLE PRECISION :: fx(1+Order*D,1)
   DOUBLE PRECISION :: rx(Ns,1)
   DOUBLE PRECISION :: res(1,1)

   INTEGER :: ii,jj,fdim
   fdim = 1 + Order*D
   
   CALL eye(I,Ns)
   CALL invertR(R,I,Rinv,Ns)
   CALL construct_beta(beta,F,Rinv,Y,Order,D,Ns)

   YmFb = Y
   ! YmFb = Y-F*beta
   call DGEMM('n','n',ns,1,fdim,(-1.0d0),F,Ns,beta,fdim,1.0d0,YmFb,Ns)

   ! RinvYmFb = Rinv * YmFb
   !  (ns,ns) x (ns,1) = (ns,1)
   CALL DGEMM('n','n',ns,1,ns,1.0d0,Rinv,Ns,YmFb,ns,0.0d0,RinvYmFb,ns)
   sigma2 = get_sigma2(YmFb,Rinv,Order,D,Ns)
   
   DO ii=1,NsNew
      ! construct f(x)
      call construct_f(fx(:,1),XNEW(ii,:),Order,D,Ns)
      
      ! construct r(x)
      DO jj=1,ns
         rx(jj,1) = get_rxy(theta,xnew(ii,:),x(jj,:),D,Pc)
      END DO
      
      ! YNEW(ii,1) = fx'*beta + rx'*Rinv*YmFb
      !  first res = fx'*beta
      CALL DGEMM('t','n',1,1,fdim,1.0d0,fx,fdim,beta,fdim,0.0d0,res,1)

      !  then res = rx'*RinvYmFb + res
      CALL DGEMM('t','n',1,1,ns,1.0d0,rx,Ns,RinvYmFb,Ns,1.0d0,res,1)
      YNEW(ii,1) = res(1,1)
   END DO 


END SUBROUTINE
