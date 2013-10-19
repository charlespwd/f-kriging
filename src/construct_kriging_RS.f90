! this file contains:
!  construct_kriging_rs
!  get_mse

! construct_kriging_rs
!  this routine makes a response surface assuming the theta is fixed.
!  it also calculates the MSE error at XNEW locations.
SUBROUTINE construct_kriging_RS(YNEW,XNEW,MSE,Y,X,F,R,theta,Order,D,Ns,NsNew)
   USE PARAMS, ONLY:Pc
   USE matrix, only: eye
   use correlation, only: invertr,get_rxy
   use regression, only: construct_f
   use mle, only : construct_beta, get_sigma2
   IMPLICIT NONE

   ! ARGUMENTS
   INTEGER, INTENT(IN) :: D, Ns, NsNew, Order
   DOUBLE PRECISION, INTENT(IN) :: X(Ns,D), XNEW(NsNew,D)
   DOUBLE PRECISION, INTENT(IN) :: Y(Ns,1)
   DOUBLE PRECISION, INTENT(IN) :: F(NS,1+Order*D)
   DOUBLE PRECISION, INTENT(INOUT) :: R(Ns,Ns)
   DOUBLE PRECISION, INTENT(IN) :: theta(D)
   DOUBLE PRECISION, INTENT(OUT) :: YNEW(NsNew,1)
   DOUBLE PRECISION, INTENT(OUT) :: MSE(nsnew)

   ! Functions
   DOUBLE PRECISION :: martinmse, get_mse

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

      MSE(ii) = martinmse(sigma2, fx, rx, F, R, ns, fdim)
   END DO 
END SUBROUTINE

!! 
! get_mse
! computed as abs(sigma2*(
! 1 - (r**T * Rinv * r) 
! + (1 - || F**T * Rinv * r ||_2 ^2) / || F**T * Rinv * F||_2 ))
double precision function get_mse(sigma2,F,Rinv,rx,D,Ns,fdim)
   use matrix, only : spectral_norm
   implicit none
   
   ! inputs
   integer, intent(in) :: D, Ns
   integer, intent(in) :: fdim
   double precision, intent(in) :: sigma2
   double precision, intent(in) :: F(Ns,fdim)
   double precision, intent(in) :: Rinv(Ns,Ns)
   double precision, intent(in) :: rx(Ns,1)
   
   ! work variables
   double precision :: rtRinv(1,Ns), rtRinvr(1,1), FtRinvR(fdim,1)
   double precision :: FtRinv(fdim,Ns), FtRinvF(fdim,fdim)
   double precision :: term1, term2, term3, term4
   integer :: foo=0

   !get RtRinvr
   ! r**T * Rinv
   ! [1,ns] = [ns,1]**T * [ns,ns]
   call dgemm('t','n',1,ns,ns,1.0d0,rx,ns,Rinv,ns,0.0d0,rtRinv,1)

   ! rtRinv * rx
   ! [1,1] = [1,ns] [ns,1]
   call dgemm('n','n',1,1,ns,1.0d0,rtRinv,1,rx,ns,0.0d0,rtRinvr,1)
   
   !get F'RinvF
   ! F**T * Rinv
   ! [fdim, Ns] = [ns,fdim] ** t * [ns,ns]
   call dgemm('t','n',fdim,ns,ns,1.0d0,F,ns,Rinv,ns,0.0d0,FtRinv,fdim)

   ! FtRinv * F
   ! [fdim,fdim] = [fdim,ns] [ns,fdim]
   call dgemm('n','n',fdim,fdim,ns,1.0d0,FtRinv,fdim,F,ns,0.0d0,FtRinvF,fdim)

   ! FtRinv * rx
   ! [fdim,1] = [fdim,ns] * [ns,1]
   call dgemm('n','n',fdim,1,ns,1.0d0,FtRinv,fdim,rx,ns,0.0d0,FtRinvR,fdim)

   ! CALCULATE MSE..........
   term1 = 1.0d0
   term2 = rtRinvr(1,1)
   term3 = (1.0d0 - spectral_norm(FtRinvR,fdim,1)) ** 2
   term4 = 1.0d0 / spectral_norm(FtRinvF,fdim,fdim)
   term3 = term3 * term4
   get_mse = abs(sigma2 * (term1 - term2 + term3))
end function 

 ! sigma2 * (1 - [ft(x); rt(x)]**T * [ 0 Ft ; F R] * [f(x); r(x)]) 
double precision function martinmse(sigma2, fx, rx, F, R, ns, fdim)
   implicit none
   integer, intent(in) :: fdim, ns
   double precision, intent(in) :: sigma2
   double precision, intent(in) :: fx(fdim,1), rx(ns,1)
   double precision, intent(in) :: F(ns,fdim), R(ns,ns)
   ! work variables
   double precision :: Ft(fdim,ns)
   double precision :: fxrx(fdim+ns,1)
   double precision :: ctrblock(fdim+ns,fdim+ns)
   double precision :: tmp(1,fdim+ns)
   double precision :: res(1,1)
   integer :: ii, jj
   fxrx(1:fdim,1) = fx(1:fdim,1)
   fxrx(fdim+1:ns,1) = rx(1:ns,1)
   Ft = transpose(F)

   ! construct [0 Ft; F R]
   !  0
   ctrblock(1:fdim,1:fdim) = 0.d0
   !  F
   ctrblock(fdim+1:fdim+ns, 1:fdim) = F(1:ns,1:fdim)
   !  Ft 
   ctrblock(1:fdim, fdim+1:fdim+ns) = Ft(1:fdim,1:ns)
   !  R
   ctrblock(fdim+1:fdim+ns, fdim+1:fdim+ns) = R(1:ns,1:ns)

   ! tmp = fxrx**T * cntrblock 
   ! [1,fdim+ns] = [fdim+ns,1] ** T * [fdim+ns,fdim+ns]
   call dgemm('t','n',1,fdim+ns,fdim+ns,1.0d0, fxrx,fdim+ns, &
      ctrblock,fdim+ns, 0.0d0,tmp,1)

   ! res = (-1.0d0) * tmp * fxrx + 1.0d0 * res
   res(1,1) = 1.0d0
   call dgemm('n','n',1,1,fdim+ns,(-1.0d0), tmp,1, fxrx,fdim+ns, &
      1.0d0, res,1)
   
   martinmse = sigma2 * res(1,1)
end function
