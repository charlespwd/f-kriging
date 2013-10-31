! this file contains:
!  construct_kriging_rs
!  get_mse

! construct_kriging_rs
!  this routine makes a response surface assuming the theta is fixed.
!  it also calculates the MSE error at XNEW locations.
subroutine construct_kriging_RS(YNEW,XNEW,MSE,Y,X,F,R,theta,Order,D,Ns,NsNew)
   use PARAMS, only:Pc
   use matrix, only: eye
   use correlation, only: invertr,get_rxy
   use regression, only: construct_f
   use mle, only : construct_beta, get_sigma2
   use error
   implicit none

   ! ARGUMENTS
   integer, intent(in) :: D, Ns, NsNew, Order
   double precision, intent(in) :: X(Ns,D), XNEW(NsNew,D)
   double precision, intent(in) :: Y(Ns,1)
   double precision, intent(in) :: F(NS,1+Order*D)
   double precision, intent(INOUT) :: R(Ns,Ns)
   double precision, intent(in) :: theta(D)
   double precision, intent(out) :: YNEW(NsNew,1)
   double precision, intent(out) :: MSE(nsnew)

   ! Work variables
   double precision :: Rinv(Ns,Ns), I(Ns,Ns)
   double precision :: beta(1+Order*D,1) 
   double precision :: YmFb(Ns,1)
   double precision :: RinvYmFb(Ns,1)
   double precision :: sigma2
   double precision :: fx(1+Order*D,1)
   double precision :: rx(Ns,1)
   double precision :: res(1,1)

   integer :: ii,jj,fdim
   fdim = 1 + Order*D
   
   call eye(I,Ns)
   call invertR(R,I,Rinv,Ns)
   call construct_beta(beta,F,Rinv,Y,Order,D,Ns)

   YmFb = Y
   ! YmFb = Y-F*beta
   call DGEMM('n','n',ns,1,fdim,(-1.0d0),F,Ns,beta,fdim,1.0d0,YmFb,Ns)

   ! RinvYmFb = Rinv * YmFb
   !  (ns,ns) x (ns,1) = (ns,1)
   call DGEMM('n','n',ns,1,ns,1.0d0,Rinv,Ns,YmFb,ns,0.0d0,RinvYmFb,ns)
   sigma2 = get_sigma2(YmFb,Rinv,Order,D,Ns)
   
   do ii=1,NsNew
      ! construct f(x)
      call construct_f(fx(:,1),XNEW(ii,:),Order,D,Ns)
      
      ! construct r(x)
      do jj=1,ns
         rx(jj,1) = get_rxy(theta,xnew(ii,:),x(jj,:),D,Pc)
      end do
      
      ! YNEW(ii,1) = fx'*beta + rx'*Rinv*YmFb
      !  first res = fx'*beta
      call DGEMM('t','n',1,1,fdim,1.0d0,fx,fdim,beta,fdim,0.0d0,res,1)

      !  then res = rx'*RinvYmFb + res
      call DGEMM('t','n',1,1,ns,1.0d0,rx,Ns,RinvYmFb,Ns,1.0d0,res,1)
      YNEW(ii,1) = res(1,1)

      MSE(ii) = lophaven_mse(sigma2, fx, rx, F, R, Rinv, ns, fdim)
   end do 
end subroutine

