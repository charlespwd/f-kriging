! !! KRIGING(XNEW,YNEW,theta,MSE,XOLD,Y,D,Ns,NewNs)
! This subroutine builds a response surface using kriging and an iterative
!  MLE to obtain theta. 

! WARNING: FOR BETTER PERFORMANCE, XNEW and X SHOULD ALWAYS BE PRECONDITIONED
! VIA THE PRECONDITION PROCEDURE IN THIS MODULE.
!
! Arguments:
!  XNEW : the new set of locations where the function should be evaluated
!  YNEW : the kriging-estimated values
!  theta: the correlation hyper-parameter, most of the time is consumed finding
!         the correct value of this parameter. With him you can do everything
!  MSE:   the built-in error estimator
!  XOLD:  the set of original snapshot positions
!  Y:     the set of original snapshot values
!  Order: the order of the polynomial regression used. 
!  D:     the # of dimesions in the domain
!  Ns:    the # of snapshots
!  NewNs: the # of new snapshot locations
!
!  Kriging boils down to : y(x) = P(x) + Z(x) 
!  P(x) is a polynomial regression
!  Z(x) is a statistical measure
!  
!  Here, we calculate y(x) for all xNEW and store this in YNEW
!  P(x) = beta * f(x), 
!     where f(x) is a vector function corresponding to the chosen
!        regression Order. (This parameter is set in the params module), see
!        construct_f for more details about the implementation
!     where beta is calculated from (F'*Rinv*F) \ (F'*Rinv*Y)
!     where F corresponds to a matrix built with f(XNEW) 
!     where R is a correlation matrix, Rij = exp(-SUM(theta .* abs(x_i-x_j)))
!  Z(x) = r'(x)*Rinv*(Y-F*beta) 
!     where r_i(x) = exp(-SUM(theta .* abs(x-x_i)))
!
!  And that's all there is to it, getting the theta is a bit complicated but you
!   are invited to read the paper by Lappo for the details on the iterative MLE
!   solver. See optimize_theta_mle for implementation. 
module kriging_module
implicit none
   contains

      subroutine precondition(x,xn,xold,xnew,D,ns,nsnew)
         use matrix, only: rescale
         integer, intent(in) :: D, ns, nsnew
         double precision, intent(in) :: XOLD(ns,d), XNEW(nsnew,D)
         double precision, intent(out) :: X(Ns,D), XN(nsnew,D) ! scaled set of snapshots
         double precision :: xmin(D), xmax(D)
         XMIN = minval(XOLD,1)
         XMAX = maxval(XOLD,1)
         X = XOLD
         XN = XNEW
         call rescale(X,D,Ns,XMIN,XMAX)  ! Precondition X by "normalizing" it
         call rescale(XN,D,Nsnew,XMIN,XMAX) ! " "
      end subroutine

      subroutine KRIGING(XNEW,YNEW,theta,MSE,X,Y,Order,D,Ns,NewNs)
         use PARAMS, only: Pc
         use regression, only: construct_fmat
         use correlation, only : construct_r
         use mle, only : optimize_theta_mle, init_theta
         implicit none
         
         ! ARGUMENTS
         integer, intent(in) :: D, Ns, NewNs, Order
         double precision, intent(in) :: X(Ns,D)      ! unscaled set of snapshots
         double precision, intent(in) :: XNEW(NewNs,D)   ! unscaled set of untried loc
         double precision, intent(in) :: Y(Ns,1)         ! set of snapshot values
         double precision, intent(INOUT) :: theta(D)     ! correlation parameter
         double precision, intent(out) :: YNEW(NewNs,1)  ! set of untried values
         double precision, intent(out) :: MSE(NewNs)     ! error estimation

         ! Work variables
         double precision :: F(Ns,1+Order*D) ! regression matrix
         double precision :: R(Ns,Ns), Rinv(Ns,Ns) ! correlation matrix and inverse
         double precision :: bounds(D,2) ! bounds on theta 
         integer :: I_MIN=1, I_MAX=2

         ! initialize theta and bounds
         call init_theta(theta,bounds,X,D,Ns)

         ! optimize theta with iterative MLE
         call optimize_theta_mle(theta,bounds,X,Y,Order,D,Ns)

         ! get the response surface and the MSE
         call kriging_rs_and_mse(YNEW,XNEW,MSE,Y,X,theta,Order,D,Ns,NewNs)
      end subroutine

      ! kriging_rs_and_mse 
      !  Calculates YNEW and the MSE of all these points provided an optimized
      !  theta
      subroutine kriging_rs_and_mse(YNEW,XNEW,MSE,Y,X,theta,Order,D,Ns,NsNew)
         use PARAMS, only:Pc
         use matrix, only: eye
         use correlation, only: invertr, construct_R
         use regression, only: construct_fmat
         use mle, only : construct_beta, get_sigma2
         implicit none

         ! ARGUMENTS
         integer, intent(in) :: D, Ns, NsNew, Order
         double precision, intent(in) :: X(Ns,D), XNEW(NsNew,D)
         double precision, intent(in) :: Y(Ns,1)
         double precision, intent(in) :: theta(D)
         double precision, intent(out) :: YNEW(NsNew,1)
         double precision, intent(out) :: MSE(nsnew)

         ! Work variables
         double precision :: F(NS,1+Order*D)
         double precision :: R(Ns,Ns)
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

          ! construct F, R
         call construct_fmat(F,X,Order,D,Ns)
         call construct_R(R,theta,X,D,Ns,Pc)
        
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

         ! get YNEW
         call fast_kriging_function(YNEW, XNEW, X, theta, beta, RinvYmFb, sigma2, Order, D, Ns, Nsnew)

         ! get MSE
         call get_mse(mse, XNEW, X, theta, sigma2, F, R, Rinv, Order, D, ns, nsnew, fdim)
      end subroutine

      ! fast_kriging_function
      !  this function assumes you inverted R already, did some matrix
      !  multiplications that are constant and independent of the new snapshot
      !  location and calculates y^(x^)  
      subroutine fast_kriging_function(YNEW, XNEW, X, theta, beta, &
            RinvYmFb, sigma2, Order, D, Ns, Nsnew)
         use PARAMS, only:Pc
         use correlation, only: get_rxy
         use regression, only: construct_f
         implicit none

         ! ARGUMENTS
         integer, intent(in) :: D, Ns, NsNew, Order
         double precision, intent(out) :: YNEW(NsNew,1)
         double precision, intent(in) :: XNEW(NsNew,D)
         double precision, intent(in) :: X(Ns,D)
         double precision, intent(in) :: theta(D)
         double precision, intent(in) :: beta(1+Order*D,1) 
         double precision, intent(in) :: RinvYmFb(Ns,1)
         double precision, intent(in) :: sigma2

         ! Work variables
         double precision :: fx(1+Order*D,1)
         double precision :: rx(Ns,1)
         double precision :: res(1,1)

         integer :: ii,jj,fdim
         fdim = 1 + Order*D
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
         end do 
      end subroutine

      ! get_mse
      !  this function calculates the mse of all xnew's
      subroutine get_mse(mse, XNEW, X, theta, sigma2, F, R, Rinv, Order, D, ns, nsnew, fdim)
         use PARAMS, only:Pc
         use correlation, only: get_rxy
         use regression, only: construct_f
         use error, only: lophaven_mse, martin_mse
         integer, intent(in) :: D, Ns, NsNew, Order, fdim
         double precision, intent(in) :: XNEW(NsNew,D)
         double precision, intent(out) :: MSE(nsnew)
         double precision, intent(in) :: X(Ns,D)
         double precision, intent(in) :: theta(D)
         double precision, intent(in) :: sigma2
         double precision, intent(in) :: F(Ns,fdim), R(ns,ns), Rinv(ns,ns)

         ! Work variables
         double precision :: fx(1+Order*D,1)
         double precision :: rx(Ns,1)

         integer :: ii, jj
         do ii=1,NsNew
            ! construct f(x)
            call construct_f(fx(:,1),XNEW(ii,:),Order,D,Ns)
            
            ! construct r(x)
            do jj=1,ns
               rx(jj,1) = get_rxy(theta,xnew(ii,:),x(jj,:),D,Pc)
            end do

            MSE(ii) = martin_mse(sigma2, fx, rx, F, R, ns, fdim)
!            MSE(ii) = lophaven_mse(sigma2, fx, rx, F, R, Rinv, ns, fdim)
         end do 
      end subroutine
        
end module

