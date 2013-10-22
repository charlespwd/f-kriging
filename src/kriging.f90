module krigingmodule
   implicit none

   contains
      ! !! KRIGING(XNEW,YNEW,theta,MSE,XOLD,Y,D,Ns,NewNs)
      ! This subroutine builds a response surface using kriging and an iterative
      !  MLE to obtain theta. 
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
      SUBROUTINE KRIGING(XNEW,YNEW,theta,MSE,XOLD,Y,Order,D,Ns,NewNs)
         USE PARAMS, ONLY: Pc
         use matrix, only: rescale
         use regression, only: construct_fmat
         use correlation, only : construct_r
         use mle, only : optimize_theta_mle, init_theta
         IMPLICIT NONE
         
         ! ARGUMENTS
         INTEGER, INTENT(IN) :: D, Ns, NewNs, Order
         DOUBLE PRECISION, INTENT(IN) :: XOLD(Ns,D)      ! unscaled set of snapshots
         DOUBLE PRECISION, INTENT(IN) :: XNEW(NewNs,D)   ! unscaled set of untried loc
         DOUBLE PRECISION, INTENT(IN) :: Y(Ns,1)         ! set of snapshot values
         DOUBLE PRECISION, INTENT(INOUT) :: theta(D)     ! correlation parameter
         DOUBLE PRECISION, INTENT(OUT) :: YNEW(NewNs,1)  ! set of untried values
         DOUBLE PRECISION, INTENT(OUT) :: MSE(NewNs)     ! error estimation

         ! Work variables
         DOUBLE PRECISION :: X(Ns,D), XN(NewNs,D) ! scaled set of snapshots
         DOUBLE PRECISION :: XMAX(D), XMIN(D) ! min's and max for scaling
         DOUBLE PRECISION :: F(Ns,1+Order*D) ! regression matrix
         DOUBLE PRECISION :: R(Ns,Ns), Rinv(Ns,Ns) ! correlation matrix and inverse
         DOUBLE PRECISION :: bounds(D,2) ! bounds on theta 
         INTEGER :: I_MIN=1, I_MAX=2

         XMIN = minval(XOLD,1)
         XMAX = maxval(XOLD,1)
         
         X = XOLD
         XN = XNEW
         
         CALL rescale(X,D,Ns,XMIN,XMAX)  ! Precondition X by "normalizing" it
         CALL rescale(XN,D,NewNs,XMIN,XMAX) ! " "

         ! initialize theta and bounds
         CALL init_theta(theta,bounds,X,D,Ns)

         ! optimize theta with iterative MLE
         CALL optimize_theta_mle(theta,bounds,X,Y,Order,D,Ns)

         ! construct F, R
         CALL construct_fmat(F,X,Order,D,Ns)
         CALL construct_R(R,theta,X,D,Ns,Pc)

         CALL construct_kriging_RS(YNEW,XN,MSE,Y,X,F,R,theta,Order,D,Ns,NewNs)
      end subroutine

      subroutine construct_kriging_RS(YNEW,XNEW,MSE,Y,X,F,R,theta,Order,D,Ns,NsNew)
         use PARAMS, only:Pc, MSEMODE
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
         double precision, intent(inout) :: R(Ns,Ns)
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
            end DO
            
            ! YNEW(ii,1) = fx'*beta + rx'*Rinv*YmFb
            !  first res = fx'*beta
            call DGEMM('t','n',1,1,fdim,1.0d0,fx,fdim,beta,fdim,0.0d0,res,1)

            !  then res = rx'*RinvYmFb + res
            call DGEMM('t','n',1,1,ns,1.0d0,rx,Ns,RinvYmFb,Ns,1.0d0,res,1)
            YNEW(ii,1) = res(1,1)

            select case (MSEMODE)
               case (MODE_ARTHUR)
                  MSE(ii) = get_mse(sigma2,F,Rinv,rx,D,Ns,fdim)
               case (MODE_MARTIN)
                  MSE(ii) = martinmse(sigma2, fx, rx, F, R, ns, fdim)
               case (MODE_LOPHAVEN)
                  MSE(ii) = lophaven_mse(sigma2, fx, rx, F, R, Rinv, ns, fdim)
               case default
                  print*, 'mode ', msemode, 'not supported. Terminating program.'
                  stop
            end select 
         end do 
      end subroutine
end module
