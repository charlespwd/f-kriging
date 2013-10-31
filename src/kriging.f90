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
module kriging_module
   contains
subroutine kRIGING(XNEW,YNEW,theta,MSE,XOLD,Y,Order,D,Ns,NewNs)
   use PARAMS, only: Pc
   use matrix, only: rescale
   use regression, only: construct_fmat
   use correlation, only : construct_r
   use mle, only : optimize_theta_mle, init_theta
   implicit none
   
   ! ARGUMENTS
   integer, intent(in) :: D, Ns, NewNs, Order
   double precision, intent(in) :: XOLD(Ns,D)      ! unscaled set of snapshots
   double precision, intent(in) :: XNEW(NewNs,D)   ! unscaled set of untried loc
   double precision, intent(in) :: Y(Ns,1)         ! set of snapshot values
   double precision, intent(INOUT) :: theta(D)     ! correlation parameter
   double precision, intent(out) :: YNEW(NewNs,1)  ! set of untried values
   double precision, intent(out) :: MSE(NewNs)     ! error estimation

   ! Work variables
   double precision :: X(Ns,D), XN(NewNs,D) ! scaled set of snapshots
   double precision :: XMAX(D), XMIN(D) ! min's and max for scaling
   double precision :: F(Ns,1+Order*D) ! regression matrix
   double precision :: R(Ns,Ns), Rinv(Ns,Ns) ! correlation matrix and inverse
   double precision :: bounds(D,2) ! bounds on theta 
   integer :: I_MIN=1, I_MAX=2

   XMIN = minval(XOLD,1)
   XMAX = maxval(XOLD,1)
   
   X = XOLD
   XN = XNEW
   
   call rescale(X,D,Ns,XMIN,XMAX)  ! Precondition X by "normalizing" it
   call rescale(XN,D,NewNs,XMIN,XMAX) ! " "

   ! initialize theta and bounds
   call init_theta(theta,bounds,X,D,Ns)

   ! optimize theta with iterative MLE
   call optimize_theta_mle(theta,bounds,X,Y,Order,D,Ns)

   ! construct F, R
   call construct_fmat(F,X,Order,D,Ns)
   call construct_R(R,theta,X,D,Ns,Pc)

   call construct_kriging_RS(YNEW,XN,MSE,Y,X,F,R,theta,Order,D,Ns,NewNs)
end subroutine
