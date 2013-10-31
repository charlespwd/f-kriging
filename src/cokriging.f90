! subroutine COKRIGING(XNEW,YNEW,theta,MSE,XOLD,YOLD,Grad,Raug,D,Ns,NsNew)
!  This routine is a wrapper to kriging, but augments the data set with 
!  information about the gradient around the snapshots
!  
! Arguments:
!  XNEW (in): 'new snapshots' location, [NsNew,D]
!  YNEW (out): 'new snapshots' values, [NsNew,1]
!  theta(inout), enters as theta initial, exits as optimized theta
!  MSE (out): error
!  XOLD (in): set of snapshot locations [Ns,D]
!  YOLD (in): set of snapshot values [Ns,1]
!  Grad (in): the gradient of Y at the snapshot locations in d'th direction [Ns,D]
!  Raug (in): the 'radius' of the augmentation as fraction of range(X), typically 1% 
!  Order(in): the order of the polynomial regression 
!  D (in): # of dimensions
!  Ns (in): # of snapshots
!  NsNew: # of new snapshots
!
!  Cokriging works by augmenting the initial data set as follows:
!  For each snasphots, an increment Raug*range(X,d) is added and subtracted
!     i.e. for each dimension X = XOLD +/- Raug*range(X)
!  For each snapshot locations, the snapshot value is derived from a taylor
!     expansion, i.e. for a corresponding X, Y = Yold + Grad * Raug * range(X)
!  After that, kriging is performed on this new, augmented, data set

! module is here to enable optional arguments
module cokrigingmodule
   CONTAINS
subroutine COKRIGING(XNEW,YNEW,theta,MSE,XOLD,YOLD,Grad,Raug,Order,D,Ns,NsNew,S)
   use sensitivity, only: construct_sensitivity
   use kriging_module, only: kriging, precondition
   implicit none

   ! arguments
   integer, intent(in) :: D, Ns, NsNew, Order
   double precision, intent(in) :: YOLD(Ns,1) 
   double precision, intent(in) :: XOLD(Ns,D)
   double precision, intent(in) :: XNEW(NsNew,D)
   double precision, intent(in) :: Grad(Ns,D)
   double precision, intent(in) :: Raug
   double precision, intent(INOUT) :: theta(D)
   double precision, intent(out) :: YNEW(NsNew,1)
   double precision, intent(out) :: MSE(NsNew)
   double precision, intent(out), optional, target :: S(:,:)
  
   ! work variables
   integer :: Naug
   integer :: dd, nn, kk
   double precision,ALLOCATABLE :: X(:,:), Y(:,:)
   double precision :: delta(D)
   double precision :: XMAX(D),XMIN(D)
   double precision,allocatable :: x_pre(:,:)
   double precision :: xnew_pre(nsnew,D)
   
   Naug=(2*D+1)*Ns
   allocate(X(Naug,D))
   allocate(Y(Naug,1))
   allocate(x_pre(Naug,D))
   
   XMAX = MAXVAL(XOLD,1)
   XMIN = MINVAL(XOLD,1)

   do dd=1,D
      delta(dd) = (XMAX(dd) - XMIN(dd)) * Raug / (ns ** (1.0d0 / D))
   end do

   X(1:Ns,:) = XOLD(1:Ns,:)
   Y(1:Ns,:) = YOLD(1:Ns,:)

   kk = Ns + 1
   do nn = 1,Ns
      do dd =1,D
         X(kk,:) = X(nn,:)
         X(kk,dd) = X(nn,dd) + delta(dd)
         Y(kk,1) = Y(nn,1) + delta(dd) * Grad(nn,dd)
         
         kk = kk + 1

         X(kk,:) = X(nn,:)
         X(kk,dd) = X(nn,dd) - delta(dd)
         Y(kk,1) = Y(nn,1) - delta(dd) * Grad(nn,dd)

         kk = kk + 1
      end do
   end do

   ! precondition for better kriging! Seriously, it performs better in a range
   ! of length 1.
   x_pre = x
   xnew_pre = xnew
   call precondition(x_pre,xnew_pre,x,xnew,D,naug,nsnew)
   
   ! perform kriging normally with new set of X and Y 
   call KRIGING(XNEW_pre,YNEW,theta,MSE,x_pre,Y,Order,D,Naug,NsNew)

   if ( present(S) ) then
      call construct_sensitivity(S,XNEW,YOLD,XOLD,Grad,theta,Order,D,Ns,NsNew)
   end if
end subroutine
end module
