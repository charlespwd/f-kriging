module surrogate_model
   implicit none
   integer :: ns, order, d
   double precision, allocatable :: x(:,:) ! is rescaled
   double precision, allocatable :: theta(:)
   double precision, allocatable :: beta(:,:)
   double precision, allocatable :: RinvYmFb(:,:)
   double precision, allocatable :: XMIN(:), XMAX(:) ! for scaling/rescaling
   double precision :: sigma2 

   contains

   function kriging_function(xnew, D)
      use kriging_module, only : fast_kriging_function
      integer, intent(in) :: D
      double precision, intent(in) :: xnew(D,1)
      double precision :: ynew(1,1)
      double precision :: kriging_function
      double precision :: xtnew(1,D)
      xtnew = transpose(xnew)
   
      call fast_kriging_function(ynew, xnew, x, theta, beta, rinvymfb, &
      sigma2, order, d, ns, nsnew=1)
      kriging_function = ynew(1,1)
   end function

   subroutine init_kriging_constants(in_xold, in_yold, in_theta, in_order, in_d, in_ns) 
      use params, only : Pc
      use matrix, only : rescale, unscale, eye
      use mle, only : construct_beta, get_sigma2 
      use correlation, only: invertr, construct_R
      use regression, only: construct_fmat

      integer, intent(in) :: in_d, in_ns, in_order
      double precision, intent(in) :: in_xold(in_ns, in_d) ! IS NOT RESCALED.
      double precision, intent(in) :: in_yold(in_ns, 1) 
      double precision, intent(in) :: in_theta(in_d)

      double precision :: Y(in_ns, 1)
      double precision :: F(in_ns,1+in_order*in_d)
      double precision :: R(in_ns,in_ns)
      double precision :: Rinv(in_ns,in_ns), I(in_ns,in_ns)
      double precision :: YmFb(in_ns,1)
      double precision :: res(1,1)
      integer :: ii,jj,fdim

      if (allocated(x)) then
         deallocate(x)
         deallocate(theta)
         deallocate(beta)
         deallocate(rinvymfb)
         deallocate(xmin)
         deallocate(xmax)
      end if
      
      ns = in_ns
      d = in_d
      order = in_order
      fdim = 1 + Order*D
      allocate(x(ns,d))
      allocate(theta(d))
      allocate(beta(1+order*D, 1))
      allocate(rinvymfb(ns,1))
      allocate(xmin(d))
      allocate(xmax(d))

      x = in_xold
      xmin = minval(x,1)
      xmax = maxval(x,1)  
      ! the kriging is performed on a rescaled domain.
      call rescale(x, d, ns, xmin, xmax)

      y = in_yold
      
      theta = in_theta

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
   end subroutine

end module
