! THis module permits to call functions with a string command
! the beauty of this is that you don't need twenty interfaces 
! to achieve something beautiful
module ANALYTICAL_FUNCTIONS
   implicit none
   double precision, parameter :: PI = 4.0D0*DATAN(1.0D0)

   CONTAINS

   function Y_GRADIENT(X,D,Ns,func_name)
      ! arguments
      integer, intent(in) :: D,Ns
      double precision, intent(in) :: X(Ns,D)
      character(len=20), intent(in) :: func_name
      double precision, dimension(Ns,1+D) :: Y_GRADIENT
      
      integer :: nn
      double precision :: tmp(1,D)

      select case (adjustl(func_name))
         case("--drag","-d")
            do nn=1,Ns
               tmp = fdrag(X(nn,:),D)
               Y_GRADIENT(nn,1:(D+1)) = tmp(1,1:(D+1))
            end do
         case("--branin", "-b")
            do nn=1,Ns
               tmp = fbranin(X(nn,:),D)
               Y_GRADIENT(nn,1:(D+1)) = tmp(1,1:(D+1))
            end do
         case("--cosine", "-c")
            do nn=1,Ns
               tmp = fcosine(X(nn,:),D)
               Y_GRADIENT(nn,1:(D+1)) = tmp(1,1:(D+1))
            end do
         case("--rosenbrock", "-r");
            do nn=1,Ns
               tmp = frosenbrock(X(nn,:),D)
               Y_GRADIENT(nn,1:(D+1)) = tmp(1,1:(D+1))
            end do
         case DEFAULT
            print * , 'option ',adjustl(func_name), ' not supported'
            stop
      end select
   end function

   function fdrag(xv,D)
      use PARAMS, only : PI
      implicit none
      integer, intent(in) :: D
      double precision, intent(in) :: xv(1,D)
      double precision, dimension(1,1+2) :: fdrag
      double precision :: x,y
      ! D should be two
      if (D /= 2) then
         print *, 'D should be 2'
         stop
      end if

      x = xv(1,1)
      y = xv(1,2)

      ! y
      fdrag(1,1) =  0.001d0*(1+y**2 /(1.08d0-y)) *x;
      ! gradx
      fdrag(1,2) =  0.001d0 *(1+y**2 /(1.08d0-y));
      ! grady
      fdrag(1,3) =  0.001d0 * x  * (2 * y  / (1.08d0 - y) + y**2  / ((1.08d0-y)**2));
   end function

   ! fbranin(x,D) 
   ! returns the function y and the gradient of y @ x
   !  arguments:
   !  x: loc vector, [1,D]
   !  D: # of dimensions of domain, 
   !  fbranin: [y,gradx,grady], [D+1,1]
   function fbranin(xv,D)
      use PARAMS, only : PI
      implicit none
      integer, intent(in) :: D
      double precision, intent(in) :: xv(1,D)
      double precision, dimension(1,1+2) :: fbranin
      double precision :: x,y
      ! D should be two
      if (D /= 2) then
         print *, 'D should be 2'
         stop
      end if

      x = xv(1,1)
      y = xv(1,2)

      ! y
      fbranin(1,1) =  (y - 5.1d0/(4*PI**2) * x**2 + 5.d0/PI * x - 6)**2 + 10 * (1 - 1.0d0/(8*PI))*cos(x) + 10
      ! gradx
      fbranin(1,2) =  2 * (5.d0/PI - 0.258369d0*x) * (-6 + 5*x/PI - 0.129185d0*x**2 + y) - 10 * (1-1.0d0/(8*PI)) * sin(x)
      ! grady
      fbranin(1,3) = 2*(-6 + 5.d0*x/pi - 0.129185d0*x**2 + y)
   end function   
   
   function fcosine(xv,D)
      use PARAMS, only : PI
      implicit none
      integer, intent(in) :: D
      double precision, intent(in) :: xv(1,D)
      double precision, dimension(1,1+2) :: fcosine
      double precision :: x,y
      ! D should be two
      if (D /= 2) then
         print *, 'D should be 2'
         stop
      end if

      x = xv(1,1)
      y = xv(1,2)

      ! y
      fcosine(1,1) =  cos(10*x) + sin(10*y) + x*y
      ! gradx
      fcosine(1,2) =  -10*sin(10*x) + y
      ! grady
      fcosine(1,3) =  10*cos(10*y) + x;
   end function

   function frosenbrock(xv,D)
      use PARAMS, only : PI
      implicit none
      integer, intent(in) :: D
      double precision, intent(in) :: xv(1,D)
      double precision, dimension(1,1+2) :: frosenbrock
      double precision :: x,y
      ! D should be two
      if (D /= 2) then
         print *, 'D should be 2'
         stop
      end if

      x = xv(1,1)
      y = xv(1,2)

      ! y
      frosenbrock(1,1) = (1 - x) ** 2 + 100 * (y - x ** 2) ** 2
      ! gradx
      frosenbrock(1,2) = 2 * (200 * x ** 3 - 200 * x * y + x - 1) 
      ! grady
      frosenbrock(1,3) = 200 * (y - x ** 2) 
   end function
end module

!program p
!   use ANALYTICAL_FUNCTIONS
!   double precision :: X(3,2)
!   double precision :: Y(3,3) 
!   integer :: i
!   character(len=20) :: func_name
!   func_name = '-d'
!   X(:,1) = (/ 1,2,3 /)
!   X(:,2) = (/ 1,2,3 /)
!   Y = YGRAD(X,2,3,func_name)
!   do i =1,3
!      print * , Y(i,:)
!   end do
!end program
