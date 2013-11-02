! THis module permits to call functions with a string command
! the beauty of this is that you don't need twenty interfaces 
! to achieve something beautiful
module ANALYTICAL_FUNCTIONS
   implicit none
   double precision, parameter :: PI = 4.0D0*DATAN(1.0D0)

   CONTAINS

   function Y_GRADIENT(x,D,Ns,func_name)
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
               tmp = ftdrag(X(nn,:),D)
               Y_GRADIENT(nn,1:(D+1)) = tmp(1,1:(D+1))
            end do
         case("--branin", "-b")
            do nn=1,Ns
               tmp = ftbranin(X(nn,:),D)
               Y_GRADIENT(nn,1:(D+1)) = tmp(1,1:(D+1))
            end do
         case("--cosine", "-c")
            do nn=1,Ns
               tmp = ftcosine(X(nn,:),D)
               Y_GRADIENT(nn,1:(D+1)) = tmp(1,1:(D+1))
            end do
         case("--rosenbrock", "-r");
            do nn=1,Ns
               tmp = ftrosenbrock(X(nn,:),D)
               Y_GRADIENT(nn,1:(D+1)) = tmp(1,1:(D+1))
            end do
         case DEFAULT
            print * , 'option ',adjustl(func_name), ' not supported'
            stop
      end select
   end function

   function ftdrag(xv,D)
      integer, intent(in) :: D
      double precision, intent(in) :: xv(1,D)
      double precision, dimension(1,1+2) :: ftdrag
      double precision :: x,y
      ! D should be two
      if (D /= 2) then
         print *, 'D should be 2'
         stop
      end if

      x = xv(1,1)
      y = xv(1,2)

      ! y
      ftdrag(1,1) =  0.001d0*(1+y**2 /(1.08d0-y)) *x;
      ! gradx
      ftdrag(1,2) =  0.001d0 *(1+y**2 /(1.08d0-y));
      ! grady
      ftdrag(1,3) =  0.001d0 * x  * (2 * y  / (1.08d0 - y) + y**2  / ((1.08d0-y)**2));
   end function

   function fdrag(xv,D)
      integer, intent(in) :: D
      double precision, intent(in) :: xv(D,1)
      double precision :: fdrag
      double precision :: x,y
      ! D should be two
      if (D /= 2) then
         print *, 'D should be 2'
         stop
      end if

      x = xv(1,1)
      y = xv(2,1)

      fdrag =  0.001d0*(1+y**2 /(1.08d0-y)) *x;
   end function 

   function gdrag(xv,D)
      integer, intent(in) :: D
      double precision, intent(in) :: xv(D,1)
      double precision :: gdrag(D,1)
      double precision :: x,y
      ! D should be two
      if (D /= 2) then
         print *, 'D should be 2'
         stop
      end if

      x = xv(1,1)
      y = xv(2,1)
      
      ! gradx
      gdrag(1,1) =  0.001d0 *(1+y**2 /(1.08d0-y));
      ! grady
      gdrag(2,1) =  0.001d0 * x  * (2 * y  / (1.08d0 - y) + y**2  / ((1.08d0-y)**2));
   end function

   ! ftbranin(xv,D) 
   ! returns the function y and the gradient of y @ x
   !  arguments:
   !  x: loc vector, [1,D]
   !  D: # of dimensions of domain, 
   !  ftbranin: [y,gradx,grady], [D+1,1]
   function ftbranin(xv,D)
      use PARAMS, only : PI
      implicit none
      integer, intent(in) :: D
      double precision, intent(in) :: xv(1,D)
      double precision, dimension(1,1+2) :: ftbranin
      double precision :: x,y
      ! D should be two
      if (D /= 2) then
         print *, 'D should be 2'
         stop
      end if

      x = xv(1,1)
      y = xv(1,2)

      ! y
      ftbranin(1,1) =  (y - 5.1d0/(4*PI**2) * x**2 + 5.d0/PI * x - 6)**2 + 10 * (1 - 1.0d0/(8*PI))*cos(x) + 10
      ! gradx
      ftbranin(1,2) =  2 * (5.d0/PI - 0.258369d0*x) * (-6 + 5*x/PI - 0.129185d0*x**2 + y) - 10 * (1-1.0d0/(8*PI)) * sin(x)
      ! grady
      ftbranin(1,3) = 2*(-6 + 5.d0*x/pi - 0.129185d0*x**2 + y)
   end function   
   
   function fbranin(xv,D)
      integer, intent(in) :: D
      double precision, intent(in) :: xv(D,1)
      double precision :: fbranin
      double precision :: x,y
      ! D should be two
      if (D /= 2) then
         print *, 'D should be 2'
         stop
      end if

      x = xv(1,1)
      y = xv(2,1)

      fbranin =  (y - 5.1d0/(4*PI**2) * x**2 + 5.d0/PI * x - 6)**2 + 10 * (1 - 1.0d0/(8*PI))*cos(x) + 10
   end function 

   function gbranin(xv,D)
      integer, intent(in) :: D
      double precision, intent(in) :: xv(D,1)
      double precision :: gbranin(D,1)
      double precision :: x,y
      ! D should be two
      if (D /= 2) then
         print *, 'D should be 2'
         stop
      end if

      x = xv(1,1)
      y = xv(2,1)
      
      ! gradx
      gbranin(1,1) =   2 * (5.d0/PI - 0.258369d0*x) * (-6 + 5*x/PI - 0.129185d0*x**2 + y) - 10 * (1-1.0d0/(8*PI)) * sin(x)
      ! grady
      gbranin(2,1) =  2*(-6 + 5.d0*x/pi - 0.129185d0*x**2 + y)
   end function

   function ftcosine(xv,D)
      use PARAMS, only : PI
      implicit none
      integer, intent(in) :: D
      double precision, intent(in) :: xv(1,D)
      double precision, dimension(1,1+2) :: ftcosine
      double precision :: x,y
      ! D should be two
      if (D /= 2) then
         print *, 'D should be 2'
         stop
      end if

      x = xv(1,1)
      y = xv(1,2)

      ! y
      ftcosine(1,1) =  cos(10*x) + sin(10*y) + x*y
      ! gradx
      ftcosine(1,2) =  -10*sin(10*x) + y
      ! grady
      ftcosine(1,3) =  10*cos(10*y) + x;
   end function

   function fcosine(xv,D)
      integer, intent(in) :: D
      double precision, intent(in) :: xv(D,1)
      double precision :: fcosine
      double precision :: x,y
      ! D should be two
      if (D /= 2) then
         print *, 'D should be 2'
         stop
      end if

      x = xv(1,1)
      y = xv(2,1)

      fcosine =  cos(10*x) + sin(10*y) + x*y
   end function 

   function gcosine(xv,D)
      integer, intent(in) :: D
      double precision, intent(in) :: xv(D,1)
      double precision :: gcosine(D,1)
      double precision :: x,y
      ! D should be two
      if (D /= 2) then
         print *, 'D should be 2'
         stop
      end if

      x = xv(1,1)
      y = xv(2,1)
      
      ! gradx
      gcosine(1,1) =   -10*sin(10*x) + y
      ! grady
      gcosine(2,1) =  10*cos(10*y) + x;
   end function

   function ftrosenbrock(xv,D)
      use PARAMS, only : PI
      implicit none
      integer, intent(in) :: D
      double precision, intent(in) :: xv(1,D)
      double precision, dimension(1,1+2) :: ftrosenbrock
      double precision :: x,y
      ! D should be two
      if (D /= 2) then
         print *, 'D should be 2'
         stop
      end if

      x = xv(1,1)
      y = xv(1,2)

      ! y
      ftrosenbrock(1,1) = (1 - x) ** 2 + 100 * (y - x ** 2) ** 2
      ! gradx
      ftrosenbrock(1,2) = 2 * (200 * x ** 3 - 200 * x * y + x - 1) 
      ! grady
      ftrosenbrock(1,3) = 200 * (y - x ** 2) 
   end function

   function frosenbrock(xv,D)
      integer, intent(in) :: D
      double precision, intent(in) :: xv(D,1)
      double precision :: frosenbrock
      double precision :: x,y
      ! D should be two
      if (D /= 2) then
         print *, 'D should be 2'
         stop
      end if

      x = xv(1,1)
      y = xv(2,1)

      frosenbrock =  (1 - x) ** 2 + 100 * (y - x ** 2) ** 2
   end function 

   function grosenbrock(xv,D)
      integer, intent(in) :: D
      double precision, intent(in) :: xv(D,1)
      double precision :: grosenbrock(D,1)
      double precision :: x,y
      ! D should be two
      if (D /= 2) then
         print *, 'D should be 2'
         stop
      end if

      x = xv(1,1)
      y = xv(2,1)
      
      ! gradx
      grosenbrock(1,1) = 2 * (200 * x ** 3 - 200 * x * y + x - 1) 
      ! grady
      grosenbrock(2,1) =  200 * (y - x ** 2) 
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
!   Y = YGRAD(xv,2,3,func_name)
!   do i =1,3
!      print * , Y(i,:)
!   end do
!end program
