program test
   use analytical_functions
   use sequential_sampling 
   implicit none

   integer :: D=2
   double precision :: xminmax(1,2)
   xminmax(1:D,1) = -1.d0
   xminmax(1:D,2) = 1.d0

   ngridstart = 3
   nfinal = 60
   mode = MODE_SENSITIVITY
   order = 1
   call adaptive_ssck(fcosine, gcosine, xminmax, 2)
   
   contains 
      function func(x,D)
         integer,intent(in) :: D
         double precision, intent(in) :: x(D,1)
         double precision :: func
         func = (x(1,1) - 1.0d0) ** 2
      end function

      function grad(x,D)
         integer, intent(in) :: D
         double precision, intent(in) :: x(D,1)
         double precision :: grad(D,1)
         grad = 2 * (x(1,1) - 1) 
      end function

      function f2(x,D)
         integer,intent(in) :: D
         double precision, intent(in) :: x(D,1)
         double precision :: f2
         f2 = (x(1,1) - 1.0d0) ** 2 + (x(2,1) - 1.5d0) ** 2
      end function

      function g2(x,D)
         integer, intent(in) :: D
         double precision, intent(in) :: x(D,1)
         double precision :: g2(D,1)
         g2(1,1) = 2 * (x(1,1) - 1.d0) 
         g2(2,1) = 2 * (x(2,1) - 1.5d0) 
      end function
end program

