program test
   use sequential_sampling, only : adaptive_ssck, mode, &
      mode_sensitivity, order
   implicit none

   integer :: D=1
   double precision :: xminmax(1,2)
   xminmax(1,1) = 0.d0
   xminmax(1,2) = 4.d0

   mode = MODE_SENSITIVITY
   order = 0
   call adaptive_ssck(func, grad, xminmax, D)
   
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
end program

