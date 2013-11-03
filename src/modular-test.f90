program test
   use analytical_functions
   use sequential_sampling 
   use linesearch_module
   implicit none

   integer :: D=2
   character(len=20) :: fname

   tol = 1.d-12
   ngridstart = 3
   nfinal = 50
!   mode = MODE_SENSITIVITY
   order = 2
   fname = "rosenbrock"
   call adaptive_ssck(frosenbrock, grosenbrock, get_range(fname, d), 2, optimize=1)
   
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
      
      function get_range(fname,D)
         integer :: D
         character(len=20) :: fname
         double precision :: get_range(D,2)
         select case(fname)
            case ("drag","cosine")
               get_range(1:D,1) = 0.d0
               get_range(1:D,2) = 1.d0
            case ("branin")
               get_range(1:D,1) = (/-5.d0, 0.d0/)
               get_range(1:D,2) = (/5.d0, 15.d0/)
            case ("rosenbrock") 
               get_range(1:D,1) = -2.d0
               get_range(:,2) = 2.d0
            case default
               get_range(:,1) = 0.d0
               get_range(:,2) = 1.d0
         end select 
      end function


end program

