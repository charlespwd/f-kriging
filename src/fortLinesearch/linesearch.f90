! requires BLAS library (uses dgemm)
module linesearch_module
   implicit none
   integer :: iimax = 100
   double precision :: alpha_0 = 10.d0 
   double precision :: alpha_multiplier = 0.75d0
   double precision :: c_1 = 1.d-4 ! 10^-4
   double precision :: tol = 1.d-8 ! sum of grads
   double precision :: h = 1.d-10 ! cfde step size

   contains
      subroutine linesearch(x_star, x_initial, f, D)
         ! arguments 
         integer, intent(in) :: D
         double precision, intent(out) :: x_star(D,1)
         double precision, intent(in) :: x_initial(D,1)
         interface func
            double precision pure function f(x,D)
               integer, intent(in) :: D
               double precision, dimension(D,1),intent(in) :: x
            end function
         end interface func

         ! work variables
         double precision :: f_k
         double precision :: g_k(D,1)
         double precision :: p_k(D,1)
         double precision :: alpha_k
         double precision :: x_k(D,1)
         double precision :: x_kp(D,1)
         double precision :: condition
         integer :: ii,jj

         ii = 0
         x_kp = x_initial
         g_k = grad_cfde(f,x_kp,D)
         condition = sum(abs(g_k)) 
         do while (condition > tol .and. ii < iimax)
            ! update x_k
            x_k = x_kp

            ! get p_k
            f_k = f(x_k,D)
            p_k = -1.d0 * g_k

            ! get alpha_k
            alpha_k = get_alpha_backtrack(f, x_k, f_k, g_k, p_k, D)

            ! update x and gradient
            x_kp = x_k + alpha_k * p_k
            ii = ii + 1
            g_k = grad_cfde(f,x_kp,D)
            condition = sum(abs(g_k))
         enddo
         x_star = x_kp
      end subroutine

      ! grad_cfde
      ! calculate the gradient of input function f by central finite difference 
      function grad_cfde(f,x_k,D)
         ! arguments
         integer, intent(in) :: D
         double precision, intent(in) :: x_k(D,1)
         double precision :: grad_cfde(D,1)
         interface fnc
            double precision pure function f(x,D)
               integer, intent(in) :: D
               double precision, intent(in) :: x(D,1)
            end function
         end interface fnc
         
         ! work
         double precision :: ei(D,1)
         integer :: ii, jj

         do ii = 1,D
            do jj = 1,D
               ei(jj,1) = 0.d0
            enddo
            ei(ii,1) = 1.d0
            grad_cfde(ii,1) = (0.5d0 / h) * ( f(x_k + h*ei, D) - f(x_k - h*ei,D) ) 
         enddo
      end function

      double precision function get_alpha_backtrack(f, x_k, f_k, g_k, p_k, D)
         ! arguments
         integer, intent(in) :: D
         double precision, intent(in) :: x_k(D,1)
         double precision, intent(in) :: f_k
         double precision, intent(in) :: g_k(D,1)
         double precision, intent(in) :: p_k(D,1)
         interface fnc
            double precision pure function f(x,D)
               integer, intent(in) :: D
               double precision, intent(in) :: x(D,1)
            end function
         end interface fnc

         ! work
         double precision :: alpha_j

         alpha_j = alpha_0
         do while (phi(alpha_j,D) > wolfe_1(alpha_j,D)) 
            alpha_j = alpha_multiplier * alpha_j
         end do
         get_alpha_backtrack = alpha_j

         ! function contains helper functions, wolfe condition, etc. 
         contains 
            double precision function phi(alpha,D)
               integer, intent(in) :: D
               double precision, intent(in) ::  alpha
               phi = f(x_k + alpha * p_k, D)
            end function

            double precision function wolfe_1(alpha, D)
               integer, intent(in) :: D
               double precision, intent(in) :: alpha 
               double precision :: res(1,1)
               res(1,1) = f_k
               ! wolfe_1 = f_k + c_1 * alpha * g_k**T * p_k
               call dgemm('t','n', 1,1,D, c_1 * alpha, g_k,D, p_k,D, 1.0d0,res,1)
               wolfe_1 = res(1,1)
            end function
      end function
end module
