! requires BLAS library (uses dgemm)
module linesearch_module
   implicit none
   integer :: iimax = 100
   double precision :: alpha_0 = 100.d0 
   double precision :: alpha_multiplier = 0.75d0
   double precision :: c_1 = 1.d-4 
   double precision :: c_2 = 0.9d0 
   double precision :: tol = 1.d-8 ! sum of grads
   double precision :: h = 1.d-8 ! cfde step size
   double precision :: tooclose = 1.d-5

   contains
      subroutine linesearch(x_star, x_initial, f, D, y_star)
         ! arguments 
         integer, intent(in) :: D
         double precision, intent(out) :: x_star(D,1)
         double precision, intent(in) :: x_initial(D,1)
         double precision, intent(out), optional :: y_star
         interface func
            double precision  function f(x,D)
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
         if (present(y_star)) then
            y_star = f(x_star,D)
         end if
      end subroutine

      subroutine quasi_newton_bfgs(x_star, x_initial, f, D, y_star)
         use matrix, only : eye
         ! arguments 
         integer, intent(in) :: D
         double precision, intent(out) :: x_star(D,1)
         double precision, intent(in) :: x_initial(D,1)
         double precision, intent(out), optional :: y_star
         interface func
            double precision  function f(x,D)
               integer, intent(in) :: D
               double precision, dimension(D,1),intent(in) :: x
            end function
         end interface func

         ! work variables
         double precision :: f_k, f_kp
         double precision :: g_k(D,1)
         double precision :: g_kp(D,1)
         double precision :: p_k(D,1)
         double precision :: H_k(D,D)
         double precision :: H_kp(D,D)
         double precision :: y_k(D,1)
         double precision :: s_k(D,1)
         double precision :: alpha_k
         double precision :: x_k(D,1)
         double precision :: x_kp(D,1)
         double precision :: condition
         integer :: ii,jj

         ! IMPORTANT - reset the alpha_0 global to 1 as suggested by Nocedal +
         ! Wright. A Quasi-Newton algorithm should have alpha_0 = 1.
         alpha_0 = 1.d0

         ii = 0
         x_kp = x_initial
         f_kp = f(x_kp, D)
         g_kp = grad_cfde(f,x_kp,D)
         call eye(H_kp,D)
         condition = sum(abs(g_kp)) 
         do while (condition > tol .and. ii < iimax)
            x_k = x_kp

            f_k = f_kp
            g_k = g_kp
            H_k = H_kp
            ! pk = -1 * H * gk
            call dgemm('n', 'n', D, 1, D, (-1.0d0), H_k, D, g_k, D, 0.d0, p_k, D) 

            alpha_k = get_strong_alpha(f, x_k, f_k, g_k, p_k, D)

            x_kp = x_k + alpha_k * p_k
            s_k = x_kp - x_k
            
            f_kp = f(x_kp,D)
            
            g_kp = grad_cfde(f, x_kp, D)
            y_k = g_kp - g_k

            call bfgs_update(H_kp, H_k, s_k, y_k, D) 
            ii = ii + 1
            condition = sum(abs(g_kp))
         enddo
         x_star = x_kp
         if (present(y_star)) then
            y_star = f(x_star,D)
         end if

         contains
            subroutine bfgs_update(H_kp, H_k, s_k, y_k, D)
               use matrix, only : eye
               integer, intent(in) :: D
               double precision, intent(in) :: H_k(D,D)
               double precision, intent(in) :: s_k(D,1)
               double precision, intent(in) :: y_k(D,1)
               double precision, intent(out) :: H_kp(D,D)
               double precision :: ImR(D,D)
               double precision :: ImL(D,D)
               double precision :: rho(1,1)
               double precision :: ImLH(D,D)

               ! rho = 1 / (y_k' * s_k)
               call dgemm('t','n', 1, 1, D, 1.0d0, y_k, D, s_k, D, 0.d0, rho, 1)
               rho(1,1) = 1.0d0 / rho(1,1)
               
               ! ImL = I - rho * s_k * y_k'
               call eye(ImL,D)
               call dgemm('n', 't', D, D, 1, -rho(1,1), s_k, D, y_k, D, 1.0d0, &
                  ImL, D)

               ! ImR
               call eye(ImR,D)
               call dgemm('n', 't', D, D, 1, -rho(1,1), y_k, D, s_k, D, 1.0d0, &
                  ImR, D)
   
               ! ImL H_k
               call dgemm('n', 'n', D, D, D, 1.0d0, ImL, D, H_k, D, 0.d0, ImLH, D)
               ! H_kp = ImL H_K ImR
               call dgemm('n', 'n', D, D, D, 1.0d0, ImLH, D, ImR, D, 0.0d0, H_kp, D)
               ! H_kp = H_kp + rho*sk*sk**T
               call dgemm('n', 't', D, D, 1, rho(1,1), s_k, D, s_k, D, 1.d0, H_kp, D)         
            end subroutine

            subroutine dfp_update(H_kp, H_k, s_k, y_k, D)
               integer, intent(in) :: D
               double precision, intent(in) :: H_k(D,D)
               double precision, intent(in) :: s_k(D,1)
               double precision, intent(in) :: y_k(D,1)
               double precision, intent(out) :: H_kp(D,D)
               double precision :: ytH(1,D), ytHy(1,1), HyytH(D,D)
               double precision :: sst(D,D), yts(1,1), yytH(D,D)

               call dgemm('t', 'n', 1, D, D, 1.0d0, y_k, D, H_k, D, 0.0d0, ytH, 1)
               
               call dgemm('n', 'n', 1, 1, D, 1.0d0, yth, 1, y_k, D, 0.0d0, &
                  ythy, 1)

               call dgemm('n', 'n', D, D, 1, 1.0d0, y_k, D, ytH, 1, 0.d0, &
                 yytH, D)

               call dgemm('n', 'n', D, D, D, (1.0d0/ythy(1,1)), H_k, D, yyth, &
                 D, 0.0d0, HyytH, D)

               call dgemm('t', 'n', 1, 1, D, 1.0d0, y_k, D, s_k, D, 0.0d0, &
                  yts, 1)

               call dgemm('n', 't', D, D, 1, (1.0d0/yts), s_k, D, s_k, D, &
                  0.0d0, sst, D)

               H_kp = H_k - HyytH + sst
            end subroutine
      end subroutine
      ! grad_cfde
      ! calculate the gradient of input function f by central finite difference 
      function grad_cfde(f,x_k,D)
         ! arguments
         integer, intent(in) :: D
         double precision, intent(in) :: x_k(D,1)
         double precision :: grad_cfde(D,1)
         interface fnc
            double precision  function f(x,D)
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

      double precision function get_strong_alpha(f, x_k, f_k, g_k, p_k, D)
         integer, intent(in) :: D
         double precision, intent(in) :: x_k(D,1)
         double precision, intent(in) :: p_k(D,1)
         double precision, intent(in) :: g_k(D,1)
         double precision, intent(In) :: f_k
         interface func
            function f(x, D)
               integer, intent(in) :: D
               double precision, intent(in) :: x(D,1)
               double precision :: f
            end function
         end interface func

         double precision :: a_j, a_jm, a_max, multiplier
         double precision :: phi_0, dphi_0
         double precision :: phi_j, dphi_j, phi_jm, dphi_jm
         double precision :: wolfe_2  
         double precision :: res(1,1)
         integer :: i

         a_jm = 0.d0
         a_j = alpha_0
         a_max = 100
         multiplier = 1.75 

         phi_0 = f_k
         call dgemm('t', 'n', 1, 1, D, 1.0d0, g_k, D, p_k, D, 1.0d0, res, 1)
         dphi_0 = res(1,1)

         wolfe_2 = (-1.0d0) * c_2 * dphi_0

         i = 1
         phi_j = phi_0
         dphi_j = dphi_0
         do while (a_j < a_max)
            phi_jm = phi_j
            phi_j  = phi(a_j)
            dphi_jm = dphi_j
            dphi_j  = dphi(a_j)

            if (phi_j > wolfe_1(a_j) .or. (phi_j > phi_jm .and. i > 1)) then
               get_strong_alpha = zoom(a_jm, phi_jm, dphi_jm, a_j, phi_j, dphi_j) 
               return
            end if   
            
            if (abs(dphi_j) <= wolfe_2) then
               get_strong_alpha = a_j
               return
            end if

            if (dphi_j >= 0) then
               get_strong_alpha = zoom(a_j, phi_j, dphi_j, a_jm, phi_jm, dphi_jm)
               return
            end if
            
            a_jm = a_j
            a_j  = a_j * multiplier
            i = i + 1
         end do

         get_strong_alpha = a_j

         contains
            double precision function phi(alpha)
               double precision, intent(in) ::  alpha
               phi = f(x_k + alpha * p_k, D)
            end function

            double precision function dphi(alpha_j)
               double precision, intent(in) :: alpha_j
               double precision :: g_j(D,1)
               double precision :: res(1,1)
               g_j = grad_cfde(f, x_k + alpha_j * p_k, D)
               call dgemm('t', 'n', 1, 1, D, 1.0d0, g_j, D, p_k, D, 1.0d0, res, 1)
               dphi = res(1,1)
            end function

            double precision function wolfe_1(alpha)
               double precision, intent(in) :: alpha 
               wolfe_1 = phi_0 + c_1 * alpha * dphi_0
            end function

            double precision function zoom(a_lo, phi_lo, dphi_lo, a_hi, phi_hi, dphi_hi)
               double precision, intent(in) :: a_lo, phi_lo, dphi_lo
               double precision, intent(in) :: a_hi, phi_hi, dphi_hi

               double precision :: a_l, phi_l, dphi_l
               double precision :: a_h, phi_h, dphi_h
               double precision :: a_j, phi_j, dphi_j
               double precision :: qa, qb, qc
               integer :: condition
               integer :: i 

               a_l = a_lo
               phi_l = phi_lo
               dphi_l = dphi_lo

               a_h = a_hi
               phi_h = phi_hi
               dphi_h = dphi_hi

               i = 1
               do while (i < 10)
                  a_j = c_interpolate(a_l, phi_l, dphi_l, a_h, phi_h, dphi_h, &
                     phi)
                  phi_j = phi(a_j) 
                  dphi_j = dphi(a_j)
                  if (abs(a_j - a_l) < tooclose) then
                     a_h = a_h / 2.0d0
                     phi_h = phi(a_h)
                     dphi_h = dphi(a_h)
                     i = i + 1
                  elseif (phi_j > wolfe_1(a_j) .or. phi_j >= phi_l) then
                     a_h    = a_j
                     phi_h  = phi_j
                     dphi_h = dphi_j
                     i = i + 1
                  else
                     if (abs(dphi_j) <= wolfe_2) then
                        zoom = a_j 
                        return
                     end if

                     if (dphi_j * (a_h - a_l) >= 0) then
                        a_h    = a_l
                        phi_h  = phi_l
                        dphi_h = dphi_l
                     end if

                     a_l    = a_j
                     phi_l  = dphi_j
                     dphi_l = dphi_j

                     i = i + 1
                  end if
               end do
               zoom = a_j
            end function

      end function

      double precision function c_interpolate(a_i, phi_i, dphi_i, &
            a_im, phi_im, dphi_im, phi)
         double precision, intent(in) :: a_i, phi_i, dphi_i
         double precision, intent(in) :: a_im, phi_im, dphi_im
         double precision :: amin, valmin
         double precision :: d1, d2
         interface func
            double precision function phi(a)
               double precision, intent(in) :: a
            end function
         end interface
         d1 = dphi_im + dphi_i - 3 * (phi_im - phi_i) / (a_im - a_i)
         d2 = d1**2 - dphi_im * dphi_i 
         if (d2 < 0) then
            print*, 'neg d2'
            stop
         end if
         d2 = sqrt(d2)
         amin = a_i - (a_i - a_im) * &
            (dphi_i + d2 - d1) / (dphi_i - dphi_im + 2 * d2) 
         valmin = phi(amin)
         if (phi_i < valmin) then
            amin = a_i
            valmin = phi_i
         end if
         if (phi_im < valmin) then
            amin = a_im
            valmin = phi_im
         end if
         c_interpolate = amin
      end function
      
      double precision function get_alpha_backtrack(f, x_k, f_k, g_k, p_k, D)
         ! arguments
         integer, intent(in) :: D
         double precision, intent(in) :: x_k(D,1)
         double precision, intent(in) :: f_k
         double precision, intent(in) :: g_k(D,1)
         double precision, intent(in) :: p_k(D,1)
         interface fnc
            double precision  function f(x,D)
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

!program test
!   use analytical_functions, only: frosenbrock
!   use linesearch_module, only: quasi_newton_bfgs, iimax
!   implicit none
!   double precision :: x_initial(2,1)
!   double precision :: x_star(2,1)
!   x_initial(:,1) = 0.0d0 
!   
!   call quasi_newton_bfgs(x_star, x_initial, frosenbrock, 2)
!   print*, 'xstar', x_star
!   contains 
!      function f(x,D)
!         integer, intent(in) :: D
!         double precision, intent(in) :: x(D,1)
!         double precision :: f
!         integer :: i
!         f = 0.d0
!         do i = 1,D
!            f = f + (x(i,1) - i)**2 
!         enddo
!      end function
!end program
