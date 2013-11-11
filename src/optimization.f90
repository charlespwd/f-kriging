module optimization
   contains
      function carpet_bombed_min(y,x,D,Ns,y_star)
         ! arguments
         integer, intent(in) :: D, Ns
         double precision, intent(in) :: y(ns,1), x(ns,D)
         double precision, intent(out), optional :: y_star(1,1)
         double precision :: carpet_bombed_min(1,D)

         ! work
         integer :: i, imin
         double precision :: currentMin

         currentMin = huge(currentMin)
         do i=1,ns
            if (y(i,1) < currentMin) then
               currentMin = y(i,1)
               imin = i
            endif
         end do
         
         carpet_bombed_min(1,1:D) = x(imin,1:D)
         y_star(1,1) = currentMin
      end function

      function kriging_linesearched_min(x_guess, xold, yold, theta, order, d, ns, y_star, info)
         use surrogate_model, only : kriging_function, init_kriging_constants, &
            xmin, xmax
         use linesearch_module, only : quasi_newton_bfgs
         use matrix, only : rescale, unscale
         integer, intent(in) :: order, d, ns
         integer, intent(out), optional :: info ! 1 if out of search domain 
         double precision, intent(in) :: x_guess(1,D)
         double precision, intent(in) :: xold(ns,D)
         double precision, intent(in) :: yold(ns,1)
         double precision, intent(in) :: theta(D)
         double precision, intent(out) :: y_star(1,1)
         double precision :: kriging_linesearched_min(1,D)
         double precision :: x_star(D,1), xt_guess(D,1)
         double precision :: x_guess_scaled(1,D)
         integer :: i

         call init_kriging_constants(xold, yold, theta, order, d, ns)

         ! kriging works on a rescaled domain. Therefore, the guess provided to
         ! the line search algorithm should be rescaled as well.
         x_guess_scaled = x_guess
         call rescale(x_guess_scaled, D, ns=1, xmin=xmin, xmax=xmax) 

         ! the line search algorithm works with column vectors. 
         xt_guess = transpose(x_guess)

         call quasi_newton_bfgs(x_star, xt_guess, kriging_function, D, y_star(1,1))

         kriging_linesearched_min = transpose(x_star)
         call unscale(kriging_linesearched_min, D, ns=1, xmin=xmin, xmax=xmax)
         
         if (present(info)) then
            info = 0
         end if
         do i=1,D
            if (kriging_linesearched_min(1,i) < xmin(i) &
                .or. &
                kriging_linesearched_min(1,i) > xmax(i) ) then
               if(present(info)) then 
                  print *,  "min out of domain", kriging_linesearched_min(i,1),&
                     xmin(i), xmax(i)
                  info = 1
               else
                  print *,  "minimum out of domain"
                  stop
               end if
            end if
         end do
         
      end function
end module

