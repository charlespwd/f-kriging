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
         
         print*, imin
         print*, D
         print*, x(imin,1:D)
         carpet_bombed_min(1,1:D) = x(imin,1:D)
         y_star(1,1) = currentMin
      end function
end module

