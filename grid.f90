MODULE grid
   implicit none

   contains
      ! provide linspaces in all dimensions and produce a vector 
      ! mesh 
      ! e.g. :
      ! X1, y1, [z11]
      ! X1, y2, [z12]
      ! X1, y3, [z13]
      ! X2, y1, [z21]
      ! ...
      ! X3, y3, [z33]      
      SUBROUTINE vector_grid(X,LINSPACES,D,ngrid)
         integer, intent(in) :: D, ngrid
         double precision, intent(in) :: linspaces(ngrid,D)
         double precision, intent(out) :: X(ngrid**D,D)
         
         integer :: counts(D), i, row
         counts = 1
         row = 1 
         DO WHILE(counts(1) .le. ngrid)
            do i = 1,D
               X(row,i) = LINSPACES(counts(i),i)
            end do
            row = row + 1
            counts = countplusplus(counts,D,ngrid)
         end do
      end subroutine
      
      recursive function countplusplus(counts,D,ngrid) result(results)
         integer :: D,ngrid
         integer, intent(in) :: counts(d)
         integer, allocatable :: r(:)
         integer,dimension(D) :: results

         if (D>1) then
            allocate(r(D-1))
         end if
         results = counts
         results(D) = results(D) + 1
         if (results(D) .gt. ngrid) then
            if (D .gt. 1) then
               results(D) = 1
               r = countplusplus(results(1:D-1),D-1,ngrid)
               results(1:D-1) = r(1:D-1)
            else
               ! do nothing
            end if
         end if 
      end function
END MODULE
