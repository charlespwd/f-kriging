module utils
   implicit none
   contains
      subroutine printer(x,y,ntotal,nrow,D,filename)
         integer,intent(in) :: ntotal, nrow, D
         double precision, intent(in) :: x(ntotal,D), y(ntotal,1)
         double precision :: xt(D,ntotal), yt(1,ntotal)
         character(len=20), intent(in) :: filename
         integer :: ii,j, fileid = 56
         xt = transpose(x)
         yt = transpose(y)
         open(fileid,FILE=adjustl(filename),status='replace')
         j = 1
         do ii=1,ntotal
            write(fileid,*) xt(1,ii)," ",xt(2,ii)," ",yt(1,ii)
            j = j+1
            if (j > nrow) then
               write(fileid,*) " " ! a blank line
               j = 1
            end if
         end do
         close(fileid)
      end subroutine
end module
