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

      subroutine process_command_input(funcname,Ns,Ngrid)
         integer :: ns, nsnew, ngrid
         character(len=20), intent(out) :: funcname
         integer :: narg, i
         character(len=20) :: command
         narg = command_argument_count()
         ! DEFAULT VALUES
         funcname = '-d' ! by default
         Ns = 16
         Ngrid = 36
         if (narg > 0) then
            i = 1
            do while (i <= narg)
               call get_command_argument(i,command)
               select case(adjustl(command))
                  case ("-d","--drag")
                     funcname="-d"
                  case ("-b","--branin")
                     funcname="-b"
                  case ("--ns")
                     i = i+1
                     call get_command_argument(i,command)
                     read(command,'(I10)') ns
                  case ("--ngrid")
                     i = i+1
                     call get_command_argument(i,command)
                     read(command,'(I10)') ngrid
               end select
               i = i+1
            end do
         end if
      end subroutine
end module
