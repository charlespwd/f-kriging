! This module contains 
!  printer(x,y,ntotal,nrow,D,filename) - prints 2D data to a file (ready for
!     gnuplot)
!  process_command_input(funcname,ns,ngrid,xmin,xmax) - process command line
!     args
!  dumpmat,vec,tensor - dumps content of arrays to STDOUT
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

      subroutine process_command_input(funcname,Ns,Ngrid,xmin,xmax)
         integer :: ns, nsnew, ngrid
         character(len=20), intent(out) :: funcname
         double precision, dimension(2) :: xmin,xmax
         integer :: narg, i
         character(len=20) :: command
         narg = command_argument_count()
         ! DEFAULT VALUES
         funcname = '-d' ! by default
         Ns = 16
         Ngrid = 36
         xmin = (/0.d0,0.d0/)
         xmax = (/1.d0,1.d0/)
         if (narg > 0) then
            i = 1
            do while (i <= narg)
               call get_command_argument(i,command)
               select case(adjustl(command))
                  case ("-d","--drag")
                     funcname="-d"
                     xmin = (/0.d0,0.d0/)
                     xmax = (/1.d0,1.d0/)
                  case ("-b","--branin")
                     funcname="-b"
                     xmin = (/-5.d0,0.d0/)
                     xmax = (/5.d0,15.d0/)
                  case ("-c","--cosine")
                     funcname="-c"
                     xmin = (/0.d0,0.d0/)
                     xmax = (/1.d0,1.d0/)
                  case ("--ns")
                     i = i+1
                     call get_command_argument(i,command)
                     read(command,'(I10)') ns
                  case ("--ngrid")
                     i = i+1
                     call get_command_argument(i,command)
                     read(command,'(I10)') ngrid
                  case default
                     write(*,*) 'option "',trim(command),'" not supported' 
                     STOP
               end select
               i = i+1
            end do
         end if
      end subroutine

      SUBROUTINE DUMPMAT(A,dim1,dim2)
         INTEGER :: dim1, dim2 
         DOUBLE PRECISION :: A(dim1,dim2)
         INTEGER :: ii,jj
         900 format (I3,a1,i3,a2,F33.16)
         DO ii=1,dim1
         DO jj=1,dim2
            write(*,900) ii,",",jj,": ",A(ii,jj) 
         END DO
         END DO
      END SUBROUTINE

      SUBROUTINE DUMPVEC(A,DIMS)
         INTEGER :: DIMS
         DOUBLE PRECISION :: A(DIMS)
         INTEGER :: ii
         900 format (I3,a2,F8.5)
         DO ii=1,DIMS
            write(*,900) ii,": ",A(ii)
         end do
      END SUBROUTINE

      SUBROUTINE DUMPTENsOR(A,dim1,dim2,dim3)
         INTEGER :: dim1,dim2,dim3
         DOUBLE PRECISION :: A(dim1,dim2,dim3)
         INTEGER :: ii,jj,kk
         900 format (I3,a1,i3,a1,i3,a2,F8.5)
         DO ii=1,dim1
         DO jj=1,dim2
         DO kk=1,dim3
            write(*,900) ii,",",jj,",",kk,": ",A(ii,jj,kk) 
         END DO
         END DO
         END DO
      END SUBROUTINE

end module
