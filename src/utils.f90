! This module contains 
!  printer(x,y,ntotal,nrow,D,filename) - prints 2D data to a file (ready for
!     gnuplot)
!  process_command_input(funcname,ns,ngrid,xmin,xmax) - process command line
!     args
!  dumpmat,vec,tensor - dumps content of arrays to STDOUT
module utils
   use sequential_sampling, only : MODE_MSE, MODE_SENSITIVITY
   implicit none
   contains
      ! assumes D = 2
      subroutine printer(x,y,ntotal,nrow,D,filename,datafolder,iterate)
         integer, intent(in) :: ntotal, nrow, D
         integer, intent(in),optional :: iterate
         character(len=20),intent(in), optional :: datafolder
         double precision, intent(in) :: x(ntotal,D), y(ntotal,1)
         double precision :: xt(D,ntotal), yt(1,ntotal)
         character(len=20), intent(in) :: filename
         character(len=20) :: wfilename
         integer :: ii,j, fileid = 56
         xt = transpose(x)
         yt = transpose(y)
         if (present(datafolder)) then
            wfilename = trim(datafolder)//'/'//trim(filename)
         else
            wfilename = adjustl(filename)
         endif
         ! append number to filename
         if (present(iterate)) then
            write(wfilename,'(a,i3.3)') trim(wfilename), iterate 
         endif

         open(fileid,FILE=adjustl(wfilename),status='replace')
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

      subroutine process_command_input(funcname,Ns,Ngrid,xmin,xmax,nfinal, &
            mode,deltans, order)
         integer :: ns, nsnew, ngrid
         integer, intent(inout), optional :: nfinal
         integer, intent(inout), optional :: mode !! mse=0, sensitivy=1
         integer, intent(inout), optional :: deltans
         integer, intent(inout), optional :: order
         character(len=20), intent(out) :: funcname
         double precision, dimension(2) :: xmin,xmax
         integer :: narg, i
         character(len=20) :: command
         narg = command_argument_count()
         ! DEFAULT VALUES
         funcname = '-d' ! by default
         Ns = 3 
         Ngrid = 40
         Nfinal = 30
         xmin = (/0.d0,0.d0/)
         xmax = (/1.d0,1.d0/)
         ! mode default value
         if (present(mode)) then
            mode = MODE_MSE 
         endif
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
                  case ("-r","--rosenbrock")
                     funcname="-r"
                     xmin = (/-2.d0,-2.d0/)
                     xmax = (/2.d0,2.d0/)
                  case ("--ns")
                     i = i+1
                     call get_command_argument(i,command)
                     read(command,'(I10)') ns
                  case ("--ngrid")
                     i = i+1
                     call get_command_argument(i,command)
                     read(command,'(I10)') ngrid
                  case ("--nfinal")
                     i = i+1
                     call get_command_argument(i,command)
                     if (present(nfinal)) then
                        read(command,'(I10)') nfinal
                     else 
                        print*, 'nfinal is not supported'
                        stop 
                     endif
                  case ("-s","--sensitivity")
                     if (present(mode)) then
                        mode = mode_sensitivity 
                     else 
                        print*, 'mode is not supported'
                        stop
                     endif
                  case ("--deltans")
                     if (present(deltans)) then
                        i = i+1
                        call get_command_argument(i,command)
                        read(command,'(i10)') deltans
                     else 
                        print*, 'deltans is not supported'
                        stop
                     endif
                  case ("--order","-o")
                     if (present(order)) then
                        i = i+1
                        call get_command_argument(i,command)
                        read(command,'(i10)') order
                     else 
                        print*, 'order was not passed'
                        stop
                     endif
                  case ("-m","--mse")
                     if (present(mode)) then
                        mode = MODE_MSE 
                     else
                        print*, 'mode is not supported'
                        stop
                     endif
                  case ("-h","--help") 
                     print*, "SYNOPSIS"
                     print*, "     ./bar [options]"
                     print*, ""
                     print*, "DESCRIPTION"
                     print*, "      -b, --branin"
                     print*, "            use the branin function"
                     print*, "      -c, --cosine"
                     print*, "            use the cosine function"
                     print*, "      -d, --drag"
                     print*, "            use the drag function, by default"
                     print*, "      --deltans"
                     print*, "            set the number of points to be added &
                        between every sampling iterations."
                     print*, "      -m, --mse"
                     print*, "            set the sampling criterion to be &
                        strictly mse, by default"
                     print*, "      --ngrid"
                     print*, "            set the number of points per &
                        dimension in the grid, default 40"
                     print*, "      --nfinal"
                     print*, "            set the maximum number of snapshots &
                        to be included in the set."
                     print*, "      --ns"
                     print*, "            set the initial number of snapshots per &
                        dimension, default 3 (for 3x3 initial grid)"
                     print*, "      -o, --order"
                     print*, "            set the order of the regression, default 0"
                     print*, "      -s, --sensitivity"
                     print*, "            set the sampling criterion to be &
                       using the sensitivity analysis"
                     print*, ""
                     print*, "AUTHOR"
                     print*, "      Written by Charles-Philippe Clermont"  
                     stop
                  case default
                     write(*,*) 'option "',trim(command),'" not supported' 
                     stop
               end select
               i = i+1
            end do
         end if
      end subroutine

      double precision function l1error(ytrue,ynew,nsnew)
         integer, intent(in) :: nsnew
         double precision, intent(in) :: ytrue(nsnew,1), ynew(nsnew,1)
         integer :: ii
         double precision :: meanl1
         MeanL1 = 0.d0
         do ii=1,NsNew
            MeanL1 = MeanL1 + abs(ytrue(ii,1) - ynew(ii,1))
         end do
         l1error = meanl1 
      end function

      subroutine DUMPMAT(A,dim1,dim2)
         integer :: dim1, dim2 
         double precision :: A(dim1,dim2)
         integer :: ii,jj
         900 format (I3,a1,i3,a2,F33.16)
         do ii=1,dim1
         do jj=1,dim2
            write(*,900) ii,",",jj,": ",A(ii,jj) 
         end do
         end do
      end subroutine

      subroutine DUMPVEC(A,DIMS)
         integer :: DIMS
         double precision :: A(DIMS)
         integer :: ii
         900 format (I3,a2,F8.5)
         do ii=1,DIMS
            write(*,900) ii,": ",A(ii)
         end do
      end subroutine

      subroutine DUMPTENsOR(A,dim1,dim2,dim3)
         integer :: dim1,dim2,dim3
         double precision :: A(dim1,dim2,dim3)
         integer :: ii,jj,kk
         900 format (I3,a1,i3,a1,i3,a2,F8.5)
         do ii=1,dim1
         do jj=1,dim2
         do kk=1,dim3
            write(*,900) ii,",",jj,",",kk,": ",A(ii,jj,kk) 
         end do
         end do
         end do
      end subroutine

end module
