
! This module contains (at a high level)
!  vector_grid(X,LINSPACES,D,NGrid) - makes a long grid vector, (e.g. snap_pos)
!  lhs(Xmin,Xmax,D,Ns) - makes an LHS of the domain 
module grid
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
      subroutine vector_grid(X,LINSPACES,D,ngrid)
         integer, intent(in) :: D, ngrid
         double precision, intent(in) :: linspaces(ngrid,D)
         double precision, intent(out) :: X(ngrid**D,D)
         
         integer :: counts(D), i, row
         counts = 1
         row = 1 
         do while(counts(1) .le. ngrid)
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

      function LHS(XMIN,XMAX,D,Ns)
         implicit none
         integer, intent(in) :: D, Ns
         double precision, intent(in) :: xmin(D), xmax(D)
         double precision, dimension(NS,D) :: LHS

         integer :: idx(Ns,1)
         double precision :: ran(NS,D), P(Ns,1) 
         integer :: ii

         call init_random_seed()
         call RANDOM_NUMBER(ran)

         do ii=1,D
            call rperm(Ns,idx(1:Ns,1))
            P(1:Ns,1) = (idx(1:Ns,1)-ran(1:Ns,ii))/Ns
            LHS(1:Ns,ii) = xmin(ii) + P(1:Ns,1) * (XMAX(ii) - XMIN(ii))
         end do
            
      end function

      ! creates a linearly spaced vector spanning xmin to xmax
      function linspace(xmin,xmax,n)
         integer, intent(in) :: n
         double precision, intent(in) :: xmin,xmax
         double precision, dimension(n,1) :: linspace
         integer :: ii
         linspace(1:n,1) = (/(xmin + (ii-1) * ((xmax - xmin)/(n-1)),ii=1,n)/)
      end function

      ! makes a grid for xnew (linear lhs)
      subroutine columngrid(xnew,xmin,xmax,d,ngrid)
         implicit none
         integer, intent(in) :: d, ngrid
         double precision, intent(in) :: xmin(d), xmax(d)
         double precision, intent(out) :: xnew(ngrid ** D)
         double precision :: linspaces(ngrid,d)
         double precision :: tmpspace(ngrid,1)
         integer :: ii
         
         do ii=1,D
            tmpspace = linspace(xmin(ii),xmax(ii),ngrid)
            linspaces(1:ngrid,ii) = tmpspace(1:ngrid,1)
         enddo
         call vector_grid(xnew,linspaces,D,ngrid)
      end subroutine

      ! from the web
      subroutine init_random_seed()
         implicit none
         integer, allocatable :: seed(:)
         integer :: i, n, un, istat, dt(8), pid, t(2), s
         integer(8) :: count, tms
       
         call random_seed(size = n)
         allocate(seed(n))
         ! First try if the OS provides a random number generator
         open(newunit=un, file="/dev/urandom", access="stream", &
              form="unformatted", action="read", status="old", iostat=istat)
         if (istat == 0) then
            read(un) seed
            close(un)
         else
            ! Fallback to XOR:ing the current time and pid. The PID is
            ! useful in case one launches multiple instances of the same
            ! program in parallel.
            call system_clock(count)
            if (count /= 0) then
               t = transfer(count, t)
            else
               call date_and_time(values=dt)
               tms = (dt(1) - 1970) * 365_8 * 24 * 60 * 60 * 1000 &
                    + dt(2) * 31_8 * 24 * 60 * 60 * 1000 &
                    + dt(3) * 24 * 60 * 60 * 60 * 1000 &
                    + dt(5) * 60 * 60 * 1000 &
                    + dt(6) * 60 * 1000 + dt(7) * 1000 &
                    + dt(8)
               t = transfer(tms, t)
            end if
            s = ieor(t(1), t(2))
            pid = getpid() + 1099279 ! Add a prime
            s = ieor(s, pid)
            if (n >= 3) then
               seed(1) = t(1) + 36269
               seed(2) = t(2) + 72551
               seed(3) = pid
               if (n > 3) then
                  seed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
               end if
            else
               seed = s + 37 * (/ (i, i = 0, n - 1 ) /)
            end if
         end if
         call random_seed(put=seed)
      end subroutine init_random_seed

      ! Knuth shuffle, from the web
      subroutine RPERM(N, P)

         integer, intent(in) :: N
         integer, dimension(:), intent(out) :: P

         integer :: I
         integer :: K, J, IPJ, ITEMP, M
         double precision, dimension(100) :: U

         P = (/ (I, I=1,N) /)

         ! GENERATE UP TO 100 U(0,1) NUMBERS AT A TIME.
         do I=1,N,100
         M = MIN(N-I+1, 100)
         call RANDOM_NUMBER(U)
         do J=1,M
         IPJ = I+J-1
         K = INT(U(J)*(N-IPJ+1)) + IPJ
         ITEMP = P(IPJ)
         P(IPJ) = P(K)
         P(K) = ITEMP
         end do
         end do
         RETURN

      end subroutine RPERM

		subroutine grow2d(array,deltarow) 
			integer :: deltarow
			double precision, allocatable, intent(inout) :: array(:,:)
			double precision, allocatable :: tmparray(:,:)
		
			allocate(tmparray(size(array,1)+deltarow,size(array,2)))
			tmparray(1:size(array,1),1:size(array,2)) = array
			deallocate(array)
			allocate(array(size(tmparray,1),size(tmparray,2)))
			array = tmparray
		end subroutine

end module
