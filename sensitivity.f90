module sensitivity   
   implicit none

   contains

      ! constructs the density function 
      !  defined as 1 - exp(mindist to a snapshot) 
      !  so that this density function is 0 at snapshot locations
      !  and more than that everywhere else
      !  the result is 'normalized'
      subroutine construct_density_function(Psi,Xnew,X,XMAX,XMIN,D,nsnew,ns)
         use matrix, only : normalize
         ! arguments
         integer, intent(in) :: nsnew, ns, d
         double precision, intent(in) :: Xnew(nsnew,D), X(ns,D)
         double precision, intent(in) :: XMAX(D), XMIN(D)
         double precision, intent(out) :: Psi(nsnew,1)
         
         ! work variables
         double precision :: mindist,dist
         integer :: ii,jj,pp
         
         do ii=1,nsnew
            ! initialize at maximum
            mindist = maxval(xmax)
            do jj=1,ns
               dist=0.0d0
               do pp=1,D
                  dist = dist + (xnew(ii,pp) - x(jj,pp)) ** 2
               enddo
               dist = dist ** (0.5d0)
               if ( dist < mindist ) then
                  mindist = dist
               endif
            enddo
            Psi(ii,1) = 1.0d0 - exp(-mindist)
         enddo
         call normalize(psi(:,1),nsnew)
      end subroutine

end module

!program test
!   use sensitivity, only : construct_density_function
!   double precision :: xmax(2), xmin(2)
!   double precision :: xnew(5,2), xold(4,2)
!   double precision :: psi(5) 
!   integer :: ns, d, nsnew
!   ns = 4
!   nsnew = 5
!   d = 2
!   xmax = (/1,1/)
!   xmin = (/0,0/)
!   xnew(:,1) = (/0,0,0,0,0/)
!   xnew(:,2) = (/0.0d0,0.2d0,0.4d0,0.6d0,0.8d0/)
!   xold(:,1) = (/0,0,1,1/)
!   xold(:,2) = (/0,1,0,1/)
!   call construct_density_function(psi, xnew, xold,xmax,xmin,d,nsnew,ns)
!   print * , psi
!   
!end program
