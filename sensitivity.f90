! dp get_sampling_radius(xmin,xmax,D,ns)
! construct_density_function(Psi,Xnew,X,XMAX,XMIN,D,nsnew,ns)
! construct_sensitivity(S,XNEW,Y,X,Grad,F,R,theta,D,Ns,NsNew)

module sensitivity   
   use la_precision, only:wp=>dp
   use f95_lapack, only:LA_GESV
   use matrix, only : eye, normalize
   USE PARAMS, ONLY:Order,Pc
   use correlation, only:get_rxy,invertr
   use regression, only:construct_f

   implicit none

   private
   public construct_density_function, get_sampling_radius, &
      construct_sensitivity

   contains

      ! constructs the density function 
      !  defined as 1 - exp(mindist to a snapshot) 
      !  so that this density function is 0 at snapshot locations
      !  and more than that everywhere else
      !  the result is 'normalized'
      subroutine construct_density_function(Psi,Xnew,X,XMAX,XMIN,D,nsnew,ns)
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

      ! get the sampling radius as per defined in Arthur's paper
      double precision function get_sampling_radius(xmin,xmax,D,ns)
         ! arguments
         integer, intent(in) :: D,ns
         double precision, intent(in) :: xmax(D),xmin(D)
         
         ! work
         double precision :: radius=0
         integer :: ii

         do ii=1,D
           Radius = Radius + (xmax(ii) - xmin(ii)) ** 2
         enddo
         Radius = 0.25d0 * (Radius ** 0.5d0) / (ns ** (1.0d0/D))
         get_sampling_radius = Radius
      end function

      ! construct the sensitivity vectort
      SUBROUTINE construct_sensitivity(S,XNEW,Y,X,Grad,F,R,theta,D,Ns,NsNew)
         IMPLICIT NONE

         ! ARGUMENTS
         INTEGER, INTENT(IN) :: D,Ns,NsNew
         DOUBLE PRECISION, INTENT(IN) :: X(Ns,D), XNEW(NsNew,D)
         DOUBLE PRECISION, INTENT(IN) :: Y(Ns,1)
         double precision, intent(in) :: Grad(Ns,D)
         DOUBLE PRECISION, INTENT(IN) :: F(NS,1+Order*D)
         DOUBLE PRECISION, INTENT(INOUT) :: R(Ns,Ns)
         DOUBLE PRECISION, INTENT(IN) :: theta(D)
         DOUBLE PRECISION, INTENT(OUT) :: S(NsNew,1)

         ! Work variables
         DOUBLE PRECISION :: Rinv(Ns,Ns), I(Ns,Ns)
         DOUBLE PRECISION :: fx(1,1+Order*D)
         DOUBLE PRECISION :: rx(Ns,1)
         DOUBLE PRECISION :: res(1,1)
         double precision :: psi(ns,1) 
         double precision :: gradsum(ns,1)
         double precision :: zeta(1+Order*D,1)
         double precision :: DRHS(ns,1)


         INTEGER :: ii,jj,dd,fdim
         fdim = 1 + Order*D
         
         CALL eye(I,Ns)
         CALL invertR(R,I,Rinv,Ns)

         ! construct gradsum
         do jj=1,ns
            gradsum(jj,1) =0.0d0
            do dd=1,D
               gradsum(jj,1) = gradsum(jj,1) + abs(grad(jj,dd))
            enddo
         enddo
         call normalize(gradsum(:,1),ns)

         DO ii=1,NsNew
            ! construct f(x)
            call construct_f(fx,XNEW(ii,:),Order,D,Ns)
            
            ! construct r(x)
            DO jj=1,ns
               rx(jj,1) = get_rxy(theta,xnew(ii,:),x(jj,:),D,Pc)
            END DO

            ! construct PSI
            psi=0
            S(ii,1) = 0.0d0
            do jj=1,ns
               psi(jj,1) = gradsum(jj,1)
               call construct_zeta(zeta,F,Rinv,Psi,D,Ns,fdim)
               call construct_DRHS(DRHS,F,Rinv,Psi,zeta,D,Ns,fdim)
               S(ii,1) = S(ii,1) + get_S_j(fx,zeta,rx,DRHS,D,Ns,fdim)
               psi(jj,1) = 0.0d0
            enddo
         END DO 
         call normalize(S(:,1),nsnew)

      END SUBROUTINE

      ! -----------------------------------------------------------------
      ! private functions to compute the sentitivity 
      ! -----------------------------------------------------------------
         ! zeta=(F'*Rinv*F)\(F'*Rinv*Psi);
         subroutine construct_zeta(zeta,F,Rinv,Psi,D,Ns,fdim)
            integer, intent(in) :: d, ns, fdim
            double precision, intent(in) :: F(ns,fdim), Rinv(ns,ns)
            double precision, intent(in) :: Psi(ns,1)
            double precision, intent(out) :: zeta(fdim,1)
            
            ! work vars
            double precision :: ftRinv(fdim,ns), FtRinvF(fdim,fdim)
            double precision :: FtRinvPsi(fdim,1)

            ! F**T * Rinv
            call dgemm('t','n',fdim,ns,ns,1.0d0,F,ns,Rinv,ns,0.0d0,FtRinv,fdim)

            ! F**T * Rinv * F
            call dgemm('n','n',fdim,fdim,ns,1.0d0,FtRinv,fdim,F,ns,0.0d0,FtRinvF,fdim)

            ! F**T * Rinv * Psi
            ! fdim,ns * ns,1 = fdim,1
            call dgemm('n','n',fdim,1,ns,1.0d0,FtRinv,fdim,psi,ns,0.0d0,FtRinvPsi,fdim)

            ! FTRinvF zeta = FtRinvPsi
            zeta = FtRinvPsi
            
            call LA_GESV(FtRinvF,zeta) 
         end subroutine

         ! DRHS=Rinv*(Psi-F*zeta);
         subroutine construct_DRHS(DRHS,F,Rinv,Psi,zeta,D,Ns,fdim)
            integer, intent(in) :: D, Ns, fdim
            double precision, intent(in) :: F(ns,fdim), Rinv(ns,ns)
            double precision, intent(in) :: Psi(ns,1), zeta(fdim,1)
            double precision, intent(out) :: DRHS(ns,1)
            !work
            double precision :: PsiMFzeta(ns,1)
            
            ! PsiMFzeta <- (-1.0d0) * F * zeta + 1.0d0 Psi
            PsiMFzeta = Psi
            call dgemm('n','n',ns,1,fdim,(-1.0d0),F,ns,zeta,fdim,(1.0d0),PsiMFzeta,ns)

            ! DRHS = Rinv * PsiMFzeta
            call dgemm('n','n',ns,1,ns,1.0d0,Rinv,ns,PsiMFzeta,ns,0.0d0,DRHS,ns)
         end subroutine

         ! abs(fx*zeta+r'*DRHS);
         double precision function get_S_j(fx,zeta,rx,DRHS,D,Ns,fdim)
            integer, intent(in) :: D, Ns, fdim
            double precision, intent(in) :: fx(fdim,1),zeta(fdim,1),rx(ns,1)
            double precision, intent(in) :: DRHS(ns,1)
            ! work
            double precision :: fxTzeta(1,1), res(1,1)

            ! fx**T * zeta
            call dgemm('t','n',1,1,fdim,1.0d0,fx,fdim,zeta,fdim,0.0d0,fxtzeta,1)

            ! res = abs(fx*zeta+r'*DRHS);
            res = fxtzeta
            call dgemm('t','n',1,1,ns,1.0d0,rx,ns,DRHS,ns,1.0d0,res,1)
            
            get_S_j = abs(res(1,1))
         end function
      ! -----------------------------------------------------------------
      ! end - private functions to compute the sentitivity 
      ! -----------------------------------------------------------------
        
end module

!program testsensitivity
!   use sensitivity, only: construct_sensitivity
!   implicit none
!   integer :: D, Ns, NsNew, order=2
!   double precision,allocatable :: S(:,:),X(:,:),Y(:,:),Grad(:,:),F(:,:),R(:,:)
!   double precision,allocatable :: theta(:), XNEW(:,:) 
!   integer :: ii,fdim
!
!   ! D
!   open(unit=15,file='csvd.dat',status='old')
!      read(15,*) D
!      print*, 'd= ', D
!   close(15)
!
!   ! ns
!   open(unit=15,file='csvns.dat',status='old')
!      read(15,*) ns
!      print*, 'ns= ', ns
!   close(15)
!
!   ! nsnew
!   open(unit=15,file='csvnsnew.dat',status='old')
!      read(15,*) nsnew
!      print*, 'nsnew =', nsnew
!   close(15)
!
!   fdim = 1+Order*D
!   print*, 'fdim  =', fdim 
!
!   allocate(X(ns,D))
!   allocate(XNEW(nsnew,D))
!   allocate(Y(ns,1))
!   allocate(Grad(ns,D))
!   allocate(F(ns,fdim))
!   allocate(R(ns,ns))
!   allocate(theta(D))
!   allocate(S(nsnew,1))
!
!   ! x
!   open(unit=15,file='csvx.dat',status='old')
!   do ii=1,D
!      read(15,*) X(:,ii)
!   enddo
!   close(15)
!
!   ! xnew
!   open(unit=15,file='csvxnew.dat',status='old')
!   do ii=1,D
!      read(15,*) XNEW(:,ii)
!   enddo
!   close(15)
!
!   ! Y
!   open(unit=15,file='csvy.dat',status='old')
!      read(15,*) Y(:,1)
!   close(15)
!
!    ! GRAD
!   open(unit=15,file='csvgrad.dat',status='old')
!   do ii=1,D
!      read(15,*) GRAD(:,ii)
!   enddo
!   close(15)
!
!   ! F
!   open(unit=15,file='csvf.dat',status='old')
!   do ii=1,fdim
!      read(15,*) F(:,ii)
!   enddo
!   close(15)
!
!   ! R
!   open(unit=15,file='csvr.dat',status='old')
!   do ii=1,ns
!      read(15,*) R(:,ii)
!   enddo
!   close(15)
!
!   ! theat
!   open(unit=15,file='csvtheta.dat',status='old')
!      read(15,*) theta
!      print*, 'theta =', theta
!   close(15)
!
!   call construct_sensitivity(S,XNEW,Y,X,Grad,F,R,theta,D,Ns,NsNew)
!   do ii=1,15
!    print*, 'S(1)= ', S(ii,1) 
!   enddo
!end program

!program testsampling
!   use sensitivity, only : get_sampling_radius
!   double precision :: xmax(2), xmin(2)
!   xmax = (/4,3/)
!   xmin = (/0,0/)
!   print *, get_sampling_radius(xmin,xmax,2,4)
!end program

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
