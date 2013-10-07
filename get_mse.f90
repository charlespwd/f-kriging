!! 
! get_mse
! computed as abs(sigma2*(
! 1 - (r**T * Rinv * r) 
! + (1 - || F**T * Rinv * r ||_2 ^2) / || F**T * Rinv * F||_2 ))
double precision function get_mse(sigma2,F,Rinv,rx,D,Ns,fdim)
   use matrixmath, only : spectral_norm
   implicit none
   
   ! inputs
   integer, intent(in) :: D, Ns
   integer, intent(in) :: fdim
   double precision, intent(in) :: sigma2
   double precision, intent(in) :: F(Ns,fdim)
   double precision, intent(in) :: Rinv(Ns,Ns)
   double precision, intent(in) :: rx(Ns,1)
   
   ! work variables
   double precision :: rtRinv(1,Ns), rtRinvr(1,1), FtRinvR(fdim,1)
   double precision :: FtRinv(fdim,Ns), FtRinvF(fdim,fdim)
   double precision :: term1, term2, term3, term4
   integer :: foo=0

   !get RtRinvr
   ! r**T * Rinv
   ! [1,ns] = [ns,1]**T * [ns,ns]
   call dgemm('t','n',1,ns,ns,1.0d0,rx,ns,Rinv,ns,0.0d0,rtRinv,1)

   ! rtRinv * rx
   ! [1,1] = [1,ns] [ns,1]
   call dgemm('n','n',1,1,ns,1.0d0,rtRinv,1,rx,ns,0.0d0,rtRinvr,1)
   
   !get F'RinvF
   ! F**T * Rinv
   ! [fdim, Ns] = [ns,fdim] ** t * [ns,ns]
   call dgemm('t','n',fdim,ns,ns,1.0d0,F,ns,Rinv,ns,0.0d0,FtRinv,fdim)

   ! FtRinv * F
   ! [fdim,fdim] = [fdim,ns] [ns,fdim]
   call dgemm('n','n',fdim,fdim,ns,1.0d0,FtRinv,fdim,F,ns,0.0d0,FtRinvF,fdim)

   ! FtRinv * rx
   ! [fdim,1] = [fdim,ns] * [ns,1]
   call dgemm('n','n',fdim,1,ns,1.0d0,FtRinv,fdim,rx,ns,0.0d0,FtRinvR,fdim)

   ! CALCULATE MSE..........
   term1 = 1.0d0
   term2 = rtRinvr(1,1)
   term3 = (1.0d0 - spectral_norm(FtRinvR,fdim,1)) ** 2
   term4 = 1.0d0 / spectral_norm(FtRinvF,fdim,fdim)
   term3 = term3 * term4
   get_mse = abs(sigma2 * (term1 - term2 + term3))
end function 

! unit test.
!program test
!   integer :: ns=189, fdim=6.
!   double precision :: sigma2
!   double precision,allocatable :: rx(:,:), F(:,:), Rinv(:,:)
!   double precision :: amse, mse
!   double precision :: get_mse
!   integer :: i
!   allocate(rx(ns,2))
!   allocate(F(ns,fdim))
!   allocate(Rinv(ns,ns))
!
!   open(unit=115,file= 'csvsigma2.dat', status='old') 
!   read(115,*)sigma2
!   print*, sigma2
!   close(unit=115)
!   open(unit=116, file='csvrx.dat', status='old')
!   read(116,*) rx
!   close(116)
!   open(unit=117,file='csvF.dat',status='old')
!   do i=1,fdim
!      read(117,*) F(:,i)
!   enddo
!   close(117)
!   open(unit=118,file='csvRinv.dat',status='old')
!   do i=1,ns
!      read(118,*) Rinv(:,i)
!   enddo
!   close(unit=118)
!   print *, get_mse(sigma2,F,Rinv,rx,3,Ns,fdim)
!end program
