module error
   implicit none
   contains
!! 
! get_mse
! computed as abs(sigma2*(
! 1 - (r**T * Rinv * r) 
! + (1 - || F**T * Rinv * r ||_2 ^2) / || F**T * Rinv * F||_2 ))
double precision function get_mse(sigma2,F,Rinv,rx,D,Ns,fdim)
   use matrix, only : spectral_norm
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

 ! MARTIN HAS A TYPO IN IT, it says ^-1 in Sacks' 
 ! sigma2 * (1 - [ft(x); rt(x)]**T * [ 0 Ft ; F R]^-1 * [f(x); r(x)]) 
double precision function martin_mse(sigma2, fx, rx, F, R, ns, fdim)
   use LA_PRECISION, only: WP=>DP
   use F95_LAPACK, only: LA_GESV
   use matrix, only: eye
   implicit none
   integer, intent(in) :: fdim, ns
   double precision, intent(in) :: sigma2
   double precision, intent(in) :: fx(fdim,1), rx(ns,1)
   double precision, intent(in) :: F(ns,fdim), R(ns,ns)
   ! work variables
   double precision :: Ft(fdim,ns)
   double precision :: fxrx(fdim+ns,1)
   double precision :: ctrblock(fdim+ns,fdim+ns)
   double precision :: I(fdim+ns, fdim+ns)
   double precision :: tmp(1,fdim+ns)
   double precision :: res(1,1)
   integer :: ii, jj
   fxrx(1:fdim,1) = fx(1:fdim,1)
   fxrx(fdim+1:ns,1) = rx(1:ns,1)
   Ft = transpose(F)

   ! construct [0 Ft; F R]
   !  0
   ctrblock(1:fdim,1:fdim) = 0.d0
   !  F
   ctrblock(fdim+1:fdim+ns, 1:fdim) = F(1:ns,1:fdim)
   !  Ft 
   ctrblock(1:fdim, fdim+1:fdim+ns) = Ft(1:fdim,1:ns)
   !  R
   ctrblock(fdim+1:fdim+ns, fdim+1:fdim+ns) = R(1:ns,1:ns)

   call eye(I,fdim+ns)
   call LA_GESV(ctrblock,I)

   ! tmp = fxrx**T * cntrblock^-1
   ! [1,fdim+ns] = [fdim+ns,1] ** T * [fdim+ns,fdim+ns]
   call dgemm('t','n',1,fdim+ns,fdim+ns,1.0d0, fxrx,fdim+ns, &
      I,fdim+ns, 0.0d0,tmp,1)

   ! res = (-1.0d0) * tmp * fxrx + 1.0d0 * res
   res(1,1) = 1.0d0
   call dgemm('n','n',1,1,fdim+ns,(-1.0d0), tmp,1, fxrx,fdim+ns, &
      1.0d0, res,1)
   
   martin_mse = sigma2 * res(1,1)
end function

! sigma2 * (1 + u**T * (F**T * Rinv * F)inv * u - r**T *Rinv * r)
double precision function lophaven_mse(sigma2, fx, rx, F, R, Rinv, ns, fdim)
   use LA_PRECISION, only:WP=>DP
   use F95_LAPACK, only: LA_GESVX
   use matrix, only: eye
   implicit none
   integer, intent(in) :: ns, fdim
   double precision, intent(in) :: sigma2
   double precision, intent(in) :: fx(fdim, 1), rx(ns,1)
   double precision, intent(in) :: F(ns,fdim), R(ns,ns)
   double precision, intent(in) :: Rinv(ns,ns)
   ! work variables
   integer :: info
   double precision :: u(fdim,1), FtRinv(fdim,ns), FtRinvF(fdim,fdim)
   double precision :: FtRinvFinv(fdim,fdim), utFtRinvFinv(1,fdim)
   double precision :: rtRinv(1,ns)
   double precision :: feye(fdim,fdim)
   double precision :: res1(1,1), res2(1,1)

   call eye(feye,fdim) 

   !! u = F**T * Rinv * rx - fx
   !  suboperation F**T * Rinv
   call dgemm('t','n',fdim,ns,ns,1.0d0,F,ns,Rinv,ns,0.0d0,FtRinv,fdim)

   ! u =FtRinv * rx - fx
   u = fx
   call dgemm('n','n',fdim,1,ns,1.0d0,FtRinv,fdim,rx,ns,(-1.0d0),u,fdim)

   !! FTRinv*F inv
   !  suboperation FtRinv * F
   call dgemm('n','n',fdim,fdim,ns,1.0d0,FtRinv,fdim,F,ns,0.0d0,FtRinvF,fdim)

   ! invert FtRinvF
   call LA_GESVX(FtRinvF,feye,FtRinvFinv,INFO=info) 

   !! u**T * FtRinvFinv * u
   ! suboperation : u**T * FtRinvFinv
   call dgemm('t','n',1,fdim,fdim,1.0d0,u,fdim,FtRinvFinv,fdim, & 
      0.0d0, utFtRinvFinv,1)

   ! utFtRinvFinv * u
   call dgemm('n','n',1,1,fdim,1.0d0,utftrinvfinv,1,u,fdim,0.0d0,res1,1)

   !! rx**T * Rinv * rx
   ! suboperation rx ** T * Rinv
   call dgemm('t','n', 1, ns, ns, 1.0d0, rx, ns, Rinv, ns, 0.0d0, rtRinv, 1)

   ! rtRinv*rx
   call dgemm('n','n', 1, 1, ns, 1.0d0, rtRinv, 1, rx, ns, 0.0d0, res2, 1)

   lophaven_mse = sigma2 * (1 + res1(1,1) - res2(1,1))
end function
end module


