! construct_sensitivity(S,XNEW,Y,X,Grad,theta,D,Ns,NsNew)
module sensitivity   
   use la_precision, only:wp=>dp
   use f95_lapack, only:LA_GESV
   use matrix, only : eye, normalize
   use PARAMS, only:Pc
   use correlation, only:get_rxy,invertr,construct_R
   use regression, only:construct_f,construct_fmat

   implicit none

   contains
      ! construct the sensitivity vectort
      SUBROUTINE construct_sensitivity(S,XNEW,Y,X,Grad,theta,Order,D,Ns,NsNew)
         implicit none

         ! ARGUMENTS
         integer, intent(in) :: D, Ns, NsNew, Order
         double precision, intent(in) :: X(Ns,D), XNEW(NsNew,D)
         double precision, intent(in) :: Y(Ns,1)
         double precision, intent(in) :: Grad(Ns,D)
         double precision, intent(in) :: theta(D)
         double precision, intent(out) :: S(NsNew,1)

         ! Work variables
         double precision :: F(NS,1+Order*D)
         double precision :: R(Ns,Ns)
         double precision :: Rinv(Ns,Ns), I(Ns,Ns)
         double precision :: fx(1+Order*D,1)
         double precision :: rx(Ns,1)
         double precision :: res(1,1)
         double precision :: psi(ns,1) 
         double precision :: gradsum(ns,1)
         double precision :: zeta(1+Order*D,1)
         double precision :: DRHS(ns,1)

         integer :: ii,jj,dd,fdim
         fdim = 1 + Order*D
         
         call construct_fmat(F,X,Order,D,Ns)
         call construct_R(R,theta,X,D,Ns,Pc)
         call eye(I,Ns)
         call invertR(R,I,Rinv,Ns)

         ! construct gradsum
         do jj=1,ns
            gradsum(jj,1) =0.0d0
            do dd=1,D
               gradsum(jj,1) = gradsum(jj,1) + abs(grad(jj,dd))
            enddo
         enddo
         call normalize(gradsum(:,1),ns)

         do ii=1,NsNew
            ! construct f(x)
            call construct_f(fx(:,1),XNEW(ii,:),Order,D,Ns)
            
            ! construct r(x)
            do jj=1,ns
               rx(jj,1) = get_rxy(theta,xnew(ii,:),x(jj,:),D,Pc)
            end do

            ! construct S 
            psi=0
            S(ii,1) = 0.0d0
            do jj=1,ns
               ! psi_jj = [0,...,0,sum(grad)_jj,0,...,0]
               psi(jj,1) = gradsum(jj,1)
               ! zeta=(F'*Rinv*F)\(F'*Rinv*Psi);
               call construct_zeta(zeta,F,Rinv,Psi,D,Ns,fdim)
               ! DRHS = Rinv*(Psi-F*zeta);
               call construct_DRHS(DRHS,F,Rinv,Psi,zeta,D,Ns,fdim)
               ! Si = sum (abs(fx*zeta+r'*DRHS))
               S(ii,1) = S(ii,1) + get_S_j(fx,zeta,rx,DRHS,D,Ns,fdim)
               psi(jj,1) = 0.0d0
            enddo
         end do 
         call normalize(S(:,1),nsnew)
      end SUBROUTINE

      ! -----------------------------------------------------------------
      ! private functions to compute the sentitivity 
      ! -----------------------------------------------------------------
         ! zeta=(F'*Rinv*F)\(F'*Rinv*Psi);
         SUBROUTINE construct_zeta(zeta,F,Rinv,Psi,D,Ns,fdim)
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
         end SUBROUTINE

         ! DRHS=Rinv*(Psi-F*zeta);
         SUBROUTINE construct_DRHS(DRHS,F,Rinv,Psi,zeta,D,Ns,fdim)
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
         end SUBROUTINE

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

