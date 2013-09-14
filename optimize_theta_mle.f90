!!
! Optimizes the theta extra-parameter with Maximum Likelihood Estimation
! The method used is iterative as described in Arthur's thesis.
!
! Arguments:
!  Theta(inout) : in : the initial guesses at theta
!               : out: the optimized set of extraparameters
!  X(in) : the set of snapshots
!  Y(in) : the set of snapshot values
!  D(in) : # of dimesnions in domain
!  Ns(in) : # of snapshots
SUBROUTINE optimize_theta_mle(theta,X,Y,D,Ns)
   USE PARAMS, ONLY: ORDER, Pc, MAXCOUNT
   USE LA_PRECISION, ONLY:WP=>DP
   USE F95_LAPACK, ONLY:LA_GESV 
   IMPLICIT NONE
   
   ! Arguments
   INTEGER, INTENT(IN) :: D, Ns
   DOUBLE PRECISION, INTENT(IN) :: X(Ns,Ns), Y(Ns,1)
   DOUBLE PRECISION, INTENT(INOUT) :: theta(D)

   ! Functions
   DOUBLE PRECISION :: trace, get_sigma2

   ! Work variables
   DOUBLE PRECISION :: F(Ns,1+D*ORDER) ! Regression Matrix
   DOUBLE PRECISION :: R(Ns,Ns), DR(Ns,Ns,D) ! Correlation and derivative
   DOUBLE PRECISION :: Rinv(Ns,Ns), RinvDRij(Ns,Ns) 
   DOUBLE PRECISION :: Rinv_DR(Ns,Ns,D) ! Rinv * R
   DOUBLE PRECISION :: YmFb(Ns,1) ! Y-F*beta
   DOUBLE PRECISION :: beta(1+Order*D,1) 
   DOUBLE PRECISION :: delta(D,1) ! stores del ln / del theta_old
   DOUBLE PRECISION :: trace_RDR(D,1) ! the trace of (Rinv*R)
   DOUBLE PRECISION :: B(D,D), Binv(D,D) ! as per MLE
   DOUBLE PRECISION :: ID(D,D) ! DxD identity matrix
   DOUBLE PRECISION :: I(Ns,Ns) ! Ns x Ns identity matrix
   DOUBLE PRECISION :: sigma2 ! sigma^2 <BS>
   DOUBLE PRECISION :: RCOND ! condition # reciprocal
   INTEGER :: ii,jj,dd,kk
   INTEGER :: INFO,fdim

   fdim = 1+Order*D
   
   ! make identity matrices
   call eye(I,ns)
   call eye(ID,D)

   !DO ii=1,MAXCOUNT
      call construct_FMAT(F,X,Order,D,Ns)
      call construct_R(R,theta,X,D,Ns,Pc)
      call construct_DR(DR,R,X,D,Ns)
      call invertR(R,I,Rinv,Ns)
      call construct_beta(beta,F,Rinv,Y,Order,D,Ns)

      YmFb = Y
      ! YmFb = Y-F*beta
      call DGEMM('n','n',ns,1,fdim,(-1.0d0),F,Ns,beta,fdim,1.0d0,YmFb,Ns)

      sigma2 = get_sigma2(YmFb,Rinv,Order,D,Ns)
      print *, 'sigma2 = ', sigma2
      ! ------------------------------------------------
      ! Maximum likelihood as per Lappo
      ! Intention : theta_new = theta_old + Binv * delta
      
      DO dd=1,D
         ! Rinv_DR(:,:,d) := Rinv*DR(:,:,d);
         call DGEMM('n','n',ns,ns,ns,1.0D0,Rinv,Ns,DR(:,:,dd),Ns,0.0d0,Rinv_DR(:,:,dd),Ns)
         trace_RDR(dd,1) = trace(Rinv_DR(:,:,dd),ns)
      END DO

      ! construct the B_matrix
      DO jj=1,D ! j before i because fortran is column major.
      DO ii=1,D
         ! RinvDRij := Rinv_DR(:,:,i) * Rinv_DR(:,:,j)
         call DGEMM('n','n',Ns,Ns,Ns,1.0D0,Rinv_DR(:,:,ii),Ns,Rinv_DR(:,:,jj),Ns,0.0D0,RinvDRij,Ns)
         B(ii,jj) = 0.5d0 * trace(RinvDRij,Ns)
      END DO
      END DO
      call dumpmat(B,D,D)
      
      Binv = ID! to prevent overwrite on LAPACK routine, although we don't really
               ! need it, it's more for readability. If you don't have room for a
               ! 2x2 matrix, you have a serious RAM problem...     
     
      ! Binv = B\ID
      print *, 'bb' 
      call LA_GESV(B,Binv) 
      ! careful if you use B, it's no longer the same matrix. 
      
      !DO dd=1,D
         !delta(1,d) = -0.5*(trace_RDR(d) - (1.0 / sigma2) *...
			!	((Y-F*beta)'*Rinv_DR(:,:,d)*Rinv*(Y-F*beta)));


END SUBROUTINE
