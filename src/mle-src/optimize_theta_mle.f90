!!
! Optimizes the theta extra-parameter with Maximum Likelihood Estimation
! The method used is iterative as described in Arthur's thesis.
!
! Arguments:
!  Theta(inout) : in : the initial guesses at theta
!               : out: the optimized set of extraparameters
!  bounds(in) : the set of lower and upper bounds for theta.
!  X(in) : the set of snapshots
!  Y(in) : the set of snapshot values
!  D(in) : # of dimesnions in domain
!  Ns(in) : # of snapshots
SUBROUTINE optimize_theta_mle(theta,bounds,X,Y,D,Ns)
   USE PARAMS, ONLY: ORDER, Pc, MAXCOUNT, TOL
   USE LA_PRECISION, ONLY:WP=>DP
   USE F95_LAPACK, ONLY:LA_GESV, LA_GESVX
   use matrix, only: eye, get_trace
   use regression, only: construct_fmat
   use correlation, only: construct_rt, invertr
   IMPLICIT NONE
   
   ! Arguments
   INTEGER, INTENT(IN) :: D, Ns
   DOUBLE PRECISION, INTENT(IN) :: X(Ns,D), Y(Ns,1)
   DOUBLE PRECISION, INTENT(IN) :: bounds(D,2)
   DOUBLE PRECISION, INTENT(INOUT) :: theta(D)

   ! Work variables
   DOUBLE PRECISION :: F(Ns,1+D*ORDER) ! Regression Matrix
   DOUBLE PRECISION :: R(Ns,Ns), DR(Ns,Ns,D) ! Correlation and derivative
   DOUBLE PRECISION :: Rinv(Ns,Ns), RinvDRij(Ns,Ns) 
   DOUBLE PRECISION :: Rinv_DR(Ns,Ns,D) ! Rinv * R
   DOUBLE PRECISION :: YmFb(Ns,1) ! Y-F*beta
   DOUBLE PRECISION :: beta(1+Order*D,1) 
   DOUBLE PRECISION :: delta(D,1) ! stores del ln / del theta_old
   DOUBLE PRECISION :: increment(D,1)
   DOUBLE PRECISION :: trace_RDR(D,1) ! the trace of (Rinv*R)
   DOUBLE PRECISION :: B(D,D), Binv(D,D) ! as per MLE
   DOUBLE PRECISION :: ID(D,D) ! DxD identity matrix
   DOUBLE PRECISION :: I(Ns,Ns) ! Ns x Ns identity matrix
   DOUBLE PRECISION :: sigma2 ! sigma^2 <BS>
   DOUBLE PRECISION :: RCOND ! condition # reciprocal
   DOUBLE PRECISION :: NRM
   DOUBLE PRECISION :: XT(D,Ns)
   INTEGER :: I_MIN=1, I_MAX=2
   INTEGER :: ii,jj,dd,kk
   INTEGER :: INFO,fdim

   ! Performance hack, construct_RT makes the same job but with loops
   ! organized differently (column major). It reduced by 25-30% the cost
   ! of constructing R. 
   XT = TRANSPOSE(X)

   fdim = 1+Order*D
   
   ! make identity matrices
   call eye(I,ns)
   call eye(ID,D)

   DO kk=1,MAXCOUNT
      call construct_FMAT(F,X,Order,D,Ns)
      call construct_RT(R,theta,XT,D,Ns,Pc)
      call construct_DR(DR,R,X,D,Ns)
      call invertR(R,I,Rinv,Ns)
      call construct_beta(beta,F,Rinv,Y,Order,D,Ns)

      YmFb = Y
      ! YmFb = Y-F*beta
      call DGEMM('n','n',ns,1,fdim,(-1.0d0),F,Ns,beta,fdim,1.0d0,YmFb,Ns)

      sigma2 = get_sigma2(YmFb,Rinv,Order,D,Ns)
      ! ------------------------------------------------
      ! Maximum likelihood as per Lappo
      ! Intention : theta_new = theta_old + Binv * delta
      
      ! construct Rinv_DR and trace_RDR
      DO dd=1,D
         ! Rinv_DR(:,:,d) := Rinv*DR(:,:,d);
         call DGEMM('n','n',ns,ns,ns,1.0D0,Rinv,Ns,DR(:,:,dd),Ns,0.0d0,Rinv_DR(:,:,dd),Ns)
         trace_RDR(dd,1) = get_trace(Rinv_DR(:,:,dd),ns)
      END DO

      ! construct the B_matrix
      DO jj=1,D ! j before i because fortran is column major.
      DO ii=1,D
         ! RinvDRij := Rinv_DR(:,:,i) * Rinv_DR(:,:,j)
         call DGEMM('n','n',Ns,Ns,Ns,1.0D0,Rinv_DR(:,:,ii),Ns,Rinv_DR(:,:,jj),Ns,0.0D0,RinvDRij,Ns)
         B(ii,jj) = 0.5d0 * get_trace(RinvDRij,Ns)
      END DO
      END DO
      
      Binv = ID ! Binv = B\ID
      call LA_GESVX(B,ID,Binv,RCOND=RCOND,INFO=INFO)

      IF (INFO > 0) THEN
         print *, 'B is singular, stopped at iteration #',kk
         EXIT
      END IF 

      ! construct the delta matrix
      call construct_delta(delta,Rinv,Rinv_DR,trace_RDR,YmFb,sigma2,Order,D,Ns)

      ! Increment = Binv*delta
      ! [D,D] x [D,1]
      call DGEMM('n','n',D,1,D,1.0d0,Binv,D,delta,D,0.0d0,increment,D)
      
      ! Very important stuff. 
      DO ii=1,D
         DO WHILE (((theta(ii) + increment(ii,1)) < bounds(ii,I_MIN)) .OR. ( (theta(ii) + increment(ii,1)) > bounds(ii,I_MAX)))
            increment(ii,1) = increment(ii,1) / 2.0d0
         END DO
      END DO

      ! theta_new = theta_old + Binv*Delta
      DO ii=1,D
         theta(ii) = theta(ii) + increment(ii,1)
      END DO
     
      ! Calculate the norm of the increment,  
      NRM = 0.0d0
      DO ii=1,D
         NRM=NRM+increment(ii,1)**2
      END DO
      NRM = sqrt(NRM)

      IF (NRM < TOL) THEN
         PRINT *, 'Increment converged, count= ',kk
         EXIT
      END IF
   END DO

END SUBROUTINE
