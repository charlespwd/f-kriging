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
   IMPLICIT NONE
   
   ! Arguments
   INTEGER, INTENT(IN) :: D, Ns
   DOUBLE PRECISION, INTENT(IN) :: X(Ns,Ns), Y(Ns,1)
   DOUBLE PRECISION, INTENT(INOUT) :: theta(D)

   ! Work variables
   DOUBLE PRECISION :: F(Ns,1+D*ORDER) ! Regression Matrix
   DOUBLE PRECISION :: R(Ns,Ns), DR(Ns,Ns,D) ! Correlation and derivative
   DOUBLE PRECISION :: Rinv_DR(Ns,Ns,D) ! Rinv * R
   DOUBLE PRECISION :: trace_RDR(D,1) ! the trace of (Rinv*R)
   DOUBLE PRECISION :: b_mat(D,D) ! as per MLE
   DOUBLE PRECISION :: ID(D,D) ! DxD identity matrix
   DOUBLE PRECISION :: I(Ns,Ns) ! Ns x Ns identity matrix
   DOUBLE PRECISION :: sigma2 ! sigma^2 <BS>
   DOUBLE PRECISION :: RCOND ! condition # reciprocal
   INTEGER :: ii,dd,kk
   INTEGER :: INFO
   
   ! make identity matrices
   call eye(I,nsnap)
   call eye(ID,D)

   DO ii=1,MAXCOUNT
      call construct_FMAT(F,X,Order,D,Ns)
      call construct_R(R,theta,X,D,Ns,Pc)
      call construct_DR(DR,R,X,D,Ns)
      call invertR(R,I,Rinv,Ns)

   

END SUBROUTINE
