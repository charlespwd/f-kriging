! SUBROUTINE COKRIGING(XNEW,YNEW,theta,MSE,XOLD,YOLD,Grad,Raug,D,Ns,NsNew)
!  This routine is a wrapper to kriging, but augments the data set with 
!  information about the gradient around the snapshots
!  
! Arguments:
!  XNEW (in): 'new snapshots' location, [NsNew,D]
!  YNEW (out): 'new snapshots' values, [NsNew,1]
!  theta(inout), enters as theta initial, exits as optimized theta
!  MSE (out): error
!  XOLD (in): set of snapshot locations [Ns,D]
!  YOLD (in): set of snapshot values [Ns,1]
!  Grad (in): the gradient of Y at the snapshot locations in d'th direction [Ns,D]
!  Raug (in): the 'radius' of the augmentation as fraction of range(X), typically 1% 
!  D (in): # of dimensions
!  Ns (in): # of snapshots
!  NsNew: # of new snapshots
!
!  Cokriging works by augmenting the initial data set as follows:
!  For each snasphots, an increment Raug*range(X,d) is added and subtracted
!     i.e. for each dimension X = XOLD +/- Raug*range(X)
!  For each snapshot locations, the snapshot value is derived from a taylor
!     expansion, i.e. for a corresponding X, Y = Yold + Grad * Raug * range(X)
!  After that, kriging is performed on this new, augmented, data set
SUBROUTINE COKRIGING(XNEW,YNEW,theta,MSE,XOLD,YOLD,Grad,Raug,D,Ns,NsNew)
   IMPLICIT NONE

   ! arguments
   INTEGER, INTENT(IN) :: D,Ns,NsNew
   DOUBLE PRECISION, INTENT(IN) :: XOLD(Ns,D),YOLD(Ns,1) 
   DOUBLE PRECISION, INTENT(IN) :: XNEW(NsNew,D)
   DOUBLE PRECISION, INTENT(IN) :: Grad(Ns,D)
   DOUBLE PRECISION, INTENT(IN) :: Raug
   DOUBLE PRECISION, INTENT(INOUT) :: theta(D)
   DOUBLE PRECISION, INTENT(OUT) :: YNEW(NsNew,1)
   DOUBLE PRECISION, INTENT(OUT) :: MSE(NsNew)
  
   ! work variables
   INTEGER :: Naug
   INTEGER :: dd, nn, kk
   DOUBLE PRECISION,ALLOCATABLE :: X(:,:), Y(:,:)
   DOUBLE PRECISION :: delta(D)
   DOUBLE PRECISION :: XMAX(D),XMIN(D)
   
   Naug=(2*D+1)*Ns
   ALLOCATE(X(Naug,D))
   ALLOCATE(Y(Naug,1))
   
   XMAX = MAXVAL(XOLD,1)
   XMIN = MINVAL(XOLD,1)

   DO dd=1,D
      delta(dd) = (XMAX(dd) - XMIN(dd)) * Raug
   END DO

   X(1:Ns,:) = XOLD(1:Ns,:)
   Y(1:Ns,:) = YOLD(1:Ns,:)

   kk = Ns + 1
   DO nn = 1,Ns
      DO dd =1,D
         X(kk,:) = X(nn,:)
         X(kk,dd) = X(nn,dd) + delta(dd)
         Y(kk,1) = Y(nn,1) + delta(dd) * Grad(nn,dd)
         
         kk = kk + 1

         X(kk,:) = X(nn,:)
         X(kk,dd) = X(nn,dd) - delta(dd)
         Y(kk,1) = Y(nn,1) - delta(dd) * Grad(nn,dd)

         kk = kk + 1
      END DO
   END DO
   
   ! perform kriging normally with new set of X and Y 
   CALL KRIGING(XNEW,YNEW,theta,MSE,X,Y,D,Naug,NsNew)

END SUBROUTINE
