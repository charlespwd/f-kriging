module analytical_solver
   contains

   SUBROUTINE solver(XNEW,YNEW,theta,MSE,XMIN,XMAX,X,Y,GRAD,Order,D,Ns,NsNew,&
         func_name, S, deltans)
      USE PARAMS, ONLY: Raug
      USE ANALYTICAL_FUNCTIONS, ONLY: Y_GRADIENT
      use cokrigingmodule, only:cokriging
      IMPLICIT NONE
      ! arguments
      INTEGER, INTENT(IN) :: D,Ns,NsNew,Order
      integer, intent(in), optional :: deltans
      DOUBLE PRECISION, INTENT(IN) :: XNEW(NSNEW,D)
      DOUBLE PRECISION, INTENT(IN) :: XMIN(D), XMAX(D)
      DOUBLE PRECISION, INTENT(INOUT) :: theta(D)
      DOUBLE PRECISION, INTENT(OUT) :: MSE(NsNew)
      DOUBLE PRECISION, INTENT(OUT) :: YNEW(NsNew,1)
      DOUBLE PRECISION, intent(inout) :: GRAD(Ns,D)
      DOUBLE PRECISION, INTENT(OUT), optional, target :: S(nsnew,1)
      CHARACTER(len=20),INTENT(IN) :: func_name

      ! work variables
      DOUBLE PRECISION :: YGRAD(Ns,D+1)
      double precision,allocatable :: tmpgrad(:,:)
      DOUBLE PRECISION, intent(inout) :: Y(Ns,1)
      DOUBLE PRECISION, intent(in) :: X(Ns,D)

      ! only calculate grad and y @ new locations
      if (present(deltans)) then
         allocate(tmpgrad(deltans,D+1))
         ! copy old stuff
         YGRAD(1:ns-deltans,1) = Y(1:ns,1)
         YGRAD(1:ns-deltans,2:D+1) = Grad(1:ns,1:D)
         ! calculate new stuff
         tmpgrad = Y_gradient(X(ns-deltans+1:ns,:),D,deltans,func_name)
         ! put it in ygrad
         YGRAD(ns-deltans+1:ns,1:D+1) = tmpgrad(1:deltans,1:D+1)
      else      
         YGRAD = Y_GRADIENT(X,D,Ns,func_name)
      endif

      Y(1:NS,1) = YGRAD(1:NS,1)
      GRAD(1:NS,1:D) = YGRAD(1:NS,2:(D+1))

      if (present(S)) then
         CALL COKRIGING(XNEW,YNEW,theta,MSE,X,Y,GRAD,Raug,Order,D,Ns,NsNew,S)
      else 
         CALL COKRIGING(XNEW,YNEW,theta,MSE,X,Y,GRAD,Raug,Order,D,Ns,NsNew)
      endif

   END SUBROUTINE
end module
