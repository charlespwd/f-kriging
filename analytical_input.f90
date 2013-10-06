SUBROUTINE analytical_solver(XNEW,YNEW,theta,MSE,XMIN,XMAX,X,Y,D,Ns,NsNew,func_name)
   USE PARAMS, ONLY: Raug
   USE LHSU, ONLY:LHS
   USE ANALYTICAL_FUNCTIONS, ONLY: Y_GRADIENT
   IMPLICIT NONE
   ! arguments
   INTEGER, INTENT(IN) :: D,Ns,NsNew
   DOUBLE PRECISION, INTENT(IN) :: XNEW(NSNEW,D)
   DOUBLE PRECISION, INTENT(IN) :: XMIN(D), XMAX(D)
   DOUBLE PRECISION, INTENT(INOUT) :: theta(D)
   DOUBLE PRECISION, INTENT(OUT) :: MSE
   DOUBLE PRECISION, INTENT(OUT) :: YNEW(NsNew,1)
   CHARACTER(len=20),INTENT(IN) :: func_name

   ! work variables
   DOUBLE PRECISION :: YGRAD(Ns,D+1)
   DOUBLE PRECISION :: GRAD(Ns,D)
   DOUBLE PRECISION, intent(out) :: Y(Ns,1), X(Ns,D)
   
   X(1:Ns-4,:) = LHS(XMIN,XMAX,D,Ns-4);
   ! Strong corners (makes prettier graphs)
   X(Ns-3,:) = (/XMIN(1),XMIN(2)/)
   X(NS-2,:) = (/XMIN(1),XMAX(2)/)
   X(NS-1,:) = (/XMAX(1),XMIN(2)/)
   X(NS,:) = (/XMAX(1),XMAX(2)/) 

   YGRAD = Y_GRADIENT(X,D,Ns,func_name)

   Y(1:NS,1) = YGRAD(1:NS,1)
   GRAD(1:NS,1:D) = YGRAD(1:NS,2:(D+1))

   CALL COKRIGING(XNEW,YNEW,theta,MSE,X,Y,GRAD,Raug,D,Ns,NsNew)
END SUBROUTINE
