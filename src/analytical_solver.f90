module analytical_solver
   contains

   subroutine solver(XNEW,YNEW,theta,MSE,XMIN,XMAX,X,Y,GRAD,Order,D,Ns,NsNew,&
         func_name, S, deltans)
      use PARAMS, only: Raug
      use ANALYTICAL_FUNCTIONS, only: Y_GRADIENT
      use cokrigingmodule, only:cokriging
      implicit none
      ! arguments
      integer, intent(in) :: D,Ns,NsNew,Order
      integer, intent(in), optional :: deltans
      double precision, intent(in) :: XNEW(NSNEW,D)
      double precision, intent(in) :: XMIN(D), XMAX(D)
      double precision, intent(INOUT) :: theta(D)
      double precision, intent(out) :: MSE(NsNew)
      double precision, intent(out) :: YNEW(NsNew,1)
      double precision, intent(inout) :: GRAD(Ns,D)
      double precision, intent(out), optional, target :: S(nsnew,1)
      character(len=20),intent(in) :: func_name

      ! work variables
      double precision :: YGRAD(Ns,D+1)
      double precision,allocatable :: tmpgrad(:,:)
      double precision, intent(inout) :: Y(Ns,1)
      double precision, intent(in) :: X(Ns,D)

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
         call COKRIGING(XNEW,YNEW,theta,MSE,X,Y,GRAD,Raug,Order,D,Ns,NsNew,S)
      else 
         call COKRIGING(XNEW,YNEW,theta,MSE,X,Y,GRAD,Raug,Order,D,Ns,NsNew)
      endif

   end subroutine
end module
