PROGRAM f
   USE ANALYTICAL_FUNCTIONS, ONLY: Y_GRADIENT
   IMPLICIT NONE
   integer :: D=2, Ns=20, NsNew=5
   double precision,allocatable :: xnew(:,:),ynew(:,:)
   double precision,allocatable :: xmin(:), xmax(:)
   double precision,allocatable :: ygrad(:,:)
   double precision,allocatable :: analytical_solution(:,:)
   double precision,allocatable :: theta(:)
   double precision :: MSE
   character(len=20) :: func_name
   integer :: ii

   allocate(xnew(nsnew,D))
   allocate(ynew(nsnew,1))
   allocate(ygrad(nsnew,D+1))
   allocate(analytical_solution(nsnew,1))
   allocate(theta(d))
   allocate(xmin(d))
   allocate(xmax(d))
   
   func_name = "--drag"
   xnew(1:NsNew,1) = 0.0d0
   xnew(1:nsNew,2) = (/(ii,ii=1,nsnew)/)

   ! so to get [-5:5]x[0,15]
   xmin = (/-5,0/) 
   xmax = (/5,15/)
   
   ygrad = Y_GRADIENT(xnew,D,nsnew,func_name)
   analytical_solution(1:NsNew,1) = ygrad(1:nsnew,1)

   theta = (/-1,-1/)
   call analytical_solver(xnew,ynew,theta,mse,xmin,xmax,D,Ns,NsNew,func_name)
   print *, 'analytical vs RS'
   do ii=1,NsNew
      print *, analytical_solution(ii,1),' : ', ynew(ii,1)
   end do
end program

