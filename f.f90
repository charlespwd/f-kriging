PROGRAM f
   USE ANALYTICAL_FUNCTIONS, ONLY: Y_GRADIENT
   USE grid, only:vector_grid
   USE utils, only: printer, process_command_input
   IMPLICIT NONE
   integer :: D=2, Ns, NsNew, ngrid
   double precision,allocatable :: xnew(:,:),ynew(:,:)
   double precision,allocatable :: x(:,:),y(:,:)
   double precision,allocatable :: linspace(:,:)
   double precision,allocatable :: xmin(:), xmax(:)
   double precision,allocatable :: ygrad(:,:)
   double precision,allocatable :: ytrue(:,:)
   double precision,allocatable :: theta(:)
   double precision :: MSE
   double precision :: MeanL1
   character(len=20) :: rsfile='d_rs.dat', truefile='d_true.dat', dotsfile='d_dots.dat'
   character(len=20) :: func_name
   integer :: ii

   ! set default values or get arguments from command line
   allocate(xmin(d))
   allocate(xmax(d))
   call process_command_input(func_name,Ns,ngrid,xmin,xmax) 
   nsnew = ngrid ** D
   allocate(linspace(ngrid,D))
   allocate(xnew(nsnew,D))
   allocate(ynew(nsnew,1))
   allocate(ygrad(nsnew,D+1))
   allocate(ytrue(nsnew,1))
   allocate(x(ns,D))
   allocate(y(ns,1))
   allocate(theta(d))

   linspace(1:ngrid,1) = (/(xmin(1) + (ii-1) * ((xmax(1)-xmin(1))/(ngrid-1)), ii=1, ngrid)/) 
   linspace(1:ngrid,2) = (/(xmin(2) + (ii-1) * ((xmax(2)-xmin(2))/(ngrid-1)), ii=1, ngrid)/)
   call vector_grid(xnew,linspace,D,ngrid)

   ygrad = Y_GRADIENT(xnew,D,nsnew,func_name)
   ytrue(1:NsNew,1) = ygrad(1:nsnew,1)

   theta = (/-1,-1/)
   call analytical_solver(xnew,ynew,theta,mse,xmin,xmax,x,y,D,Ns,NsNew,func_name)
   
   ! make fancy graphs
   call printer(x,y,ns,1,D,dotsfile)
   call printer(xnew,ynew,nsnew,ngrid,D,rsfile)
   call printer(xnew,ytrue,nsnew,ngrid,D,truefile)
   MeanL1 = 0.d0
   do ii=1,NsNew
      MeanL1 = MeanL1 + abs(ytrue(ii,1) - ynew(ii,1))
      ! print *, '%err: ',(ytrue(ii,1) - ynew(ii,1))/ytrue(ii,1)*100
   end do
   print *, 'L1 error: ', MeanL1
   print *, 'Mean L1 Error: ', MeanL1/NsNew;

END PROGRAM


