PROGRAM f
   USE ANALYTICAL_FUNCTIONS, ONLY: Y_GRADIENT
   USE grid, only:vector_grid, columngrid, LHS
   USE utils, only: printer, process_command_input
   IMPLICIT NONE
   integer :: D=2, Ns, NsNew, ngrid
   double precision,allocatable :: xnew(:,:),ynew(:,:)
   double precision,allocatable :: x(:,:),y(:,:)
   double precision,allocatable :: xmin(:), xmax(:)
   double precision,allocatable :: ygrad(:,:)
   double precision,allocatable :: ytrue(:,:)
   double precision,allocatable :: theta(:)
   double precision,allocatable :: MSE(:)
   double precision :: MeanL1
   character(len=20) :: rsfile='d_rs.dat', truefile='d_true.dat', dotsfile='d_dots.dat'
   character(len=20) :: func_name
   character(len=20) :: datadir
   integer :: ii

   datadir = 'data'
   ! set default values or get arguments from command line
   allocate(xmin(d))
   allocate(xmax(d))
   call process_command_input(func_name,Ns,ngrid,xmin,xmax) 
   nsnew = ngrid ** D
   allocate(xnew(nsnew,D))
   allocate(ynew(nsnew,1))
   allocate(ygrad(nsnew,D+1))
   allocate(ytrue(nsnew,1))
   allocate(x(ns,D))
   allocate(y(ns,1))
   allocate(theta(d))
   allocate(mse(nsnew))

   call columngrid(xnew,xmin,xmax,d,ngrid)

   ygrad = Y_GRADIENT(xnew,D,nsnew,func_name)
   ytrue(1:NsNew,1) = ygrad(1:nsnew,1)
   
   ! make initial grid
   X(1:Ns-4,:) = LHS(XMIN,XMAX,D,Ns-4);
   ! Strong corners (makes prettier graphs)
   X(Ns-3,:) = (/XMIN(1),XMIN(2)/)
   X(NS-2,:) = (/XMIN(1),XMAX(2)/)
   X(NS-1,:) = (/XMAX(1),XMIN(2)/)
   X(NS,:) = (/XMAX(1),XMAX(2)/) 

   theta = (/-1,-1/)
   call analytical_solver(xnew,ynew,theta,mse,xmin,xmax,x,y,D,Ns,NsNew,func_name)
   
   ! make fancy graphs
   call printer(x,y,ns,1,D,dotsfile,datadir)
   call printer(xnew,ynew,nsnew,ngrid,D,rsfile,datadir)
   call printer(xnew,ytrue,nsnew,ngrid,D,truefile,datadir)
   MeanL1 = 0.d0
   do ii=1,NsNew
      MeanL1 = MeanL1 + abs(ytrue(ii,1) - ynew(ii,1))
      ! print *, '%err: ',(ytrue(ii,1) - ynew(ii,1))/ytrue(ii,1)*100
   end do
   print *, 'L1 error: ', MeanL1
   print *, 'Mean L1 Error: ', MeanL1/NsNew/(maxval(ytrue(:,1))-minval(ytrue(:,1)));

END PROGRAM


