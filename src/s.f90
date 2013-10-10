PROGRAM sprogram
   use analytical_functions, only: y_gradient
   use analytical_solver, only:solver
   USE grid, only:vector_grid, columngrid, LHS, grow2d
   USE utils, only: printer, process_command_input
   use sequential_sampling, only : get_sampling_radius, &
      construct_density_function, samplingcriterion, sampler
   IMPLICIT NONE
   integer :: D=2, Ns, NsNew, ngrid, deltans, nfinal, loopcount
   integer,allocatable :: newsamples(:)
   double precision,allocatable :: xnew(:,:),ynew(:,:)
   double precision,allocatable :: x(:,:),y(:,:)
   double precision,allocatable :: grad(:,:)
   double precision,allocatable :: xmin(:), xmax(:)
   double precision,allocatable :: ygrad(:,:)
   double precision,allocatable :: ytrue(:,:)
   double precision,allocatable :: theta(:)
   double precision,allocatable :: psi(:,:)
   double precision,allocatable :: MSE(:,:)
   double precision,allocatable :: S(:,:)
   double precision,allocatable :: Eest(:,:)
   double precision :: MeanL1
   double precision :: samplingRadius 
   character(len=20) :: rsfile='d_rs.dat', truefile='d_true.dat', dotsfile='d_dots.dat'
   character(len=20) :: func_name, errfile
   character(len=20) :: datadir
   integer :: ii

   errfile='./data/e.dat'
   deltans=3
   nfinal=30
   loopcount = 0
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
   allocate(x(ns**D,D))
   allocate(y(ns**D,1))
   allocate(grad(ns**D,D))
   allocate(theta(d))
   allocate(mse(nsnew,1))
   allocate(eest(nsnew,1))
   allocate(psi(nsnew,1))
   allocate(S(nsnew,1))
   allocate(newsamples(deltans))

   call columngrid(xnew,xmin,xmax,d,ngrid)

   ygrad = Y_GRADIENT(xnew,D,nsnew,func_name)
   ytrue(1:NsNew,1) = ygrad(1:nsnew,1)

   ! make initial grid
   call columngrid(x,xmin,xmax,d,ns)
   theta = (/-1,-1/)
   ns = ns ** D
   ! perform cokriging, sensitivity analysis, etc. 
   call solver(xnew,ynew,theta,mse,xmin,xmax,x,y,Grad,D,Ns,NsNew,func_name,&
      S=S) 

    ! make fancy graphs
   call printer(x,y,ns,1,D,dotsfile,datadir,loopcount)
   call printer(xnew,ynew,nsnew,ngrid,D,rsfile,datadir,loopcount)
   call printer(xnew,ytrue,nsnew,ngrid,D,truefile,datadir,loopcount)
   open(67,file=adjustl(errfile),status='replace')
   MeanL1 = 0.d0
   do ii=1,NsNew
      MeanL1 = MeanL1 + abs(ytrue(ii,1) - ynew(ii,1))
      ! print *, '%err: ',(ytrue(ii,1) - ynew(ii,1))/ytrue(ii,1)*100
   end do
   write(67,*) 'ns= ',ns,', ngrid=',ngrid,'x',ngrid
   write(67,*) '  L1 error: ', MeanL1
   write(67,*) '  Mean L1 Error: ', MeanL1/NsNew/(maxval(ytrue(:,1))-minval(ytrue(:,1)));
   write(67,*) ''

   loopcount = 0
   ! sequential sampling loop
   do while (ns<nfinal)
      print*, 'ns = ', ns 

      ! reinit theta
      theta = (/-1,-1/)
           
      ! update the sampling radius because ns changed
      samplingradius = get_sampling_radius(xmin,xmax,D,ns)

      ! construct the new density function
      call construct_density_function(psi,xnew,x,xmax,xmin,D,nsnew,ns)

      ! construct the new eest stack
      call samplingcriterion(eest,mse,nsnew,psi,S=S)

      ! find the points to add to the set
      call sampler(newsamples,Eest,XNEW,X,samplingRadius,deltaNs,D,nsnew,ns)

      ! add the points to the existing set
      call grow2d(x,deltans)
      call grow2d(y,deltans)
      call grow2d(grad,deltans)
      do ii=1,deltans
         !! add newsamples to the set
         x(ns+ii,1:D) = xnew(newsamples(ii),1:D)    
      enddo
      
      ! increment number of snapshots
      ns = ns + deltans

      ! perform cokriging, sensitivity analysis, etc. 
      print*, func_name
      call solver(xnew,ynew,theta,mse,xmin,xmax,x,y,grad,D,Ns,NsNew,func_name, &
         S=S,DELTANS=deltans)
      
      loopcount = loopcount+1
      ! make fancy graphs
      call printer(x,y,ns,1,D,dotsfile,datadir,loopcount)
      call printer(xnew,ynew,nsnew,ngrid,D,rsfile,datadir,loopcount)
      MeanL1 = 0.d0
      do ii=1,NsNew
         MeanL1 = MeanL1 + abs(ytrue(ii,1) - ynew(ii,1))
         ! print *, '%err: ',(ytrue(ii,1) - ynew(ii,1))/ytrue(ii,1)*100
      end do
      write(67,*) 'ns= ',ns,', ngrid=',ngrid,'x',ngrid
      write(67,*) '  L1 error: ', MeanL1
      write(67,*) '  Mean L1 Error: ', MeanL1/NsNew/(maxval(ytrue(:,1))-minval(ytrue(:,1)));
      write(67,*) ''

   enddo

   close(67)

   MeanL1 = 0.d0
   do ii=1,NsNew
      MeanL1 = MeanL1 + abs(ytrue(ii,1) - ynew(ii,1))
      ! print *, '%err: ',(ytrue(ii,1) - ynew(ii,1))/ytrue(ii,1)*100
   end do
   print *, 'L1 error: ', MeanL1
   print *, 'Mean L1 Error: ', MeanL1/NsNew/(maxval(ytrue(:,1))-minval(ytrue(:,1)));

END PROGRAM


