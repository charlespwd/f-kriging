program sprogram
   use analytical_functions, only: y_gradient
   use analytical_solver, only:solver
   use grid, only:vector_grid, columngrid, LHS, grow2d
   use utils, only: printer, process_command_input, l1error
   use matrix, only: vector_range
   use sequential_sampling, only : get_sampling_radius, &
      construct_density_function, samplingcriterion, sampler
   implicit none
   integer :: D=2, Ns, NsNew, ngrid, deltans, nfinal, loopcount, mode, Order
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
   double precision :: maxerror
   real :: start, finish
   character(len=20) :: rsfile='d_rs.dat', truefile='d_true.dat', dotsfile='d_dots.dat'
   character(len=20) :: cL1ErrorFile='./data/c_l1err.dat', cLinfErrorFile='./data/c_linferr.dat'
   character(len=20) :: func_name, errfile
   character(len=50) :: datadir
   integer :: ii, cl1, clinf
   logical :: exists

   cl1 = 70
   clinf = 71
   inquire(file=adjustl(cl1errorfile), exist=exists)
   if (exists) then
      open(cl1,file=adjustl(cl1errorfile), status="old", position="append", action="write")
      open(clinf,file=adjustl(clinferrorfile), status="old", position="append", action="write")
   else
      open(cl1,file=adjustl(cl1errorfile), status="new")
      open(clinf,file=adjustl(clinferrorfile), status="new")
   end if

   errfile='./data/e.dat'
   deltans=4
   nfinal=30
   loopcount = 0
   datadir = 'data'
   Order=0 !default value
   ! set default values or get arguments from command line
   allocate(xmin(d))
   allocate(xmax(d))
   call process_command_input(func_name,Ns,ngrid,xmin,xmax, & 
      nfinal=nfinal, mode=mode,deltans=deltans, Order=Order) 
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
   call cpu_time(start)
   call solver(xnew,ynew,theta,mse,xmin,xmax,x,y,Grad,Order,D,Ns,NsNew,func_name,&
      S=S) 

    ! make fancy graphs
   call printer(x,y,ns,1,D,dotsfile,datadir,loopcount)
   call printer(xnew,ynew,nsnew,ngrid,D,rsfile,datadir,loopcount)
   call printer(xnew,ytrue,nsnew,ngrid,D,truefile,datadir,loopcount)
   open(67,file=adjustl(errfile),status='replace')
   ! print error to file
   maxerror = l1error(ytrue,ynew,nsnew) 
   call cpu_time(finish)
   write(67, *) ns, l1error(ytrue,ynew,nsnew), &
      l1error(ytrue,ynew,nsnew) / maxerror, finish-start
   
   loopcount = 0
   ! sequential sampling loop
   do while (ns<=nfinal)
      print*, 'ns = ', ns 
           
      ! update the sampling radius because ns changed
      samplingradius = get_sampling_radius(xmin,xmax,D,ns)

      ! construct the new density function
      call construct_density_function(psi,xnew,x,xmax,xmin,D,nsnew,ns)

      ! construct the new eest stack
      call samplingcriterion(eest,mse,nsnew,psi,mode,S=S)

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
      call solver(xnew,ynew,theta,mse,xmin,xmax,x,y,grad,Order,D,Ns,NsNew,func_name, &
         S=S,DELTANS=deltans)
      
      loopcount = loopcount+1
      ! make fancy graphs
      call printer(x,y,ns,1,D,dotsfile,datadir,loopcount)
      call printer(xnew,ynew,nsnew,ngrid,D,rsfile,datadir,loopcount)
      call cpu_time(finish)
      write(67, *) ns, l1error(ytrue,ynew,nsnew), &
         l1error(ytrue,ynew,nsnew) / maxerror, finish-start

   enddo

   write(cl1, *) ngrid, l1error(ytrue,ynew,nsnew) / nsnew, nsnew
   write(clinf, *) ngrid, maxval(abs(ytrue-ynew)), nsnew 

   close(67)
   open(unit=67,file='loopcount.dat',status='replace')
   write(67,*) loopcount
   close(67)

end program


