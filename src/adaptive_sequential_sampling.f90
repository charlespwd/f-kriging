module adaptive_sequential_sampling
   use sequential_sampling, only : MODE_MSE, MODE_SENSITIVITY
   implicit none

   ! module globals ( editable parameters ) & defaults
   integer :: nGridStart = 3     ! for a nGridStart**D, evenly spaced, initial
                                 !  sample set
   integer :: delta_ns = 4       ! number of points added to sample @ new iter
   integer :: nFinal = 40        ! maximum number of points to kriging model
   integer :: nGridRS = 50       ! for a nGridRs**D response surface grid
   integer :: mode = MODE_MSE    ! the adaptive scheme's mode (MSE || MSE+S)
   integer :: order = 2          ! regression order 

   contains
      
      ! adaptive_ssck - adaptive sequential sampling cokriging routine 
      subroutine adaptive_ssck(func, grad, xminxmax, D, iterate, datadirectory, &
         fileprefix)
         use analytical_solver, only:solver
         use grid, only:vector_grid, columnGridRs, LHS, grow2d
         use utils, only: printer, process_command_input, l1error
         use matrix, only: vector_range
         use sequential_sampling, only : get_sampling_radius, &
            construct_density_function, samplingcriterion, sampler

         ! arguments
         integer, intent(in) :: D
         integer, intent(in), optional :: iterate
         double precision, intent(in) :: xminmax(D,2) 
         character(len=50), intent(in), optional :: datadirectory
         character(len=50), intent(in), optional :: fileprefix
         interface true_function
            function func(x,D) 
               integer, intent(in) :: D
               double precision, intent(in) :: x(D,1)
               double precision, intent(out) :: func
            end function
         end interface true_function
         interface true_gradient
            function grad(x,D)
               integer, intent(in) :: D
               double precision, intent(in) :: x(D,1)
               double precision, intent(out) :: grad(D,1)
            end function
         end interface true_gradient

         ! work variables
         integer :: Ns, NsNew, loopcount
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
         character(len=20) :: rsfile='d_rs.dat', truefile='d_true.dat', dotsfile='d_dots.dat'
         character(len=20) :: func_name, errfile
         character(len=20) :: datadir
         integer :: ii

         xmin(1:D) = xminmax(1:D,1)
         xmax(1:D) = xminmax(1:D,2)

         datadir = './data'
         if (present(datadirectory)) then
            datadir = datadirectory
         end if
         errfile= trim(datadirectory)//'/'//'e.dat'
         loopcount = 0
         ! set default values or get arguments from command line
         allocate(xmin(d))
         allocate(xmax(d))

         ns = nStartGrid
         nsnew = nGridRs ** D
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
         allocate(newsamples(delta_ns))

         call columnGridRs(xnew,xmin,xmax,d,nGridRs)

         ygrad = Y_GRADIENT(xnew,D,nsnew,func_name)
         ytrue(1:NsNew,1) = ygrad(1:nsnew,1)

         ! make initial grid
         call columnGridRs(x,xmin,xmax,d,ns)
         theta = (/-1,-1/)
         ns = ns ** D
         ! perform cokriging, sensitivity analysis, etc. 
         call solver(xnew,ynew,theta,mse,xmin,xmax,x,y,Grad,Order,D,Ns,NsNew,func_name,&
            S=S) 

          ! make fancy graphs
         call printer(x,y,ns,1,D,dotsfile,datadir,loopcount)
         call printer(xnew,ynew,nsnew,nGridRs,D,rsfile,datadir,loopcount)
         call printer(xnew,ytrue,nsnew,nGridRs,D,truefile,datadir,loopcount)
         open(67,file=adjustl(errfile),status='replace')
         ! print error to file
         maxerror = l1error(ytrue,ynew,nsnew) 
         write(67, *) ns, l1error(ytrue,ynew,nsnew), l1error(ytrue,ynew,nsnew) / maxerror
         
         loopcount = 0
         ! sequential sampling loop
         do while (ns<=nFinal)
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
            call grow2d(x,delta_ns)
            call grow2d(y,delta_ns)
            call grow2d(grad,delta_ns)
            do ii=1,delta_ns
               !! add newsamples to the set
               x(ns+ii,1:D) = xnew(newsamples(ii),1:D)    
            enddo
            
            ! increment number of snapshots
            ns = ns + delta_ns
            
            ! perform cokriging, sensitivity analysis, etc. 
            call solver(xnew,ynew,theta,mse,xmin,xmax,x,y,grad,Order,D,Ns,NsNew,func_name, &
               S=S,DELTANS=delta_ns)
            
            loopcount = loopcount+1
            ! make fancy graphs
            call printer(x,y,ns,1,D,dotsfile,datadir,loopcount)
            call printer(xnew,ynew,nsnew,nGridRs,D,rsfile,datadir,loopcount)
            write(67, *) ns, l1error(ytrue,ynew,nsnew), l1error(ytrue,ynew,nsnew) / maxerror

         enddo

         close(67)
         open(unit=67,file='loopcount.dat',status='replace')
         write(67,*) loopcount
         close(67)

         contains 
            function Y_GRADIENT(x, D, ns)
               ! arguments
               integer, intent(in) :: D, ns
               double precision, intent(in) :: x(Ns,D)
               double precision, intent(out) :: Y_GRADIENT(NS,1+D)
               
               ! work
               integer :: ii
               double precision :: xt(D,ns)
               double precision :: yt(1+D, ns)
               xt = transpose(x)
               do ii=1,ns
                  Yt(1,ii) = func(xt(:,ii), D)
                  Yt(2:(D+1),ii) = grad(xt(:,ii),D)
               enddo
               Y_GRADIENT = transpose(yt)
            end function 
      end subroutine
end module

program test

   integer :: D=1
   double precision :: xminmax(1,2)
   xminmax(1,1) = 0
   xminmax(1,2) = 4

   call adaptive_ssck(func, grad, xminxmax, D)
   
   contains 
      function func(x,D)
         integer :: D
         double precision :: x(D,1)
         double precision :: func
         func = (x(1,1) - 1.0d0) ** 2
      end function

      function grad(x,D)
         integer :: D
         double precision :: x(D,1)
         double precision :: grad(D,1)
         grad = 2 * (x(1,1) - 1) 
      end function
end program

