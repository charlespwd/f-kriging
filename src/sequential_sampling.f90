! this module contains everything related to sequential sampling
!  get_sampling_radius(xmin,xmax,D,ns)
!  construct_density_function(Psi,Xnew,X,XMAX,XMIN,D,nsnew,ns)
!  samplingcriterion(Eest,MSE,nsnew,PSI,S,CV)
!  construct_insertion_stack(eest,ns) 
module sequential_sampling
   use fort_arrange, only: sort
   use matrix, only: normalize, distance
   implicit none

   ! module globals ( editable parameters ) & defaults
   integer, parameter :: MODE_MSE = 0
   integer, parameter :: MODE_SENSITIVITY = 1
   integer :: nGridStart = 3     ! for a nGridStart**D, evenly spaced, initial
                                 !  sample set
   integer :: delta_ns = 4       ! number of points added to sample @ new iter
   integer :: nFinal = 40        ! maximum number of points to kriging model
   integer :: nGridRS = 50       ! for a nGridRs**D response surface grid
   integer :: mode = MODE_MSE    ! the adaptive scheme's mode (MSE || MSE+S)
   integer :: order = 2          ! regression order 

   contains
     ! adaptive_ssck - adaptive sequential sampling cokriging routine 
      subroutine adaptive_ssck(func, fgrad, xminmax, D, iterate, datadirectory, &
            optimize)
         use grid, only:vector_grid, columnGrid, LHS, grow2d
         use utils, only: printer, l1error
         use matrix, only: vector_range
         use optimization, only: carpet_bombed_min

         ! arguments
         integer, intent(in) :: D
         integer, intent(in), optional :: iterate
         integer, intent(in), optional :: optimize
         double precision, intent(in) :: xminmax(D,2) 
         character(len=50), intent(in), optional :: datadirectory
         interface true_function
            function func(x,D) 
               integer, intent(in) :: D
               double precision, intent(in) :: x(D,1)
               double precision :: func
            end function
         end interface true_function
         interface true_gradient
            function fgrad(x,D)
               integer, intent(in) :: D
               double precision, intent(in) :: x(D,1)
               double precision :: fgrad(D,1)
            end function
         end interface true_gradient

         ! work variables
         integer :: fileid_err, fileid_optim
         integer :: Ns, NsNew, loopcount
         integer,allocatable :: newsamples(:)
         double precision,allocatable :: xnew(:,:),ynew(:,:)
         double precision,allocatable :: x(:,:),y(:,:)
         double precision,allocatable :: x_star(:,:)
         double precision,allocatable :: grad(:,:)
         double precision,allocatable :: xmin(:), xmax(:)
         double precision,allocatable :: ygrad(:,:)
         double precision,allocatable :: ytrue(:,:)
         double precision,allocatable :: theta(:)
         double precision,allocatable :: psi(:,:)
         double precision,allocatable :: MSE(:,:)
         double precision,allocatable :: S(:,:)
         double precision,allocatable :: Eest(:,:)
         double precision :: y_star(1,1)
         double precision :: MeanL1
         double precision :: samplingRadius 
         double precision :: maxerror
         character(len=20) :: rsfile='d_rs.dat', truefile='d_true.dat', dotsfile='d_dots.dat'
         character(len=20) :: minfile='d_min.dat', optimfile='d_op.dat' 
         character(len=20) :: func_name, errfile, prefix
         character(len=50) :: datadir
         integer :: ii

         datadir = './data' ! by default
         if (present(datadirectory)) then
            datadir = datadirectory
         end if
         errfile= trim(datadir)//'/'//'e.dat'
         optimfile= trim(datadir)//'/'//trim(optimfile)

         ns = nGridStart ** D
         nsnew = nGridRs ** D

         allocate(xmin(d))
         allocate(xmax(d))
         allocate(xnew(nsnew,D))
         allocate(ynew(nsnew,1))
         allocate(ygrad(nsnew,D+1))
         allocate(ytrue(nsnew,1))
         allocate(x(ns,D))
         allocate(y(ns,1))
         allocate(x_star(1,D))
         allocate(grad(ns,D))
         allocate(theta(d))
         allocate(mse(nsnew,1))
         allocate(eest(nsnew,1))
         allocate(psi(nsnew,1))
         allocate(S(nsnew,1))

         allocate(newsamples(delta_ns))

         xmin(1:D) = xminmax(1:D,1)
         xmax(1:D) = xminmax(1:D,2)

         ! make grids
         call columnGrid(x,xmin,xmax,d,nGridStart)
         call columnGrid(xnew,xmin,xmax,d,nGridRs)

         ! true solution (ouch)
         ! TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO 
         ! TODO make sure to turn that off if the function isnt analytical
         ! TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO 
         ygrad = Y_GRADIENT(xnew,D,nsnew)
         ytrue(1:NsNew,1) = ygrad(1:nsnew,1)

         ! initialize theta with defaults 
         theta = -1 

         ! perform cokriging, sensitivity analysis, etc. 
         loopcount = 0
         call c_solver(xnew,ynew,theta,mse,xmin,xmax,x,y,Grad,Order,D,Ns,NsNew,Y_GRADIENT,&
            S=S) 

          ! make fancy graphs
         call printer(x,y,ns,1,D,dotsfile,datadir,loopcount)
         call printer(xnew,ynew,nsnew,nGridRs,D,rsfile,datadir,loopcount)
         call printer(xnew,ytrue,nsnew,nGridRs,D,truefile,datadir,loopcount)

         if (present(optimize)) then
            x_star = carpet_bombed_min(ynew, xnew, D , nsnew , y_star)
            print*, x_star
            call printer(x_star, y_star, 1, 1, D, minfile, datadir, loopcount)
            fileid_optim = 68
            open(fileid_optim, file=adjustl(optimfile), status='replace')
            write(fileid_optim, *) x_star, y_star
         end if

         fileid_err = 67
         open(fileid_err,file=adjustl(errfile),status='replace')
         ! print error to file
         maxerror = l1error(ytrue,ynew,nsnew) 
         write(fileid_err, *) ns, l1error(ytrue,ynew,nsnew), l1error(ytrue,ynew,nsnew) / maxerror

         loopcount = 0
         ! sequential sampling loop
         do while (ns<=nFinal)
            print*, 'ns = ', ns 
                 
            ! update the sampling radius because ns changed
            samplingradius = get_sampling_radius(xmin,xmax,D,ns)

            ! construct the new eest stack
            call construct_density_function(psi,xnew,x,xmax,xmin,D,nsnew,ns)
            call samplingcriterion(eest,mse,nsnew,psi,mode,S=S)

            ! find the points to add to the set
            call sampler(newsamples,Eest,XNEW,X,samplingRadius,delta_ns,D,nsnew,ns)

            ! add the points to the existing set
            call grow2d(x,delta_ns)
            call grow2d(y,delta_ns)
            call grow2d(grad,delta_ns)
            do ii=1,delta_ns
               x(ns+ii,1:D) = xnew(newsamples(ii),1:D)    
            enddo
            
            ! increment number of snapshots
            ns = ns + delta_ns
            
            ! perform cokriging, sensitivity analysis, etc. 
            call c_solver(xnew,ynew,theta,mse,xmin,xmax,x,y,grad,Order,D,Ns,NsNew,Y_GRADIENT, &
               S=S,DELTANS=delta_ns)
            
            loopcount = loopcount+1

            ! make fancy graphs
            call printer(x,y,ns,1,D,dotsfile,datadir,loopcount)
            call printer(xnew,ynew,nsnew,nGridRs,D,rsfile,datadir,loopcount)
            write(fileid_err, *) ns, l1error(ytrue,ynew,nsnew), l1error(ytrue,ynew,nsnew) / maxerror
            if (present(optimize)) then
               x_star = carpet_bombed_min(ynew, xnew, D, Nsnew, y_star)
               call printer(x_star, y_star, 1, 1, D, minfile, datadir, loopcount)
               write(fileid_optim, *) x_star, y_star
            end if
         enddo

         close(fileid_err)
         if(present(optimize)) then 
            close(fileid_optim)
         end if
         open(unit=67,file='loopcount.dat',status='replace')
         write(67,*) loopcount
         close(67)

         contains 
            function Y_GRADIENT(x, D, ns)
               ! arguments
               integer, intent(in) :: D, ns
               double precision, intent(in) :: x(Ns,D)
               double precision :: Y_GRADIENT(NS,1+D)
               
               ! work
               integer :: ii
               double precision :: xt(D,ns)
               double precision :: yt(1+D, ns)
               xt = transpose(x)
               do ii=1,ns
                  Yt(1,ii) = func(xt(:,ii), D)
                  Yt(2:(D+1),ii:ii) = fgrad(xt(:,ii),D)
               enddo
               Y_GRADIENT = transpose(yt)
            end function 
      end subroutine

      ! cokriging solver 
      subroutine c_solver(XNEW,YNEW,theta,MSE,XMIN,XMAX,X,Y,GRAD,Order,D,Ns,NsNew,&
            Y_GRADIENT, S, deltans)
         use PARAMS, only: Raug
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
         interface true_function
            function Y_GRADIENT(x,D,ns)
               integer, intent(in) :: D, ns
               double precision, intent(in) :: x(ns,D)
               double precision :: Y_GRADIENT(ns,1+D)
            end function
         end interface true_function

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
            tmpgrad = Y_GRADIENT(X(ns-deltans+1:ns,:),D,deltans)
            ! put it in ygrad
            YGRAD(ns-deltans+1:ns,1:D+1) = tmpgrad(1:deltans,1:D+1)
         else      
            YGRAD = Y_GRADIENT(X,D,Ns)
         endif

         Y(1:NS,1) = YGRAD(1:NS,1)
         GRAD(1:NS,1:D) = YGRAD(1:NS,2:(D+1))

         if (present(S)) then
            call COKRIGING(XNEW,YNEW,theta,MSE,X,Y,GRAD,Raug,Order,D,Ns,NsNew,S)
         else 
            call COKRIGING(XNEW,YNEW,theta,MSE,X,Y,GRAD,Raug,Order,D,Ns,NsNew)
         endif
      end subroutine

      ! sampler - finds the indexes of the points with the largest error
      !   and that is not within radius of the old points or the newly added
      !   points
      subroutine sampler(newsamples,Eest,XNEW,XOLD,Radius,deltaNs,D,nsnew,ns)
         integer, intent(in) :: deltaNs,D,nsnew,ns
         integer, intent(out) :: newsamples(deltaNs)
         double precision, intent(in) :: Eest(nsnew,1),Xnew(nsnew,D),Xold(ns,D)
         double precision, intent(in) :: radius
         !work variables
         integer :: ii, jj, inserted_count, toadd, ierr
         double precision :: errorstack(nsnew,2)

         ! sort errors in ascending order, and have the index map in the 
         ! second column
         
         errorstack = construct_insertion_stack(eest,nsnew)

         ! check the errorstack from then end to the beginning
         inserted_count = 0
         do ii=nsnew,1,-1
            toadd = 1
            ierr = errorstack(ii,2)
            ! check if new point is close to old snapshots
            do jj=1,ns
               if (distance(XNEW(ierr,1:D),XOLD(jj,1:D),D) <= Radius) then
                  toadd = 0
               endif
            enddo
            
            ! check if new point is close to newly added point
            if (inserted_count > 0) then
               do jj=1,inserted_count
                  if (distance(XNEW(ierr,1:D),XNEW(newsamples(jj),1:D),D) <= Radius) then
                     toadd = 0
                  endif
               enddo
            endif

            ! check if still ok to add the point
            if (toadd == 1) then
               inserted_count = inserted_count + 1
               newsamples(inserted_count) = ierr
               if (inserted_count == deltaNs) then
                  ! we added enough points
                  exit
               endif
            endif
          
         enddo      
      end subroutine

      ! get the sampling radius as per defined in Arthur's paper
      double precision function get_sampling_radius(xmin,xmax,D,ns)
         ! arguments
         integer, intent(in) :: D,ns
         double precision, intent(in) :: xmax(D),xmin(D)
         
         ! work
         double precision :: radius=0
         integer :: ii

         do ii=1,D
           Radius = Radius + (xmax(ii) - xmin(ii)) ** 2
         enddo
         Radius = 0.25d0 * (Radius ** 0.5d0) / (ns ** (1.0d0/D))
         get_sampling_radius = Radius
      end function

      ! constructs the density function 
      !  defined as 1 - exp(mindist to a snapshot) 
      !  so that this density function is 0 at snapshot locations
      !  and more than that everywhere else
      !  the result is 'normalized'
      subroutine construct_density_function(Psi,Xnew,X,XMAX,XMIN,D,nsnew,ns)
         ! arguments
         integer, intent(in) :: nsnew, ns, d
         double precision, intent(in) :: Xnew(nsnew,D), X(ns,D)
         double precision, intent(in) :: XMAX(D), XMIN(D)
         double precision, intent(out) :: Psi(nsnew,1)
         
         ! work variables
         double precision :: mindist,dist
         integer :: ii,jj,pp
         
         do ii=1,nsnew
            ! initialize at maximum
            mindist = maxval(xmax)
            do jj=1,ns
               dist=0.0d0
               do pp=1,D
                  dist = dist + (xnew(ii,pp) - x(jj,pp)) ** 2
               enddo
               dist = dist ** (0.5d0)
               if ( dist < mindist ) then
                  mindist = dist
               endif
            enddo
            Psi(ii,1) = 1.0d0 - exp(-mindist)
         enddo
         call normalize(psi(:,1),nsnew)
      end subroutine

      ! calculates the sampling criterion 
      ! S and CV are optional parameters
      ! S takes priority over CV. So if you want CV, make sure to 
      ! call samplingcriterion(Eest,MSE,nsnew,Psi,CV=CV)
      ! psi is the density function as called here.
      subroutine samplingcriterion(Eest,MSE,nsnew,PSI,mode,S)
         integer,intent(in) :: nsnew
         integer, intent(in) :: mode
         double precision, intent(in) :: MSE(nsnew,1)
         double precision, intent(in) :: psi(nsnew,1)
         double precision, intent(in), optional, target :: S(nsnew,1)
         double precision, intent(out) :: Eest(nsnew,1)
         double precision :: wmse(nsnew,1), ws(nsnew,1)
         ! work
         integer :: ii

         wmse = mse
         call normalize(wmse(:,1),nsnew)

         if(present(s)) then
            ws = s
            call normalize(ws(:,1),nsnew)
         endif
         
         do ii=1,nsnew         
         select case (mode) 
               case(MODE_MSE) 
                  Eest(ii,1) = WMSE(ii,1) 
               case(MODE_SENSITIVITY)
                  if ( .not. present(S)) then
                     print*, 'in sampling criterion, S not present'
                     stop
                  endif
                  Eest(ii,1) = (wMSE(ii,1)+wS(ii,1)) * psi(ii,1)
               case default
                  print*, 'mode not recognized'
                  stop
            end select
         enddo
         call normalize(eest,nsnew)
      end subroutine

      ! construct_insertion_stack
      !  this function constructs an array where the first column
      !  represents the error at the second column's index in the 
      !  xnew grid.
      !  this array is sorted with the minimum error at the first index 
      function construct_insertion_stack(eest,ns) 
         integer, intent(in) :: ns
         double precision, intent(in) :: eest(ns,1)
         double precision,dimension(ns,2) :: construct_insertion_stack
         integer :: i
         construct_insertion_stack(1:ns,1) = eest(1:ns,1)
         construct_insertion_stack(1:ns,2) = (/(i,i=1,ns)/)
         call sort(construct_insertion_stack)
      end function
end module

