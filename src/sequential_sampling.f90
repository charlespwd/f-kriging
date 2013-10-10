! this module contains everything related to sequential sampling
!  get_sampling_radius(xmin,xmax,D,ns)
!  construct_density_function(Psi,Xnew,X,XMAX,XMIN,D,nsnew,ns)
!  samplingcriterion(Eest,MSE,nsnew,PSI,S,CV)
!  construct_insertion_stack(eest,ns) 
module sequential_sampling
   use fort_arrange, only: sort
   use matrix, only: normalize, distance

   implicit none

   contains

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
      subroutine samplingcriterion(Eest,MSE,nsnew,PSI,S,CV)
         integer,intent(in) :: nsnew
         double precision, intent(in) :: MSE(nsnew,1)
         double precision, intent(in) :: psi(nsnew,1)
         double precision, intent(in), optional, target :: S(nsnew,1)
         double precision, intent(in), optional, target :: CV(nsnew,1)
         double precision, intent(out) :: Eest(nsnew,1)
         ! work
         integer :: ii
         
         do ii=1,nsnew         
            if ( present(S) ) then
               Eest(ii,1) = (MSE(ii,1)+S(ii,1)) * psi(ii,1)
            elseif ( present(CV) ) then
               Eest(ii,1) = (MSE(ii,1)+CV(ii,1)) * psi(ii,1)
            else 
               Eest(ii,1) = MSE(ii,1) 
            endif
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
         double precision :: construct_insertion_stack(ns,2)
         integer :: i
         construct_insertion_stack(1:ns,1) = eest(1:ns,1)
         construct_insertion_stack(1:ns,2) = (/(i,i=1,ns)/)
         call sort(construct_insertion_stack,1)
      end function
end module

