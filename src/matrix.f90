! this module contains
!  spectral_norm(A,m,n)
!  rescale(X,D,Ns,XMIN,XMAX)
!  normalize(X,m) - x is X(m)
!  eye(I,m) 
module matrix
   use F95_LAPACK,only:DGESVD_F95
   implicit none

   contains
      function spectral_norm(A,m,n)
         integer, intent(in) :: m,n
         double precision, intent(in) :: A(m,n)
         double precision :: spectral_norm
         double precision :: wA(m,n)
         double precision,allocatable :: SingularValues(:)
         integer :: INFO
         wA = A
         
         allocate(SingularValues(min(m,n)))
         
         ! DGESVD destroys the content of A. 
         call DGESVD_F95(wA,SingularValues)
         spectral_norm = SingularValues(1)
      end function

      !!
      ! This routine rescales the vector x in a range 
      ! of distance 1. As suggested in kriging
      !
      ! ARGUMENTS:
      !  X (inout) : the vector to be scaled
      !  D (in) : # of dimensions
      !  Ns (in) : # of snapshots
      !  Xmax (in) : array of maximums per column of X (i.e. per dimension)
      !  Xmin (in) : array of min per column of X (i.e. per dimension)
      !  Xout (out), optional : if present, x isn't overwritten, the scaled
      subroutine rescale(X, D, Ns, XMIN, XMAX)
         integer, intent(in) :: Ns, D
         double precision, intent(in) :: XMAX(D), XMIN(D)
         double precision, intent(INOUT) :: X(Ns,D)
         integer :: dd, ii

         double precision :: mymax

         do dd=1,D
            X(:,dd) = X(:,dd) / (xmax(dd) - xmin(dd))  
         end do
      end subroutine

      ! does the opposite of rescale
      subroutine unscale(X, D, Ns, XMIN, XMAX)
         integer, intent(in) :: Ns, D
         double precision, intent(in) :: XMAX(D), XMIN(D)
         double precision, intent(INOUT) :: X(Ns,D)
         integer :: dd, ii

         double precision :: mymax

         do dd=1,D
            X(:,dd) = X(:,dd) * (xmax(dd) - xmin(dd))  
         end do
      end subroutine

      ! normalizes positive vector between 0 and 1
      subroutine normalize(X,m)
         integer, intent(in) :: m
         double precision, intent(inout) :: X(m)
         double precision :: xmax,xmin
         integer :: ii, signs
         X = ABS(X)
         xmax = maxval(X)
         X = X / xmax
      end subroutine

      double precision function vector_range(v, D)
         integer,  intent(in) :: D
         double precision, intent(in) :: v(D)
         double precision :: vmax, vmin
         vmax = maxval(v)
         vmin = minval(v)
         vector_range = vmax - vmin
      end function

      !! make an identiy matrix, 
      subroutine eye(I,N)
         integer :: ii,N
         double precision :: I(N,N)
         I(:,:) = 0.0D0
         do ii=1,N
            I(ii,ii) = 1
         end do
      end subroutine

      ! computes the trace of a square matrix
      ! Arguments
      !  A: square matrix, mxm
      !  m: size of matrix
      double precision function get_trace(A,m)
         implicit none
         integer, intent(in) :: m
         double precision, intent(in) :: A(m,m)
         integer :: ii
         get_trace = 0
         do ii=1,m
            get_trace = get_trace + A(ii,ii)
         end do
      end function

      ! calculates the distance betwen two vectors
      double precision function distance(X,Y,D)
         integer, intent(in) :: D
         double precision, intent(in) :: X(1,D), Y(1,D)
         ! work
         integer :: ii
         double precision :: tmp
         tmp = 0.0d0
         do ii = 1,D
           tmp = tmp + (X(1,ii) - Y(1,ii)) ** 2.0d0
         enddo
         distance = tmp ** (0.5d0) 
      end function
end module
