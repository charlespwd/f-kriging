module matrix
   USE F95_LAPACK,ONLY:DGESVD_F95
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
         
         ! LA_GESVD destroys the content of A. 
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
      SUBROUTINE rescale(X, D, Ns, XMIN, XMAX)
         INTEGER, INTENT(IN) :: Ns, D
         DOUBLE PRECISION, INTENT(IN) :: XMAX(D), XMIN(D)
         DOUBLE PRECISION, INTENT(INOUT) :: X(Ns,D)
         INTEGER :: dd

         DO dd=1,D
            X(:,dd) = X(:,dd) / (XMAX(dd) - XMIN(dd)) 
         END DO
      END SUBROUTINE
end module



