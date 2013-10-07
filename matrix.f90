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
end module



