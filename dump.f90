SUBROUTINE DUMPMAT(A,dim1,dim2)
   INTEGER :: dim1, dim2 
   DOUBLE PRECISION :: A(dim1,dim2)
   INTEGER :: ii,jj
   900 format (I3,a1,i3,a2,F33.16)
   DO ii=1,dim1
   DO jj=1,dim2
      write(*,900) ii,",",jj,": ",A(ii,jj) 
   END DO
   END DO
END SUBROUTINE

SUBROUTINE DUMPVEC(A,DIMS)
   INTEGER :: DIMS
   DOUBLE PRECISION :: A(DIMS)
   INTEGER :: ii
   900 format (I3,a2,F8.5)
   DO ii=1,DIMS
      write(*,900) ii,": ",A(ii)
   end do
END SUBROUTINE

SUBROUTINE DUMPTENsOR(A,dim1,dim2,dim3)
   INTEGER :: dim1,dim2,dim3
   DOUBLE PRECISION :: A(dim1,dim2,dim3)
   INTEGER :: ii,jj,kk
   900 format (I3,a1,i3,a1,i3,a2,F8.5)
   DO ii=1,dim1
   DO jj=1,dim2
   DO kk=1,dim3
      write(*,900) ii,",",jj,",",kk,": ",A(ii,jj,kk) 
   END DO
   END DO
   END DO
END SUBROUTINE

!PROGRAM test 
!   double precision :: a(2,2,2)
!   a(1,1,1) = 1
!   a(1,1,2) = 2
!   a(1,2,1) = 3
!   a(1,2,2) = 4
!   a(2,1,1) = 5
!   a(2,1,2) = 6
!   a(2,2,1) = 7
!   a(2,2,2) = 8
!   call dumptensor(a,(/2,2,2/))
!   call dumpmat(a(1,:,:),(/2,2/))
!   call dumpvec(a(:,:,1),2)
!END PROGRAM
