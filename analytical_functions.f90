! THis module permits to call functions with a string command
! the beauty of this is that you don't need twenty interfaces 
! to achieve something beautiful
MODULE ANALYTICAL_FUNCTIONS
   IMPLICIT NONE
   DOUBLE PRECISION, PARAMETER :: PI = 4.0D0*DATAN(1.0D0)

   CONTAINS

   FUNCTION Y_GRADIENT(X,D,Ns,func_name)
      ! arguments
      INTEGER, INTENT(IN) :: D,Ns
      DOUBLE PRECISION, INTENT(IN) :: X(Ns,D)
      CHARACTER(len=20), INTENT(IN) :: func_name
      DOUBLE PRECISION, DIMENSION(Ns,1+D) :: Y_GRADIENT
      
      INTEGER :: nn
      DOUBLE PRECISION :: tmp(1,D)

      SELECT CASE (adjustl(func_name))
         CASE("--drag","-d")
            DO nn=1,Ns
               tmp = fdrag(X(nn,:),D)
               Y_GRADIENT(nn,1:(D+1)) = tmp(1,1:(D+1))
            END DO
         CASE DEFAULT
            PRINT * , 'option ',adjustl(func_name), ' not supported'
            STOP
      END SELECT
   END FUNCTION

   INCLUDE 'fdrag.f90'
END MODULE

!PROGRAM p
!   USE ANALYTICAL_FUNCTIONS
!   double precision :: X(3,2)
!   double precision :: Y(3,3) 
!   integer :: i
!   character(len=20) :: func_name
!   func_name = '-d'
!   X(:,1) = (/ 1,2,3 /)
!   X(:,2) = (/ 1,2,3 /)
!   Y = YGRAD(X,2,3,func_name)
!   do i =1,3
!      print * , Y(i,:)
!   end do
!end program
