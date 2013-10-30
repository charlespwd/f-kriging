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
         CASE("--branin", "-b")
            DO nn=1,Ns
               tmp = fbranin(X(nn,:),D)
               Y_GRADIENT(nn,1:(D+1)) = tmp(1,1:(D+1))
            END DO
         CASE("--cosine", "-c")
            DO nn=1,Ns
               tmp = fcosine(X(nn,:),D)
               Y_GRADIENT(nn,1:(D+1)) = tmp(1,1:(D+1))
            END DO
         case("--rosenbrock", "-r");
            DO nn=1,Ns
               tmp = frosenbrock(X(nn,:),D)
               Y_GRADIENT(nn,1:(D+1)) = tmp(1,1:(D+1))
            END DO
         CASE DEFAULT
            PRINT * , 'option ',adjustl(func_name), ' not supported'
            STOP
      END SELECT
   END FUNCTION

   FUNCTION fdrag(xv,D)
      USE PARAMS, ONLY : PI
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: D
      DOUBLE PRECISION, INTENT(IN) :: xv(1,D)
      DOUBLE PRECISION, DIMENSION(1,1+2) :: fdrag
      DOUBLE PRECISION :: x,y
      ! D should be two
      IF (D /= 2) THEN
         print *, 'D should be 2'
         STOP
      END IF

      x = xv(1,1)
      y = xv(1,2)

      ! y
      fdrag(1,1) =  0.001d0*(1+y**2 /(1.08d0-y)) *x;
      ! gradx
      fdrag(1,2) =  0.001d0 *(1+y**2 /(1.08d0-y));
      ! grady
      fdrag(1,3) =  0.001d0 * x  * (2 * y  / (1.08d0 - y) + y**2  / ((1.08d0-y)**2));
   END FUNCTION

   ! fbranin(x,D) 
   ! returns the function y and the gradient of y @ x
   !  arguments:
   !  x: loc vector, [1,D]
   !  D: # of dimensions of domain, 
   !  fbranin: [y,gradx,grady], [D+1,1]
   FUNCTION fbranin(xv,D)
      USE PARAMS, ONLY : PI
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: D
      DOUBLE PRECISION, INTENT(IN) :: xv(1,D)
      DOUBLE PRECISION, DIMENSION(1,1+2) :: fbranin
      DOUBLE PRECISION :: x,y
      ! D should be two
      IF (D /= 2) THEN
         print *, 'D should be 2'
         STOP
      END IF

      x = xv(1,1)
      y = xv(1,2)

      ! y
      fbranin(1,1) =  (y - 5.1d0/(4*PI**2) * x**2 + 5.d0/PI * x - 6)**2 + 10 * (1 - 1.0d0/(8*PI))*cos(x) + 10
      ! gradx
      fbranin(1,2) =  2 * (5.d0/PI - 0.258369d0*x) * (-6 + 5*x/PI - 0.129185d0*x**2 + y) - 10 * (1-1.0d0/(8*PI)) * sin(x)
      ! grady
      fbranin(1,3) = 2*(-6 + 5.d0*x/pi - 0.129185d0*x**2 + y)
   END FUNCTION   
   
   FUNCTION fcosine(xv,D)
      USE PARAMS, ONLY : PI
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: D
      DOUBLE PRECISION, INTENT(IN) :: xv(1,D)
      DOUBLE PRECISION, DIMENSION(1,1+2) :: fcosine
      DOUBLE PRECISION :: x,y
      ! D should be two
      IF (D /= 2) THEN
         print *, 'D should be 2'
         STOP
      END IF

      x = xv(1,1)
      y = xv(1,2)

      ! y
      fcosine(1,1) =  cos(10*x) + sin(10*y) + x*y
      ! gradx
      fcosine(1,2) =  -10*sin(10*x) + y
      ! grady
      fcosine(1,3) =  10*cos(10*y) + x;
   END FUNCTION

   FUNCTION frosenbrock(xv,D)
      USE PARAMS, ONLY : PI
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: D
      DOUBLE PRECISION, INTENT(IN) :: xv(1,D)
      DOUBLE PRECISION, DIMENSION(1,1+2) :: frosenbrock
      DOUBLE PRECISION :: x,y
      ! D should be two
      IF (D /= 2) THEN
         print *, 'D should be 2'
         STOP
      END IF

      x = xv(1,1)
      y = xv(1,2)

      ! y
      frosenbrock(1,1) = (1 - x) ** 2 + 100 * (y - x ** 2) ** 2
      ! gradx
      frosenbrock(1,2) = 2 * (200 * x ** 3 - 200 * x * y + x - 1) 
      ! grady
      frosenbrock(1,3) = 200 * (y - x ** 2) 
   END FUNCTION
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
