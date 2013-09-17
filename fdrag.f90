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
