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
