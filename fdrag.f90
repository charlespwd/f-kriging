! fdrag(x,D) 
! returns the function y and the gradient of y @ x
!  arguments:
!  x: loc vector, [1,D]
!  D: # of dimensions of domain, 
!  fdrag: [y,gradx,grady], [D+1,1]
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
   fdrag(1,1) =  (y - 5.1d0/(4*PI**2) * x**2 + 5.d0/PI * x - 6)**2 + 10 * (1 - 1.0d0/(8*PI))*cos(x) + 10
   ! gradx
   fdrag(1,2) =  2 * (5.d0/PI - 0.258369d0*x) * (-6 + 5*x/PI - 0.129185d0*x**2 + y) - 10 * (1-1.0d0/(8*PI)) * sin(x)
   ! grady
   fdrag(1,3) = 2*(-6 + 5.d0*x/pi - 0.129185d0*x**2 + y)
END FUNCTION
