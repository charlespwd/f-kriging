!!
! This function computes the spatial correlation between x and y
!
! inputs:
!	theta: the correlation factor;
! 	x: the first vector; 1xD
! 	y: the second vector; 1xD
!
! outputs: 
! 	r(x,y);	 
REAL FUNCTION rxy(theta,x,y,D,Pc)
   IMPLICIT NONE
   INTEGER :: D, Pc
   REAL :: r
   REAL :: theta(D), x(D), y(D)
   INTEGER :: dd
   REAL :: tmp
   tmp = 0
   DO dd=1,D
      tmp = tmp & 
         - theta(dd) * (abs(x(dd) - y(dd)) ** Pc)
   END DO
   rxy = exp(tmp)
   RETURN
END FUNCTION
      
   
