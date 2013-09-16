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
!
! THIS FUNCTION HAS BEEN INLINED IN CONSTRUCT_R.F90, YOU HAVE TO UPDATE
! IT THERE AS WELL IF YOU MODIFY ANYTHING THAT IS IN HERE. 
DOUBLE PRECISION FUNCTION get_rxy(theta,x,y,D,Pc)
   IMPLICIT NONE
   INTEGER :: D, Pc
   DOUBLE PRECISION :: r
   DOUBLE PRECISION :: theta(D), x(D), y(D)
   INTEGER :: dd
   DOUBLE PRECISION :: tmp
   tmp = 0
   DO dd=1,D
      tmp = tmp & 
         - theta(dd) * (abs(x(dd) - y(dd)) ** Pc)
   END DO
   get_rxy = exp(tmp)
   RETURN
END FUNCTION
      
   
