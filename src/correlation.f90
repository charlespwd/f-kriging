module correlation
   implicit none

   contains
      !!
      ! This function computes the spatial correlation matrix R vectorially
      ! 
      ! inputs:
      !  R: the correlation matrix, NxN
      !	theta: the correlation factor;
      ! 	X: the vector of snapshot positions
      !  D: # of dimensions
      !  Ns: # of snapshots
      !  Pc: Power of correlation 
      !
      ! outputs: 
      ! 	R	 
      SUBROUTINE construct_R(R,theta,X,D,Ns,Pc)
         IMPLICIT NONE
         DOUBLE PRECISION :: get_rxy
         INTEGER :: D, Ns, Pc
         INTEGER :: ii, jj, dd
         DOUBLE PRECISION,intent(out) :: R(Ns,Ns)
         DOUBLE PRECISION,intent(in) :: X(Ns,D),theta(D)
         DOUBLE PRECISION :: tmp
         DO ii=1,Ns
            DO jj=ii,Ns
               ! inlined for performance, this gets called so much that the overhead
               ! isn't worth the better readability. I left the call there
               ! instead
      !         R(ii,jj) = get_rxy(theta,X(ii,:),X(jj,:),D,Pc)
               tmp = 0.0d0
               DO dd=1,D
                  tmp = tmp &
                     - theta(dd) * (abs(X(ii,dd) - X(jj,dd)) ** Pc)
               END DO
               R(ii,jj) = exp(tmp)
               R(jj,ii) = R(ii,jj)
            END DO
         END DO
      END SUBROUTINE

      SUBROUTINE construct_RT(R,theta,X,D,Ns,Pc)
         IMPLICIT NONE
         DOUBLE PRECISION :: get_rxy
         INTEGER :: D, Ns, Pc
         INTEGER :: ii, jj, dd
         DOUBLE PRECISION,intent(out) :: R(Ns,Ns)
         DOUBLE PRECISION,intent(in) :: X(D,Ns),theta(D)
         DOUBLE PRECISION :: tmp
         DO jj=1,Ns
            DO ii=jj,Ns
               ! inlined for performance, this gets called so much that the overhead
               ! isn't worth the better readability. I left the call there
               ! instead
      !         R(ii,jj) = get_rxy(theta,X(ii,:),X(jj,:),D,Pc)
               tmp = 0.0d0
               DO dd=1,D
                  tmp = tmp &
                     - theta(dd) * (abs(X(dd,ii) - X(dd,jj)) ** Pc)
               END DO
               R(ii,jj) = exp(tmp)
               R(jj,ii) = R(ii,jj)
            END DO
         END DO
      END SUBROUTINE   

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
            
      !!
      ! A routine to inverse the correlation matrix
      ! it is semi intelligent as it corrects it by calling nuggetcorrect
      ! in the event the condition number is too large. 
      ! 
      ! ARGUMENTS:
      !  R (inout) : the correlation matrix, gets modified if nuggetcorrected.
      !  I (in) : NsxNs, the identity matrix
      !  Rinv (out) : the inverse of the correlation matrix, ns x ns
      !  Ns : number of snapshots, i.e. dimension of matrix;
      SUBROUTINE invertR(R,I,Rinv,Ns)
         USE LA_PRECISION,ONLY:WP=>DP
         USE F95_LAPACK,ONLY:LA_GESVX
         IMPLICIT NONE 
         INTEGER, INTENT(IN) :: Ns
         DOUBLE PRECISION, INTENT(INOUT) :: I(Ns,Ns)
         DOUBLE PRECISION, INTENT(INOUT) :: R(Ns,Ns)
         DOUBLE PRECISION, INTENT(OUT) :: Rinv(Ns,Ns)
         DOUBLE PRECISION :: RCOND
         INTEGER :: INFO

         CALL LA_GESVX(R,I,Rinv,RCOND=RCOND,INFO=INFO)
         CALL NUGGETCORRECT(R,I,Rinv,Ns,RCOND ** (-1))
      END SUBROUTINE

      !! 
      ! nugget correct, this function is CALLed by inverse R
      ! When a matrix is close to singular, a small epsilon is added to the 
      ! diagonal of the correlation matrix to permit it to solve
      SUBROUTINE nuggetcorrect(R,I,Rinv,Ns,cond)
         USE PARAMS,ONLY:NUGGET,CONDTOL
         USE LA_PRECISION,ONLY:WP=>DP
         USE F95_LAPACK,ONLY:LA_GESVX
         IMPLICIT NONE
         INTEGER,INTENT(IN) :: Ns
         DOUBLE PRECISION,INTENT(inout) :: Rinv(Ns,Ns), R(Ns,Ns) 
         DOUBLE PRECISION,INTENT(inout) :: I(Ns,Ns) 
         INTEGER :: ii,info
         DOUBLE PRECISION :: eps, cond,rcon

      ! could possibly optimize this by giving the factored L matrix since
      ! you need to calculate LA_GESVX to get the condition number. 
      !   CALL LA_GESVX(R,I,Rinv,RCOND=rcon)
      !   cond = rcon ** (-1) 
         IF (cond > condtol) THEN 
      !      WRITE(*,'(a25,E11.5)') 'applying nugget to cond= ',cond
            eps = (10 + Ns) * 10D0 ** (-NUGGET) 
            DO ii=1,Ns
              R(ii,ii) = R(ii,ii) + eps
            END DO

            ! we do not want R or I to be modified by the LAPACK subroutines

            CALL LA_GESVX(R,I,Rinv,RCOND=rcon,INFO=info) 
      !      WRITE(*,'(a25,E11.5)') 'cond after nugget= ', rcon**(-1)
         END IF
      END SUBROUTINE
end module
