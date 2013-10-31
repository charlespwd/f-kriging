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
      subroutine construct_R(R,theta,X,D,Ns,Pc)
         implicit none
         double precision :: get_rxy
         integer :: D, Ns, Pc
         integer :: ii, jj, dd
         double precision,intent(out) :: R(Ns,Ns)
         double precision,intent(in) :: X(Ns,D),theta(D)
         double precision :: tmp
         do ii=1,Ns
            do jj=ii,Ns
               ! inlined for performance, this gets called so much that the overhead
               ! isn't worth the better readability. I left the call there
               ! instead
      !         R(ii,jj) = get_rxy(theta,X(ii,:),X(jj,:),D,Pc)
               tmp = 0.0d0
               do dd=1,D
                  tmp = tmp &
                     - theta(dd) * (abs(X(ii,dd) - X(jj,dd)) ** Pc)
               end do
               R(ii,jj) = exp(tmp)
               R(jj,ii) = R(ii,jj)
            end do
         end do
      end subroutine

      subroutine construct_RT(R,theta,X,D,Ns,Pc)
         implicit none
         double precision :: get_rxy
         integer :: D, Ns, Pc
         integer :: ii, jj, dd
         double precision,intent(out) :: R(Ns,Ns)
         double precision,intent(in) :: X(D,Ns),theta(D)
         double precision :: tmp
         do jj=1,Ns
            do ii=jj,Ns
               ! inlined for performance, this gets called so much that the overhead
               ! isn't worth the better readability. I left the call there
               ! instead
      !         R(ii,jj) = get_rxy(theta,X(ii,:),X(jj,:),D,Pc)
               tmp = 0.0d0
               do dd=1,D
                  tmp = tmp &
                     - theta(dd) * (abs(X(dd,ii) - X(dd,jj)) ** Pc)
               end do
               R(ii,jj) = exp(tmp)
               R(jj,ii) = R(ii,jj)
            end do
         end do
      end subroutine   

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
      ! THIS function HAS BEEN INLINED in CONSTRUCT_R.F90, YOU HAVE TO UPDATE
      ! IT THERE AS WELL if YOU MODIFY ANYTHING THAT IS in HERE. 
      double precision function get_rxy(theta,x,y,D,Pc)
         implicit none
         integer :: D, Pc
         double precision :: r
         double precision :: theta(D), x(D), y(D)
         integer :: dd
         double precision :: tmp
         tmp = 0
         do dd=1,D
            tmp = tmp & 
               - theta(dd) * (abs(x(dd) - y(dd)) ** Pc)
         end do
         get_rxy = exp(tmp)
         RETURN
      end function
            
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
      subroutine invertR(R,I,Rinv,Ns)
         use LA_PRECISION,only:WP=>DP
         use F95_LAPACK,only:LA_GESVX
         implicit none 
         integer, intent(in) :: Ns
         double precision, intent(INOUT) :: I(Ns,Ns)
         double precision, intent(INOUT) :: R(Ns,Ns)
         double precision, intent(out) :: Rinv(Ns,Ns)
         double precision :: RCOND
         integer :: INFO

         call LA_GESVX(R,I,Rinv,RCOND=RCOND,INFO=INFO)
         call NUGGETCORRECT(R,I,Rinv,Ns,RCOND ** (-1))
      end subroutine

      !! 
      ! nugget correct, this function is CALLed by inverse R
      ! When a matrix is close to singular, a small epsilon is added to the 
      ! diagonal of the correlation matrix to permit it to solve
      subroutine nuggetcorrect(R,I,Rinv,Ns,cond)
         use PARAMS,only:NUGGET,CONDTOL
         use LA_PRECISION,only:WP=>DP
         use F95_LAPACK,only:LA_GESVX
         implicit none
         integer,intent(in) :: Ns
         double precision,intent(inout) :: Rinv(Ns,Ns), R(Ns,Ns) 
         double precision,intent(inout) :: I(Ns,Ns) 
         integer :: ii,info
         double precision :: eps, cond,rcon

      ! could possibly optimize this by giving the factored L matrix since
      ! you need to calculate LA_GESVX to get the condition number. 
      !   call LA_GESVX(R,I,Rinv,RCOND=rcon)
      !   cond = rcon ** (-1) 
         if (cond > condtol) then 
      !      write(*,'(a25,E11.5)') 'applying nugget to cond= ',cond
            eps = (10 + Ns) * 10D0 ** (-NUGGET) 
            do ii=1,Ns
              R(ii,ii) = R(ii,ii) + eps
            end do

            ! we do not want R or I to be modified by the LAPACK subroutines

            call LA_GESVX(R,I,Rinv,RCOND=rcon,INFO=info) 
      !      write(*,'(a25,E11.5)') 'cond after nugget= ', rcon**(-1)
         end if
      end subroutine
end module
