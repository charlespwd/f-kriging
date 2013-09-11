!! 
! This function builds the regression matrix F.
! F is built as [f(X_1),f(X_2),...,f(X_NSNAP)]^T
! where, in this case, 
! f(X_1) = [1,x_1^1,x_2^1,...,x_D^1,x_1^2,...,x_D^2,...,x_1^O,...,x_D^]^T
! 
! example: (D = 2, O = 3) :
! f(X_1) = [1,  x_1^1, x_2^1,  x_1^2, x_2^2,  x_1^3, x_2^3]^T
!
! inputs:
!  snap_pos: the vector of snapshot positions
!  nsnap: number of snapshots
!  D: number of dimensions
!  Order: Order of the polynomial chosen.
!
! Output:
!  F: The regression matrix without the coefficients. 
SUBROUTINE construct_fmat(F,snap_pos,Order,D,nsnap)
   IMPLICIT NONE
   INTEGER :: Order, D, nsnap
   INTEGER :: nn
   REAL :: F(nsnap,1 + Order*D)
   REAL :: snap_pos(nsnap,D)
   DO nn=1,nsnap
      call construct_f(F(nn,:),snap_pos(nn,:),Order,D,nsnap)
   END DO
END SUBROUTINE
