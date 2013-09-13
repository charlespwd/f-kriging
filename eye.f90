SUBROUTINE eye(I,N)
   INTEGER :: ii,N
   DOUBLE PRECISION :: I(N,N)
   I(:,:) = 0.0D1
   do ii=1,N
      I(ii,ii) = 1
   end do
END SUBROUTINE

