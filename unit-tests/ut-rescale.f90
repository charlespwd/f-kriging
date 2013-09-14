program k
   double precision :: x(3,2), xma(2), xmi(2)
   x(1,:) = (/0,-2/)
   x(2,:) = (/4,0/)
   x(3,:) = (/8,2/)
   call dumpmat(x,3,2)
   xma = maxval(x,1)
   xmi = minval(x,1)
   call rescale(x,2,3,xma,xmi)
   write(*,*) 'rescaled.xout is present'
   call dumpmat(x,3,2)
end program
