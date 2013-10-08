
program testsensitivity
   use sensitivity, only: construct_sensitivity
   implicit none
   integer :: D, Ns, NsNew, order=2
   double precision,allocatable :: S(:,:),X(:,:),Y(:,:),Grad(:,:),F(:,:),R(:,:)
   double precision,allocatable :: theta(:), XNEW(:,:) 
   integer :: ii,fdim

   ! D
   open(unit=15,file='./unit-tests/csvd.dat',status='old')
      read(15,*) D
      print*, 'd= ', D
   close(15)

   ! ns
   open(unit=15,file='./unit-tests/csvns.dat',status='old')
      read(15,*) ns
      print*, 'ns= ', ns
   close(15)

   ! nsnew
   open(unit=15,file='./unit-tests/csvnsnew.dat',status='old')
      read(15,*) nsnew
      print*, 'nsnew =', nsnew
   close(15)

   fdim = 1+Order*D
   print*, 'fdim  =', fdim 

   allocate(X(ns,D))
   allocate(XNEW(nsnew,D))
   allocate(Y(ns,1))
   allocate(Grad(ns,D))
   allocate(F(ns,fdim))
   allocate(R(ns,ns))
   allocate(theta(D))
   allocate(S(nsnew,1))

   ! x
   open(unit=15,file='./unit-tests/csvx.dat',status='old')
   do ii=1,D
      read(15,*) X(:,ii)
   enddo
   close(15)

   ! xnew
   open(unit=15,file='./unit-tests/csvxnew.dat',status='old')
   do ii=1,D
      read(15,*) XNEW(:,ii)
   enddo
   close(15)

   ! Y
   open(unit=15,file='./unit-tests/csvy.dat',status='old')
      read(15,*) Y(:,1)
   close(15)

    ! GRAD
   open(unit=15,file='./unit-tests/csvgrad.dat',status='old')
   do ii=1,D
      read(15,*) GRAD(:,ii)
   enddo
   close(15)

   ! F
   open(unit=15,file='./unit-tests/csvf.dat',status='old')
   do ii=1,fdim
      read(15,*) F(:,ii)
   enddo
   close(15)

   ! R
   open(unit=15,file='./unit-tests/csvr.dat',status='old')
   do ii=1,ns
      read(15,*) R(:,ii)
   enddo
   close(15)

   ! theat
   open(unit=15,file='./unit-tests/csvtheta.dat',status='old')
      read(15,*) theta
      print*, 'theta =', theta
   close(15)

   call construct_sensitivity(S,XNEW,Y,X,Grad,F,R,theta,D,Ns,NsNew)
   do ii=1,15
    print*, 'S(1)= ', S(ii,1) 
   enddo
end program

!program testsampling
!   use sensitivity, only : get_sampling_radius
!   double precision :: xmax(2), xmin(2)
!   xmax = (/4,3/)
!   xmin = (/0,0/)
!   print *, get_sampling_radius(xmin,xmax,2,4)
!end program

!program test
!   use sensitivity, only : construct_density_function
!   double precision :: xmax(2), xmin(2)
!   double precision :: xnew(5,2), xold(4,2)
!   double precision :: psi(5) 
!   integer :: ns, d, nsnew
!   ns = 4
!   nsnew = 5
!   d = 2
!   xmax = (/1,1/)
!   xmin = (/0,0/)
!   xnew(:,1) = (/0,0,0,0,0/)
!   xnew(:,2) = (/0.0d0,0.2d0,0.4d0,0.6d0,0.8d0/)
!   xold(:,1) = (/0,0,1,1/)
!   xold(:,2) = (/0,1,0,1/)
!   call construct_density_function(psi, xnew, xold,xmax,xmin,d,nsnew,ns)
!   print * , psi
!   
!end program
