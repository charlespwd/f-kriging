module MLE
   implicit none

   contains
      include "./mle-src/optimize_theta_mle.f90" 
      include "./mle-src/construct_beta.f90" 
      include "./mle-src/construct_delta.f90" 
      include "./mle-src/construct_DR.f90" 
      include "./mle-src/init_theta.f90" 
end module

