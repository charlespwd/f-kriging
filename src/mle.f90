module MLE
   implicit none

   contains
      include "optimize_theta_mle.f90" 
      include "construct_beta.f90" 
      include "construct_delta.f90" 
      include "construct_DR.f90" 
      include "init_theta.f90" 
end module

