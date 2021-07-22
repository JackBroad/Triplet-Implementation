module GP_variables
  implicit none

  double precision :: hyperParams(3)
  double precision, allocatable ::  alpha(:)
  integer ::  Perm(6,3)
  double precision, allocatable :: trainData(:,:)
  integer:: N_tp,  nArgs, N_p 
  
end module GP_variables
