module global_Flags

  logical :: textOutput = .true.

end module global_Flags
  

module energiesData_Module

  type energiesData
     
     double precision :: Utotal
     integer, allocatable :: distancesIntMat(:,:), triMat(:,:), changedTriInd(:)
     double precision, allocatable :: interatomicDistances(:,:), tripletEnergies(:)
     double precision, allocatable :: expMatrix(:,:,:), changedTriU(:)

  end type energiesData

end module EnergiesData_Module

module positionData_Module

  type positionData

    integer :: N_a, N_tri, N_distances, N_changed_triplets
    double precision, allocatable :: posArray(:,:)

  end type positionData

end module positionData_Module
