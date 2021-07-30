module global_Flags
  logical :: textOutput = .true.

end module global_Flags
  

module energiesData_Module

  type energiesData
     
     double precision:: Utotal
     integer, allocatable :: distancesIntMat(:,:)
     double precision, allocatable :: interatomicDistances(:,:), tripletEnergies(:)
     double precision, allocatable :: expMatrix(:,:,:)

  end type energiesData

end module EnergiesData_Module

