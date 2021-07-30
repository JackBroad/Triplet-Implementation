module energiesData_Module

  type energiesData
     
     double precision:: Utotal
     integer, allocatable :: distancesIntMat(:,:)
     double precision, allocatable :: interatomicDistances(:,:), tripletEnergies(:)
     double precision, allocatable :: expMatrix(:,:,:)

  end type EnergiesData

end module EnergiesData_Module

