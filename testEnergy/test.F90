program test
  use energy_test
  implicit none

  double precision :: xStar(3), U

  xStar = (/ 1d0/3d0,1d0/3d0,1d0/3d0 /)

  call load_GP_Data()
  U = PES_GP(xStar)
  print *, U

end program test
