module threeBody_mod
  use GP_variables
  use mpi_variables
  use dataStructure_variables
  use expShare_variables
  use pbcAndMic_variables
  use energiesData_Module, only: energiesData
  use positionData_Module, only: positionData
  implicit none


  public


contains


function threeBodyLRC_getTypeTwoTable(nU,nTheta,nOneTwo,expo) result(table)
  implicit none
  integer :: nU, nTheta, nGrid, nOneTwo
  integer :: counter, i, j, k
  real(dp) :: pi = 3.1415926535897932
  double precision :: rOneTwo, thetaMin, thetaMax
  double precision :: uMin, uMax, sinTheta, cosTheta
  double precision :: theta, u, delTheta, delU, delR
  double precision :: xPrime, Rc, gp, energy, expo
  double precision, allocatable :: table(:,:)

  ! Allocate memory to hold table of evaluations of Type 2 integrand
  nGrid = nOneTwo*nU*nTheta
  allocate(table(nGrid,4))

  ! Set limits of integration for theta (set uMax inside loop over
  ! theta)
  thetaMin = 0d0
  thetaMax = pi/2d0
  uMin = 0.0000000000001

  ! Set non-inv cut-off and counter
  Rc = 1d0/Rcut
  counter = 0

  ! Find amount by which to increases rOneTwo, and theta
  delTheta = thetaMax/(nTheta-1d0)
  delR = (Rc-1.5)/(nOneTwo-1d0) ! Rc-1.5 as 1.5 is minimum possible triplet dist

  do i = 1, nOneTwo
    rOneTwo = 1.5000001 + (i-1)*delR
    do j = 1, nTheta
      ! Set theta value
      theta = (j-1)*delTheta
      cosTheta = cos(theta)
      if (cosTheta .lt. 0d0) then
        cosTheta = 0d0
      end if
      sinTheta = sin(theta)
      if (sinTheta .lt. 0d0) then
        sinTheta = 0d0
      end if

      ! Set upper limit for u at this theta
      xPrime = 4d0*(Rc**2)
      xPrime = xPrime - (rOneTwo**2)*(sinTheta**2)
      xPrime = xPrime**0.5
      xPrime = xPrime - rOneTwo*cosTheta
      xPrime = xPrime/2d0
      uMax = xPrime**(-1d0*expo)

      delU = uMax/(nU-1d0)
      do k = 1, nU
        ! Set u value and update counter
        u = (k-1)*delU
        if (u .eq. 0d0) then
          u = u + uMin
        end if
        counter = counter + 1

        ! Get GP energy for current theta if all r_ij<r_cut
        gp = 1d0 ! PLACEHOLDER

        ! Get energy for this u, theta and rOneTwo
        energy = threeBodyLRC_evaluateTypeTwoIntegrand(gp,u,sinTheta,cosTheta,&
                                                       rOneTwo,xPrime,expo)

        ! Save theta, inv distance and energy in table
        table(counter,1) = rOneTwo
        table(counter,2) = theta
        table(counter,3) = u**(-1d0/expo)
        table(counter,4) = energy
      end do
    end do
  end do

return
end function threeBodyLRC_getTypeTwoTable


double precision function threeBodyLRC_evaluateTypeTwoIntegrand(gpEn,u,sinT,cosT,&
                                                                r,xPrime,expo)
  implicit none
  real(dp) :: pi = 3.1415926535897932
  double precision :: gpEn, u, sinT, cosT, r, xPrime
  double precision :: numerator, denominator, expo

  numerator = 0.5 - (cosT**2)
  numerator = numerator * (r**2)*(xPrime**2)
  numerator = numerator + (1d0/16d0)*(r**4)
  numerator = numerator + xPrime**4
  numerator = numerator**(3d0/2d0)

  denominator = 0.5 - (cosT**2)
  denominator = denominator * (r**2)*(u**(-2d0/expo))
  denominator = denominator + (1d0/16d0)*(r**4)
  denominator = denominator + u**(-4d0/expo)
  denominator = denominator**(3d0/2d0)
  denominator = denominator * expo

  threeBodyLRC_evaluateTypeTwoIntegrand = numerator/denominator
  threeBodyLRC_evaluateTypeTwoIntegrand = threeBodyLRC_evaluateTypeTwoIntegrand &
                                          * gpEn
  threeBodyLRC_evaluateTypeTwoIntegrand = threeBodyLRC_evaluateTypeTwoIntegrand &
                                          * (sinT**2)
  threeBodyLRC_evaluateTypeTwoIntegrand = threeBodyLRC_evaluateTypeTwoIntegrand &
                                          * (u**(-1-(3d0/expo)))
  threeBodyLRC_evaluateTypeTwoIntegrand = threeBodyLRC_evaluateTypeTwoIntegrand &
                                          * 4d0*pi

return
end function threeBodyLRC_evaluateTypeTwoIntegrand


end module threeBody_mod
