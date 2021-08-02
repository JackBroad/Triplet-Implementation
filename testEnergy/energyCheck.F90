module energyCheck

contains

subroutine energyCheckCalc(xStar, PES_GP)
  implicit none
  double precision :: xStar(3), alpha(3)
  double precision :: PES_GP, xTraining(3,3), xTrainingPerm(3,6,3)
  double precision :: lScale=1.010933217823705432e-01
  double precision :: expVar=5.065189615948498530e-05
  integer :: i, j, k, nTraining, nPerms, nDim, perm(3,6), l, m, n
  double precision :: kSqExpAllPerms, kSqExpJthPerm, kKernTotal

  !xStar = (/ 4.9530743964810102E-002, 2.1838616589923775E-002, 3.0229538987517568E-002/)
  perm(:,1) = (/1, 2, 3/)
  perm(:,2) = (/1, 3, 2/)
  perm(:,3) = (/2, 1, 3/)
  perm(:,4) = (/2, 3, 1/)
  perm(:,5) = (/3, 1, 2/)
  perm(:,6) = (/3, 2, 1/)

  ! Get the no of TPs and alpha values for each
  open(2, file='smallAlpha.txt', status='old')
  read(2,*) nTraining
  read(2,*) (alpha(l), l=1,nTraining)
  close(2)

  ! Read in nArgs and the distances for each TP
  open(3, file='smallTrainingSet.txt', status='old')
  read(3,*) nDim
  do m = 1, nTraining
    read(3,*) (xTraining(n,m), n=1,nDim)
  end do
  close(3)

  nPerms = 6
  do i=1,nDim
     do j=1,nPerms
        do k=1,nTraining
           xTrainingPerm(i,j,k)=xTraining(perm(i,j),k)
        end do
     end do
  end do
  do m = 1, nPerms
    print *, xTrainingPerm(:,m,3) 
  end do

  kKernTotal=0
  do i=1,nTraining
     kSqExpAllPerms=0
     do j=1,nPerms
        kSqExpJthPerm=1        
        do k=1,nDim
           kSqExpJthPerm  =  kSqExpJthPerm * &
           (exp( - (xStar(k)-xTrainingPerm(k,j,i))**2 /2.0/lScale**2))
        end do !Dimensions (k)
        kSqExpAllPerms = kSqExpAllPerms + kSqExpJthPerm
     end do !Permuations (j)
     kKernTotal = kKernTotal + alpha(i) * kSqExpAllPerms
  end do !Training points (i)
  
  PES_GP=kKernTotal * expVar
  print *, 'Non-additive E:', PES_GP
end subroutine energyCheckCalc

end module
