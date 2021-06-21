PROGRAM PES_GP
  implicit none
  double precision, dimension(3) :: xStar
  double precision, dimension(3,337) :: xTraining
  double precision, dimension(3,6,337) :: xTrainingPerm
  double precision, dimension(337) :: alpha
  double precision :: U, leng, sig_f
  integer :: Perm(6,3)
  integer :: i,j,k,nDims,p,q,o
  double precision kSqExpAllPerms, kSqExpJthPerm, kKernTotal

  leng = 1.010933217823705432e-01
  sig_f = 5.065189615948498530e-05
  xStar = (/0.5, 0.5, 0.5/)
  open(3, file='alpha.txt', status='old')
  read(3,*) nDims
  read(3,*) (alpha(o), o=1,337)
  close(3)
  open(4, file='trainingSet.txt', status='old')
  read(4,*) nDims
  do p = 1, 337
    read(4,*) (xTraining(q,p), q=1,3)
  end do
  !xTraining = TRANSPOSE(xTraining)
  close(4)
  Perm(1,:) = (/1, 2, 3/)
  Perm(2,:) = (/1, 3, 2/)
  Perm(3,:) = (/2, 1, 3/)
  Perm(4,:) = (/2, 3, 1/)
  Perm(5,:) = (/3, 1, 2/)
  Perm(6,:) = (/3, 2, 1/) 

  do i=1,3
     do j=1,6
        do k=1,337
           xTrainingPerm(i,j,k)=xTraining(Perm(j,i),k)
        end do
     end do
  end do

  !! Non-symmetric way
  !kKernTotal=0
  !do i=1,nTraining
  !   kSqExpAllPerms=1.0;
  !   do k=1,nDim
  !      kSqExpAllPerms=kSqExpAllPerms * &
  !           ( exp( - (xStar(k)-xTraining(k,i))**2 /2.0/lScale(k)**2) )
  !   end do
  !   kKernTotal = kKernTotal + alpha(i)*kSqExpAllPerms
  !end do

  !Symmetric way
  kKernTotal=0
  do i=1,337
     kSqExpAllPerms=0
     do j=1,6
        kSqExpJthPerm=1        
        do k=1,3
           kSqExpJthPerm  =  kSqExpJthPerm * &
                ( exp( - (xStar(k)-xTrainingPerm(k,j,i))**2/2.0/leng**2) )
        end do !Dimensions (k)
        kSqExpAllPerms = kSqExpAllPerms + kSqExpJthPerm
     end do !Permuations (
     kKernTotal = kKernTotal + alpha(i) * kSqExpAllPerms
  end do !Training points (i)
  
  U=kKernTotal * sig_f
  print *, "The energy of the triplet is", U
END PROGRAM PES_GP
