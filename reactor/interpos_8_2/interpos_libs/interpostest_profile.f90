program interpostest_profile
  USE prec_rkind
  USE interpos_module
  implicit none
  integer, parameter :: nout=3000
  real(rkind), allocatable :: xin(:), yin(:), yinnew(:), sigma(:)
  real(rkind) :: xout(nout), yout(nout), youtp(nout), youtpp(nout), youtint(nout)
  real(rkind) :: taus
  integer :: nin, i
  !
  print *,'% test for rho profile fitting'
  open(9,file='Teprofin.data')
  read(9,*) nin
  allocate(xin(nin+1))
  allocate(yin(nin+1))
  allocate(yinnew(nin+1))
  allocate(sigma(nin+1))
  read(9,*) (xin(i),yin(i),sigma(i),i=2,nin+1)
  ! add a point at rho=0, with a large error bar, to impose dy/drho=0
  xin(1)=0._rkind
  yin(1)=yin(2)
  sigma(1)=1000._rkind * sigma(2)
  !
  do i=1,nout
    xout(i)=real(i-1,rkind)*1._rkind/real(nout-1,rkind)
  end do
  !
  taus=-0.1_rkind
  call interpos(xin,yin,nin+1,nout,taus,xout=xout,yout=yout,youtp=youtp, &
& nbc=(/1,0/),ybc=(/0._rkind,0._rkind/),sigma=sigma)
  print *,' in=[ '
  write(*,'(1p3e25.15)') (xin(i), yin(i),sigma(i), i=1,nin+1)
  print *,'];'
  print *,' out=[ '
  write(*,'(1p3e25.15)') (xout(i), yout(i), youtp(i), i=1,nout)
  print *,'];'
  taus=0._rkind
  call interpos(xin,yin,nin+1,nout,tension=taus,xout=xout,yout=yout,youtp=youtp, &
    & nbc=(/1,0/),ybc=(/0._rkind,0._rkind/),sigma=sigma)
  print *,' out0=[ '
  write(*,'(1p3e25.15)') (xout(i), yout(i), youtp(i), i=1,nout)
  print *,'];'

  call interpos(xin,yin,nin,nout,tension=taus,xout=xout,yout=yout,youtp=youtp)
  print *,' out02=[ '
  write(*,'(1p3e25.15)') (xout(i), yout(i), youtp(i), i=1,nout)
  print *,'];'

end program interpostest_profile
