program interpostest
  USE prec_rkind
!!$  USE interpos
  implicit none
  integer, parameter :: nin=30, nout=3000
  real(rkind) :: xin(nin), yin(nin), yinspl(nin), yinsplpp(nin)
  real(rkind) :: xout(nout), yout(nout), youtp(nout), youtpp(nout), youtint(nout)
  integer :: i, ioptder, iextrapol, nbc(2), iflag, icount, icount2, irate, icmax
  real(rkind) :: taus, taus0, ybc(6), sigma(nin)
  REAL(RKIND), PARAMETER :: TWOPI=6.283185307179586476925286766559005768394_rkind
  real(rkind) :: rannumbers(nin), zepsilon
  integer :: info_time(8), i_msec, nsize_random, i_change_random_seed
  external cbsplgen0
  !
  print *,'% do same test as interpostest case 3'
  print *,'% sin(x+randn(1,length(x))*epsilon).^2;'
  ! gets a random number
  call random_number(rannumbers)
  print *,'% 1st 5 randoms: '
  write(*,'(a,1p5e12.4)') '% ',rannumbers(1:5)
  call date_and_time(VALUES=info_time)
  write(*,'(a15,8I5)') '% info_time= ',info_time
  i_msec = 1000 * info_time(7) + info_time(8)
  call random_seed(nsize_random)
  print *,'% nsize_random= ',nsize_random
  i_change_random_seed = 0
  if (i_change_random_seed.EQ.1) then
    print *,'% change random seed: ',(i*i_msec,i=1,nsize_random)
    call random_seed(PUT=(/(i*i_msec,i=1,nsize_random)/))
  end if
  call random_number(rannumbers)
  print *,'% 2nd 5 randoms: '
  write(*,'(a,1p5e12.4)') '% ',rannumbers(1:5)
  !
  ! x-mesh between 0 and 6 pi
  ! y=sin(x+random*epsilon)^2
  zepsilon=0.03_rkind;
  zepsilon=0.05_rkind;
  ! zepsilon=0.0_rkind;
  call random_number(rannumbers)
  do i=1,nin
    xin(i)=TWOPI*real(i-1,rkind)/real(nin-1,rkind)
    yin(i)=sin(xin(i))**2 + 2.*(rannumbers(i)-0.5_rkind)*zepsilon
  end do
  do i=1,nout
    xout(i)=TWOPI/2._rkind*(-1._RKIND+3.5_rkind*real(i-1,rkind)/real(nout-1,rkind))
  end do
  !
  ! prepare parameters for cbsplgen0
  ioptder = 3
  iextrapol = 1
  nbc = 0
  ybc = 0.0_rkind
  ! more complicated option
  !      nbc(1)=2
  !      ybc(1)=1._rkind
  !
  ! call spline with tension
  taus = 1.0e-02_rkind
  ! taus = 3.0e-03_rkind
  ! taus = 1.0e-06_rkind
  print *,'% taus= ',taus
  sigma = taus
  call system_clock(icount, irate, icmax)
  write(7,*) ' icount, irate, icmax= ',icount, irate, icmax
  call cbsplgen0(xin,yin,nin,xout,yout,youtp,youtpp,youtint,nout, &
    & ioptder,iextrapol,sigma,nbc,ybc,iflag)
  call system_clock(icount2, irate, icmax)
  write(7,*) ' icount2, irate, icmax, icount2-icount= ',icount2, irate, icmax, icount2-icount
  write(7,*) ' dtime = ',real(icount2-icount,rkind) / real(irate,rkind)
  !
  print *,'% iflag = ',iflag
  print *,' in=[ '
  write(*,'(1p2e25.15)') (xin(i), yin(i), i=1,nin)
  print *,'];'
  print *,' out=[ '
  write(*,'(1p5e25.15)') (xout(i), yout(i), youtp(i), youtpp(i), youtint(i), i=1,nout)
  print *,'];'

  ! test using splibnd
  ! first calculated spline function and 2nd der: yinspl and yinsplpp on xin
  call system_clock(icount)
  call cbsplgen0(xin,yin,nin,xin,yinspl,youtp,yinsplpp,youtint,nin, &
    & 2,iextrapol,sigma,nbc,ybc,iflag)
  call system_clock(icount2)
  write(7,*) ' icount2-icount, dtime = ',icount2-icount,real(icount2-icount,rkind) / real(irate,rkind)
  ! Then can calculate spline on any new x value
  call system_clock(icount)
  CALL SPLIBNDA(xin,yinspl,yinsplpp,nin,xout,yout, youtp, youtpp,youtint,nout,iextrapol,ioptder)
!!$      &      iextrapol)
!!$  do i=1,nout
!!$    CALL SPLIBND(xin,yinspl,yinsplpp,nin,xout(i),yout(i), youtp(i), youtpp(i), &
!!$      &      iextrapol)
!!$  end do
  call system_clock(icount2)
  write(7,*) ' icount2-icount, dtime = ',icount2-icount,real(icount2-icount,rkind) / real(irate,rkind)
  !
  print *,' outsplibnd=[ '
  write(*,'(1p5e25.15)') (xout(i), yout(i), youtp(i), youtpp(i), youtint(i), i=1,nout)
  print *,'];'
  !
  ! call spline with tension=0
  taus0=0._rkind
  sigma=taus0
  call system_clock(icount, irate, icmax)
  write(7,*) ' icount, irate, icmax= ',icount, irate, icmax
  call cbsplgen0(xin,yin,nin,xout,yout,youtp,youtpp,youtint,nout, &
    & ioptder,iextrapol,sigma,nbc,ybc,iflag)
  call system_clock(icount2, irate, icmax)
  write(7,*) ' icount2, irate, icmax, icount2-icount= ',icount2, irate, icmax, icount2-icount
  write(7,*) ' dtime = ',real(icount2-icount,rkind) / real(irate,rkind)
  !
  print *,'% iflag = ',iflag
  print *,' out0=[ '
  write(*,'(1p5e15.6)') (xout(i), yout(i), youtp(i), youtpp(i), youtint(i), i=1,nout)
  print *,'];'
  !
  ! test with periodic conditions
  ! call spline with tension
  call system_clock(icount, irate, icmax)
  write(7,*) ' icount, irate, icmax= ',icount, irate, icmax
  nbc(1)=-1
  ybc(1)=xin(nin)-xin(1)
  sigma=taus
  call cbsplgenp0(xin(1:nin-1),yin(1:nin-1),nin-1,xout,yout,youtp,youtpp,youtint,nout, &
    & ioptder,sigma,nbc,ybc,iflag)
  call system_clock(icount2, irate, icmax)
  write(7,*) ' icount2, irate, icmax, icount2-icount= ',icount2, irate, icmax, icount2-icount
  write(7,*) ' dtime = ',real(icount2-icount,rkind) / real(irate,rkind)
  !
  print *,'% iflag = ',iflag
  print *,' outper=[ '
  write(*,'(1p5e15.6)') (xout(i), yout(i), youtp(i), youtpp(i), youtint(i), i=1,nout)
  print *,'];'

end program interpostest
