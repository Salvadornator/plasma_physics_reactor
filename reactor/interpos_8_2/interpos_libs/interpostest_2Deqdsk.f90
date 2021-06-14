program interpostest_2Deqdsk
  USE prec_rkind
  USE interpos_module
  implicit none
  integer, parameter :: nrout=257, nzout=257, nbout=512
  REAL(RKIND), PARAMETER :: TWOPI=6.283185307179586476925286766559005768394_rkind
  real(rkind), allocatable :: psimesh(:), dummyprof(:), Rplas(:), Zplas(:), psi(:,:), rmesh(:), zmesh(:), &
    & rwall(:), zwall(:), psitmp(:,:), pprime(:), ffprime(:), thetain(:), rhoin(:), theta_sorted(:), rho_sorted(:)
  real(rkind) :: psimeshout(nrout), dummyprofout(nrout), Rplasout(nbout), Zplasout(nbout), &
    & psiout(nrout,nzout), rmeshout(nrout), zmeshout(nzout), &
    & thetaout(nbout), rhoout(nbout), rhooutp(nbout), rhooutpp(nbout),  &
    & rhoout0(nbout), rhooutp0(nbout), rhooutpp0(nbout),  &
    & pprime_rz(nrout,nzout), ffprime_rz(nrout,nzout), gradshaf(nrout,nzout), rjphi(nrout,nzout), &
    & dpsidr(nrout,nzout), d2psidr2(nrout,nzout), d2psidz2(nrout,nzout)
  real(rkind) :: tension1D, rboxlen, zboxlen, r0exp, rboxlft, zboxmid, &
    & raxis, zaxis, psiaxis, psiedge, b0exp, currt, zdummy, mu0
  integer :: nin, i, ir, iz, nrin, nzin, nbound, nwall, nbound_eff
  integer, allocatable :: i_sorted(:)
  character(LEN=48) :: dumheader
  !
  print *,'% test for 2D interpolation and periodic plasma boundary smoothing'
  !  open(9,file='EQDSK.29866t1.5003')
  !  open(9,file='eqdsksigns.36151t0.5000')
  open(9,file='eqdsksigns.31837t1.0000')
  read(9,'(A48,3I4)') dumheader,i,nrin,nzin
  allocate(psimesh(nrin))
  allocate(dummyprof(nrin))
  allocate(ffprime(nrin))
  allocate(pprime(nrin))
  allocate(psi(nrin,nzin))
  allocate(rmesh(nrin))
  allocate(zmesh(nzin))
9381 FORMAT(1P,5E16.9)
  write(0,*) 'dumheader,i,nrin,nzin= ',dumheader,i,nrin,nzin
  read(9,9381) rboxlen,zboxlen,r0exp,rboxlft,zboxmid
  read(9,9381) raxis,zaxis,psiaxis,psiedge,b0exp
  read(9,9381) currt
  read(9,9381) zdummy
  write(0,*) 'dumheader,i,nrin,nzin= ',dumheader,i,nrin,nzin
  !     f
  read(9,9381) (dummyprof(i), i=1,nrin)
  !     p
  read(9,9381) (dummyprof(i), i=1,nrin)
  !     ff'
  write(0,*) 'ffprime,i,nrin,nzin= ',dumheader,i,nrin,nzin
  read(9,9381) (ffprime(i), i=1,nrin)
  !     p'
  write(0,*) 'pprime,i,nrin,nzin= ',dumheader,i,nrin,nzin
  read(9,9381) (pprime(i), i=1,nrin)
  !     psi(r,z)
  write(0,*) 'psi,i,nrin,nzin= ',dumheader,i,nrin,nzin
  read(9,9381) ((psi(ir,iz),ir=1,nrin),iz=1,nzin)
  !     q
  write(0,*) 'dumheader,i,nrin,nzin= ',dumheader,i,nrin,nzin
  read(9,*) (dummyprof(i),i=1,nrin)
  !     plasma and wall boundary
  read(9,*) nbound, nwall
  write(0,*) 'nbound,nwall= ',nbound,nwall
  allocate(rplas(nbound))
  allocate(zplas(nbound))
  allocate(rwall(nwall))
  allocate(zwall(nwall))
  read(9,9381) (rplas(i),  zplas(i), i=1,nbound)
  read(9,9381) (rwall(i),  rwall(i), i=1,nwall)

  !     construct input and output meshes
  zdummy=rboxlen/real(nrin-1,rkind)
  rmesh=(/ (rboxlft + real((i-1),rkind)*zdummy,i=1,nrin) /)
  zdummy=zboxlen/real(nzin-1,rkind)
  zmesh=(/ (zboxmid-0.5_rkind*zboxlen+real((i-1),rkind)*zdummy,i=1,nzin) /)
  zdummy=rboxlen/real(nrout-1,rkind)
  rmeshout=(/ (rboxlft + real((i-1),rkind)*zdummy,i=1,nrout) /)
  zdummy=zboxlen/real(nzout-1,rkind)
  zmeshout=(/ (zboxmid-0.5_rkind*zboxlen+real((i-1),rkind)*zdummy,i=1,nzout) /)
  zdummy=twopi/real(nbout-1,rkind)
  thetaout=(/ (real((i-1),rkind)*zdummy,i=1,nbout) /)

  zdummy=(psiedge-psiaxis)/real(nrin-1,rkind)
  psimesh=(/ (psiaxis + real((i-1),rkind)*zdummy,i=1,nrin) /)
  !
  ! 1. Example just to fit 2D psi(r,z) on psiout(rout,zout) with 2 times 1D spline
  ! Start from finer mesh since expect most important (although probably does not matter)
  !
  ! default tension:
  ! tension1D=0.5.*(((rmesh(2)-rmesh(1))+(zmesh(2)-zmesh(1)))./2).^3
  tension1D = -0.5_rkind
  if ((rmesh(2)-rmesh(1)) .le. (zmesh(2)-zmesh(1))) then
    allocate(psitmp(nrout,nzin))
    do iZ=1,nzin
      call interpos(rmesh,psi(:,iZ),nrin,tension=tension1D,xout=rmeshout,yout=psitmp(:,iZ))
    end do
    do iR=1,nrout
      call interpos(zmesh,psitmp(iR,:),nzin,tension=tension1D,xout=zmeshout,yout=psiout(iR,:));
    end do
  else
    allocate(psitmp(nrin,nzout))
    do iR=1,nrin
      call interpos(zmesh,psi(iR,:),nzin,tension=tension1D,xout=zmeshout,yout=psitmp(iR,:))
    end do
    do iZ=1,nzout
      call interpos(rmesh,psitmp(:,iZ),nrin,tension=tension1D,xout=rmeshout,yout=psiout(:,iZ));
    end do
  end if
  !
  ! 2. Calculate 1st and 2nd derivatives and then check Grad-Shafranov equation. 
  ! Since psiout known, it is already good std spline (thus use tension=0.)
  ! (could calculate derivatives at the same time as calculating psiout)
  !
  do iZ=1,nzout
    call interpos(rmeshout,psiout(:,iZ),nrout,tension=0._rkind,youtp=dpsidr(:,iZ),youtpp=d2psidr2(:,iZ));
  end do
  ! to compute d2/dz2 on full (nrout,nzout) needs extra call, can use final psi with 0 tension
  do iR=1,nrout
    call interpos(zmeshout,psiout(iR,:),nzout,tension=0._rkind,youtpp=d2psidz2(iR,:));
  end do

  ! need pprime, ffprime and R on full (rout,zout) mesh
  mu0=2.e-07_rkind*twopi
  do iR=1,nrout
    call interpos(psimesh,pprime,nrin,tension=1.e-10_rkind,xout=psiout(iR,:),yout=pprime_rz(iR,:),option=-63) ! 0 outside
    call interpos(psimesh,ffprime,nrin,tension=1.e-10_rkind,xout=psiout(iR,:),yout=ffprime_rz(iR,:),option=-63) ! 0 outside
    gradshaf(iR,:) = d2psidr2(iR,:) -1._rkind/rmeshout(iR)*dpsidr(iR,:) + d2psidz2(iR,:)
    rjphi(iR,:) = -rmeshout(iR)**2.*mu0*pprime_rz(iR,:)-ffprime_rz(iR,:)
  end do
!!$  iR=120
!!$  write(0,*) 'd2psidr2(iR,1)= ',d2psidr2(iR,1)
!!$  write(0,*) 'gradshaf(iR,1)= ',gradshaf(iR,1)
!!$  !write(0,*) 'gradshaf(:,1)= ',(iR,gradshaf(iR,1),iR=1,nrout)
!!$  write(0,*) 'dpsidr(iR,1)= ',dpsidr(iR,1)
!!$  write(0,*) 'd2psidz2(iR,1)= ',d2psidz2(iR,1)
!!$  write(0,*) 'pprime_rz(iR,1)= ',pprime_rz(iR,1)
!!$  write(0,*) 'ffprime_rz(iR,1)= ',ffprime_rz(iR,1)
!!$  write(0,*) 'rmesh(iR)= ',rmesh(iR)
!!$  stop
  !
  ! 3. Example to smooth (experimental) plasma boundary
  !
  allocate(thetain(nbound))
  allocate(rhoin(nbound))
  do i=1,nbound
    thetain(i) = atan2((zplas(i)-zaxis),(rplas(i)-raxis))
    if (thetain(i) .lt. 0._rkind) thetain(i) = thetain(i) + twopi
    rhoin(i)=sqrt((rplas(i)-raxis)**2 + (zplas(i)-zaxis)**2)
  end do
  if (abs(thetain(nbound)-thetain(1)) .lt. 1e-07_rkind) then
    write(0,*) 'end theta points given twice, remove it'
    nbound_eff=nbound - 1
  end if
  ! sort rho mesh
  allocate(theta_sorted(nbound_eff))
  allocate(rho_sorted(nbound_eff))
  allocate(i_sorted(nbound_eff))
  call indexx(nbound_eff,thetain(1:nbound_eff),i_sorted)
  ! write(*,'(2i5,1p2e15.4)') (i,i_sorted(i),thetain(i),thetain(i_sorted(i)),i=1,nbound_eff)
  theta_sorted = thetain(i_sorted(1:nbound_eff))
  rho_sorted = rhoin(i_sorted(1:nbound_eff))
  ! write(*,'(i5,1p2e15.4)') (i,theta_sorted(i),rho_sorted(i),i=1,nbound_eff)
  call interpos(theta_sorted,rho_sorted,nbound_eff,nout=nbout,tension=-0.1_rkind, &
    & xout=thetaout,yout=rhoout,youtp=rhooutp,youtpp=rhooutpp,nbc=-1,ybc=twopi)
  do i=1,nbout
    Rplasout(i) = raxis + cos(thetaout(i))*rhoout(i)
    Zplasout(i) = zaxis + sin(thetaout(i))*rhoout(i)
  end do
  ! standard spline for comparison
  call interpos(theta_sorted,rho_sorted,nbound_eff,nout=nbout,tension=0._rkind, &
    & xout=thetaout,yout=rhoout0,youtp=rhooutp0,youtpp=rhooutpp0,nbc=-1,ybc=twopi)

  ! write data ready for matlab
  print *,' rmesh=[ '
  write(*,'(1p1e25.15)') (rmesh(i), i=1,nrin)
  print *,'];'
  print *,' zmesh=[ '
  write(*,'(1p1e25.15)') (zmesh(i), i=1,nzin)
  print *,'];'
  print *,' psi=[ '
  write(*,'(1p1e25.15)') ((psi(ir,iz), ir=1,nrin),iz=1,nzin)
  print *,'];'

  print *,' rmeshout=[ '
  write(*,'(1p1e25.15)') (rmeshout(i), i=1,nrout)
  print *,'];'
  print *,' zmeshout=[ '
  write(*,'(1p1e25.15)') (zmeshout(i), i=1,nzout)
  print *,'];'
  print *,' psiout=[ '
  write(*,'(1p1e25.15)') ((psiout(ir,iz), ir=1,nrout),iz=1,nzout)
  print *,'];'

  print *,' gradshaf=[ '
  write(*,'(1p1e25.15)') ((gradshaf(ir,iz), ir=1,nrout),iz=1,nzout)
  print *,'];'
  print *,' rjphi=[ '
  write(*,'(1p1e25.15)') ((rjphi(ir,iz), ir=1,nrout),iz=1,nzout)
  print *,'];'

  print *,' rzplas=[ '
  write(*,'(1p2e25.15)') (rplas(i), zplas(i), i=1,nbound)
  print *,'];'

  print *,' rzplasout=[ '
  write(*,'(1p2e25.15)')  (rplasout(i), zplasout(i), i=1,nbout)
  print *,'];'

  do i=1,nbout
    Rplasout(i) = raxis + cos(thetaout(i))*rhoout0(i)
    Zplasout(i) = zaxis + sin(thetaout(i))*rhoout0(i)
  end do
  print *,' rzplasout0=[ '
  write(*,'(1p2e25.15)')  (rplasout(i), zplasout(i), i=1,nbout)
  print *,'];'

  print *,' thetarhoin=[ '
  write(*,'(1p2e25.15)') (theta_sorted(i), rho_sorted(i), i=1,nbound_eff)
  print *,'];'

  print *,' thetarhopppout=[ '
  write(*,'(1p4e25.15)') (thetaout(i), rhoout(i), rhooutp(i), rhooutpp(i), i=1,nbout)
  print *,'];'

  print *,' thetarhopppout0=[ '
  write(*,'(1p4e25.15)') (thetaout(i), rhoout0(i), rhooutp0(i), rhooutpp0(i), i=1,nbout)
  print *,'];'


end program interpostest_2Deqdsk
