SUBROUTINE cbsplgen0(XIN,YIN,YNEW,YINPP,nin,xout,yout,youtp,youtpp,youtint, &
  &  nout,ioptder,iextrapo,psig,nbc,ybc,IFLAG)
  !
  USE prec_rkind
  implicit none
  integer LENREAL
  parameter(LENREAL=8)
  integer nin,nout, nbc(2)
  REAL(RKIND), dimension(:) ::  psig(:)
  REAL(RKIND) ::  xin(nin), yin(nin), ynew(NIN), yinpp(NIN), ybc(6)
  REAL(RKIND) ::  xout(nout),yout(nout), youtp(nout), youtpp(nout), youtint(nout)
  INTEGER ioptder, NBCLFT,NBCRGT, iextrapo
  REAL(RKIND) :: XBCLFT,XBCRGT, YBCLFT,YBCRGT
  REAL(RKIND) :: PXEXP0,PXEXPDL
  REAL(RKIND), DIMENSION(:), ALLOCATABLE :: pamat
  REAL(RKIND), DIMENSION(:), ALLOCATABLE :: pwork
  INTEGER, DIMENSION(:), ALLOCATABLE :: kwork
  !
  integer iflag, idamat, mdamat, nbytes
  !
  !
  !%OS      NBCLFT = 1
  !%OS      NBCRGT = 1
  !%OS      YBCLFT = 1.E32_RKIND
  !%OS      YBCRGT = 1.E32_RKIND
  !%OS      print *,'xin,yin= ',xin,yin
  print *,'size(psig)= ',size(psig)
  NBCLFT = nbc(1)
  NBCRGT = nbc(2)
  YBCLFT = ybc(1)
  YBCRGT = ybc(2)
  if (NBCLFT .ge. 10) XBCLFT = ybc(3)
  if (NBCRGT .ge. 10) XBCRGT = ybc(4)
  !
  IDAMAT = 3
  IF (PSIG(1) .EQ. 0._RKIND) IDAMAT = 2
  IF (NBCLFT.GE.10 .OR. NBCRGT.GE.10 .OR. NBCLFT.EQ.2 &
    &  .OR. NBCRGT.EQ.2) IDAMAT = 3*(IDAMAT-1)+1
  IF (NBCLFT .GE. 10) IDAMAT = IDAMAT + 1
  IF (NBCRGT .GE. 10) IDAMAT = IDAMAT + 2
  mdamat = IDAMAT*nin
  ALLOCATE(pamat(mdamat))
  ALLOCATE(pwork(2*nin))
  ALLOCATE(kwork(5*nin))
  PXEXP0 = ybc(5)
  PXEXPDL= ybc(6)
  CALL CBSPLGEN(XIN,YIN,YNEW,YINPP,Nin,XOUT,YOUT,YOUTP,YOUTPP,YOUTINT, &
    &    nout,ioptder,PSIG,KWORK,PWORK,PAMAT,mdamat,NBCLFT, &
    &    NBCRGT,XBCLFT,XBCRGT,YBCLFT,YBCRGT,iextrapo,PXEXP0,PXEXPDL &
    &    ,IFLAG)
  !
  deallocate(pamat)
  deallocate(pwork)
  deallocate(kwork)
END SUBROUTINE cbsplgen0
