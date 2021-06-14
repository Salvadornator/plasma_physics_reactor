MODULE interpos_module

  implicit none
  public :: interpos
  interface interpos
     ! simple cbsplgen0_s, and full cases so far: cbsplgen0_f
     module procedure cbsplgen0_def, cbsplgen0_sfull, cbsplgen0_vfull, cbsplgen0_ssimple, cbsplgen0_vsimple
  end interface

contains

!!$  SUBROUTINE cbsplgen0_def(XIN,YIN,nin,xout,yout,nout,taus)
!!$
!!$  SUBROUTINE cbsplgen0_defp(XIN,YIN,nin,xout,yout,youtp,nout,taus)
!!$  !
!!$  SUBROUTINE cbsplgen0_defpp(XIN,YIN,nin,xout,yout,youtp,youtpp,nout,taus)
!!$  !
!!$  SUBROUTINE cbsplgen0_defppint(XIN,YIN,nin,xout,yout,youtp,youtpp,youtint,nout,taus)
!!$  !
!!$  SUBROUTINE cbsplgen0_def2(XIN,YIN,YNEW,YINPP,nin,xout,yout,nout,taus)
!!$
!!$  SUBROUTINE cbsplgen0_def2p(XIN,YIN,YNEW,YINPP,nin,xout,yout,youtp,nout,taus)
!!$  !
!!$  SUBROUTINE cbsplgen0_def2pp(XIN,YIN,YNEW,YINPP,nin,xout,yout,youtp,youtpp,nout,taus)
!!$  !
!!$  SUBROUTINE cbsplgen0_def2ppint(XIN,YIN,YNEW,YINPP,nin,xout,yout,youtp,youtpp,youtint,nout,taus)
!!$  !
!!$    SUBROUTINE cbsplgen0_sfull(XIN,YIN,YNEW,YINPP,nin,xout,yout,youtp,youtpp,youtint, &
!!$         &  nout,ioptder,iextrapo,taus,nbc,ybc,IFLAG)
!!$  !
!!$  SUBROUTINE cbsplgen0_sfull(XIN,YIN,YNEW,YINPP,nin,xout,yout,youtp,youtpp,youtint, &
!!$  &  nout,ioptder,iextrapo,taus,nbc,ybc,IFLAG)
!!$  !
!!$  SUBROUTINE cbsplgen0_sfull(XIN,YIN,YNEW,YINPP,nin,xout,yout,youtp,youtpp,youtint, &
!!$  &  nout,ioptder,iextrapo,taus,nbc,ybc,IFLAG)
!!$  !
!!$  SUBROUTINE cbsplgen0_sfull(XIN,YIN,YNEW,YINPP,nin,xout,yout,youtp,youtpp,youtint, &
!!$  &  nout,ioptder,iextrapo,taus,nbc,ybc,IFLAG)
!!$  !
!!$  SUBROUTINE cbsplgen0_sfull(XIN,YIN,YNEW,YINPP,nin,xout,yout,youtp,youtpp,youtint, &
!!$  &  nout,ioptder,iextrapo,taus,nbc,ybc,IFLAG)
!!$  !
!!$  SUBROUTINE cbsplgen0_def(XIN,YIN,nin,xout,yout,nout,taus)

  SUBROUTINE cbsplgen0_def(XIN,YIN,nin,xout,yout,nout,taus)
    USE prec_rkind
    implicit none
    integer :: nin,nout
    REAL(RKIND) :: xin(nin), yin(nin), xout(nout)
    REAL(RKIND) :: yout(nout)
    REAL(RKIND), optional ::  taus
    !
    REAL(RKIND) ::  ynew(NIN), yinpp(NIN), ybc(6),psig(nin)
    REAL(RKIND) ::  youtp(nout), youtpp(nout), youtint(nout)
    INTEGER ioptder, NBCLFT,NBCRGT, iextrapo, nbc(2)
    REAL(RKIND) :: XBCLFT,XBCRGT, YBCLFT,YBCRGT
    REAL(RKIND) :: PXEXP0,PXEXPDL,zz
    REAL(RKIND), DIMENSION(:), ALLOCATABLE :: pamat
    REAL(RKIND), DIMENSION(:), ALLOCATABLE :: pwork
    INTEGER, DIMENSION(:), ALLOCATABLE :: kwork
    !
    integer iflag, idamat, mdamat, nbytes
    !
    ! default values
    !
    nbc=0
    ybc=0._rkind
    ioptder=0
    iextrapo=32
    !
    if ((present(taus)) .and. (taus .ge. 0._rkind)) then
      print *,'% in def: taus= ',taus
      psig=taus
    else
      zz=minval(xin(2:nin)-xin(1:nin-1))
      psig=zz**3
      print *,'% in def: taus defined as = ',psig(1)
    end if
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

  end SUBROUTINE cbsplgen0_def

  SUBROUTINE cbsplgen0_sfull(XIN,YIN,YNEW,YINPP,nin,xout,yout,youtp,youtpp,youtint, &
  &  nout,ioptder,iextrapo,taus,nbc,ybc,IFLAG)
  !
  USE prec_rkind
  implicit none
  integer nin,nout, nbc(2)
  REAL(RKIND) ::  taus
  REAL(RKIND) ::  xin(nin), yin(nin), ynew(NIN), yinpp(NIN), ybc(6),psig(nin)
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
  psig=taus
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
END SUBROUTINE cbsplgen0_sfull

  SUBROUTINE cbsplgen0_vfull(XIN,YIN,YNEW,YINPP,nin,xout,yout,youtp,youtpp,youtint, &
  &  nout,ioptder,iextrapo,psig,nbc,ybc,IFLAG)
  !
  USE prec_rkind
  implicit none
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
END SUBROUTINE cbsplgen0_vfull

  SUBROUTINE cbsplgen0_ssimple(XIN,YIN,nin,xout,yout,youtp,youtpp,youtint, &
  &  nout,ioptder,iextrapo,taus,nbc,ybc,IFLAG)
  !
  USE prec_rkind
  implicit none
  integer nin,nout, nbc(2)
  REAL(RKIND) ::  taus
  REAL(RKIND) ::  xin(nin), yin(nin), ynew(NIN), yinpp(NIN), ybc(6),psig(nin)
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
  psig=taus
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
END SUBROUTINE cbsplgen0_ssimple

  SUBROUTINE cbsplgen0_vsimple(XIN,YIN,nin,xout,yout,youtp,youtpp,youtint, &
  &  nout,ioptder,iextrapo,psig,nbc,ybc,IFLAG)
  !
  USE prec_rkind
  implicit none
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
END SUBROUTINE cbsplgen0_vsimple

end MODULE interpos_module
