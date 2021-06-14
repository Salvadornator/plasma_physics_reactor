SUBROUTINE cbsplgenp0(XIN,YIN,PYNEW,PYINPP,nin,xout,yout,youtp,youtpp,youtint, &
  &  nout,ioptder,psig,nbc,ybc,IFLAG)
  !
  USE prec_rkind
  implicit none
  integer nin,nout, nbc(2)
  REAL(RKIND) ::  psig(nin)
  REAL(RKIND) ::  xin(nin), yin(nin), pynew(NIN), pyinpp(NIN), ybc(6)
  REAL(RKIND) ::  xout(nout),yout(nout), youtp(nout), youtpp(nout), youtint(nout)
  INTEGER ioptder, NBCLFT,NBCRGT, iextrapo
  REAL(RKIND) :: XBCLFT,XBCRGT, YBCLFT,YBCRGT
  REAL(RKIND) :: PXEXP0,PXEXPDL
  REAL(RKIND), ALLOCATABLE :: pamat(:)
  REAL(RKIND), ALLOCATABLE :: pwork(:)
  INTEGER, ALLOCATABLE :: kwork(:)
  REAL(RKIND) ::  PERIOD
  INTEGER :: NP1_INCLD
  !
  integer iflag, idamat, mdamat, nbytes
  !
  ! Default for periodic boundary conditions: NBCLFT=-1, PERIOD=YBCLFT
  NBCLFT=-1
  NBCLFT = nbc(1)
  IF (NBCLFT .GE. 0) THEN
    PRINT *,'SHOULD CALL STANDARD CBSPLGEN0'
    STOP 'CBSPLGENP0 1'
  ELSEIF (NBCLFT .EQ. -1) THEN
    ! PERIOD=YBCLFT
    YBCLFT = ybc(1)
    PERIOD=YBCLFT
    NP1_INCLD=0
    IF (ABS(XIN(NIN)-XIN(1)-PERIOD) .LT. 1.0e-10_RKIND*(XIN(2)-XIN(1))) THEN
      NP1_INCLD=1
      PRINT *,'IT SEEMS PERIODIC POINT X(n)=X(1)+PERIOD IS INCLUDED IN INPUT'
    END IF
  ELSEIF (NBCLFT .EQ. -2) THEN
    PERIOD=YBCLFT
  ELSE
    PRINT *,'OPTION NBCLFT= ',NBCLFT,' NOT YET AVAILABLE'
    STOP 'CBSPLGENP0 2'
  END IF
  !
  IDAMAT = 13
  IF (PSIG(1) .EQ. 0._RKIND) IDAMAT = 7
  mdamat = IDAMAT*nin
  ALLOCATE(pamat(mdamat))
  ALLOCATE(pwork(3*nin))
  ALLOCATE(kwork(6*nin))
  PXEXP0 = ybc(5)
  PXEXPDL= ybc(6)
  CALL CBSPLGNP(XIN,YIN,PYNEW,PYINPP,Nin,XOUT,YOUT,YOUTP,YOUTPP,YOUTINT, &
    &    nout,ioptder,PSIG,PERIOD,KWORK,PWORK,PAMAT,mdamat,PXEXP0,PXEXPDL,IFLAG)
  DEALLOCATE(pamat)
  DEALLOCATE(pwork)
  DEALLOCATE(kwork)
  !
END SUBROUTINE cbsplgenp0
