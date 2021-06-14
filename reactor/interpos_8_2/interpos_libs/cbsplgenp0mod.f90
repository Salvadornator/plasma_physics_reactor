MODULE cbsplgenp0mod

  implicit none
  public :: cbsplgenp0
  interface cbsplgenp0
     module procedure cbsplgenp0_def, cbsplgenp0_sfull, cbsplgenp0_vfull, cbsplgenp0_ssimple, cbsplgenp0_vsimple
  end interface

contains

  SUBROUTINE cbsplgenp0_def(XIN,YIN,nin,xout,yout,nout,taus,period)
    !
    USE prec_rkind
    implicit none
    !
    integer :: nin,nout
    REAL(RKIND) :: xin(nin), yin(nin), xout(nout)
    REAL(RKIND) :: yout(nout)
    REAL(RKIND), optional ::  taus
    REAL(RKIND), optional ::  PERIOD
    !
    integer nbc(2)
    REAL(RKIND) ::  psig(nin)
    REAL(RKIND) ::  pynew(NIN), pyinpp(NIN), ybc(6)
    REAL(RKIND) ::  youtp(nout), youtpp(nout), youtint(nout)
    INTEGER ioptder, NBCLFT,NBCRGT, iextrapo
    REAL(RKIND) :: XBCLFT,XBCRGT, YBCLFT,YBCRGT
    REAL(RKIND) :: PXEXP0,PXEXPDL, period_eff, zz
    INTEGER :: NP1_INCLD, nin_eff
    !
    integer iflag
    !
    ! Define defaults from this short inputs version:
    !
    nbc=-1
    ioptder=0
    nin_eff = nin
    !
    if (present(period)) then
      ybc=period
    else
      ! Assumes last x point corresponds to x(1)+period and is redundant
      ybc=xin(nin_eff)-xin(1)
      nin_eff = nin-1
    end if
    period_eff = ybc(1)
    !
    if ((present(taus)) .and. (taus .ge. 0._rkind)) then
      print *,' in def: taus= ',taus
      psig=taus
    else
      zz=minval(xin(2:nin_eff)-xin(1:nin_eff-1))
      psig=zz**3
      print *,' in def: taus defined as = ',psig(1)
    end if

    NP1_INCLD=0
    IF (ABS(XIN(NIN_EFF)-XIN(1)-PERIOD_eff) .LT. 1.0e-10_RKIND*(XIN(2)-XIN(1))) THEN
      NP1_INCLD=1
      PRINT *,'IT SEEMS PERIODIC POINT X(n)=X(1)+PERIOD IS INCLUDED IN INPUT'
      nin_eff = nin_eff-1
    END IF
    !
    PXEXP0 = -1._rkind
    PXEXPDL= 1000._rkind
    CALL CBSPLGNP(XIN,YIN,PYNEW,PYINPP,Nin_eff,XOUT,YOUT,YOUTP,YOUTPP,YOUTINT, &
      &    nout,ioptder,PSIG,PERIOD,PXEXP0,PXEXPDL,IFLAG)
    !

  end SUBROUTINE cbsplgenp0_def

  SUBROUTINE cbsplgenp0_sfull(XIN,YIN,PYNEW,PYINPP,nin,xout,yout,youtp,youtpp,youtint, &
    &  nout,ioptder,taus,nbc,ybc,IFLAG)
    !
    USE prec_rkind
    implicit none
    integer :: nin,nout, nbc(2)
    REAL(RKIND) ::  taus
    REAL(RKIND) ::  xin(nin), yin(nin), ybc(6)
    REAL(RKIND) ::  pynew(NIN), pyinpp(NIN)
    REAL(RKIND) ::  xout(nout),yout(nout), youtp(nout), youtpp(nout), youtint(nout)
    INTEGER :: ioptder
    REAL(RKIND) :: XBCLFT,XBCRGT, psig(nin)
    REAL(RKIND) :: PXEXP0,PXEXPDL
    REAL(RKIND) ::  PERIOD, zz
    INTEGER :: NP1_INCLD, NBCLFT,NBCRGT, nin_eff
    !
    integer iflag
    !
    nin_eff = nin
    if (taus .lt. 0._rkind) then
      zz=minval(xin(2:nin_eff)-xin(1:nin_eff-1))
      taus=zz**3
      print *,' in def: taus defined as = ',taus
    end if
    psig = taus
    !
    ! Default for periodic boundary conditions: NBCLFT=-1, PERIOD=YBCLFT
    NBCLFT = nbc(1)
    IF (NBCLFT .GE. 0) THEN
      PRINT *,'SHOULD CALL STANDARD CBSPLGEN0'
      STOP 'CBSPLGENP0 1'
    ELSE
      ! PERIOD=YBCLFT
      PERIOD = ybc(1)
      NP1_INCLD=0
      IF (ABS(XIN(NIN_EFF)-XIN(1)-PERIOD) .LT. 1.0e-10_RKIND*(XIN(2)-XIN(1))) THEN
        NP1_INCLD=1
        PRINT *,'IT SEEMS PERIODIC POINT X(n)=X(1)+PERIOD IS INCLUDED IN INPUT'
        IF (NBCLFT .EQ. -1) THEN
          ! ! Use last x point corresponding to x(1)+period to define period
          period = XIN(NIN_EFF)-XIN(1)
          nin_eff = nin_eff - 1
        else
          ! Refuse to assume that last x point is redundant and just returns
          print *,'use nbc(1)=-1 if you want to assume that last point corresponds to x(1)+period'
          return
        end IF
      END IF
    end IF
    !
    PXEXP0 = ybc(5)
    PXEXPDL= ybc(6)
    CALL CBSPLGNP(XIN,YIN,PYNEW,PYINPP,Nin_eff,XOUT,YOUT,YOUTP,YOUTPP,YOUTINT, &
      &    nout,ioptder,PSIG,PERIOD,PXEXP0,PXEXPDL,IFLAG)
    !
  END SUBROUTINE cbsplgenp0_sfull

  SUBROUTINE cbsplgenp0_vfull(XIN,YIN,PYNEW,PYINPP,nin,xout,yout,youtp,youtpp,youtint, &
    &  nout,ioptder,psig,nbc,ybc,IFLAG)
    !
    USE prec_rkind
    implicit none
    integer :: nin,nout, nbc(2)
    REAL(RKIND) ::  psig(nin)
    REAL(RKIND) ::  xin(nin), yin(nin), ybc(6)
    REAL(RKIND) ::  pynew(NIN), pyinpp(NIN)
    REAL(RKIND) ::  xout(nout),yout(nout), youtp(nout), youtpp(nout), youtint(nout)
    INTEGER :: ioptder
    REAL(RKIND) :: XBCLFT,XBCRGT
    REAL(RKIND) :: PXEXP0,PXEXPDL, zz
    REAL(RKIND) ::  PERIOD
    INTEGER :: NP1_INCLD, NBCLFT,NBCRGT, nin_eff
    !
    integer iflag
    !
    nin_eff = nin
    if (psig(1) .lt. 0._rkind) then
      zz=minval(xin(2:nin_eff)-xin(1:nin_eff-1))
      psig=zz**3
      print *,' in def: taus defined as = ',psig(1)
    end if
    !
    ! Default for periodic boundary conditions: NBCLFT=-1, PERIOD=YBCLFT
    NBCLFT = nbc(1)
    IF (NBCLFT .GE. 0) THEN
      PRINT *,'SHOULD CALL STANDARD CBSPLGEN0'
      STOP 'CBSPLGENP0 1'
    ELSE
      ! PERIOD=YBCLFT
      PERIOD = ybc(1)
      NP1_INCLD=0
      IF (ABS(XIN(NIN_EFF)-XIN(1)-PERIOD) .LT. 1.0e-10_RKIND*(XIN(2)-XIN(1))) THEN
        NP1_INCLD=1
        PRINT *,'IT SEEMS PERIODIC POINT X(n)=X(1)+PERIOD IS INCLUDED IN INPUT'
        IF (NBCLFT .EQ. -1) THEN
          ! ! Use last x point corresponding to x(1)+period to define period
          period = XIN(NIN_EFF)-XIN(1)
          nin_eff = nin_eff - 1
        else
          ! Refuse to assume that last x point is redundant and just returns
          print *,'use nbc(1)=-1 if you want to assume that last point corresponds to x(1)+period'
          return
        end IF
      END IF
    end IF
    !
    PXEXP0 = ybc(5)
    PXEXPDL= ybc(6)
    CALL CBSPLGNP(XIN,YIN,PYNEW,PYINPP,Nin_eff,XOUT,YOUT,YOUTP,YOUTPP,YOUTINT, &
      &    nout,ioptder,PSIG,PERIOD,PXEXP0,PXEXPDL,IFLAG)
    !
  END SUBROUTINE cbsplgenp0_vfull

  SUBROUTINE cbsplgenp0_ssimple(XIN,YIN,nin,xout,yout,youtp,youtpp,youtint, &
    &  nout,ioptder,taus,nbc,ybc,IFLAG)
    !
    USE prec_rkind
    implicit none
    integer :: nin,nout, nbc(2)
    REAL(RKIND) ::  taus
    REAL(RKIND) ::  xin(nin), yin(nin), ybc(6)
    REAL(RKIND) ::  xout(nout),yout(nout), youtp(nout), youtpp(nout), youtint(nout)
    INTEGER :: ioptder
    REAL(RKIND) ::  pynew(NIN), pyinpp(NIN)
    REAL(RKIND) :: XBCLFT,XBCRGT, psig(nin)
    REAL(RKIND) :: PXEXP0,PXEXPDL
    REAL(RKIND) ::  PERIOD, zz
    INTEGER :: NP1_INCLD, NBCLFT,NBCRGT, nin_eff
    !
    integer iflag
    !
    nin_eff = nin
    if (taus .lt. 0._rkind) then
      zz=minval(xin(2:nin_eff)-xin(1:nin_eff-1))
      taus=zz**3
      print *,' in def: taus defined as = ',taus
    end if
    psig = taus
    !
    ! Default for periodic boundary conditions: NBCLFT=-1, PERIOD=YBCLFT
    NBCLFT = nbc(1)
    IF (NBCLFT .GE. 0) THEN
      PRINT *,'SHOULD CALL STANDARD CBSPLGEN0'
      STOP 'CBSPLGENP0 1'
    ELSE
      ! PERIOD=YBCLFT
      PERIOD = ybc(1)
      NP1_INCLD=0
      IF (ABS(XIN(NIN_EFF)-XIN(1)-PERIOD) .LT. 1.0e-10_RKIND*(XIN(2)-XIN(1))) THEN
        NP1_INCLD=1
        PRINT *,'IT SEEMS PERIODIC POINT X(n)=X(1)+PERIOD IS INCLUDED IN INPUT'
        IF (NBCLFT .EQ. -1) THEN
          ! ! Use last x point corresponding to x(1)+period to define period
          period = XIN(NIN_EFF)-XIN(1)
          nin_eff = nin_eff - 1
        else
          ! Refuse to assume that last x point is redundant and just returns
          print *,'use nbc(1)=-1 if you want to assume that last point corresponds to x(1)+period'
          return
        end IF
      END IF
    end IF
    !
    PXEXP0 = ybc(5)
    PXEXPDL= ybc(6)
    CALL CBSPLGNP(XIN,YIN,PYNEW,PYINPP,Nin_eff,XOUT,YOUT,YOUTP,YOUTPP,YOUTINT, &
      &    nout,ioptder,PSIG,PERIOD,PXEXP0,PXEXPDL,IFLAG)
    !
  END SUBROUTINE cbsplgenp0_ssimple

  SUBROUTINE cbsplgenp0_vsimple(XIN,YIN,nin,xout,yout,youtp,youtpp,youtint, &
    &  nout,ioptder,psig,nbc,ybc,IFLAG)
    !
    USE prec_rkind
    implicit none
    integer :: nin,nout, nbc(2)
    REAL(RKIND) ::  psig(nin)
    REAL(RKIND) ::  xin(nin), yin(nin), ybc(6)
    REAL(RKIND) ::  xout(nout),yout(nout), youtp(nout), youtpp(nout), youtint(nout)
    INTEGER :: ioptder
    REAL(RKIND) ::  pynew(NIN), pyinpp(NIN)
    REAL(RKIND) :: XBCLFT,XBCRGT
    REAL(RKIND) :: PXEXP0,PXEXPDL
    REAL(RKIND) ::  PERIOD, zz
    INTEGER :: NP1_INCLD, NBCLFT,NBCRGT, nin_eff
    !
    integer iflag
    !
    nin_eff = nin
    if (psig(1) .lt. 0._rkind) then
      zz=minval(xin(2:nin_eff)-xin(1:nin_eff-1))
      psig=zz**3
      print *,' in def: taus defined as = ',psig(1)
    end if
    !
    ! Default for periodic boundary conditions: NBCLFT=-1, PERIOD=YBCLFT
    NBCLFT = nbc(1)
    IF (NBCLFT .GE. 0) THEN
      PRINT *,'SHOULD CALL STANDARD CBSPLGEN0'
      STOP 'CBSPLGENP0 1'
    ELSE
      ! PERIOD=YBCLFT
      PERIOD = ybc(1)
      NP1_INCLD=0
      IF (ABS(XIN(NIN_EFF)-XIN(1)-PERIOD) .LT. 1.0e-10_RKIND*(XIN(2)-XIN(1))) THEN
        NP1_INCLD=1
        PRINT *,'IT SEEMS PERIODIC POINT X(n)=X(1)+PERIOD IS INCLUDED IN INPUT'
        IF (NBCLFT .EQ. -1) THEN
          ! ! Use last x point corresponding to x(1)+period to define period
          period = XIN(NIN_EFF)-XIN(1)
          nin_eff = nin_eff - 1
        else
          ! Refuse to assume that last x point is redundant and just returns
          print *,'use nbc(1)=-1 if you want to assume that last point corresponds to x(1)+period'
          return
        end IF
      END IF
    end IF
    !
    PXEXP0 = ybc(5)
    PXEXPDL= ybc(6)
    CALL CBSPLGNP(XIN,YIN,PYNEW,PYINPP,Nin_eff,XOUT,YOUT,YOUTP,YOUTPP,YOUTINT, &
      &    nout,ioptder,PSIG,PERIOD,PXEXP0,PXEXPDL,IFLAG)
    !
  END SUBROUTINE cbsplgenp0_vsimple

end MODULE cbsplgenp0mod
