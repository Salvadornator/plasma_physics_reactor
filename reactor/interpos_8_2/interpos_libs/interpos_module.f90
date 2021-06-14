MODULE interpos_module
  !
  ! Examples of calls to this generic routine:
  !
  ! call interpos(xin,yin,nin,nout,xout=xout,yout=yout)
  ! call interpos(xin,yin,nin,nout,tension,xout,yout)
  ! call interpos(xin,yin,nin,nout,tension,xout,yout,youtp)
  ! call interpos(xin,yin,nin,nout,tension,xout=xout,yout=yout,nbc=bc_type,ybc=bc_values)
  ! call interpos(xin,yin,nin,nout,tension,xout,yout,nbc=bc_type,ybc=bc_values,option=32)
  !
  ! The full possibilities can be expressed as follows:
  !
  ! call interpos(xin,yin[[,yinnew],yinpp,]nin[,nout,tension],xout,yout,[youtp,youtpp,youtint,nbc,ybc,sigma,option,info])
  ! with nout optional if equals one (xout, yout, etc are scalars)
  !
  ! These are grouped into 3 main routines:
  !
  ! main:
  ! call interpos(xin,yin,nin,[xout,yout,nout,tension,youtp,youtpp,youtint,nbc,ybc,sigma,option])
  ! These two to compute 1st yinpp and then to just evaluate for different xout values with second calls
  ! call interpos(xin,yin,yinnew,yinpp,nin,[tension,nbc,ybc,sigma])                 ! only computes yinpp
  ! call interpos(xin,yinnew,yinpp,nin,nout,xout,yout,youtp,youtpp,youtint,option])    ! only computes yout's
  !
  ! Each of these allow for scalar xout, yout, etc which means that nout is also optional in this case
  !
  !
  implicit none
  public :: interpos
  interface interpos
     module procedure interpos_def, interpos_defper, interpos_yinpp, interpos_interp &
       & ,interpos_defxscal, interpos_defperxscal, interpos_interpxscal
  end interface

contains

  SUBROUTINE interpos_def(XIN,YIN,NIN,NOUT,tension,xout,yout,youtp,youtpp,youtint,nbc,ybc,sigma,option,info)
    !
    USE prec_rkind
    implicit none
    !
    INTERFACE
       SUBROUTINE INTLINEAR(PXIN,PYIN,KNIN,PXOUT,PYOUT,PYOUTP,PYOUTPP,PYOUTINT,KNOUT,KOPTDER,KEXTRAPO)
         USE PREC_RKIND
         IMPLICIT NONE
         ! arguments
         INTEGER :: KNIN, KNOUT, KOPTDER, KEXTRAPO
         REAL(RKIND) :: PXIN(KNIN), PYIN(KNIN), PXOUT(KNOUT)
         REAL(RKIND) :: PYOUT(KNOUT), PYOUTP(KNOUT), PYOUTPP(KNOUT), PYOUTINT(KNOUT)
       END SUBROUTINE INTLINEAR
       SUBROUTINE INTQUADRATIC(PXIN,PYIN,KNIN,PXOUT,PYOUT,PYOUTP,PYOUTPP,PYOUTINT,KNOUT,KOPTDER,KOPTXPOL,NBC)
         USE PREC_RKIND
         IMPLICIT NONE
         ! arguments
         INTEGER :: KNIN, KNOUT, KOPTDER, KOPTXPOL
         INTEGER, OPTIONAL ::  NBC
         REAL(RKIND) :: PXIN(KNIN), PYIN(KNIN), PXOUT(KNOUT)
         REAL(RKIND):: PYOUT(KNOUT), PYOUTP(KNOUT), PYOUTPP(KNOUT), PYOUTINT(KNOUT)
         !
       END SUBROUTINE INTQUADRATIC
       SUBROUTINE CBSPLGEN(PXIN,PYIN,PYINNEW,PYINPP,KNIN,PXOUT,PYOUT, &
            &  PYOUTP,PYOUTPP,PYOUTINT,KNOUT,KOPT,PSIG,NBCLFT &
            &  ,NBCRGT,XBCLFT,XBCRGT,YBCLFT,YBCRGT,KEXTRPOL,PXEXP0,PXEXPDL,KFLAG)
         USE prec_rkind
         implicit none
         REAL(RKIND) :: EPSILON
         PARAMETER(EPSILON = 1.0E-10_RKIND)
         ! arguments
         integer :: knin, knout, kopt, nbclft, nbcrgt, kextrpol
         integer:: kflag
         real(rkind) :: PXIN(KNIN), PYIN(KNIN), PXOUT(KNOUT), PSIG(KNIN), &
              & xbclft, xbcrgt, ybclft, ybcrgt, pxexp0, pxexpdl
         real(rkind):: PYINNEW(KNIN), PYINPP(KNIN), PYOUT(KNOUT), PYOUTP(KNOUT), &
              & PYOUTPP(KNOUT), PYOUTINT(KNOUT)
       END SUBROUTINE CBSPLGEN
    END INTERFACE
    !
    !
    INTEGER :: NIN
    INTEGER,optional :: NOUT
    REAL(RKIND) ::  XIN(NIN), YIN(NIN)
    REAL(RKIND), optional ::  tension
    REAL(RKIND), optional ::  xout(:)
    REAL(RKIND), optional ::  yout(:), youtp(:), youtpp(:), youtint(:)
    REAL(RKIND), optional ::  ybc(:), sigma(:)
    INTEGER, optional :: nbc(2), option
    INTEGER, optional :: info
    ! need variables computed even if not given as argument
    REAL(RKIND), allocatable ::  xout2(:), yout2(:), youtp2(:), youtpp2(:), youtint2(:), xin_eff(:), yin_eff(:)
    REAL(RKIND) ::  sigma2(NIN), zdx
    integer :: info2, nout2, nin_eff, i, j
    !
    INTEGER ioptder, NBCLFT,NBCRGT, iextrapo, inttype, optabs
    REAL(RKIND) :: XBCLFT,XBCRGT, YBCLFT,YBCRGT
    REAL(RKIND) ::  ynew(NIN), yinpp(NIN), sigmamin
    REAL(RKIND) :: PXEXP0,PXEXPDL, zz
    !
    info2 = -1
    !
    ! 1. Deals with various optional arguments
    !
    ! check that there are not 2 points too close  in rho: now 1e-5 from average dx
    zdx=(xin(nin)-xin(1)) / real(nin,rkind)
    nin_eff = nin
    allocate(xin_eff(nin_eff))
    allocate(yin_eff(nin_eff))
    i=1
    xin_eff(1)=xin(1)
    yin_eff(1)=yin(1)
    do j=2,NIN_eff
      if (abs(xin(j)-xin_eff(i)) .gt. 1e-5_rkind*zdx) then
        i=i+1;
        xin_eff(i) = xin(j)
        yin_eff(i) = yin(j)
      end if
    end do
    nin_eff = i;
    if (nin_eff .ne. nin) write(0,*) 'There were multiple input x points, changed nin=',nin, &
      & ' to nin_eff=',nin_eff,' in interpos_def to avoid this'
    !
    ! xout
    if (present(xout)) then
      if (present(nout)) then
        nout2=nout
      else
        ! print *,'% nout not present, use size(xout) instead'
        nout2 = size(xout)
      end if
      allocate(xout2(nout2))
      xout2 = xout(1:nout2)
    else
      nout2 = nin
!      if (present(nout)) nout = nin
      allocate(xout2(nout2))
      xout2 = xin_eff(1:nout2)
    end if
    allocate(yout2(nout2))
    allocate(youtp2(nout2))
    allocate(youtpp2(nout2))
    allocate(youtint2(nout2))
    ! sigma
    if (present(sigma)) then
      ! normalize to minimum value, so that taus=tension at this position
      sigmamin = minval(sigma)
      if (sigmamin .le. 0._rkind) then
        print *,' ERROR: min(sigma)=',sigmamin,' is .le. 0.'
        if (present(info)) info = 1
        return
      end if
      sigma2 = sigma/sigmamin
    else
      sigma2 = 1._rkind
    end if
    ! tension
    if (present(tension)) then
       if (tension .ge. 0._rkind) then
          sigma2=tension*sigma2
       else
          ! tension=-1 or negative, uses default value |tension| * mean(Delta_x)**3
         ! zz=minval(xin_eff(2:nin_eff)-xin_eff(1:nin_eff-1)) (minval not good if there is one very small dx like rho mesh from (R,Z))
          zz=sum(xin_eff(2:nin_eff)-xin_eff(1:nin_eff-1)) / real(nin_eff-1,rkind)
          ! tension=abs(tension)*zz**3 ! if intent(inout) but does not allow call: tension=-1.
          sigma2=abs(tension)*zz**3*sigma2
       end if
    else
       ! zz=minval(xin_eff(2:nin_eff)-xin_eff(1:nin_eff-1))
       zz=sum(xin_eff(2:nin_eff)-xin_eff(1:nin_eff-1)) / real(nin_eff-1,rkind)
       sigma2=zz**3*sigma2
    end if
    !
    ! boundary conditions (if nbc given, ybc should also be given)
    if (present(nbc)) then
      if (present(ybc)) then
        NBCLFT = nbc(1)
        NBCRGT = nbc(2)
        YBCLFT = ybc(1)
        YBCRGT = ybc(2)
        if (size(ybc) .GE. 4) then
          if (NBCLFT .ge. 10) XBCLFT = ybc(3)
          if (NBCRGT .ge. 10) XBCRGT = ybc(4)
          if (size(ybc) .GE. 6) then
            PXEXP0 = ybc(5)
            PXEXPDL= ybc(6)
          end if
        end if
      else
        print *,' problems in interpos, nbc is given but not ybc'
        if (present(info)) info = 2
        return
      end if
    else
      ! default is 2nd derivative=0 at each end
      ! For quadratic interpolation, default is NBCLFT=0 is used for quadratic defined with points and discontinuous 1st derivative
      NBCLFT = 0
      NBCRGT = 0
      YBCLFT = 0._rkind
      YBCRGT = 0._rkind
      ! not used but sets value to make sure:
      XBCLFT = xin_eff(1)
      XBCRGT = xin_eff(nin_eff)
      PXEXP0 = xin_eff(1)-xin_eff(nin_eff)
      PXEXPDL= 1._rkind
    end if
    ! compute y (0), also yp (1), also ypp (2) or also yint(3)
    ioptder = 0
    if (present(youtp)) ioptder=1
    if (present(youtpp)) ioptder=2
    if (present(youtint)) ioptder=3
    ! option: interpolation and extrapolation choice
    if (.NOT. present(option)) then
      inttype=3 ! cubic
      iextrapo = 32 !extrapolation cubic on 1st dx outside and then quadratic
    else
      optabs = abs(option)
      select case (optabs)
      case(1, 11, 21)
        ! linear interpolation
        inttype = 1
        ! extrapolation
        if (optabs .eq. 1)  iextrapo = 0 ! no extrapolation
        if (optabs .eq. 11) iextrapo = sign(1,option) ! linear extrapolation with edge-deriv for 11; edge-val for -11
        if (optabs .eq. 21) iextrapo = sign(10,option) ! constant outside: y=yedge for option=21, y=0. for -21
      case(2,12,22,32,42)
        ! quadratic interpolation
        inttype = 2
        ! extrapolation
        if (optabs .eq. 2)  iextrapo = 0 ! no extrapolation
        if (optabs .eq. 12) iextrapo = sign(21,option) ! quadratic on dx then linear extrapol. with edge-deriv for 12; edge-val for -12
        if (optabs .eq. 22) iextrapo = sign(1,option) ! linear extrapol. with edge-deriv for 22; edge-val for -22
        if (optabs .eq. 32) iextrapo = sign(2,option) ! quadratic extrapol. with edge-deriv for 32; edge-val for -32
        if (optabs .eq. 42) iextrapo = sign(10,option) ! constant outside: y=yedge for option=42, y=0. for -42
      case(3,13,23,33,43,53,63)
        ! cubic interpolation
        inttype = 3
        ! extrapolation
        if (optabs .eq. 3)  iextrapo = 0 ! no extrapolation
        if (optabs .eq. 13) iextrapo = sign(32,option) ! cubic on dx then quadratic extrapol. with edge-deriv for 13; edge-val for -13
        if (optabs .eq. 23) iextrapo = sign(1,option) ! linear extrapol. with edge-deriv for 22; edge-val for -22
        if (optabs .eq. 33) iextrapo = sign(2,option) ! quadratic extrapol. with edge-deriv for 33; edge-val for -33
        if (optabs .eq. 43) iextrapo = sign(3,option) ! cubic extrapol. with edge-deriv for 43; edge-val for -43
        if (optabs .eq. 53) iextrapo = sign(31,option) ! cubic on dx then linear extrapol. with edge-deriv for 53; edge-val for -53
        if (optabs .eq. 63) iextrapo = sign(10,option) ! constant outside: y=yedge for option=63, y=0. for -63
      case default
        inttype=mod(optabs,10)
        select case (inttype)
        case (1)
          iextrapo = 1
        case (2)
          iextrapo = 21
        case (3)
          iextrapo = 32
        case default
          print *,'option= ',option, &
            & ' is not valid. The last digit should be 1, 2 or 3 for linear, quadratic or cubic interpolation'
          if (present(info)) info = 3
          return
        end select
      end select
    end if
    !
    ! 2. call relevant interpolation routine
    !
    select case (inttype)
    case (1)
      ! linear
      call intlinear(XIN_EFF,YIN_EFF,NIN_EFF,xout2,yout2,youtp2,youtpp2,youtint2,NOUT2,ioptder,iextrapo)
    case (2)
      ! quadratic
      call intquadratic(XIN_EFF,YIN_EFF,NIN_EFF,xout2,yout2,youtp2,youtpp2,youtint2,NOUT2,ioptder,iextrapo,NBCLFT)
    case (3)
      ! cubic
      ! general cubic spline routine with non-periodic boundary conditions

      CALL CBSPLGEN(XIN_EFF,YIN_EFF,YNEW,YINPP,Nin_eff,XOUT2,YOUT2,YOUTP2,YOUTPP2,YOUTINT2, &
        &    nout2,ioptder,sigma2,NBCLFT, &
        &    NBCRGT,XBCLFT,XBCRGT,YBCLFT,YBCRGT,iextrapo,PXEXP0,PXEXPDL,info2)
      !
    case default
      print *,' inttype=',inttype,' ; error in this routine, this should not be possible'
      if (present(info)) info = 4
      return
    end select
    !
    ! 3. Fill in output values
    !
    ! yout, p, pp, int
    if (present(yout)) yout(1:nout2) = yout2
    if (present(youtp)) youtp(1:nout2) = youtp2
    if (present(youtpp)) youtpp(1:nout2) = youtpp2
    if (present(youtint)) youtint(1:nout2) = youtint2
    !
    if (present(info)) info = 0
    !
    deallocate(xout2)
    deallocate(yout2)
    deallocate(youtp2)
    deallocate(youtpp2)
    deallocate(youtint2)
    return
  END SUBROUTINE interpos_def

  SUBROUTINE interpos_defper(XIN,YIN,NIN,NOUT,tension,xout,yout,youtp,youtpp,youtint,NBC,YBC,sigma,info)
    ! for periodic boundary conditions, nbc=-1 and ybc=period should be given
    ! option not needed since no extrapolation and only valid for cubic spline here
    USE prec_rkind
    implicit none
    interface
       SUBROUTINE CBSPLGNP(PXIN,PYIN,PYINNEW,PYINPP,KNIN,PXOUT,PYOUT, &
            &  PYOUTP,PYOUTPP,PYOUTINT,KNOUT,KOPT,PSIG,PERIOD,PXEXP0,PXEXPDL,KFLAG)
         USE prec_rkind
         implicit none
         !
         REAL(RKIND) :: PXIN(KNIN), PYIN(KNIN) &
              &  ,PXOUT(KNOUT), PSIG(KNIN)
         REAL(RKIND) :: PYINNEW(KNIN), PYINPP(KNIN) &
              &  ,PYOUT(KNOUT), PYOUTP(KNOUT), PYOUTPP(KNOUT), PYOUTINT(KNOUT)
         !
         REAL(RKIND) :: PERIOD, PXEXP0, PXEXPDL
         INTEGER :: KOPT, KNIN, KNOUT
         INTEGER :: KFLAG
       END SUBROUTINE CBSPLGNP
    end interface
    INTEGER :: NIN, NBC
    INTEGER,optional :: NOUT
    REAL(RKIND) ::  XIN(NIN), YIN(NIN)
    REAL(RKIND), optional ::  tension
    REAL(RKIND), optional ::  xout(:)
    REAL(RKIND), optional ::  yout(:), youtp(:), youtpp(:), youtint(:)
    REAL(RKIND), optional ::   ybc, sigma(:)
    INTEGER, optional :: info
    ! need variables computed even if not given as argument
    REAL(RKIND), allocatable ::  xout2(:), yout2(:), youtp2(:), youtpp2(:), youtint2(:), xin_eff(:), yin_eff(:)
    REAL(RKIND) ::  sigma2(NIN), zdx, tension_eff
    integer :: info2, nout2, nin_eff, i, j
    !
    INTEGER ioptder
    REAL(RKIND) :: period
    REAL(RKIND) ::  ynew(NIN), yinpp(NIN), sigmamin
    REAL(RKIND) :: PXEXP0,PXEXPDL, zz
    !
    info2 = -1
    !
    ! 1. Deals with various optional arguments
    !
    ! period
    nin_eff = nin
    if (present(ybc)) then
      period = ybc
      ! check if extra point is given, if yes skip it with warning:
      if (abs(xin(nin_eff)-xin(1)-period) .lt. 1.0e-10_rkind*(xin(2)-xin(1))) then
        print *,'% it seems periodic point x(n)=x(1)+period is included in input, use nin-1 points'
        nin_eff = nin_eff-1
      end if
    else
      ! Assumes last x point corresponds to x(1)+period and is redundant
      period=xin(nin_eff)-xin(1)
      nin_eff = nin_eff-1
    end if
    !
    ! check that there are not 2 points too close  in rho: now 1e-5 from average dx
    zdx=period / real(nin_eff,rkind)
    allocate(xin_eff(nin_eff))
    allocate(yin_eff(nin_eff))
    i=1
    xin_eff(1)=xin(1)
    yin_eff(1)=yin(1)
    do j=2,NIN_eff
      if (abs(xin(j)-xin_eff(i)) .gt. 1e-5_rkind*zdx) then
        i=i+1;
        xin_eff(i) = xin(j)
        yin_eff(i) = yin(j)
      end if
    end do
    nin_eff = i;
    if (nin_eff .ne. nin) write(0,*) 'There were multiple input x points, changed nin=',nin, &
      & ' to nin_eff=',nin_eff,' in interpos_defper to avoid this'
    !
    ! xout
    if (present(xout)) then
      if (present(nout)) then
        nout2=nout
      else
        ! print *,'% nout not present, use size(xout) instead'
        nout2 = size(xout)
      end if
      allocate(xout2(nout2))
      xout2 = xout(1:nout2)
    else
      nout2 = nin_eff
!      if (present(nout)) nout = nin_eff
      allocate(xout2(nout2))
      xout2 = xin(1:nout2)
    end if
    allocate(yout2(nout2))
    allocate(youtp2(nout2))
    allocate(youtpp2(nout2))
    allocate(youtint2(nout2))
    ! sigma
    if (present(sigma)) then
      ! normalize to minimum value, so that taus=tension at this position
      sigmamin = minval(sigma)
      if (sigmamin .le. 0._rkind) then
        print *,' ERROR: min(sigma)=',sigmamin,' is .le. 0.'
        if (present(info)) info = 1
        return
      end if
      sigma2 = sigma/sigmamin
    else
      sigma2 = 1._rkind
    end if
    ! tension
    if (present(tension)) then
      ! For imposing periodic boundary conditions, one needs a finite tension. 
      ! If tension=0, issue a warning and use 1e-3*default as default
      tension_eff = tension
      if (tension_eff .eq. 0._rkind) then
        write(0,*) 'Warning, tension=0 with periodic interpos. Needs a finite tension, uses -1e-3'
        tension_eff = -0.001_RKIND
      end if
      if (tension_eff .gt. 0._rkind) then
        sigma2=tension_eff*sigma2
      else
        ! tension_eff=-1 or negative, uses default value |tension_eff| * Delta_x**3
        !zz=minval(xin_eff(2:nin_eff)-xin_eff(1:nin_eff-1)) (minval not good if there is a very small dx as rho mesh from (R,Z))
        zz=sum(xin_eff(2:nin_eff)-xin_eff(1:nin_eff-1)) / real(nin_eff-1,rkind)
        ! tension_eff=abs(tension_eff)*zz**3 ! if intent(inout) but does not allow call: tension_eff=-1.
        sigma2=abs(tension_eff)*zz**3*sigma2
      end if
    else
      ! zz=minval(xin_eff(2:nin_eff)-xin_eff(1:nin_eff-1))
      zz=sum(xin_eff(2:nin_eff)-xin_eff(1:nin_eff-1)) / real(nin_eff-1,rkind)
      sigma2=zz**3*sigma2
    end if
    !
    ! boundary conditions
    if (nbc .ne. -1) then
      print *,'nbc when scalar should be -1 for periodic boundary conditions, otherwise vector of size 2'
      if (present(info)) info = 2
      return
    end if

    ! compute y (0), also yp (1), also ypp (2) or also yint(3)
    ioptder = 0
    if (present(youtp)) ioptder=1
    if (present(youtpp)) ioptder=2
    if (present(youtint)) ioptder=3
    !
    ! 2. call relevant interpolation routine
    !
    ! cubic
    ! general cubic spline routine with periodic boundary conditions
    CALL CBSPLGNP(XIN_EFF,YIN_EFF,YNEW,YINPP,Nin_eff,XOUT2,YOUT2,YOUTP2,YOUTPP2,YOUTINT2, &
      &    nout2,ioptder,sigma2,period,PXEXP0,PXEXPDL,info2)
      !
    !
    ! 3. Fill in output values
    !
    ! yout, p, pp, int
    ! if xout = xin, could add extra periodic point if nin_eff was set to nin-1
    if (present(yout)) yout(1:nout2) = yout2
    if (present(youtp)) youtp(1:nout2) = youtp2
    if (present(youtpp)) youtpp(1:nout2) = youtpp2
    if (present(youtint)) youtint(1:nout2) = youtint2
    !
    if (present(info)) info = 0
    !
    deallocate(xin_eff)
    deallocate(yin_eff)
    deallocate(xout2)
    deallocate(yout2)
    deallocate(youtp2)
    deallocate(youtpp2)
    deallocate(youtint2)
    return
  END SUBROUTINE interpos_defper

  SUBROUTINE interpos_yinpp(XIN,YIN,YINNEW,YINPP,NIN,tension,nbc,ybc,sigma,info)
    !
    USE prec_rkind
    implicit none
    INTEGER :: NIN
    REAL(RKIND) ::  XIN(NIN), YIN(NIN)
    REAL(RKIND) ::   YINNEW(NIN), YINPP(NIN)
    REAL(RKIND), optional ::  tension
    REAL(RKIND), optional ::  ybc(:), sigma(:)
    INTEGER, optional :: nbc(2)
    INTEGER, optional :: info
    ! need variables computed even if not given as argument
    REAL(RKIND) ::  sigma2(NIN)
    integer :: info2
    !
    INTEGER NBCLFT,NBCRGT
    REAL(RKIND) :: XBCLFT,XBCRGT, YBCLFT,YBCRGT
    REAL(RKIND) ::  sigmamin
    REAL(RKIND) :: PXEXP0,PXEXPDL, zz
    !
    ! 1. Deals with various optional arguments
    !
    write(0,*) ' in interpos_yinpp'
    !
    ! sigma
    if (present(sigma)) then
      ! normalize to minimum value, so that taus=tension at this position
      sigmamin = minval(sigma)
      if (sigmamin .le. 0._rkind) then
        print *,' ERROR: min(sigma)=',sigmamin,' is .le. 0.'
        if (present(info)) info = 1
        return
      end if
      sigma2 = sigma/sigmamin
    else
      sigma2 = 1._rkind
    end if
    ! tension
    if (present(tension)) then
       if (tension .ge. 0._rkind) then
          sigma2=tension*sigma2
       else
          ! tension=-1 or negative, uses default value |tension| * Delta_x**3
         ! zz=minval(xin(2:nin)-xin(1:nin-1)) (minval not good if there is a very small dx as rho mesh from (R,Z))
          zz=sum(xin(2:nin)-xin(1:nin-1)) / real(nin-1,rkind)
          ! tension=abs(tension)*zz**3 ! if intent(inout) but does not allow call: tension=-1.
          sigma2=abs(tension)*zz**3*sigma2
       end if
    else
       ! zz=minval(xin(2:nin)-xin(1:nin-1))
       zz=sum(xin(2:nin)-xin(1:nin-1)) / real(nin-1,rkind)
       sigma2=zz**3*sigma2
    end if
    !
    ! boundary conditions (if nbc given, ybc should also be given)
    if (present(nbc)) then
      if (present(ybc)) then
        NBCLFT = nbc(1)
        NBCRGT = nbc(2)
        YBCLFT = ybc(1)
        YBCRGT = ybc(2)
        if (size(ybc) .GE. 4) then
          if (NBCLFT .ge. 10) XBCLFT = ybc(3)
          if (NBCRGT .ge. 10) XBCRGT = ybc(4)
          if (size(ybc) .GE. 6) then
            PXEXP0 = ybc(5)
            PXEXPDL= ybc(6)
          end if
        end if
      else
        print *,' problems in interpos, nbc is given but not ybc'
        if (present(info)) info = 2
        return
      end if
    else
      ! default is 2nd derivative=0 at each end
      ! For quadratic interpolation, default is NBCLFT=0 is used for quadratic defined with points and discontinuous 1st derivative
      NBCLFT = 0
      NBCRGT = 0
      YBCLFT = 0._rkind
      YBCRGT = 0._rkind
      ! not used but sets value to make sure:
      XBCLFT = xin(1)
      XBCRGT = xin(nin)
      PXEXP0 = xin(1)-xin(nin)
      PXEXPDL= 1._rkind
    end if
    !
    ! 2. call relevant interpolation routine
    !
    CALL CBFITBND(XIN,YIN,YINNEW,NIN,YINPP,SIGMA2,NBCLFT,NBCRGT, &
      &  XBCLFT,XBCRGT,YBCLFT,YBCRGT,PXEXP0,PXEXPDL)
    !
    if (present(info)) info = 0
    !
    return
  END SUBROUTINE interpos_yinpp

  SUBROUTINE interpos_interp(XIN,YIN,YPP,NIN,NOUT,xout,yout,youtp,youtpp,youtint,option,info)
    ! Note: calls arg ypp instead of yinpp otherwise there is ambiguity with interpos_yinpp
    USE prec_rkind
    implicit none
    INTEGER :: NIN
    INTEGER,optional :: NOUT
    REAL(RKIND) ::  XIN(NIN), YIN(NIN), YPP(NIN)
    REAL(RKIND), optional ::  xout(:)
    REAL(RKIND), optional ::  yout(:), youtp(:), youtpp(:), youtint(:)
    INTEGER, optional :: option
    INTEGER, optional :: info
    ! need variables computed even if not given as argument
    REAL(RKIND), allocatable ::  xout2(:), yout2(:), youtp2(:), youtpp2(:), youtint2(:)
    integer :: info2, nout2, optabs
    !
    INTEGER :: ioptder, iextrapo, inttype
    !
    ! 1. Deals with various optional arguments
    !
    ! xout
    if (present(xout)) then
      if (present(nout)) then
        nout2=nout
      else
        ! print *,'% nout not present, use size(xout) instead'
        nout2 = size(xout)
      end if
      allocate(xout2(nout2))
      xout2 = xout(1:nout2)
    else
      nout2 = nin
!      if (present(nout)) nout = nin
      allocate(xout2(nout2))
      xout2 = xin(1:nout2)
    end if
    allocate(yout2(nout2))
    allocate(youtp2(nout2))
    allocate(youtpp2(nout2))
    allocate(youtint2(nout2))
    ! compute y (0), also yp (1), also ypp (2) or also yint(3)
    ioptder = 0
    if (present(youtp)) ioptder=1
    if (present(youtpp)) ioptder=2
    if (present(youtint)) ioptder=3
    ! option: interpolation and extrapolation choice
    if (.NOT. present(option)) then
      iextrapo = 32 !extrapolation cubic on 1st dx outside and then quadratic
    else
      optabs = abs(option)
      select case (optabs)
      case(1, 11, 21)
        ! linear interpolation
        print *,' should call it only for cubic since uses ypp already calculated by cubic interpolation assumption'
        if (present(info)) info = 3
        return
      case(2,12,22,32,42)
        ! quadratic interpolation
        print *,' should call it only for cubic since uses ypp already calculated by cubic interpolation assumption'
        if (present(info)) info = 4
        return
      case(3,13,23,33,43,53,63)
        ! cubic interpolation
        inttype = 3
        ! extrapolation
        if (optabs .eq. 3)  iextrapo = 0 ! no extrapolation
        if (optabs .eq. 13) iextrapo = sign(32,option) ! cubic on dx then quadratic extrapol. with edge-deriv for 13; edge-val for -13
        if (optabs .eq. 23) iextrapo = sign(1,option) ! linear extrapol. with edge-deriv for 22; edge-val for -22
        if (optabs .eq. 33) iextrapo = sign(2,option) ! quadratic extrapol. with edge-deriv for 33; edge-val for -33
        if (optabs .eq. 43) iextrapo = sign(3,option) ! cubic extrapol. with edge-deriv for 43; edge-val for -43
        if (optabs .eq. 53) iextrapo = sign(31,option) ! cubic on dx then linear extrapol. with edge-deriv for 53; edge-val for -53
        if (optabs .eq. 63) iextrapo = sign(10,option) ! constant outside: y=yedge for option=63, y=0. for -63
      case default
        inttype=mod(optabs,10)
        select case (inttype)
        case (1)
          print *,' should call it only for cubic since uses ypp already calculated by cubic interpolation assumption'
          if (present(info)) info = 3
          return
        case (2)
          print *,' should call it only for cubic since uses ypp already calculated by cubic interpolation assumption'
          if (present(info)) info = 4
          return
        case (3)
          iextrapo = 32
        case default
          print *,'option= ',option, &
            & ' is not valid. The last digit should be 1, 2 or 3 for linear, quadratic or cubic interpolation'
          if (present(info)) info = 3
          return
        end select
      end select
    end if
    !
    ! 2. call relevant interpolation routine
    !
    ! cubic interpolation
    CALL SPLIBNDA(XIN,YIN,YPP,NIN,XOUT2,YOUT2,YOUTP2,YOUTPP2,YOUTINT2,NOUT2, &
      & iextrapo, ioptder)
    !
    ! 3. Fill in output values
    !
    ! yout, p, pp, int
    if (present(yout)) yout(1:nout2) = yout2
    if (present(youtp)) youtp(1:nout2) = youtp2
    if (present(youtpp)) youtpp(1:nout2) = youtpp2
    if (present(youtint)) youtint(1:nout2) = youtint2
    !
    if (present(info)) info = 0
    !
    deallocate(xout2)
    deallocate(yout2)
    deallocate(youtp2)
    deallocate(youtpp2)
    deallocate(youtint2)
    return
  END SUBROUTINE interpos_interp

  SUBROUTINE interpos_defxscal(X,Y,N,xscal,tension,yscal,yscalp,yscalpp,yscalint,nbcscal,ybcscal,sigma,option,info)
    !
    USE prec_rkind
    implicit none
    !
    INTERFACE
       SUBROUTINE INTLINEAR(PXIN,PYIN,KNIN,PXOUT,PYOUT,PYOUTP,PYOUTPP,PYOUTINT,KNOUT,KOPTDER,KEXTRAPO)
         USE PREC_RKIND
         IMPLICIT NONE
         ! arguments
         INTEGER :: KNIN, KNOUT, KOPTDER, KEXTRAPO
         REAL(RKIND) :: PXIN(KNIN), PYIN(KNIN), PXOUT(KNOUT)
         REAL(RKIND) :: PYOUT(KNOUT), PYOUTP(KNOUT), PYOUTPP(KNOUT), PYOUTINT(KNOUT)
       END SUBROUTINE INTLINEAR
       SUBROUTINE INTQUADRATIC(PXIN,PYIN,KNIN,PXOUT,PYOUT,PYOUTP,PYOUTPP,PYOUTINT,KNOUT,KOPTDER,KOPTXPOL,NBC)
         USE PREC_RKIND
         IMPLICIT NONE
         ! arguments
         INTEGER :: KNIN, KNOUT, KOPTDER, KOPTXPOL
         INTEGER, OPTIONAL ::  NBC
         REAL(RKIND) :: PXIN(KNIN), PYIN(KNIN), PXOUT(KNOUT)
         REAL(RKIND):: PYOUT(KNOUT), PYOUTP(KNOUT), PYOUTPP(KNOUT), PYOUTINT(KNOUT)
         !
       END SUBROUTINE INTQUADRATIC
       SUBROUTINE CBSPLGEN(PXIN,PYIN,PYINNEW,PYINPP,KNIN,PXOUT,PYOUT, &
            &  PYOUTP,PYOUTPP,PYOUTINT,KNOUT,KOPT,PSIG,NBCLFT &
            &  ,NBCRGT,XBCLFT,XBCRGT,YBCLFT,YBCRGT,KEXTRPOL,PXEXP0,PXEXPDL,KFLAG)
         USE prec_rkind
         implicit none
         REAL(RKIND) :: EPSILON
         PARAMETER(EPSILON = 1.0E-10_RKIND)
         ! arguments
         integer :: knin, knout, kopt, nbclft, nbcrgt, kextrpol
         integer:: kflag
         real(rkind) :: PXIN(KNIN), PYIN(KNIN), PXOUT(KNOUT), PSIG(KNIN), &
              & xbclft, xbcrgt, ybclft, ybcrgt, pxexp0, pxexpdl
         real(rkind):: PYINNEW(KNIN), PYINPP(KNIN), PYOUT(KNOUT), PYOUTP(KNOUT), &
              & PYOUTPP(KNOUT), PYOUTINT(KNOUT)
       END SUBROUTINE CBSPLGEN
    END INTERFACE
    !
    INTEGER :: N
    REAL(RKIND) ::  X(N), Y(N)
    REAL(RKIND), optional ::  tension
    REAL(RKIND) ::  xscal
    REAL(RKIND), optional ::  yscal, yscalp, yscalpp, yscalint
    REAL(RKIND), optional ::  ybcscal(:), sigma(:)
    INTEGER, optional :: nbcscal(:), option
    INTEGER, optional :: info
    ! need variables computed even if not given as argument
    REAL(RKIND) ::  xout2(1), yout2(1), youtp2(1), youtpp2(1), youtint2(1)
    REAL(RKIND) ::  sigma2(N), xin_eff(N), yin_eff(N)
    integer :: info2, nout2, i, j
    !
    INTEGER :: ioptder, NBCLFT,NBCRGT, iextrapo, inttype, optabs, nin_eff
    REAL(RKIND) :: XBCLFT,XBCRGT, YBCLFT,YBCRGT
    REAL(RKIND) ::  ynew(N), ypp(N), sigmamin
    REAL(RKIND) :: PXEXP0,PXEXPDL, zz, zdx
    !
    info2 = -1
    !
    ! 1. Deals with various optional arguments
    !
    ! check that there are not 2 points too close  in rho: now 1e-5 from average dx
    zdx=(x(n)-x(1)) / real(n,rkind)
    nin_eff = n
    i=1
    xin_eff(1)=x(1)
    yin_eff(1)=y(1)
    do j=2,NIN_eff
      if (abs(x(j)-xin_eff(i)) .gt. 1e-5_rkind*zdx) then
        i=i+1;
        xin_eff(i) = x(j)
        yin_eff(i) = y(j)
      end if
    end do
    nin_eff = i;
    if (nin_eff .ne. N) write(0,*) 'There were multiple input x points, changed N=',N, &
      & ' to nin_eff=',nin_eff,' in interpos_defxscal to avoid this'
    !
    ! xscal
    nout2 = 1
    xout2(1) = xscal
    ! sigma
    if (present(sigma)) then
      ! normalize to minimum value, so that tension=tension at this position
      sigmamin = minval(sigma)
      if (sigmamin .le. 0._rkind) then
        print *,' ERROR: min(sigma)=',sigmamin,' is .le. 0.'
        if (present(info)) info = 1
        return
      end if
      sigma2 = sigma/sigmamin
    else
      sigma2 = 1._rkind
    end if
    ! tension
    if (present(tension)) then
      if (tension .ge. 0._rkind) then
        sigma2=tension*sigma2
      else
        ! tension=-1 or negative, uses default value |tension| * Delta_x**3
        ! zz=minval(xin_eff(2:n)-xin_eff(1:nin_eff-1)) (minval not good if there is a very small dx as rho mesh from (R,Z))
        zz=sum(xin_eff(2:nin_eff)-xin_eff(1:nin_eff-1)) / real(nin_eff-1,rkind)
        ! tension=abs(tension)*zz**3 ! if intent(inout) but does not allow call: tension=-1.
        sigma2=abs(tension)*zz**3*sigma2
      end if
    else
      ! zz=minval(xin_eff(2:nin_eff)-xin_eff(1:nin_eff-1))
      zz=sum(xin_eff(2:nin_eff)-xin_eff(1:nin_eff-1)) / real(nin_eff-1,rkind)
      sigma2=zz**3*sigma2
    end if
    !
    ! boundary conditions (if nbcscal given, ybcscal should also be given)
    if (present(nbcscal)) then
      if (present(ybcscal)) then
        NBCLFT = nbcscal(1)
        NBCRGT = nbcscal(2)
        YBCLFT = ybcscal(1)
        YBCRGT = ybcscal(2)
        if (size(ybcscal) .GE. 4) then
          if (NBCLFT .ge. 10) XBCLFT = ybcscal(3)
          if (NBCRGT .ge. 10) XBCRGT = ybcscal(4)
          if (size(ybcscal) .GE. 6) then
            PXEXP0 = ybcscal(5)
            PXEXPDL= ybcscal(6)
          end if
        end if
      else
        print *,' problems in interpos, nbcscal is given but not ybcscal'
        if (present(info)) info = 2
        return
      end if
    else
      ! default is 2nd derivative=0 at each end
      ! For quadratic interpolation, default is NBCLFT=0 is used for quadratic defined with points and discontinuous 1st derivative
      NBCLFT = 0
      NBCRGT = 0
      YBCLFT = 0._rkind
      YBCRGT = 0._rkind
      ! not used but sets value to make sure:
      XBCLFT = xin_eff(1)
      XBCRGT = xin_eff(nin_eff)
      PXEXP0 = xin_eff(1)-xin_eff(nin_eff)
      PXEXPDL= 1._rkind
    end if
    ! compute y (0), also yp (1), also ypp (2) or also yint(3)
    ioptder = 0
    if (present(yscalp)) ioptder=1
    if (present(yscalpp)) ioptder=2
    if (present(yscalint)) ioptder=3
    ! option: interpolation and extrapolation choice
    if (.NOT. present(option)) then
      inttype=3 ! cubic
      iextrapo = 32 !extrapolation cubic on 1st dx outside and then quadratic
    else
      optabs = abs(option)
      select case (optabs)
      case(1, 11, 21)
        ! linear interpolation
        inttype = 1
        ! extrapolation
        if (optabs .eq. 1)  iextrapo = 0 ! no extrapolation
        if (optabs .eq. 11) iextrapo = sign(1,option) ! linear extrapolation with edge-deriv for 11; edge-val for -11
        if (optabs .eq. 21) iextrapo = sign(10,option) ! constant outside: y=yedge for option=21, y=0. for -21
      case(2,12,22,32,42)
        ! quadratic interpolation
        inttype = 2
        ! extrapolation
        if (optabs .eq. 2)  iextrapo = 0 ! no extrapolation
        if (optabs .eq. 12) iextrapo = sign(21,option) ! quadratic on dx then linear extrapol. with edge-deriv for 12; edge-val for -12
        if (optabs .eq. 22) iextrapo = sign(1,option) ! linear extrapol. with edge-deriv for 22; edge-val for -22
        if (optabs .eq. 32) iextrapo = sign(2,option) ! quadratic extrapol. with edge-deriv for 32; edge-val for -32
        if (optabs .eq. 42) iextrapo = sign(10,option) ! constant outside: y=yedge for option=42, y=0. for -42
      case(3,13,23,33,43,53,63)
        ! cubic interpolation
        inttype = 3
        ! extrapolation
        if (optabs .eq. 3)  iextrapo = 0 ! no extrapolation
        if (optabs .eq. 13) iextrapo = sign(32,option) ! cubic on dx then quadratic extrapol. with edge-deriv for 13; edge-val for -13
        if (optabs .eq. 23) iextrapo = sign(1,option) ! linear extrapol. with edge-deriv for 22; edge-val for -22
        if (optabs .eq. 33) iextrapo = sign(2,option) ! quadratic extrapol. with edge-deriv for 33; edge-val for -33
        if (optabs .eq. 43) iextrapo = sign(3,option) ! cubic extrapol. with edge-deriv for 43; edge-val for -43
        if (optabs .eq. 53) iextrapo = sign(31,option) ! cubic on dx then linear extrapol. with edge-deriv for 53; edge-val for -53
        if (optabs .eq. 63) iextrapo = sign(10,option) ! constant outside: y=yedge for option=63, y=0. for -63
      case default
        inttype=mod(optabs,10)
        select case (inttype)
        case (1)
          iextrapo = 1
        case (2)
          iextrapo = 21
        case (3)
          iextrapo = 32
        case default
          print *,'option= ',option, &
            & ' is not valid. The last digit should be 1, 2 or 3 for linear, quadratic or cubic interpolation'
          if (present(info)) info = 3
          return
        end select
      end select
    end if
    !
    ! 2. call relevant interpolation routine
    !
    select case (inttype)
    case (1)
      ! linear
      call intlinear(Xin_eff,Yin_eff,Nin_eff,xout2,yout2,youtp2,youtpp2,youtint2,NOUT2,ioptder,iextrapo)
    case (2)
      ! quadratic
      call intquadratic(Xin_eff,Yin_eff,Nin_eff,xout2,yout2,youtp2,youtpp2,youtint2,NOUT2,ioptder,iextrapo,NBCLFT)
    case (3)
      ! cubic
      ! general cubic spline routine with non-periodic boundary conditions
      CALL CBSPLGEN(Xin_eff,Yin_eff,YNEW,YPP,Nin_eff,XOUT2,YOUT2,YOUTP2,YOUTPP2,YOUTINT2, &
        &    nout2,ioptder,sigma2,NBCLFT, &
        &    NBCRGT,XBCLFT,XBCRGT,YBCLFT,YBCRGT,iextrapo,PXEXP0,PXEXPDL,info2)
      !
    case default
      print *,' inttype=',inttype,' ; error in this routine, this should not be possible'
      if (present(info)) info = 4
      return
    end select
    !
    ! 3. Fill in output values
    !
    ! yout, p, pp, int
    if (present(yscal)) yscal = yout2(1)
    if (present(yscalp)) yscalp = youtp2(1)
    if (present(yscalpp)) yscalpp = youtpp2(1)
    if (present(yscalint)) yscalint = youtint2(1)
    !
    if (present(info)) info = 0
    !
    return
  END SUBROUTINE interpos_defxscal

  SUBROUTINE interpos_defperxscal(X,Y,N,tension,xscal,yscal,yscalp,yscalpp,yscalint,NBCSCAL,YBCSCAL,sigma,info)
    ! for periodic boundary conditions, nbcscal=-1 and ybcscal=period should be given
    ! option not needed since no extrapolation and only valid for cubic spline here
    USE prec_rkind
    implicit none
    interface
       SUBROUTINE CBSPLGNP(PXIN,PYIN,PYINNEW,PYINPP,KNIN,PXOUT,PYOUT, &
            &  PYOUTP,PYOUTPP,PYOUTINT,KNOUT,KOPT,PSIG,PERIOD,PXEXP0,PXEXPDL,KFLAG)
         USE prec_rkind
         implicit none
         !
         REAL(RKIND) :: PXIN(KNIN), PYIN(KNIN) &
              &  ,PXOUT(KNOUT), PSIG(KNIN)
         REAL(RKIND) :: PYINNEW(KNIN), PYINPP(KNIN) &
              &  ,PYOUT(KNOUT), PYOUTP(KNOUT), PYOUTPP(KNOUT), PYOUTINT(KNOUT)
         !
         REAL(RKIND) :: PERIOD, PXEXP0, PXEXPDL
         INTEGER :: KOPT, KNIN, KNOUT
         INTEGER :: KFLAG
       END SUBROUTINE CBSPLGNP
    end interface
    INTEGER :: N, NBCSCAL
    REAL(RKIND) ::  X(N), Y(N)
    REAL(RKIND), optional ::  tension
    REAL(RKIND) ::  xscal
    REAL(RKIND), optional ::  yscal, yscalp, yscalpp, yscalint
    REAL(RKIND), optional ::   ybcscal, sigma(:)
    INTEGER, optional :: info
    ! need variables computed even if not given as argument
    REAL(RKIND) ::  xout2(1), yout2(1), youtp2(1), youtpp2(1), youtint2(1)
    REAL(RKIND) ::  sigma2(N), zdx, tension_eff
    REAL(RKIND), allocatable ::  xin_eff(:), yin_eff(:)
    integer :: info2, nout2, nin_eff, i, j
    !
    INTEGER ioptder
    REAL(RKIND) :: period
    REAL(RKIND) ::  ynew(N), ypp(N), sigmamin
    REAL(RKIND) :: PXEXP0,PXEXPDL, zz
    !
    info2 = -1
    !
    ! 1. Deals with various optional arguments
    !
    ! period
    nin_eff = n
    if (present(ybcscal)) then
      period = ybcscal
      ! check if extra point is given, if yes skip it with warning:
      if (abs(x(nin_eff)-x(1)-period) .lt. 1.0e-10_rkind*(x(2)-x(1))) then
        print *,'% it seems periodic point x(n)=x(1)+period is included in input, use n-1 points'
        nin_eff = nin_eff-1
      end if
    else
      ! Assumes last x point corresponds to x(1)+period and is redundant
      period=x(nin_eff)-x(1)
      nin_eff = nin_eff-1
    end if
    !
    ! check that there are not 2 points too close  in rho: now 1e-5 from average dx
    zdx=period / real(nin_eff,rkind)
    allocate(xin_eff(nin_eff))
    allocate(yin_eff(nin_eff))
    i=1
    xin_eff(1)=x(1)
    yin_eff(1)=y(1)
    do j=2,NIN_eff
      if (abs(x(j)-xin_eff(i)) .gt. 1e-5_rkind*zdx) then
        i=i+1;
        xin_eff(i) = x(j)
        yin_eff(i) = y(j)
      end if
    end do
    nin_eff = i;
    if (nin_eff .ne. N) write(0,*) 'There were multiple input x points, changed N=',N, &
      & ' to nin_eff=',nin_eff,' in interpos_defperxscal to avoid this'
    !
    ! xscal
    nout2 = 1
    xout2(1) = xscal
    ! sigma
    if (present(sigma)) then
      ! normalize to minimum value, so that tension=tension at this position
      sigmamin = minval(sigma)
      if (sigmamin .le. 0._rkind) then
        print *,' ERROR: min(sigma)=',sigmamin,' is .le. 0.'
        if (present(info)) info = 1
        return
      end if
      sigma2 = sigma/sigmamin
    else
      sigma2 = 1._rkind
    end if
    ! tension
    if (present(tension)) then
      ! For imposing periodic boundary conditions, one needs a finite tension. 
      ! If tension=0, issue a warning and use 1e-3*default as default
      tension_eff = tension
      if (tension_eff .eq. 0._rkind) then
        write(0,*) 'Warning, tension=0 with periodic interpos. Needs a finite tension, uses -1e-3'
        tension_eff = -0.001_RKIND
      end if
      if (tension_eff .ge. 0._rkind) then
        sigma2=tension_eff*sigma2
      else
        ! tension_eff=-1 or negative, uses default value |tension_eff| * Delta_x**3
        ! zz=minval(xin_eff(2:nin_eff)-xin_eff(1:nin_eff-1)) (minval not good if there is a very small dx as rho mesh from (R,Z)
        zz=sum(xin_eff(2:nin_eff)-xin_eff(1:nin_eff-1)) / real(nin_eff-1,rkind)
        ! tension_eff=abs(tension_eff)*zz**3 ! if intent(inout) but does not allow call: tension_eff=-1.
        sigma2=abs(tension_eff)*zz**3*sigma2
      end if
    else
      ! zz=minval(xin_eff(2:nin_eff)-xin_eff(1:nin_eff-1))
      zz=sum(xin_eff(2:nin_eff)-xin_eff(1:nin_eff-1)) / real(nin_eff-1,rkind)
      sigma2=zz**3*sigma2
    end if
    !
    ! boundary conditions
    if (nbcscal .ne. -1) then
      print *,'nbcscal when scalar should be -1 for periodic boundary conditions, otherwise vector of size 2'
      if (present(info)) info = 2
      return
    end if

    ! compute y (0), also yp (1), also ypp (2) or also yint(3)
    ioptder = 0
    if (present(yscalp)) ioptder=1
    if (present(yscalpp)) ioptder=2
    if (present(yscalint)) ioptder=3
    !
    ! 2. call relevant interpolation routine
    !
    ! cubic
    ! general cubic spline routine with periodic boundary conditions
    CALL CBSPLGNP(XIN_EFF,YIN_EFF,YNEW,YPP,Nin_eff,XOUT2,YOUT2,YOUTP2,YOUTPP2,YOUTINT2, &
      &    nout2,ioptder,sigma2,period,PXEXP0,PXEXPDL,info2)
      !
    !
    ! 3. Fill in output values
    !
    ! yscal, p, pp, int
    if (present(yscal)) yscal = yout2(1)
    if (present(yscalp)) yscalp = youtp2(1)
    if (present(yscalpp)) yscalpp = youtpp2(1)
    if (present(yscalint)) yscalint = youtint2(1)
    !
    if (present(info)) info = 0
    deallocate(xin_eff)
    deallocate(yin_eff)
    !
    return
  END SUBROUTINE interpos_defperxscal

  SUBROUTINE interpos_interpxscal(X,Y,YPP,N,xscal,yscal,yscalp,yscalpp,yscalint,option,info)
    ! Note: calls arg ypp instead of yinpp otherwise there is ambiguity with interpos_yinpp
    USE prec_rkind
    implicit none
    INTEGER :: N
    REAL(RKIND) ::  X(N), Y(N), YPP(N)
    REAL(RKIND) ::  XSCAL
    REAL(RKIND), optional ::  yscal, yscalp, yscalpp, yscalint
    INTEGER, optional :: option
    INTEGER, optional :: info
    ! need variables computed even if not given as argument
    REAL(RKIND) ::  xout2, yout2, youtp2, youtpp2, youtint2
    integer :: info2, nout2, optabs
    !
    INTEGER :: ioptder, iextrapo, inttype
    !
    ! 1. Deals with various optional arguments
    !
    ! xscal
    nout2 = 1
    xout2 = xscal
    ! compute y (0), also yp (1), also ypp (2) or also yint(3)
    ioptder = 0
    if (present(yscalp)) ioptder=1
    if (present(yscalpp)) ioptder=2
    if (present(yscalint)) ioptder=3
    ! option: interpolation and extrapolation choice
    if (.NOT. present(option)) then
      iextrapo = 32 !extrapolation cubic on 1st dx outside and then quadratic
    else
      optabs = abs(option)
      select case (optabs)
      case(1, 11, 21)
        ! linear interpolation
        print *,' should call it only for cubic since uses ypp already calculated by cubic interpolation assumption'
        if (present(info)) info = 3
        return
      case(2,12,22,32,42)
        ! quadratic interpolation
        print *,' should call it only for cubic since uses ypp already calculated by cubic interpolation assumption'
        if (present(info)) info = 4
        return
      case(3,13,23,33,43,53,63)
        ! cubic interpolation
        inttype = 3
        ! extrapolation
        if (optabs .eq. 3)  iextrapo = 0 ! no extrapolation
        if (optabs .eq. 13) iextrapo = sign(32,option) ! cubic on dx then quadratic extrapol. with edge-deriv for 13; edge-val for -13
        if (optabs .eq. 23) iextrapo = sign(1,option) ! linear extrapol. with edge-deriv for 22; edge-val for -22
        if (optabs .eq. 33) iextrapo = sign(2,option) ! quadratic extrapol. with edge-deriv for 33; edge-val for -33
        if (optabs .eq. 43) iextrapo = sign(3,option) ! cubic extrapol. with edge-deriv for 43; edge-val for -43
        if (optabs .eq. 53) iextrapo = sign(31,option) ! cubic on dx then linear extrapol. with edge-deriv for 53; edge-val for -53
        if (optabs .eq. 63) iextrapo = sign(10,option) ! constant outside: y=yedge for option=63, y=0. for -63
      case default
        inttype=mod(optabs,10)
        select case (inttype)
        case (1)
          print *,' should call it only for cubic since uses ypp already calculated by cubic interpolation assumption'
          if (present(info)) info = 3
          return
        case (2)
          print *,' should call it only for cubic since uses ypp already calculated by cubic interpolation assumption'
          if (present(info)) info = 4
          return
        case (3)
          iextrapo = 32
        case default
          print *,'option= ',option, &
            & ' is not valid. The last digit should be 1, 2 or 3 for linear, quadratic or cubic interpolation'
          if (present(info)) info = 3
          return
        end select
      end select
    end if
    !
    ! 2. call relevant interpolation routine
    !
    ! cubic interpolation
    CALL SPLIBNDA(X,Y,YPP,N,XOUT2,YOUT2,YOUTP2,YOUTPP2,YOUTINT2,NOUT2, &
      & iextrapo, ioptder)
    !
    ! 3. Fill in output values
    !
    ! yout, p, pp, int
    if (present(yscal)) yscal = yout2
    if (present(yscalp)) yscalp = youtp2
    if (present(yscalpp)) yscalpp = youtpp2
    if (present(yscalint)) yscalint = youtint2
    !
    if (present(info)) info = 0
    !
    return
  END SUBROUTINE interpos_interpxscal

end MODULE interpos_module
