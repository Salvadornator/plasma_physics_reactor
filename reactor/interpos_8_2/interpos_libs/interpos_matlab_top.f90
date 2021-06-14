subroutine mexfunction(nlhs, plhs, nrhs, prhs)
  !
  !   the mex function should have the name of the file. In this case interpos
  !   Then typically in matlab one calls:
  !   >> [a,b,c,..]=interpos(arg1,arg2,...)
  !
  !   nlhs, nrhs: number of input arguments on left and right-hand side respectively
  !   thus one can test and depending on number of arguments choose different options
  !
  !   plhs, prhs: pointer to the different arguments: plhs(i), i=1,nlhs
  !   (plhs and prhs are integer arrays)
  !
  !   fmex function: mxGetPr gets the pointer value in prhs, plhs
  !   .              mxCopyPtrToReal8: copy pointed values to local array
  !   .              mxCopyReal8ToPtr: copy local array values to pointed array
  !   .              mxGetM: get nb of rows in pointed array
  !   .              mxGetN: get nb of columns in pointed array
  !   .              mxCreateDoubleMatrix: creates a matrix for the return argument
  !
  !   In this example, matlab call expected:
  !
  !   >> [yout{,youtp,youtpp,youtint}] = interpos({kopt,}xin,yin{,xout}{,taus{,nbc,ybc}{,sigma}})
  !
  !   where {..} means facultative, thus one can have for example the following calls:
  !
  !   >> [yout]=interpos(xin,yin,yout)
  !   >> [yout]=interpos(kopt,xin,yin,xout)
  !   >> [yout]=interpos(xin,yin,taus)
  !   >> [yout,youtp]=interpos(xin,yin)
  !   >> [yout,youtp,youtpp]=interpos(xin,yin,xout,-1.)
  !   >> [yout,youtp,youtpp]=interpos(xin,yin,xout,taus,nbc,ybc)
  !   >> [yout,youtp,youtpp,youtint]=interpos(xin,yin,taus,nbc,ybc)
  !   >> [yout,youtp,youtpp,youtint]=interpos(xin,yin,xout,taus,nbc,ybc,sigma)
  !
  !   NOTE: If xout is a single point, taus should also be given otherwise it assumes it is a taus value:
  !          interpos(xin,yin,xxx) means taus=xxx if xxx size of 1
  !          interpos(xin,yin,xxx,0.) means xout=xxx, taus=0. if xxx size of 1
  !
  !   The arguments have the following meaning with respect to the function y(x)
  !   that one wants to interpolate/extrapolate or compute its derivatives and integrals:
  !
  !   kopt: option for interpolation and extrapolation method (default is kopt=13):
  !   .     = 1, 11 or 21 : linear interpolation
  !   .     = 2, 12 or 22 : quadratic "
  !   .     = 3, 13 or 23 : cubic spline interpolation, with tension=taus
  !   .     < 10 : send warning if need to extrapolate
  !   .     in [11,19] : extrapolate using a lower order than for interpolation
  !   .     in [21,29] : extrapolate using same order as interpolation
  !
  !   xin : array giving the input values of x
  !   yin : array giving the input values of y(x=xin(i))
  !
  !   xout: array of values of x at which the function and/or its derivatives
  !   .     are to be computed. If xout is not given, assumes xout=xin
  !   yout: interpolated values of y at x=xout(i).
  !   youtp: 1st derivative of y at x=xout(i)
  !   youtpp: 2nd derivative of y at x=xout(i)
  !   youtint: Integral of (y dx) from x=xin(1) to x=xout(i)
  !
  !   taus  : tension value for cubic spline interpolation. If not given, uses taus=0
  !           if taus < 0 (typically -1) uses a default taus value: taus=abs(taus)*default_taus
  !           (the default_taus is typically min(delta_x)**3)
  !
  !   nbc(2): [NBCLFT NBCRGT] (default: [0 0])
  !     BOUNDARY CONDITIONS, 4 TYPES DETERMINED BY THE VALUE OF (nbc(1)=NBCLFT
  !     and nbc(2) for left and right-hand side BC.):
  !
  !     0) VALUE OF SECOND DERIVATIVE AT XBCLFT OR RGT IS GIVEN (0 OR 10)
  !     1) VALUE OF 1ST        "       "   "     "  "   "   "   (1 OR 11)
  !     2) VALUE OF FUNCTION AT XBCLFT OR RGT IS GIVEN          (2 OR 12)
  !     if nbc(1)=-1, then assumes periodic boundary condition and uses ybc(1) for the period
  !
  !     The value of nbc(1 or 2) should be >= 10 if the BC is not at an end point:
  !     XBCLFT~=xin(1) or XBCRGT~=xin(end)
  !
  !     Examples:
  !       A good value for radial profiles in rho between [0,1] is to specify
  !       the first derivative = 0 at left and second derivative=0 at right:
  !       => nbc = [1 0] and ybc=[0. 0.]
  !
  !       or to obtained the first derivative at right-hand side from a
  !       lagrangian interpolation of the last points:
  !       => nbc=[1 1] and ybc=[0. 1e32]
  !
  !       Note: The B.C. can only be given anywhere but within the interval [xin(1),xin(end)]
  !       Therefore, say rho=[0.1 ... 1] and one wants to impose zero
  !       derivative at rho=0, one should add a point in the input before. One uses sigma to avoid forcing the value
  !       calling interpos:
  !             rho_eff(1)=0.;
  !             rho_eff(2:length(rho)+1) = rho;
  !       then  [yout]=interpos([0.;rho],[yin(1);yin],..,taus,[1 0],[0. 0.],[1000;ones(size(rho))]);
  !
  !   ybc(2 or 4 or 6) (default: [0 0]): [YBCLFT YBCRGT] or [YBCLFT YBCRGT XBCLFT XBCRGT] with
  !     THE VALUE IS GIVEN BY YBCLFT OR YBCRGT RESPECTIVELY.
  !
  !     FOR nbc type 1: IF (YBCLFT OR YBCRGT > 1E31 THEN DER. FROM LAGRANGIAN INTERP.
  !     FOR nbc type 1: IF (YBCLFT OR YBCRGT <-1E31 THEN DER. FROM LINEAR     INTERP.
  !
  !     IF NBCLFT OR NBCRGT IS < 10, PXIN(1) OR PXIN(KNIN) IS USED INSTEAD
  !     OF XBCLFT OR XBCRGT, RESPECTIVELY => XBCLFT OR XBCRGT NOT USED
  !
  !   ybc(6): [YBCLFT YBCRGT XBCLFT XBCRGT PXEXP0 PXEXPDL]
  !        Enables one to specify a gaussian wight to taus with PXEXP0 and PXEXPDL:
  !              taus_eff =  taus * EXP(-((XIN-PXEXP0)/PXEXPDL)**2)
  !        IF PXEXP0 NOT IN [PXIN(1),PXIN(KNIN)], EXP() IGNORED AND taus_eff=cst= taus
  !        if ybc(5:6) not given, then  PXEXP0=xin(1)-1. and pxexpdl=1. to get constant taus
  !
  !   ybc(1)=period for periodic boundary condition, if nbc(1)=-1. Defines the condition that y(x+period)=y(x)
  !          In this case xin(end) should not be equal to xin(1)+period, since it is redundant. However interpos just
  !          does not use the end point in this case.
  !        if period=-1, assumes period=xin(end)-xin(1). Useful if xin goes from 0 to 2pi for example
  !
  !    sigma: error_bar at each (x,y) point used for the fit. The effective taus value will then
  !           be taus(i)=taus .* sigma(i) ./ min(sigma). So uses the relative values of sigma.
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  USE prec_rkind
  USE interpos_module
  !
  implicit none
  !
  interface
     SUBROUTINE INTLINEAR(PXIN,PYIN,KNIN,PXOUT,PYOUT,PYOUTP,PYOUTPP,PYOUTINT,KNOUT,KOPTDER,KEXTRAPO)
       USE PREC_RKIND
       IMPLICIT NONE
       INTEGER :: KNIN, KNOUT, KOPTDER, KEXTRAPO
       REAL(RKIND) :: PXIN(KNIN), PYIN(KNIN), PXOUT(KNOUT)
       REAL(RKIND) :: PYOUT(KNOUT), PYOUTP(KNOUT), PYOUTPP(KNOUT), PYOUTINT(KNOUT)
     END SUBROUTINE INTLINEAR
     !
     SUBROUTINE INTQUADRATIC(PXIN,PYIN,KNIN,PXOUT,PYOUT,PYOUTP,PYOUTPP,PYOUTINT,KNOUT,KOPTDER,KOPTXPOL,NBC)
       USE PREC_RKIND
       IMPLICIT NONE
       REAL(RKIND) :: ZSIX, ZTHREE, ZTWO, ZONE
       PARAMETER(ZSIX=6._RKIND, ZTHREE=3._RKIND, ZTWO=2._RKIND, ZONE=1._RKIND)
       ! arguments
       INTEGER :: KNIN, KNOUT, KOPTDER, KOPTXPOL
       INTEGER, OPTIONAL ::  NBC
       REAL(RKIND) :: PXIN(KNIN), PYIN(KNIN), PXOUT(KNOUT)
       REAL(RKIND):: PYOUT(KNOUT), PYOUTP(KNOUT), PYOUTPP(KNOUT), PYOUTINT(KNOUT)
     end SUBROUTINE INTQUADRATIC

     logical function isnanos(a)
       use prec_rkind
       real(rkind) a
     end function isnanos
  end interface
  !
  integer(ITM_I8) plhs(*), prhs(*)
  integer nlhs, nrhs
  !
  integer(ITM_I8) mxGetPr
  !   pointers for matlab input arguments
  integer(ITM_I8) pkopt, pxin, pyin, pxout, ptaus, pnbc, pybc, psig
  REAL(RKIND),allocatable :: xin(:), yin(:), xout(:), sig(:)
  REAL(RKIND) :: ybc(6)
  integer kopt, nbc(2)
  REAL(RKIND) :: taus, zkopt, znbc(2)
  !   pointers for matlab output arguments
  integer(ITM_I8) pyout, pyoutp, pyoutpp, pyoutint
  REAL(RKIND),allocatable :: yout(:), youtp(:), youtpp(:), youtint(:)
  REAL(RKIND) :: XX, ZTAUEFF, ZCOFEXP, ZXEXP0, ZXEXPDL
  ! could use mxDestroyArray for rhs allocated related array
  integer(ITM_I8) mxgetm, mxgetn, mxCreateDoubleMatrix, mxDestroyArray
  integer(ITM_I8) malloc_f
  integer nbyt_tot, infomalloc
  common /memoryuse/ nbyt_tot
  REAL(RKIND) :: sigmin, fun_sigma
  integer kopt_sign, inttype, iextrapo, option, ninrow, nincol, nin, &
    &  nout, ixoutxin, i1len, i4len, i, inextrhs, idoexp_error, iybclen, &
    &  noutrow, noutcol, ioptder, iflag, i5len, i4isxout, i4arg
  integer nel
  integer(ITM_I8) ninint8, noutint8, nelint8, iybclenint8, noutrowint8, &
    &  noutcolint8, ninrowint8, nincolint8
  integer idum, nelint4, icheckNaNs, nin1koptbymistake
  !OS      integer :: inan, inan2, inan3
  !OS      real:: znan
  !OS      data inan/B'01111111100000100000000000000000'/ ! ok on lac g95
  !OS      data inan2/B'01111111100100010001001010101010'/  ! ok also (0 as first)
  !OS      data inan3/B'01111111100000000000000000000000'/ ! +Inf
  !OS  znan=transfer(inan,znan)
  !OS  print *,'znan=',znan
  !OS  znan=transfer(inan2,znan)
  !OS  print *,'znan2=',znan
  !OS  znan=transfer(inan3,znan)
  !OS  print *,'znan3=',znan
  !
  fun_sigma(XX)= ZTAUEFF*EXP(-ZCOFEXP*(XX-ZXEXP0)**2/ZXEXPDL**2)
  !
  !.......................................................................
  !
  !   1. defaults values
  !
  icheckNaNs = 0
  !
!!$  print *,' nlhs, plhs, nrhs, prhs= ',nlhs, plhs(1), nrhs, prhs(1)
  if (nrhs .eq. 0) then
    print *,'need at least 2 arrays for xin and yin as inputs'
    call flush(6)
    return
  end if
  !%OS      print *,'prhs= ',prhs(1),prhs(2),prhs(3)
  nbyt_tot=0
  taus = 0._RKIND
  !%OS      print *,' nbyt_tot= ',nbyt_tot
  !
  !   2. Input arguments
  !
  !%OS      print *,' nb of return arguments= ',nlhs
  !%OS      print *,' nb of input arguments= ',nrhs
  !
  !   2.1 kopt
  !
  inextrhs=1
  ninrowint8 = mxGetM(prhs(inextrhs))
  nincolint8 = mxGetN(prhs(inextrhs))
  i1len = int(max(ninrowint8,nincolint8))
  if ((ninrowint8.eq.0) .or. (nincolint8.eq.0)) then
    write(*,*) '1st input has length 0 => return'
    print *,'ninrowint8 = ',ninrowint8,' nincolint8 = ',nincolint8
    print *,'zero size'
    return
  end if
  !OS write(*,*) 'i1len= ',i1len
  if (i1len .eq. 1) then
    ! 1st input is probably kopt, except if nin=1
    pkopt=mxGetPr(prhs(inextrhs))
    inextrhs = inextrhs + 1
    nelint8 = int(1,8)
    call mxCopyPtrToReal8(pkopt, zkopt, nelint8)
    if (isnanos(zkopt)) then
      print *,'kopt is NaN => return'
      if (icheckNaNs .ne. 0) return
    end if
    ! print *,'abs(zkopt) - int(abs(zkopt))= ',abs(zkopt) - int(abs(zkopt))
    if (abs(abs(zkopt) - int(abs(zkopt))) <= 1.0E-12_RKIND) then
      ! Assume is kopt since integer
      kopt = int(abs(zkopt))
      kopt_sign=sign(1,int(zkopt))
      nin1koptbymistake=0
    else
      nin1koptbymistake=1
    end if
  else
    ! 1st input is xin, thus define default for kopt
    kopt = 13 ! 3 for cubic, large tenth for default extrapolation
    kopt_sign = 1
    inextrhs = 1
  end if
  !
  !%OS      print *,' kopt= ',kopt,zkopt
  !
  !   define interpolation type and extrapolation type
  !
  nel=10
  inttype = mod(kopt,10)
  iextrapo = kopt/nel
  !  print *,' kopt, iextrapo, inttype = ',kopt,iextrapo,inttype 
  !
  !   2.2 get xin, first length and allocate space to yin as well
  !
  if (inextrhs .gt. nrhs) then
    print *,'problem with number of inputs check help file'
    call flush(6)
    return
  end if
  pxin=mxGetPr(prhs(inextrhs))
  ninrowint8 = mxGetM(prhs(inextrhs))
  nincolint8 = mxGetN(prhs(inextrhs))
  !OS  print *,ninrowint8,nincolint8
  if ((ninrowint8 .eq. 0) .or. (nincolint8 .eq. 0)) then
    print *,'ninrowint8 = ',ninrowint8,' nincolint8 = ',nincolint8
    print *,'zero size for xin'
    return
  end if
  inextrhs = inextrhs + 1
  ninrow = int(ninrowint8)
  nincol = int(nincolint8)
  !%OS      print *,'nincol, ninrow= ',nincol, ninrow
  !   input 1D arrays in row or column
  nin = max(ninrow,nincol)
  !OS  print *,'nin= ',nin
  allocate(xin(nin))
  allocate(yin(nin))
  ninint8 = int(nin,8)
  call mxCopyPtrToReal8(pxin, xin, ninint8)
  !OS  print *,' nin= ',nin
  !OS  print *,' xin= ',(xin(i),i=1,nin)
  !OS  print *,'isnan=',isnan(xin(:))
  ! check xin for Nans
  if ((inttype.eq.3) .and. (nin .le. 3)) then
    if (nin .le. 2) then
      inttype = 1
      print *,'nin = ',nin,' <3 , cannot compute spline, uses linear interpolation'
    else
      inttype = 2
      print *,'nin = ',nin,' =3 , cannot compute spline, uses quadratic interpolation'
    end if
  end if
  idum = 0
  do i=1,nin
    !!$    idum=sum(transfer(isnanos(xin),idum,nin))
    if (isnanos(xin(i))) idum = idum + 1
  end do
  if (idum .gt. 0) then
    print *,' There are ',idum,' NaNs in xin'
    if (icheckNaNs .ne. 0) return
  end if
  !
  !   2.3 get yin
  !
  if (inextrhs .gt. nrhs) then
    print *,'problem with number of inputs check help file'
    call flush(6)
    return
  end if
  pyin=mxGetPr(prhs(inextrhs))     ! Get yin
  ninrowint8 = mxGetM(prhs(inextrhs))
  nincolint8 = mxGetN(prhs(inextrhs))
  if ((ninrowint8 .eq. 0) .or. (nincolint8 .eq. 0)) then
    print *,'ninrowint8 = ',ninrowint8,' nincolint8 = ',nincolint8
    print *,'zero size for yin'
    return
  end if
  inextrhs = inextrhs + 1
  call mxCopyPtrToReal8(pyin, yin, ninint8)
  !%OS      print *,' yin= ',(yin(i),i=1,nin)
  idum = 0
  do i=1,nin
    !!$  idum=sum(transfer(isnanos(yin),idum,nin))
    if (isnanos(yin(i))) idum = idum + 1
  end do
  if (idum .gt. 0) then
    print *,' There are ',idum,' NaNs in yin'
    if (icheckNaNs .ne. 0) return
  end if
  !
  !   2.4 check for 4th input
  !
  !%OS      print *,' nrhs= ',nrhs
  if (nrhs .lt. inextrhs) then
    nout = nin
    !   flag for xout=xin set to 1
    ixoutxin = 1
    !%OS        print *,' ixoutxin= ',ixoutxin
  else
    !
    ! check if next argument is an array, then it is xout.
    ! If it is a scalar and next one as well, then xout scalar and taus given, otherwise xout=xin and taus given except if
    ! linear or quadratic interpolation is called (inttype<3)
    !
    ninrowint8 = mxGetM(prhs(inextrhs))
    nincolint8 = mxGetN(prhs(inextrhs))
    if ((ninrowint8 .eq. 0) .or. (nincolint8 .eq. 0)) then
      print *,'ninrowint8 = ',ninrowint8,' nincolint8 = ',nincolint8
      print *,'zero size for xout or tension'
      return
    end if
    i4len = int(max(ninrowint8,nincolint8))
    i4arg = inextrhs
    inextrhs = inextrhs + 1
    !%OS      print *,'i4len= ',i4len
    i4isxout = 0
    if (i4len.gt.1) i4isxout = 1
    if (inttype .lt. 3) i4isxout = 1
    if (nrhs.ge.inextrhs) then
      ninrowint8 = mxGetM(prhs(inextrhs))
      nincolint8 = mxGetN(prhs(inextrhs))
      if ((ninrowint8 .eq. 0) .or. (nincolint8 .eq. 0)) then
        print *,'ninrowint8 = ',ninrowint8,' nincolint8 = ',nincolint8
        print *,'zero size for xout or tension'
        return
      end if
      i5len = int(max(ninrowint8,nincolint8))
      if (i4len.eq.1 .and. i5len.eq.1) then
        i4isxout = 1
      endif
    endif
    ! print *,'i4len, i5len, i4isxout= ',i4len, i5len, i4isxout
    if (i4isxout .eq. 0) then
      !   in this case prhs(4) gives the tau value  
      ixoutxin = 1
      ptaus = mxGetPr(prhs(i4arg))
      nelint8=int(1,8)
      call mxCopyPtrToReal8(ptaus, taus, nelint8)
      !%OS       print *,'salut11'
    else
      ixoutxin = 0
      nout = i4len
      allocate(xout(nout))
      pxout = mxGetPr(prhs(i4arg))
      noutint8 = int(nout,8)
      call mxCopyPtrToReal8(pxout, xout, noutint8)
    endif
  endif
  if (ixoutxin .eq. 1) then
    nout = nin
    noutint8 = int(nout,8)
    allocate(xout(nout))
    do i=1,nout
      xout(i) = xin(i)
    end do
  endif
  ! print *,'ixoutxin, nout= ',ixoutxin, nout
  !%OS      print *,' xout= ',(xout(i),i=1,nout)
  !
  !   2.5 5th argument is taus if xout .ne. xin
  !
  if (nrhs.ge.inextrhs .and. ixoutxin.eq.0) then
    ptaus = mxGetPr(prhs(inextrhs))
    ninrowint8 = mxGetM(prhs(inextrhs))
    nincolint8 = mxGetN(prhs(inextrhs))
    if ((ninrowint8 .eq. 0) .or. (nincolint8 .eq. 0)) then
      print *,'ninrowint8 = ',ninrowint8,' nincolint8 = ',nincolint8
      print *,'zero size for tension'
      return
    end if
    inextrhs = inextrhs + 1
    nelint8=int(1,8)
    call mxCopyPtrToReal8(ptaus, taus, nelint8)
  endif
  !%OS      print *,'taus= ',taus
  !
  !   2.6 5,6th or 6,7th arguments are nbc, ybc
  !
  nbc(1) = 0
  nbc(2) = 0
  ybc(1) = 0._RKIND
  ybc(2) = 0._RKIND
  idoexp_error=0
  if (nrhs .ge. inextrhs) then
    ninrowint8 = mxGetM(prhs(inextrhs))
    nincolint8 = mxGetN(prhs(inextrhs))
    if ((ninrowint8 .eq. 0) .or. (nincolint8 .eq. 0)) then
      print *,'ninrowint8 = ',ninrowint8,' nincolint8 = ',nincolint8
      print *,'zero size for nbc'
      return
    end if
    i4len = int(max(ninrowint8,nincolint8))
    !OS         print *,'i4len= ',i4len
    pnbc = mxGetPr(prhs(inextrhs))
    inextrhs = inextrhs + 1
    if (i4len .ne. 0) then
      if (i4len .eq. 1) then
        nelint8=int(1,8)
        call mxCopyPtrToReal8(pnbc, znbc, nelint8)
        nbc(1) = znbc(1)
        nbc(2) = znbc(1)
        ybc(1) = -1.
      else
        nelint8=int(2,8)
        call mxCopyPtrToReal8(pnbc, znbc, nelint8)
        nbc(1) = znbc(1)
        nbc(2) = znbc(2)
      end if
    end if
    !OS to fill in with NaNs when iextrapol=0 and xout outside interval, could use: (but may not work on all compiler)
    !OS         taus=-1._rkind
    !OS         taus=real(sqrt(taus),rkind)
    !
!!! can use this instead with data inan/B'01111111100100010001001010101010'/ or data inan/B'01111111100000100000000000000000'/
    !OS         znan=transfer(inan,znan)
    !OS         print *,'znan= ',znan,' inan= ',inan

    !OS         print *,'taus,nbc= ',taus,nbc,' znbc= ',znbc
    !
    if (nrhs .ge. inextrhs) then
      ninrowint8 = mxGetM(prhs(inextrhs))
      nincolint8 = mxGetN(prhs(inextrhs))
      if ((ninrowint8 .eq. 0) .or. (nincolint8 .eq. 0)) then
        print *,'ninrowint8 = ',ninrowint8,' nincolint8 = ',nincolint8
        print *,'zero size for ybc'
        return
      end if
      iybclenint8 = int(max(ninrowint8,nincolint8))
      !OS            print *,' iybclenint8= ',iybclenint8
      pybc = mxGetPr(prhs(inextrhs))
      inextrhs = inextrhs + 1
      if (iybclenint8 .eq. 0) then
        ybc(1) = 0._rkind
        ybc(2) = 0._rkind
        idoexp_error = 0
      else
        call mxCopyPtrToReal8(pybc, ybc, iybclenint8)
        iybclen=int(iybclenint8)
        if (iybclen .le. 4) then
          ybc(5) = xin(1) - xin(nin)
          ybc(6) = 1._RKIND
        else
          idoexp_error = 1
        endif
      end if
    end if
  else
    ybc(5) = xin(1) - xin(nin)
    ybc(6) = 1._RKIND
  endif
  !OS      print *,ybc
  !
  !   2.7 7th or 8th is sigmain
  allocate(sig(nin))
  if (nrhs .ge. inextrhs) then
    ninrowint8 = mxGetM(prhs(inextrhs))
    nincolint8 = mxGetN(prhs(inextrhs))
    if ((ninrowint8 .eq. 0) .or. (nincolint8 .eq. 0)) then
      print *,'ninrowint8 = ',ninrowint8,' nincolint8 = ',nincolint8
      print *,'zero size for sigma'
      return
    end if
    psig = mxGetPr(prhs(inextrhs))
    inextrhs = inextrhs + 1
    call mxCopyPtrToReal8(psig, sig, ninint8)
    !OS    sigmin=sig(1)
    !OS    do i=2,nin
    !OS      sigmin=min(sigmin,sig(i))
    !OS    end do
    !OS    do i=1,nin
    !OS      sig(i) = sig(i)/sigmin*abs(taus)
    !OS    end do
  else
    ! now multiplication by taus done in fortran interpos
    ZTAUEFF = 1._rkind
    ZXEXP0 = ybc(5)
    ZXEXPDL = ybc(6)
    ZCOFEXP = 0._RKIND
    if (idoexp_error.EQ.1) THEN
      ZCOFEXP = 1.0_RKIND
    endif
    do i=1,nin
      sig(i)=fun_sigma(xin(i))
    end do
  endif
  ! print *,'sig= ',sig
  !OS      if (taus.LT.0._RKIND) then
  !OS        sig(1)=10._RKIND*sig(1)
  !OS      endif
  !
  !   3. Allocate return arguments arrays (in same geometry as xin)
  !
  if (ninrow .eq. 1) then
    noutrow = 1
    noutcol = nout
  else
    noutrow = nout
    noutcol = 1
  endif
  !   pointers for yout, youtp and youtpp and youtint respectively
  nelint4=0
  noutrowint8 = int(noutrow,8)
  noutcolint8 = int(noutcol,8)
  !%OS      print *,' noutrowint8, noutcolint8= ',noutrowint8, noutcolint8
  !%OS      print *,' noutrow, noutcol= ',noutrow, noutcol
  plhs(1) = mxCreateDoubleMatrix(noutrowint8,noutcolint8,nelint4)
  !%OS      print *,plhs(1)
  pyout = mxGetPr(plhs(1))
  allocate(yout(nout))
  yout=0._RKIND
  ioptder = 0
  if (nlhs .ge. 2) then
    plhs(2) = mxCreateDoubleMatrix(noutrowint8,noutcolint8,nelint4)
    pyoutp = mxGetPr(plhs(2))
    allocate(youtp(nout))
    youtp=0._RKIND
    ioptder = 1
  endif
  if (nlhs .ge. 3) then
    plhs(3) = mxCreateDoubleMatrix(noutrowint8,noutcolint8,nelint4)
    pyoutpp = mxGetPr(plhs(3))
    allocate(youtpp(nout))
    youtpp=0._RKIND
    ioptder = 2
  endif
  if (nlhs .ge. 4) then
    plhs(4) = mxCreateDoubleMatrix(noutrowint8,noutcolint8,nelint4)
    pyoutint = mxGetPr(plhs(4))
    allocate(youtint(nout))
    youtint=0._RKIND
    ioptder = 3
  endif
  !
  ! At this stage we can decide results if only one point in xin
  !
  if (nin .eq. 1) then
    ! assume constant value, derivatives=0 and integral as well (from xin(1) to xin(1)) but return with message
    if (nin1koptbymistake == 0) then
      yout = yin
    else
      ! assume no kopt was given but with nin=1, xin was taken as kopt and yin as xin
      yout = xin
    end if
    print *,' Number of input points = ',nin,' Cannot decide what to do, not sure if kopt given,try to return with yout=yin=',yout
    call mxCopyReal8ToPtr(yout, pyout, noutint8)
    deallocate(yout)
    if (ioptder .ge. 1) then
      call mxCopyReal8ToPtr(youtp,pyoutp,noutint8)
      deallocate(youtp)
    endif
    if (ioptder .ge. 2) then
      call mxCopyReal8ToPtr(youtpp,pyoutpp,noutint8)
      deallocate(youtpp)
    endif
    if (ioptder .ge. 3) then
      call mxCopyReal8ToPtr(youtint,pyoutint,noutint8)
      deallocate(youtint)
    endif
    deallocate(xin)
    deallocate(yin)
    deallocate(sig)
    return
  end if
  !
  !   4. compute interpolated function and its derivatives
  !
  !  print *,' inttype, ioptder,iextrapo= ',inttype, ioptder,iextrapo
  if (inttype .eq. 1) then
    if (iextrapo .eq. 0) then
      iextrapo = 0
    else if (iextrapo .eq. 1) then
      iextrapo = 1
    else if (iextrapo .eq. 2) then
      iextrapo = 10
    else
      iextrapo = 1
    end if
    iextrapo = kopt_sign * iextrapo
    call intlinear(xin,yin,nin,xout,yout,youtp,youtpp,youtint,nout,ioptder,iextrapo)
  end if
  !
  if (inttype .eq. 2) then
    if (iextrapo .eq. 0) then
      iextrapo = 0
    else if (iextrapo .eq. 1) then
      iextrapo = 21
    else if (iextrapo .eq. 2) then
      iextrapo = 1
    else if (iextrapo .eq. 3) then
      iextrapo = 2
    else if (iextrapo .eq. 4) then
      iextrapo = 10
    else
      iextrapo = 2
    end if
    iextrapo = kopt_sign * iextrapo
    call intquadratic(xin,yin,nin,xout,yout,youtp,youtpp,youtint,nout,ioptder,iextrapo,nbc(1))
  end if
  !
  if (inttype .eq. 3) then
    if (iextrapo .eq. 0) then
      iextrapo = 0
    else if (iextrapo .eq. 1) then
      iextrapo = 32
    else if (iextrapo .eq. 2) then
      iextrapo = 1
    else if (iextrapo .eq. 3) then
      iextrapo = 2
    else if (iextrapo .eq. 4) then
      iextrapo = 3
    else if (iextrapo .eq. 5) then
      iextrapo = 31
    else if (iextrapo .eq. 6) then
      iextrapo = 10
    else
      iextrapo = 32
    endif
    !        print *,'iextrapo= ',iextrapo
    !
    !   add sign of kopt for extrapol
    !%OS        print *,' ioptder,iextrapo,inttype= ',ioptder,iextrapo,inttype
    !%OS        print *,'xin,yin,nin,xout',xin,yin,nin,xout
    if (nbc(1).lt. 0) then
      ! periodic boundary conditions
      if (ioptder .eq. 0) then
        call interpos(xin,yin,nin,nout,taus,xout,yout &
          &    ,nbc=nbc(1),ybc=ybc(1),sigma=sig,info=iflag)
      else if (ioptder .eq. 1) then
        call interpos(xin,yin,nin,nout,taus,xout,yout,youtp &
          &    ,nbc=nbc(1),ybc=ybc(1),sigma=sig,info=iflag)
      else if (ioptder .eq. 2) then
        call interpos(xin,yin,nin,nout,taus,xout,yout,youtp,youtpp &
          &    ,nbc=nbc(1),ybc=ybc(1),sigma=sig,info=iflag)
      else
        call interpos(xin,yin,nin,nout,taus,xout,yout,youtp,youtpp,youtint &
          &    ,nbc=nbc(1),ybc=ybc(1),sigma=sig,info=iflag)
      end if
    else
      if (ioptder .eq. 0) then
        call interpos(xin,yin,nin,nout,taus,xout,yout &
          &    ,nbc=nbc,ybc=ybc,sigma=sig,option=kopt*kopt_sign,info=iflag)
      else if (ioptder .eq. 1) then
      call interpos(xin,yin,nin,nout,taus,xout,yout,youtp &
          &    ,nbc=nbc,ybc=ybc,sigma=sig,option=kopt*kopt_sign,info=iflag)
      else if (ioptder .eq. 2) then
      call interpos(xin,yin,nin,nout,taus,xout,yout,youtp,youtpp &
          &    ,nbc=nbc,ybc=ybc,sigma=sig,option=kopt*kopt_sign,info=iflag)
      else
      call interpos(xin,yin,nin,nout,taus,xout,yout,youtp,youtpp,youtint &
          &    ,nbc=nbc,ybc=ybc,sigma=sig,option=kopt*kopt_sign,info=iflag)
      end if
    endif
    if (iflag .ne. 0) then
      print *
      print *,' problem in cbsplgen'
      print *,' iflag = ',iflag
      return
    end if
  end if
  !   
  !   5. copy arrays to return value pointers
  !

  !%OS      print *,'yout(1:10)',yout(1:10)
  !%OS      print *,'pyout= ',pyout
  !%OS      print *,'noutint8= ',noutint8
  call mxCopyReal8ToPtr(yout, pyout, noutint8)
  deallocate(yout)
  if (ioptder .ge. 1) then
    call mxCopyReal8ToPtr(youtp,pyoutp,noutint8)
    deallocate(youtp)
  endif
  if (ioptder .ge. 2) then
    call mxCopyReal8ToPtr(youtpp,pyoutpp,noutint8)
    deallocate(youtpp)
  endif
  if (ioptder .ge. 3) then
    call mxCopyReal8ToPtr(youtint,pyoutint,noutint8)
    deallocate(youtint)
  endif
  !
  deallocate(xin)
  deallocate(yin)
  deallocate(sig)
  !
  return
end subroutine mexfunction

logical function isnanos(a)
  use prec_rkind
  real(rkind) a
  if (a.ne.a) then
    isnanos = .true.
  else
    isnanos = .false.
  end if
  return
end function isnanos

