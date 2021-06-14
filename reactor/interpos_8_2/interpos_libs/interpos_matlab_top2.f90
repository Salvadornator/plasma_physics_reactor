subroutine mexfunction(nlhs, plhs, nrhs, prhs)
  USE prec_rkind
  USE interpos_module
  USE cbsplgenp0mod
  implicit none
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
  !   .              mxCreateFull: creates a matrix for the return argument
  !
  !   matlab call expected:
  !
  !   >> [yout{,youtp{,youtpp{,youtint}}}] = interpos(kopt,xin,yin{,xout}{,taus}{,nbc,ybc{,sigmain}})
  !
  !   where {..} means facultative, thus one can have for example the following calls:
  !
  !   >> [yout]=interpos(kopt,xin,yin)
  !   >> [yout]=interpos(kopt,xin,yin,xout)
  !   >> [yout]=interpos(kopt,xin,yin,taus)
  !   >> [yout,youtp]=interpos(kopt,xin,yin)
  !   >> [yout,youtp,youtpp]=interpos(kopt,xin,yin,xout,taus)
  !   >> [yout,youtp,youtpp]=interpos(kopt,xin,yin,xout,taus,nbc,ybc)
  !   >> [yout,youtp,youtpp]=interpos(kopt,xin,yin,taus,nbc,ybc)
  !   >> [yout,youtp,youtpp,youtint]=interpos(kopt,xin,yin,taus,nbc,ybc,sigmain)
  !
  !   NOTE: It cannot discrimenate between xout and taus if xout is a single point => xout SHOULD be an array
  !   of at least 2 dimensions, except if all 8 arguments are given. OR give both "xout, taus" as 2 scalars. 
  !   Since nbc should be an array, it will detect that both xout(1) and taus have been given. 
  !   Even for taus=0, put xout,0 in inputs.
  !
  !   where the arguments have the following meaning with respect to the function y(x)
  !   that one wants to interpolate or compute its derivatives:
  !
  !   kopt: option for interpolation and extrapolation method:
  !   .     = 1, 11 or 21 : linear interpolation
  !   .     = 2, 12 or 22 : quadratic "
  !   .     = 3, 13 or 23 : cubic spline interpolation, with tension=taus
  !   .     < 10 : send warning if need to extrapolate
  !   .     in [11,19] : extrapolate using same order as interpolation
  !   .     in [21,29] : extrapolate using a lower order than for interpolation
  !   xin : array giving the input values of x
  !   yin : array giving the input values of y(x=xin(i))
  !   xout: array of values of x at which the function and/or its derivatives
  !   .     are to be computed. If xout is not given, assumes xout=xin
  !   yout: interpolated values of y at x=xout(i).
  !   youtp: 1st derivative of y at x=xout(i)
  !   youtpp: 2nd derivative of y at x=xout(i)
  !   youtint: Integral of y at x=xout(i). Formally, primitive of cubic spline fit at x=xout
  !   taus  : tension value for cubic spline interpolation. If not given, uses taus=0
  !   nbc(2): [NBCLFT NBCRGT] (default: [0 0])
  !   yxbc(2 or 4) (default: [0 0]): [YBCLFT YBCRGT] or [YBCLFT YBCRGT XBCLFT XBCRGT] with
  !     BOUNDARY CONDITIONS, 3 TYPES DETERMINED BY THE VALUE OF (NBCLFT/RGT):
  !
  !     0) VALUE OF SECOND DERIVATIVE AT XBCLFT OR RGT IS GIVEN (0 OR 10)
  !     1) VALUE OF 1ST        "       "   "     "  "   "   "   (1 OR 11)
  !     2) VALUE OF FUNCTION AT XBCLFT OR RGT IS GIVEN          (2 OR 12)
  !
  !     THE VALUE IS GIVEN BY YBCLFT OR YBCRGT RESPECTIVELY.
  !
  !     FOR TYPE 1: IF (YBCLFT OR YBCRGT > 1E31 THEN DER. FROM LAGRANGIAN INTERP.
  !     FOR TYPE 1: IF (YBCLFT OR YBCRGT <-1E31 THEN DER. FROM LINEAR     INTERP.
  !
  !     IF NBCLFT OR NBCRGT IS < 10, PXIN(1) OR PXIN(KNIN) IS USED INSTEAD
  !     OF XBCLFT OR XBCRGT, RESPECTIVELY => XBCLFT OR XBCRGT NOT USED
  !
  !   ybc(5:6): [YBCLFT YBCRGT XBCLFT XBCRGT PXEXP0 PXEXPDL]
  !        Enables one to specify a gaussian wight to taus with PXEXP0 and PXEXPDL:
  !              taus_eff =  taus * EXP(-((XIN-PXEXP0)/PXEXPDL)**2)
  !        IF PXEXP0 NOT IN [PXIN(1),PXIN(KNIN)], EXP() IGNORED AND taus_eff=cst= taus
  !        if ybc(5:6) not given, then  PXEXP0=xin(1)-1. and pxexpdl=1. to get constant taus
  !
  !   sigmain: array of error bar at each input (xin,yin) data point (assumed=cst=1 if not given)
  !            it is normalized to minimum value so that tau*sigma(i) = tau at minimum
  !
  !
  !
  integer*8 plhs(*), prhs(*)
  integer nlhs, nrhs
  !
  integer*8 mxGetPr
  !   pointers for matlab input arguments
  integer*8 pkopt, pxin, pyin, pxout, ptaus, pnbc, pybc, psig
  REAL(RKIND),pointer :: xin(:), yin(:), xout(:), sig(:)
  REAL(RKIND) :: ybc(6)
  integer kopt, nbc(2)
  REAL(RKIND) :: taus, zkopt, znbc(2)
  !   pointers for matlab output arguments
  integer*8 pyout, pyoutp, pyoutpp, pyoutint
  REAL(RKIND),pointer :: yout(:), youtp(:), youtpp(:), youtint(:)
  REAL(RKIND) :: XX, ZTAUEFF, ZCOFEXP, ZXEXP0, ZXEXPDL
  integer*8 mxgetm, mxgetn, mxcreatefull
  integer*8 malloc_f
  integer nbyt_tot, infomalloc
  common /memoryuse/ nbyt_tot
  REAL(RKIND) :: sigmin, fun_sigma
  integer kopt_sign, inttype, iextrapo, ninrow, nincol, nin, &
    &  nout, ixoutxin, i1len, i4len, i, inextrhs, idoexp_error, iybclen, &
    &  noutrow, noutcol, ioptder, iflag, i5len, i4isxout, i4arg
  integer nel
  integer*8 ninint8, noutint8, nelint8, iybclenint8, noutrowint8, &
    &  noutcolint8, ninrowint8, nincolint8
!!$      integer :: inan
!!$      real :: znan
!!$      data inan/B'01111111100000100000000000000000'/ ! ok on lac g95
!!$!      data inan/B'01111111100100010001001010101010'/  ! ok also (0 as first)
!!$!       data inan/B'01111111100000000000000000000000'/ ! +Inf
  !
  fun_sigma(XX)= ZTAUEFF*EXP(-ZCOFEXP*(XX-ZXEXP0)**2/ZXEXPDL**2)
  !
  !.......................................................................
  !
  !   1. defaults values
  !
  !%OS      print *,' nrhs= ',nrhs
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
  i1len = int(max(mxGetM(prhs(inextrhs)),mxGetN(prhs(inextrhs))))
  if (i1len .eq. 1) then
    ! 1st input is kopt
    pkopt=mxGetPr(prhs(inextrhs))
    inextrhs = inextrhs + 1
    nelint8 = int(1,8)
    call mxCopyPtrToReal8(pkopt, zkopt, nelint8)
    kopt = int(abs(zkopt))
    kopt_sign=sign(1,int(zkopt))
  else
    ! 1st input is xin, thus define default for kopt
    kopt = 33 ! 3 for cubic, large tenth for default extrapolation
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
  print *,' kopt, iextrapo, inttype = ',kopt,iextrapo,inttype 
  !
  !   2.2 get xin, first length and allocate space to yin as well
  !
  pxin=mxGetPr(prhs(inextrhs))
  ninrowint8 = mxGetM(prhs(inextrhs))
  nincolint8 = mxGetN(prhs(inextrhs))
  inextrhs = inextrhs + 1
  ninrow = int(ninrowint8)
  nincol = int(nincolint8)
  !%OS      print *,'nincol, ninrow= ',nincol, ninrow
  !   input 1D arrays in row or column
  nin = max(ninrow,nincol)
  allocate(xin(nin))
  allocate(yin(nin))
  ninint8 = int(nin,8)
  call mxCopyPtrToReal8(pxin, xin, ninint8)
  !%OS      print *,xin(1:10)
  !%OS      print *,' nin= ',nin
  !%OS      print *,' xin= ',(xin(i),i=1,nin)
  !
  !   2.3 get yin
  !
  pyin=mxGetPr(prhs(inextrhs))     ! Get yin
  inextrhs = inextrhs + 1
  call mxCopyPtrToReal8(pyin, yin, ninint8)
  !%OS      print *,' yin= ',(yin(i),i=1,nin)
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
    !   check if next argument is an array.
    !   Assume value is taus only if length is one (scalar), cubic spline interpolation
    !   option was asked, and 8 input arguments are given. Also check if 4th and 5th argument are scalars.
    !   Otherwise the 4th argument is assumed to be xout
    !
    i4len = int(max(mxGetM(prhs(inextrhs)),mxGetN(prhs(inextrhs))))
    i4arg = inextrhs
    inextrhs = inextrhs + 1
    !%OS      print *,'i4len= ',i4len
    i4isxout = 0
    if (nrhs.ge.inextrhs) then
      i5len = int(max(mxGetM(prhs(inextrhs)),mxGetN(prhs(inextrhs))))
      if (i4len.eq.1 .and. i5len.eq.1) then
        i4isxout = 1
      endif
    endif
    if (i4len.eq.1 .and. nrhs.ne.(i4arg+4) .and. inttype.eq.3 .and.  &
      &      i4isxout.eq.0) then
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
  !%OS      print *,'ixoutxin, nout= ',ixoutxin, nout
  !%OS      print *,' xout= ',(xout(i),i=1,nout)
  !
  !   2.5 5th argument is taus if xout .ne. xin
  !
  if (nrhs.ge.inextrhs .and. ixoutxin.eq.0) then
    ptaus = mxGetPr(prhs(inextrhs))
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
    i4len = int(max(mxGetM(prhs(inextrhs)),mxGetN(prhs(inextrhs))))
!!$         print *,'i4len= ',i4len
    pnbc = mxGetPr(prhs(inextrhs))
    inextrhs = inextrhs + 1
    if (i4len .ne. 0) then
      if (i4len .eq. 1) then
        nelint8=int(1,8)
        call mxCopyPtrToReal8(pnbc, znbc, nelint8)
        nbc(1) = znbc(1)
        nbc(2) = znbc(1)
      else
        nelint8=int(2,8)
        call mxCopyPtrToReal8(pnbc, znbc, nelint8)
        nbc(1) = znbc(1)
        nbc(2) = znbc(2)
      end if
    end if
!!$ ! to fill in with NaNs when iextrapol=0 and xout outside interval, could use: (but may not work on all compiler)
!!$         taus=-1._rkind
!!$         taus=real(sqrt(taus),rkind)
    !
!!! can use this instead with data inan/B'01111111100100010001001010101010'/ or data inan/B'01111111100000100000000000000000'/
!!$         znan=transfer(inan,znan)
!!$         print *,'znan= ',znan,' inan= ',inan

!!$         print *,'taus,nbc= ',taus,nbc,' znbc= ',znbc
    !
    if (nrhs .ge. inextrhs) then
      iybclenint8 = max(mxGetM(prhs(inextrhs)),mxGetN(prhs(inextrhs)))
!!$            print *,' iybclenint8= ',iybclenint8
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
!!$      print *,ybc
  !
  !   2.7 7th or 8th is sigmain
  allocate(sig(nin))
  if (nrhs .ge. inextrhs) then
    psig = mxGetPr(prhs(inextrhs))
    inextrhs = inextrhs + 1
    call mxCopyPtrToReal8(psig, sig, ninint8)
    sigmin=sig(1)
    do i=2,nin
      sigmin=min(sigmin,sig(i))
    end do
    do i=1,nin
      sig(i) = sig(i)/sigmin*abs(taus)
    end do
  else
    ZTAUEFF = taus
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
!!$      if (taus.LT.0._RKIND) then
!!$        sig(1)=10._RKIND*sig(1)
!!$      endif
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
  nelint8=int(0,8)
  noutrowint8 = int(noutrow,8)
  noutcolint8 = int(noutcol,8)
  !%OS      print *,' noutrowint8, noutcolint8= ',noutrowint8, noutcolint8
  !%OS      print *,' noutrow, noutcol= ',noutrow, noutcol
  plhs(1) = mxCreateFull(noutrowint8,noutcolint8,nelint8)
  !%OS      print *,plhs(1)
  pyout = mxGetPr(plhs(1))
  allocate(yout(nout))
  yout=0._RKIND
  ioptder = 0
  if (nlhs .ge. 2) then
    nelint8=int(0,8)
    plhs(2) = mxCreateFull(noutrowint8,noutcolint8,nelint8)
    pyoutp = mxGetPr(plhs(2))
    allocate(youtp(nout))
    youtp=0._RKIND
    ioptder = 1
  endif
  if (nlhs .ge. 3) then
    nelint8=int(0,8)
    plhs(3) = mxCreateFull(noutrowint8,noutcolint8,nelint8)
    pyoutpp = mxGetPr(plhs(3))
    allocate(youtpp(nout))
    youtpp=0._RKIND
    ioptder = 2
  endif
  if (nlhs .ge. 4) then
    nelint8=int(0,8)
    plhs(4) = mxCreateFull(noutrowint8,noutcolint8,nelint8)
    pyoutint = mxGetPr(plhs(4))
    allocate(youtint(nout))
    youtint=0._RKIND
    ioptder = 3
  endif
  !
  !   4. compute interpolated function and its derivatives
  !
  print *,' inttype, ioptder,iextrapo= ',inttype, ioptder,iextrapo
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
      iextrapo = 1
    else if (iextrapo .eq. 2) then
      iextrapo = 2
    else if (iextrapo .eq. 3) then
      iextrapo = 21
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
      iextrapo = 3
    else if (iextrapo .eq. 2) then
      iextrapo = 31
    else if (iextrapo .eq. 3) then
      iextrapo = 32
    else if (iextrapo .eq. 4) then
      iextrapo = 1
    else if (iextrapo .eq. 5) then
      iextrapo = 2
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
      call cbsplgenp0(xin,yin,nin,xout,yout,youtp,youtpp,youtint,nout, &
        &    ioptder,sig,nbc,ybc,iflag)
    else
      call interpos(xin,yin,nin,xout,yout,youtp,youtpp,youtint,nout, &
        &    ioptder,kopt_sign*iextrapo,sig,nbc,ybc,iflag)
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
  if (ioptder .ge. 1) then
    call mxCopyReal8ToPtr(youtp,pyoutp,noutint8)
  endif
  if (ioptder .ge. 2) then
    call mxCopyReal8ToPtr(youtpp,pyoutpp,noutint8)
  endif
  if (ioptder .ge. 3) then
    call mxCopyReal8ToPtr(youtint,pyoutint,noutint8)
  endif
  !
  deallocate(xin)
  deallocate(yin)
  deallocate(sig)
  !
  return
end subroutine mexfunction
