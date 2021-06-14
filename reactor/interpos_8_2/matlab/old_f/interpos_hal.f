      subroutine mexfunction(nlhs, plhs, nrhs, prhs)
      implicit none
c   the mex function should have the name of the file. In this case interpos
c   Then typically in matlab one calls:
c   >> [a,b,c,..]=interpos(arg1,arg2,...)
c
c   nlhs, nrhs: number of input arguments on left and right-hand side respectively
c   thus one can test and depending on number of arguments choose different options
c
c   plhs, prhs: pointer to the different arguments: plhs(i), i=1,nlhs
c   (plhs and prhs are integer arrays)
c
c   fmex function: mxGetPr gets the pointer value in prhs, plhs
c   .              mxCopyPtrToReal8: copy pointed values to local array
c   .              mxCopyReal8ToPtr: copy local array values to pointed array
c   .              mxGetM: get nb of rows in pointed array
c   .              mxGetN: get nb of columns in pointed array
c   .              mxCreateFull: creates a matrix for the return argument
c
c   matlab call expected:
c
c   >> [yout{,youtp,youtpp}] = interpos(kopt,xin,yin{,xout}{,taus}{,nbc,ybc{,sigmain}})
c
c   where {..} means facultative, thus one can have for example the following calls:
c
c   >> [yout]=interpos(kopt,xin,yin)
c   >> [yout]=interpos(kopt,xin,yin,xout)
c   >> [yout]=interpos(kopt,xin,yin,taus)
c   >> [yout,youtp]=interpos(kopt,xin,yin)
c   >> [yout,youtp,youtpp]=interpos(kopt,xin,yin,xout,taus)
c   >> [yout,youtp,youtpp]=interpos(kopt,xin,yin,xout,taus,nbc,ybc)
c   >> [yout,youtp,youtpp]=interpos(kopt,xin,yin,taus,nbc,ybc)
c   >> [yout,youtp,youtpp]=interpos(kopt,xin,yin,taus,nbc,ybc,sigmain)
c
c   NOTE: It cannot discrimenate between xout and taus if xout is a single point => xout SHOULD be an array
c   of at least 2 dimensions, except if all 8 arguments are given. OR give both "xout, taus" as 2 scalars. 
c   Since nbc should be an array, it will detect that both xout(1) and taus have been given. 
c   Even for taus=0, put xout,0 in inputs.
c
c   where the arguments have the following meaning with respect to the function y(x)
c   that one wants to interpolate or compute its derivatives:
c
c   kopt: option for interpolation and extrapolation method:
c   .     = 1, 11 or 21 : linear interpolation
c   .     = 2, 12 or 22 : quadratic "
c   .     = 3, 13 or 23 : cubic spline interpolation, with tension=taus
c   .     < 10 : send warning if need to extrapolate
c   .     in [11,19] : extrapolate using same order as interpolation
c   .     in [21,29] : extrapolate using a lower order than for interpolation
c   xin : array giving the input values of x
c   yin : array giving the input values of y(x=xin(i))
c   xout: array of values of x at which the function and/or its derivatives
c   .     are to be computed. If xout is not given, assumes xout=xin
c   yout: interpolated values of y at x=xout(i).
c   youtp: 1st derivative of y at x=xout(i)
c   youtpp: 2nd derivative of y at x=xout(i)
c   taus  : tension value for cubic spline interpolation. If not given, uses taus=0
c   nbc(2): [NBCLFT NBCRGT] (default: [0 0])
c   yxbc(2 or 4) (default: [0 0]): [YBCLFT YBCRGT] or [YBCLFT YBCRGT XBCLFT XBCRGT] with
C     BOUNDARY CONDITIONS, 3 TYPES DETERMINED BY THE VALUE OF (NBCLFT/RGT):
C
C     0) VALUE OF SECOND DERIVATIVE AT XBCLFT OR RGT IS GIVEN (0 OR 10)
C     1) VALUE OF 1ST        "       "   "     "  "   "   "   (1 OR 11)
C     2) VALUE OF FUNCTION AT XBCLFT OR RGT IS GIVEN          (2 OR 12)
C
C     THE VALUE IS GIVEN BY YBCLFT OR YBCRGT RESPECTIVELY.
C
C     FOR TYPE 1: IF (YBCLFT OR YBCRGT > 1E31 THEN DER. FROM LAGRANGIAN INTERP.
C     FOR TYPE 1: IF (YBCLFT OR YBCRGT <-1E31 THEN DER. FROM LINEAR     INTERP.
C
C     IF NBCLFT OR NBCRGT IS < 10, PXIN(1) OR PXIN(KNIN) IS USED INSTEAD
C     OF XBCLFT OR XBCRGT, RESPECTIVELY => XBCLFT OR XBCRGT NOT USED
c
c   ybc(5:6): [YBCLFT YBCRGT XBCLFT XBCRGT PXEXP0 PXEXPDL]
c        Enables one to specify a gaussian wight to taus with PXEXP0 and PXEXPDL:
c              taus_eff =  taus * EXP(-((XIN-PXEXP0)/PXEXPDL)**2)
c        IF PXEXP0 NOT IN [PXIN(1),PXIN(KNIN)], EXP() IGNORED AND taus_eff=cst= taus
c        if ybc(5:6) not given, then  PXEXP0=xin(1)-1. and pxexpdl=1. to get constant taus
c
c   sigmain: array of error bar at each input (xin,yin) data point (assumed=cst=1 if not given)
c            it is normalized to minimum value so that tau*sigma(i) = tau at minimum
c
c
c
      integer lenreal
      PARAMETER(LENREAL=8)
      integer*4 plhs(*), prhs(*)
      integer nlhs, nrhs
c
      integer*4 mxGetPr
c   pointers for matlab input arguments
      integer*4 pkopt, pxin, pyin, pxout, ptaus, pnbc, pybc, psig
      pointer(iptr_xin,xin)
      pointer(iptr_yin,yin)
      pointer(iptr_xout,xout)
      pointer(iptr_sig,sig)
      real*8 xin(1), yin(1), xout(1), ybc(6), sig(1)
      integer kopt, nbc(2)
      real*8 taus, zkopt, znbc(2)
c   pointers for matlab output arguments
      integer*4 pyout, pyoutp, pyoutpp
      pointer(iptr_yout,yout)
      pointer(iptr_youtp,youtp)
      pointer(iptr_youtpp,youtpp)
      real*8 yout(1), youtp(1), youtpp(1)
      real*8 XX, ZTAUEFF, ZCOFEXP, ZXEXP0, ZXEXPDL
      integer mxgetm, mxgetn, mxcreatefull
      integer*4 malloc_f
      integer nbyt_tot, infomalloc
      common /memoryuse/ nbyt_tot
c      external mxgetm, mxgetn, mxgetm, mxgetn,
c     +  mxcreatefull
c   added statement for variables
      real*8 sigmin, fun_sigma
      integer kopt_sign, inttype, iextrapo, ninrow, nincol, nin, nbytes,
     +  nout, ixoutxin, i4len, i, inextrhs, idoexp_error, iybclen,
     +  noutrow, noutcol, nbytes8, ioptder, iflag, i5len, i4isxout
c
      fun_sigma(XX)= ZTAUEFF*EXP(-ZCOFEXP*(XX-ZXEXP0)**2/ZXEXPDL**2)
c
c.......................................................................
c
c   1. defaults values
c
      nbyt_tot=0
      taus = 0.
c
c   2. Input arguments
c
c%OS      print *,' nb of return arguments= ',nlhs
c%OS      print *,' nb of input arguments= ',nrhs
c
c   2.1 kopt
c
      pkopt=mxGetPr(prhs(1))
c   call mxCopyPtrToInteger4(pkopt, kopt, 1)  ! does not work => use real8
      call mxCopyPtrToReal8(pkopt, zkopt, 1)
      kopt = abs(zkopt)
      kopt_sign=zkopt
      kopt_sign=sign(1,kopt_sign)
c%OS      print *,' kopt= ',kopt
c
c   define interpolation type and extrapolation type
c
      inttype = mod(kopt,10)
      iextrapo = kopt/10
c
c   2.2 get xin, first length and allocate space to yin as well
c
      pxin=mxGetPr(prhs(2))
      ninrow = mxGetM(prhs(2))
      nincol = mxGetN(prhs(2))
c   input 1D arrays in row or column
      nin = max(ninrow,nincol)
      nbytes = 8*nin
c%OS      print *,nbytes
      iptr_xin = malloc_f(nbytes,infomalloc)
      iptr_yin = malloc_f(nbytes,infomalloc)
      call mxCopyPtrToReal8(pxin, xin, nin)
c%OS      print *,' nin= ',nin
c%OS      print *,' xin= ',(xin(i),i=1,nin)
c
c   2.3 get yin
c
      pyin=mxGetPr(prhs(3))     ! Get yin
      call mxCopyPtrToReal8(pyin, yin, nin)
c%OS      print *,' yin= ',(yin(i),i=1,nin)
c
c   2.4 check for 4th input
c
      if (nrhs .eq. 3) then
        nout = nin
c   flag for xout=xin set to 1
        ixoutxin = 1
      else
c
c   check if 4th argument is an array.
c   Assume value is taus only if length is one (scalar), cubic spline interpolation
c   option was asked, and 8 input arguments are given. Also check if 4th and 5th argument ara scalars.
c   Otherwise the 4th argument is assumed to be xout
c
        i4len = max(mxGetM(prhs(4)),mxGetN(prhs(4)))
        i4isxout = 0
        if (nrhs.ge.5) then
          i5len = max(mxGetM(prhs(5)),mxGetN(prhs(5)))
          if (i4len.eq.1 .and. i5len.eq.1) then
            i4isxout = 1
          endif
        endif
        if (i4len.eq.1 .and. nrhs.ne.8 .and. inttype.eq.3 .and. 
     +      i4isxout.eq.0) then
c   in this case prhs(4) gives the tau value  
          ixoutxin = 1
          ptaus = mxGetPr(prhs(4))
          call mxCopyPtrToReal8(ptaus, taus, 1)
        else
          ixoutxin = 0
          nout = i4len
          nbytes = 8*nout
          iptr_xout = malloc_f(nbytes,infomalloc)
          pxout = mxGetPr(prhs(4))
          call mxCopyPtrToReal8(pxout, xout, nout)
        endif
      endif
      if (ixoutxin .eq. 1) then
        nout = nin
        nbytes = 8*nout
        iptr_xout = malloc_f(nbytes,infomalloc)
        nbytes = LENREAL*nout
        do i=1,nout
          xout(i) = xin(i)
        end do
      endif
c%OS      print *,'ixoutxin, nout= ',ixoutxin, nout
c%OS      print *,' xout= ',(xout(i),i=1,nout)
c
c   2.5 5th argument is taus if xout .ne. xin
c
      inextrhs=5
      if (nrhs.ge.5 .and. ixoutxin.eq.0) then
        ptaus = mxGetPr(prhs(5))
        inextrhs=6
        call mxCopyPtrToReal8(ptaus, taus, 1)
      endif
c
c   2.6 5,6th or 6,7th arguments are nbc, ybc
c
      nbc(1) = 0
      nbc(2) = 0
      ybc(1) = 0.
      ybc(2) = 0.
      idoexp_error=0
      if (nrhs .ge. 6) then
        pnbc = mxGetPr(prhs(inextrhs))
        inextrhs = inextrhs + 1
        call mxCopyPtrToReal8(pnbc, znbc, 2)
        nbc(1) = znbc(1)
        nbc(2) = znbc(2)
c%OS        print *,'taus,nbc= ',taus,nbc,' znbc= ',znbc
        iybclen = max(mxGetM(prhs(inextrhs)),mxGetN(prhs(inextrhs)))
c%OS        print *,' iybclen= ',iybclen
        pybc = mxGetPr(prhs(inextrhs))
        inextrhs = inextrhs + 1
        call mxCopyPtrToReal8(pybc, ybc, iybclen)
        if (iybclen .le. 4) then
          ybc(5) = xin(1) - xin(nin)
          ybc(6) = 1.
        else
          idoexp_error = 1
        endif
      else
        ybc(5) = xin(1) - xin(nin)
        ybc(6) = 1.
      endif
c%OS        print *,ybc
c
c   2.7 7th or 8th is sigmain
      nbytes = LENREAL*nin
      iptr_sig = malloc_f(nbytes,infomalloc)
      if ((ixoutxin.eq.1 .and. nrhs.eq.7) .or.
     +  (ixoutxin.eq.0 .and. nrhs.eq.8)) then
        psig = mxGetPr(prhs(inextrhs))
        inextrhs = inextrhs + 1
        call mxCopyPtrToReal8(psig, sig, nin)
        sigmin=sig(1)
        do i=2,nin
          sigmin=min(sigmin,sig(i))
        end do
        do i=1,nin
          sig(i) = sig(i)/sigmin*abs(taus)
        end do
      else
        ZTAUEFF = abs(taus)
        ZXEXP0 = ybc(5)
        ZXEXPDL = ybc(6)
        ZCOFEXP = 0.
        if (idoexp_error.EQ.1) THEN
          ZCOFEXP = 1.0
        endif
        do i=1,nin
          sig(i)=fun_sigma(xin(i))
        end do
      endif
      if (taus.LT.0.) then
        sig(1)=10.*sig(1)
      endif
c
c   3. Allocate return arguments arrays (in same geometry as xin)
c
      if (ninrow .eq. 1) then
        noutrow = 1
        noutcol = nout
      else
        noutrow = nout
        noutcol = 1
      endif
      nbytes8 = 8*nout
      nbytes =  LENREAL*nout
c   pointers for yout, youtp and youtpp respectively
      plhs(1) = mxCreateFull(noutrow,noutcol,0)
      pyout = mxGetPr(plhs(1))
      iptr_yout = malloc_f(nbytes8,infomalloc)
      ioptder = 0
      if (nlhs .ge. 2) then
        plhs(2) = mxCreateFull(noutrow,noutcol,0)
        pyoutp = mxGetPr(plhs(2))
        iptr_youtp = malloc_f(nbytes8,infomalloc)
        ioptder = 1
      endif
      if (nlhs .ge. 3) then
        plhs(3) = mxCreateFull(noutrow,noutcol,0)
        pyoutpp = mxGetPr(plhs(3))
        iptr_youtpp = malloc_f(nbytes8,infomalloc)
        ioptder = 2
      endif
c
c   4. compute interpolated function and its derivatives
c
c%OS      print *,' ioptder,iextrapo= ',ioptder,iextrapo
      if (inttype .eq. 1) then
        call intlinear(xin,yin,nin,xout,yout,youtp,youtpp,nout
     +    ,ioptder,iextrapo)
      end if
c
      if (inttype .eq. 2) then
        call parint(xin,yin,nin,xout,yout,youtp,youtpp,
     +    nout,ioptder,iextrapo)
      end if
c
      if (inttype .eq. 3) then
        if (iextrapo .eq. 1) then
          iextrapo = 3
        else if (iextrapo .eq. 2) then
          iextrapo = 31
        else if (iextrapo .eq. 3) then
          iextrapo = 32
        else if (iextrapo .eq. 4) then
          iextrapo = 1
        else if (iextrapo .eq. 5) then
          iextrapo = 2
        else
          iextrapo = 21
        endif
C   add sign of kopt for extrapol
        call cbsplgenrid(xin,yin,nin,xout,yout,youtp,youtpp,nout,
     +    ioptder,kopt_sign*iextrapo,sig,nbc,ybc,iflag)
        if (iflag .ne. 0) then
          print *
          print *,' problem in cbsplgen'
          print *,' iflag = ',iflag
          return
        end if
      end if
c   
c   5. copy arrays to return value pointers
c
      call mxCopyReal8ToPtr(yout, pyout, nout)
      if (ioptder .ge. 1) then
        call mxCopyReal8ToPtr(youtp,pyoutp,nout)
      endif
      if (ioptder .ge. 2) then
        call mxCopyReal8ToPtr(youtpp,pyoutpp,nout)
      endif
c
      call free(iptr_xin)
      call free(iptr_yin)
c      call free(iptr_xout)
c      call free(iptr_yout)
      call free(iptr_sig)
c      call free(iptr_youtp)
c      call free(iptr_youtpp)

      return
      end
c.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
c
      subroutine intlinear(xin,yin,nin,xout,yout,youtp,youtpp,nout
     +  ,ioptder,iextrapo)
c
c   compute linear interpolation of (xin,yin) on (xout,yout).
c
c   ioptder = 0: compute only yout (youtp, youtpp not used)
c   ioptder = 1: compute also 1st der. in youtp (youtpp not used)
c   ioptder = 2: as 1 and 2nd der. in youtpp
c
c   iextrapo = 0: send message if need to extrapolate
c   iextrapo = 1: normal extrapolation
c   iextrapo = 2: y=y_edge for extrapolation => constant extrapolation
c
      implicit none
      integer nin, nout, ioptder, iextrapo, i, j
      real*8 xin(nin), yin(nin)
      real*8 xout(nout), yout(nout), youtp(nout), youtpp(nout)
c
      real*8 x1, x2, y1, y2, x, zx, flinear, flinearp
c
      flinear(x1,x2,y1,y2,x) = y1 + (x-x1)*(y2-y1)/(x2-x1)
      flinearp(x1,x2,y1,y2) = (y2-y1)/(x2-x1)
c.......................................................................
c
      do i=1,nout
        zx = xout(i)
        if (zx .lt. xin(1)) then
          if (iextrapo .eq. 0) print *,' warning point outside interval'
     +      ,' => lhs extrapolation needed'
          if (iextrapo .eq. 2) then
            yout(i) = yin(1)
            if (ioptder .ge. 1) youtp(i) = 0.0
            if (ioptder .ge. 2) youtpp(i) = 0.0
          else
            yout(i) = flinear(xin(1),xin(2),yin(1),yin(2),zx)
            if (ioptder .ge. 1) youtp(i) = flinearp(xin(1),xin(2),yin(1)
     +        ,yin(2))
            if (ioptder .ge. 2) youtpp(i) = 0.0
          endif
          go to 100
        endif
c
        do j=1,nin-1
          if (xin(j+1) .ge. zx) then
            yout(i) = flinear(xin(j),xin(j+1),yin(j),yin(j+1),zx)
            if (ioptder .ge. 1) youtp(i) = flinearp(xin(j),xin(j+1),
     +        yin(j),yin(j+1))
            if (ioptder .ge. 2) youtpp(i) = 0.
            go to 100
          endif
        end do
c   zx .gt. xin(nin)
        if (iextrapo .eq. 0) print *,' warning point outside interval'
     +    ,' => rhs extrapolation needed'
        if (iextrapo .eq. 2) then
          yout(i) = yin(nin)
          if (ioptder .ge. 1) youtp(i) = 0.0
          if (ioptder .ge. 2) youtpp(i) = 0.0
        else
          yout(i) = flinear(xin(nin-1),xin(nin),yin(nin-1),yin(nin),zx)
          if (ioptder .ge. 1) youtp(i) = flinearp(xin(nin-1),xin(nin)
     +      ,yin(nin-1),yin(nin))
          if (ioptder .ge. 2) youtpp(i) = 0.0
        endif
  100   continue
      end do
c
      return
      end
c......................................................................
c     
      SUBROUTINE PARINT(XIN,YIN,nin,XOUT,YOUT,youtp,youtpp
     +  ,nout,ioptder,iextrapo)
C   iextrapo = 0 : warn and keep last value
C   iextrapo = 1 : use quadratic for extrapolation
C   iextrapo = 2 : use linear for extrapolation
      implicit none
      INTEGER nin,nout
      INTEGER I,J
      REAL*8  XIN(nin), YIN(nin)
      REAL*8  XOUT(nout), YOUT(nout), youtp(nout), youtpp(nout)
      REAL*8 F(2)
C
      real*8 x1, x2, x3, y1, y2, y3, x, zx, fparabol, fparabolp,
     +  fparabolpp, flinear, flinearp
      integer ioptder, iextrapo
C
      fparabol(x1,x2,x3,y1,y2,y3,x) = 
     +  ((x-x1)*(x-x2)*y3)/((x3-x1)*(x3-x2))+
     +  ((x-x1)*(x-x3)*y2)/((x2-x1)*(x2-x3))+
     +  ((x-x2)*(x-x3)*y1)/((x1-x2)*(x1-x3))
      fparabolp(x1,x2,x3,y1,y2,y3,x) =
     +  ((x-x1)+(x-x2))*y3/((x3-x1)*(x3-x2))+
     +  ((x-x1)+(x-x3))*y2/((x2-x1)*(x2-x3))+
     +  ((x-x2)+(x-x3))*y1/((x1-x2)*(x1-x3))
      fparabolpp(x1,x2,x3,y1,y2,y3) =
     +  (y3/((x3-x1)*(x3-x2))+
     +  y2/((x2-x1)*(x2-x3))+
     +  y1/((x1-x2)*(x1-x3)))*2
C
      flinear(x1,x2,y1,y2,x) = y1 + (x-x1)*(y2-y1)/(x2-x1)
      flinearp(x1,x2,y1,y2) = (y2-y1)/(x2-x1)
C
      do i=1,nout
        zx = xout(i)
C   ZX <= interval
        if (zx .le. xin(1)) then
          if (iextrapo .eq. 0) then
            print *,' warning point outside interval'
     +        ,' => lhs extrapolation needed'
            yout(i) = yin(1)
            if (ioptder .ge. 1)  youtp(i) = 0.
            if (ioptder .ge. 2)  youtpp(i) = 0.
          else if (iextrapo .eq. 1) then
            yout(i) = fparabol(xin(1),xin(2),xin(3),yin(1),
     +        yin(2),yin(3),zx)
            if (ioptder .ge. 1)  youtp(i) = 
     +        fparabolp(xin(1),xin(2),xin(3),yin(1),yin(2),yin(3),zx)
            if (ioptder .ge. 2)  youtpp(i) =
     +        fparabolpp(xin(1),xin(2),xin(3),yin(1),yin(2),yin(3))
          else if (iextrapo .eq. 2) then
            yout(i) = flinear(xin(1),xin(2),yin(1),yin(2),zx)
            if (ioptder .ge. 1) youtp(i) = flinearp(xin(1),xin(2),yin(1)
     +        ,yin(2))
            if (ioptder .ge. 2) youtpp(i) = 0.0
          endif 
        else if (zx .ge. xin(nin)) then
C   ZX >= interval
          if (iextrapo .eq. 0) then
            print *,' warning point outside interval'
     +        ,' => rhs extrapolation needed'
            yout(i) = yin(nin)
            if (ioptder .ge. 1)  youtp(i) = 0.
            if (ioptder .ge. 2)  youtpp(i) = 0.
          else if (iextrapo .eq. 1) then
            yout(i) = fparabol(xin(nin-2),xin(nin-1),xin(nin),yin(nin-2)
     +        ,yin(nin-1),yin(nin),zx)
            if (ioptder .ge. 1)  youtp(i) = 
     +        fparabolp(xin(nin-2),xin(nin-1),xin(nin),yin(nin-2),
     +        yin(nin-1),yin(nin),zx)
            if (ioptder .ge. 2)  youtpp(i) =
     +        fparabolpp(xin(nin-2),xin(nin-1),xin(nin),yin(nin-2),
     +        yin(nin-1),yin(nin)) 
          else if (iextrapo .eq. 2) then
            yout(i) =flinear(xin(nin-1),xin(nin),yin(nin-1),yin(nin),zx)
            if (ioptder .ge. 1) youtp(i) = 
     +        flinearp(xin(nin-1),xin(nin),yin(nin-1),yin(nin))
            if (ioptder .ge. 2) youtpp(i) = 0.0
          end if 
C   
        else
C   ZX inside interval
          DO J = 2,nin
            IF (zx .LE. XIN(J)) THEN
              IF ( J.EQ.2 ) THEN
                yout(i) = fparabol(xin(1),xin(2),xin(3),yin(1),
     +            yin(2),yin(3),zx)
                if (ioptder .ge. 1)  youtp(i) = fparabolp(
     +            xin(1),xin(2),xin(3),yin(1),yin(2),yin(3),zx)
                if (ioptder .ge. 2)  youtpp(i) =
     +            fparabolpp(xin(1),xin(2),xin(3),yin(1),yin(2),yin(3))
              ELSE IF ( J.EQ.Nin ) THEN
                yout(i) = fparabol(xin(nin-2),xin(nin-1),xin(nin),
     +            yin(nin-2),yin(nin-1),yin(nin),zx)
                if (ioptder .ge. 1)  youtp(i) = 
     +            fparabolp(xin(nin-2),xin(nin-1),xin(nin),yin(nin-2),
     +            yin(nin-1),yin(nin),zx)
                if (ioptder .ge. 2)  youtpp(i) =
     +            fparabol(xin(nin-2),xin(nin-1),xin(nin),yin(nin-2),
     +            yin(nin-1),yin(nin),zx) 
              ELSE
                F(1) = fparabol(xin(j-2),xin(j-1),xin(j),yin(j-2),
     +            yin(j-1),yin(j),zx) 
                F(2) = fparabol(xin(j-1),xin(j),xin(j+1),yin(j-1),
     +            yin(j),yin(j+1),zx) 
                YOUT(I) = (F(1) + F(2))/2
                if (ioptder .ge. 1) then 
                  F(1) = fparabolp(xin(j-2),xin(j-1),xin(j),yin(j-2),
     +              yin(j-1),yin(j),zx) 
                  F(2) = fparabolp(xin(j-1),xin(j),xin(j+1),yin(j-1),
     +              yin(j),yin(j+1),zx) 
                  YOUTp(I) = (F(1) + F(2))/2
                END IF            
                if (ioptder .ge. 2) then 
                  F(1) = fparabolpp(xin(j-2),xin(j-1),xin(j),yin(j-2),
     +              yin(j-1),yin(j)) 
                  F(2) = fparabolpp(xin(j-1),xin(j),xin(j+1),yin(j-1),
     +              yin(j),yin(j+1)) 
                  YOUTpp(I) = (F(1) + F(2))/2
                END IF
              END IF  
C   yout defined
              GO TO 100
            END IF
          END DO
  100     CONTINUE
        ENDIF
C
      END DO   
      return 
      END
C
C.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
C
      SUBROUTINE cbsplgenrid(XIN,YIN,nin,xout,yout,youtp,youtpp,
     +  nout,ioptder,iextrapo,psig,nbc,ybc,IFLAG)
c
      implicit none
      integer LENREAL
      parameter(LENREAL=8)
      integer nin,nout, nbc(2)
      REAL*8  psig(nin)
      REAL*8  xin(nin), yin(nin), ybc(6)
      REAL*8  xout(nout),yout(nout), youtp(nout), youtpp(nout)
      INTEGER ioptder, NBCLFT,NBCRGT, iextrapo
      real*8 XBCLFT,XBCRGT, YBCLFT,YBCRGT
      REAL*8 PXEXP0,PXEXPDL
      pointer(iptr_pynew,pynew)
      real*8 pynew(1)
      pointer(iptr_pyinpp,pyinpp)
      real*8 pyinpp(1)
      pointer(iptr_pamat,pamat)
      real*8 pamat(1)
      pointer(iptr_pwork,pwork)
      real*8 pwork(1)
      pointer(iptr_kwork,kwork)
      integer kwork(1)
C
      integer iflag, idamat, mdamat, nbytes, infomalloc
      integer*4 malloc_f
C
c
c%OS      NBCLFT = 1
c%OS      NBCRGT = 1
c%OS      YBCLFT = 1.E32
c%OS      YBCRGT = 1.E32
      NBCLFT = nbc(1)
      NBCRGT = nbc(2)
      YBCLFT = ybc(1)
      YBCRGT = ybc(2)
      if (NBCLFT .ge. 10) XBCLFT = ybc(3)
      if (NBCRGT .ge. 10) XBCRGT = ybc(4)
c
      IDAMAT = 3
      IF (PSIG(1) .EQ. 0.) IDAMAT = 2
      IF (NBCLFT.GE.10 .OR. NBCRGT.GE.10 .OR. NBCLFT.EQ.2
     +  .OR. NBCRGT.EQ.2) IDAMAT = 3*(IDAMAT-1)+1
      IF (NBCLFT .GE. 10) IDAMAT = IDAMAT + 1
      IF (NBCRGT .GE. 10) IDAMAT = IDAMAT + 2
      mdamat = IDAMAT*nin
      nbytes = LENREAL*nin
      iptr_pynew = malloc_f(nbytes,infomalloc)
      iptr_pyinpp = malloc_f(nbytes,infomalloc)
      nbytes = LENREAL*mdamat
      iptr_pamat = malloc_f(nbytes,infomalloc)
      nbytes = LENREAL*(2*nin)
      iptr_pwork = malloc_f(nbytes,infomalloc)
      nbytes = LENREAL*(5*nin)
      iptr_kwork = malloc_f(nbytes,infomalloc)
c%OS      PXEXP0 = xin(1) - 1.
c%OS      PXEXPDL= 1.0
      PXEXP0 = ybc(5)
      PXEXPDL= ybc(6)
      if (infomalloc .eq. 0) then
        CALL CBSPLGEN(XIN,YIN,PYNEW,PYINPP,Nin,XOUT,YOUT,YOUTP,YOUTPP,
     +    nout,ioptder,PSIG,KWORK,PWORK,PAMAT,mdamat,NBCLFT,
     +    NBCRGT,XBCLFT,XBCRGT,YBCLFT,YBCRGT,iextrapo,PXEXP0,PXEXPDL
     +    ,IFLAG)
      endif
c
      call free(iptr_pynew)
      call free(iptr_pyinpp)
      call free(iptr_pamat)
      call free(iptr_pwork)
      call free(iptr_kwork)
c
      END 
c
c.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.     
c
      SUBROUTINE CBSPLGEN(PXIN,PYIN,PYINNEW,PYINPP,KNIN,PXOUT,PYOUT,
     +  PYOUTP,PYOUTPP,KNOUT,KOPT,PSIG,KWORK,PWORK,PAMAT,MDMATOT,NBCLFT
     +  ,NBCRGT,XBCLFT,XBCRGT,YBCLFT,YBCRGT,KEXTRPOL,PXEXP0,PXEXPDL
     +  ,KFLAG)
C     =================================================================
C
C     NOTE: THIS ROUTINE INCLUDES THE STANDARD CUBIC SPLINE IF PTAUS=0 (i.e. psig=0):
C           THEN PYINNEW IS NOT USED AND MDAMAT=3 IS SUFFICIENT
C           (=> PYINNEW(1) OR PYINNEW=PYIN IS OK)
C
C     Interpolate (pxin,pyin) on (pxout,pyout) using
C     Hirshman fitted cubic spline with ptaus value or
C     standard cubic spline if PTAUS=0 (psig=ptaus*sigma/sigmamin)
C
C     KOPT = 0: ONLY INTERPOLATE FUNCTION INTO PYOUT
C     KOPT = 1: INTERPOLATE FUNCTION INTO PYOUT AND 1ST DER. INTO PYOUTP
C     KOPT = 2: AS KOPT=1 PLUS 2ND DER. INTO PYOUTPP
C
C     MDMATOT = TOTAL DIMENSION OF PAMAT. THE REQUIRED SPACE DEPENDS IF
C     .         PTAUS IS ZERO OR NOT AND ON THE B.C. (SYMMETRIC OR NOT)
C     THUS, MDMATOT CAN VARY BETWEEN 2*KNIN AND 10*KNIN (MAX. VALUE NEEDED)
C
C     SEE COMMENTS FOR ROUTINE CBFITBND FOR MORE INFORMATION
C
C     IF LAPACK ROUTINES NOT AVAILABLE, USE spgbtrf_s.f
C     (LAPACK SOURCE COPIED FROM NETLIB.ORG)
C
C     PXIN    : INPUT ABSCISSA (GIVEN DATA)
C     PYIN    : INPUT VALUES OF FUNCTION AT PXIN(I),I=1,KNIN
C     PYINNEW : IF PTAUS.NE.0, THEN PYNEW CONTAINS ON OUTPUT THE NEW VALUES
C     .         OF THE FUNCTION AT PXIN(I) FOR THE CUBIC SPLINE FIT
C     PYINPP  : SECOND DER. OF THE CUBIC SPLINE FIT FOR (PXIN,PYIN) IF PTAUS=0 OR
C     .         ON (PXIN,PYNEW) OTHERWISE
C     KNIN    : NUMBER OF INPUT POINTS
C     PXOUT   : X VALUES AT WHICH THE FUNCTION HAS TO BE INTERPOLATED (INPUT)
C     PYOUT   : INTERPOLATED VALUES AT PXOUT(I),I=1,KNOUT (OUTPUT)
C     PYOUTP  : INTERPOLATED VALUES OF 1ST DER. OF FUNCTIONS AT PXOUT(I) (OUTPUT, IF KOPT.GE.1)
C     PYOUTPP : INTERPOLATED VALUES OF 2ND DER. OF FUNCTIONS AT PXOUT(I) (OUTPUT, IF KOPT.EQ.2)
C     KNOUT   : NUMBER OF POINTS FOR OUTPUT
C     KOPT    : SEE ABOVE
C     PSIG    : SIGMA at each point, normalized by minimum sigma, times PTAUS
C     PTAUS   : WEIGHT OF SECOND DERIVATIVE IN THE CHI**2 TO BE MINIMIZED. PTAUS=0 GIVES THE
C     .         STANDARD CUBIC SPLINE. LARGER VALUES OF PTAUS WILL SMOOTH MORE THE 2ND DER.
C     KWORK   : INTEGER WORK SPACE OF DIMENSION 5*KNIN
C     PWORK   : REAL WORK SPACE OF DIMENSION 2*KNIN
C     PAMAT   : REAL WORK SPACE FOR THE MATRIX OF DIMENSION MDMATOT
C     MDMATOT : DIMENSION OF PAMAT, WHICH VALUE SHOULD BE IN [2*KNIN,10*KNIN], DEPENDING ON
C     .         PTAUS BEING 0 OR NOT AND ON THE BOUNDARY CONDITION (SEE PARAGRAPH 1. BELOW)
C     NBCLFT  : FOR LEFT B.C., VALUE SHOULD BE 0,1,2,10,11 OR 12. (SEE ROUTINE CBFITBND BELOW)
C     NBCRGT  : FOR RIGHT B.C. (SEE ROUTINE CBFITBND BELOW)
C     XBCLFT  : FOR LEFT B.C., USED ONLY IF NBCLFT.GE.10 (SEE ROUTINE CBFITBND BELOW)
C     XBCRGT  : FOR RIGHT B.C., USED ONLY IF NBCRGT.GE.10 (SEE ROUTINE CBFITBND BELOW)
C     YBCLFT  : VALUE OF LEFT B.C.
C     YBCRGT  : VALUE OF RIGHT B.C.
C     .         STANDARD B.C. (SECOND DER. = 0) IS OBTAINED WITH:
C     .         NBCLFT = NBCRGT = 0 AND YBCLFT = YBCRGT = 0.
C     KEXTRPOL: OPTION ON HOW TO EXTRAPOLATE THE FUNCTION IF PXOUT(I) IS OUTSIDE [PXIN(1),PXIN(KNIN)]
C     .       = 0: STOP WITH ERROR MESSAGE IF OUT OF BOUND
C     .       = 1: LINEAR EXTRAPOLATION
C     .       = 2: USE QUADRATIC INTERPOLATION IF X OUT OF BOUND
C     .       = 3: USE CUBIC INTERPOLATION IF X OUT OF BOUND
C     .       = 21: USE QUADRATIC WITHIN ALFA*DELTA_X AND LINEAR FURTHER
C     .       = 31: USE CUBIC WITHIN ALFA*DELTA_X AND LINEAR    FURTHER
C     .       = 32: USE CUBIC WITHIN ALFA*DELTA_X AND QUADRATIC FURTHER
C     PXEXP0  : PTAUS IS WEIGHTED BY AN EXP(-((X-PXEXP0)/PXEXPDL)**2)
C     PXEXPDL : IF PXEXP0 NOT IN [PXIN(1),PXIN(KNIN)], EXP() IGNORED AND PTAUS=CST
C     .         (SEE ROUTINE CBFITBND BELOW)
C     .         GIVE PXEXP0=PXIN(1)-1. AND PXEXPDL=1. TO GET CONSTANT PTAUS
C     KFLAG   : ERROR FLAG: IF NOT 0, THERE IS A PROBLEM
C
C-----------------------------------------------------------------------
cc      implicit real*8 (p)
      implicit none
      real*8 EPSILON
      PARAMETER(EPSILON = 1.0E-10)
C
      integer knin, knout, kopt, mdmatot, nbclft, nbcrgt, kextrpol,
     +  kflag, i, idamat
      real*8 pxin, pyin, pyinnew, pyinpp, pxout, pyout, pyoutp, pyoutpp
     +  ,psig, pamat, pwork
      integer kwork
      DIMENSION PXIN(KNIN), PYIN(KNIN), PYINNEW(KNIN), PYINPP(KNIN),
     +  PXOUT(KNOUT), PYOUT(KNOUT), PYOUTP(KNOUT), PYOUTPP(KNOUT),
     +  PSIG(KNIN)
C   
      DIMENSION PAMAT(MDMATOT), PWORK(2*KNIN)
      DIMENSION KWORK(5*KNIN)
C
      real*8 xbclft, xbcrgt, ybclft, ybcrgt, pxexp0, pxexpdl, zy, zyp,
     +  zypp
C
C-----------------------------------------------------------------------
C     0. CHECK INPUT CONSISTENCY
C
      KFLAG = 0
      IF (PSIG(1) .EQ. 0.) THEN
        IF (NBCLFT.GE.10 .OR. NBCRGT.GE.10 .OR. NBCLFT.EQ.2
     +  .OR. NBCRGT.EQ.2) THEN
          PRINT *,' PTAUS=0, BUT NEED SMOOTHING, WHEN'
          PRINT *,'     NBCLFT.GE.10 .OR. NBCRGT.GE.10 .OR.',
     +      ' NBCLFT.EQ.2 .OR. NBCRGT.EQ.2'
          PRINT *,' NBCLFT = ',NBCLFT
          PRINT *,' NBCRGT = ',NBCRGT
c%OS          STOP 'TAU=0'
          KFLAG = 1
          return
        ENDIF
      ENDIF
C
C   PXIN in ASCENDING ORDER
C
      DO i=1,KNIN-1
        if (PXIN(i) .GE. PXIN(i+1)) then
          print *,' xin not in ascending order:'
          print *,' xin(',i,')= ',PXIN(i),'   >=   xin(',i+1,')= ',
     +      PXIN(i+1)
          KFLAG = 2
          RETURN
        endif
      END DO
C
C     1. PREPARE MATRIX DIMENSION.
C     MINIMUM REQUIRED:
C     IUP = 1 = IDOWN; IF PTAUS.NE.0 => IUP = IDOWN = 2
C     IF SYMMETRIC, USE ONLY UPPER BAND AND DIAGONAL =>IDAMAT=IUP+1
C     IF ASYMMETRIC => IDAMAT = 2*IDOWN + IUP + 1
C     IF B.C. NOT AT END OF INTERVAL => IUP = IUP + 1, AND/OR IDOWN=IDOWN+1
C
C     => ALTOGETHER, MINIMUM VALUE: IDOWN=1, IUP=1, SYM. =>IDAMAT_MAX = 2
C     => ALTOGETHER, MAXIMUM VALUE: IDOWN=3, IUP=3 =>IDAMAT_MAX = 10
C
      IDAMAT = 3
      IF (PSIG(1) .EQ. 0.) IDAMAT = 2
      IF (NBCLFT.GE.10 .OR. NBCRGT.GE.10 .OR. NBCLFT.EQ.2
     +  .OR. NBCRGT.EQ.2) IDAMAT = 3*(IDAMAT-1)+1
      IF (NBCLFT .GE. 10) IDAMAT = IDAMAT + 1
      IF (NBCRGT .GE. 10) IDAMAT = IDAMAT + 2
C
      IF (MDMATOT .LT. IDAMAT*KNIN) THEN
        PRINT *,' '
        PRINT *,' DIMENSION MDMATOT= ',MDMATOT,' IS TOO SMALL,',
     +    ' NEED IDAMAT*KNIN= ',IDAMAT,'*',KNIN,' = ',IDAMAT*KNIN
c%OS        STOP
        RETURN
      ENDIF
C
      CALL CBFITBND(PXIN,PYIN,PYINNEW,KNIN,PYINPP,PSIG,KWORK,
     +  PWORK(1),PWORK(KNIN+1),PAMAT,IDAMAT,NBCLFT,NBCRGT,
     +  XBCLFT,XBCRGT,YBCLFT,YBCRGT,PXEXP0,PXEXPDL)
C
CL    2. COMPUTE INTERPOLATED VALUE AT EACH PXOUT
C
      DO I=1,KNOUT
        IF (PSIG(1) .EQ. 0.0) THEN
          CALL SPLIBND(PXIN,PYIN   ,PYINPP,KNIN,PXOUT(I),ZY,ZYP,ZYPP,
     +      KEXTRPOL)
        ELSE
          CALL SPLIBND(PXIN,PYINNEW,PYINPP,KNIN,PXOUT(I),ZY,ZYP,ZYPP,
     +      KEXTRPOL)
        ENDIF
        PYOUT(I) = ZY
        IF (KOPT .GE. 1) PYOUTP(I) = ZYP
        IF (KOPT .EQ. 2) PYOUTPP(I) = ZYPP
      END DO
C
      RETURN
      END

C-----------------------------------------------------------------------
      SUBROUTINE CBFITBND(PXIN,PYIN,PYINNEW,KNIN,PYINPP,PSIG,
     +  KPM2,WHK,WOHK,PAMAT,MDAMAT,NBCLFT,NBCRGT,XBCLFT,XBCRGT,YBCLFT,
     +  YBCRGT,PXEXP0,PXEXPDL)
C
C     NON-PERIODIC B.C
C
C     PREPARE SECOND DERIVATIVE OF CUBIC SPLINE INTERPOLATION AND NEW
C     VALUES OF Y AT NODES YNEW FITTED SUCH THAT CHI**2 + TAUS*F''**2
C     IS MINIMIZED ACCORDING TO HIRSHMAN ET AL, PHYS. PLASMAS 1 (1994) 2280.
C     SIG = TAU*SIGMA_K/min(SIGMA_K) OF PAPER.
C
C     SETTING TAUS=0., ONE FINDS THE USUAL CUBIC SPLINE INT. WITH CHI**2=0
C     TAUS LARGE => FIT CLOSER TO STRAIGHT LINE (SECOND DERIV.=0)
C
C     IF LAPACK ROUTINES NOT AVAILABLE, USE spgbtrf_s.f
C     (LAPACK SOURCE COPIED FROM NETLIB.ORG)
C
C     IF TAUS=0, PYINNEW NOT USED => PYINNEW(1) OR PYINNEW=PYIN IS OK
C
C     BOUNDARY CONDITIONS, 3 TYPES DETERMINED BY THE VALUE OF (NBCLFT/RGT):
C
C     0) VALUE OF SECOND DERIVATIVE AT XBCLFT OR RGT IS GIVEN (0 OR 10)
C     1) VALUE OF 1ST        "       "   "     "  "   "   "   (1 OR 11)
C     2) VALUE OF FUNCTION AT XBCLFT OR RGT IS GIVEN          (2 OR 12)
C
C     THE VALUE IS GIVEN BY YBCLFT OR YBCRGT RESPECTIVELY.
C
C     FOR TYPE 1: IF (YBCLFT OR YBCRGT > 1E31 THEN DER. FROM LAGRANGIAN INTERP.
C     FOR TYPE 1: IF (YBCLFT OR YBCRGT <-1E31 THEN DER. FROM LINEAR     INTERP.
C
C     IF NBCLFT OR NBCRGT IS < 10, PXIN(1) OR PXIN(KNIN) IS USED INSTEAD
C     OF XBCLFT OR XBCRGT, RESPECTIVELY => XBCLFT OR XBCRGT NOT USED
C
C     IF END POINTS ARE USED FOR THE B.C. AND TYPE 0 OR 1, THEN USE SYMMETRY
C     OF MATRIX
C
C     IF XBCLFT OR XBCRGT ARE USED, IMPOSE B.C. ON NODE CLOSEST TO XBCLFT OR XBCRGT
C
C     TENSION TAUS(K) IS GIVEN WITH AN EXPONENTIAL FORM TO BE ABLE TO LOCALIZE
C     IT:
C     .     TAU_K = PTAUS * EXP( -COF * ((X-X0)/DX)**2)
C
C     WHERE X0 = PXEXP0 AND DX = PXEXPDL, AND:
C     COF = 1. IF PXEXP0 IN [PXIN(1),PXIN(KNIN)], 0 OTHERWISE
C     THUS SETTING PXEXP0 OUTSIDE DOMAIN GIVES A CST TAU_K VALUE
C
C-----------------------------------------------------------------------
C
      implicit none
      integer LENREAL
      parameter(LENREAL=8)
      integer knin, mdamat, nbclft, nbcrgt, kpm2, ixbc, ibctyp, iik,
     +  itauval, isym, n, i, j, k, iup, idown, idiag, ishift, ieff, ikp2
     +  ,ikp1, ikm1, ikm2, jk, jkp1, jkp2, jeff, iii, iupsofar, idwnsofa
     +  , jbc, ik, idiamik, idiapik, iklft, idima, idimrhs, irhs, info,
     +  info2, jkm1, jkm2, infomalloc
      DIMENSION PXIN(KNIN), PYIN(KNIN), PYINNEW(KNIN),
     +  PYINPP(KNIN), WHK(KNIN), WOHK(KNIN),
     +  PAMAT(MDAMAT,KNIN), ZYBC(2), PSIG(KNIN)
      DIMENSION KPM2(KNIN,-2:+2), IXBC(2), IBCTYP(2)
      pointer(iptr_ftauk,ftauk)
      dimension ftauk(1)
C
      real*8 xbclft, xbcrgt, ybclft, ybcrgt, pxexp0, pxexpdl, pxin,
     +  pyin, pyinnew, pyinpp, whk, wohk, pamat, psig
      real*8 zybc, ftauk, xtkm1, xohkm1, xohkm2, xtk, xhkm1,
     +  xohk, xtkp1, xhk, xohkp1, xykp1, xyk,
     +  xykm1, ztaueff, zcofexp, zxexp0, zxexpdl, a1, a2, a3, a4, b1,
     +  b2, b3, b4, px, fakk, fakkp1, fakkp2, frhs, zdelx, zero, zvalue,
     +  zypeff, ztohkk1, zsign, fa2, fa3, fa0, fa1,
     +  fakkm1, fakkm2, fun_ftauk, fcccc0, fcccc1, fcccc2, fcccc3
      integer*4 malloc_f
C
C
C     FUNCTIONS FOR MATRIX COEFFICIENTS
C
      real*8 zsix, zthree, ztwo, zone
      PARAMETER(zsix=6., zthree=3., ztwo=2., zone=1.)
      FAKKM2(XTKM1,XOHKM1,XOHKM2) = XTKM1*XOHKM1*XOHKM2
      FAKKM1(XTK,XTKM1,XHKM1,XOHK,XOHKM1,XOHKM2) = XHKM1/zsix
     +  - XOHKM1*(XTK*XOHK+(XTK+XTKM1)*XOHKM1 + XTKM1*XOHKM2)
      FAKK(XTKP1,XTK,XTKM1,XHK,XHKM1,XOHK,XOHKM1) = (XHK+XHKM1)/ZTHREE
     +  + XOHK*XOHK*(XTKP1+XTK)
     +  + XOHKM1*(ZTWO*XTK*XOHK+(XTK+XTKM1)*XOHKM1)
      FAKKP1(XTKP1,XTK,XHK,XOHKP1,XOHK,XOHKM1) = XHK/zsix
     +  - XOHK*(XTKP1*XOHKP1+(XTK+XTKP1)*XOHK + XTK*XOHKM1)
      FAKKP2(XTKP1,XOHKP1,XOHK) = XTKP1*XOHKP1*XOHK
C
      FRHS(XYKP1,XYK,XYKM1,XOHK,XOHKM1) = (XYKP1-XYK)*XOHK
     -  - (XYK-XYKM1)*XOHKM1
C
C     WHEN ONE WANTS AN ARRAY FOR TAU*SIGMA_K**2, THEN ONE SHOULD REPLACE
C     THE FUNCTION FTAU BY AN ARRAY
C
      fun_FTAUK(IIK)= ZTAUEFF*EXP(-ZCOFEXP*(PXIN(IIK)-ZXEXP0)**2
     +  /ZXEXPDL**2)
c%OS      FTAUK(IIK) = ZTAUEFF
C
c.......................................................................
C*COMDECK CUCCCC
C ----------------------------------------------------------------------
C --     STATEMENT FUNCTION FOR CUBIC INTERPOLATION                   --
C --                         23.04.88            AR        CRPP       --
C --                                                                  --
C -- CUBIC INTERPOLATION OF A FUNCTION F(X)                           --
C -- THE EIGHT ARGUMENTS A1,A2,A3,A4,B1,B2,B3,B4 ARE DEFINED BY:      --
C -- F(B1) = A1 , F(B2) = A2 , F(B3) = A3 , F(B4) = A4                --
C ----------------------------------------------------------------------
C
         FA3(A1,A2,A3,A4,B1,B2,B3,B4) =
     F        (A1-A2) / ((B1-B2)*(B2-B4)*(B2-B3)) +
     F        (A1-A3) / ((B4-B3)*(B3-B1)*(B3-B2)) +
     F        (A1-A4) / ((B1-B4)*(B2-B4)*(B3-B4))
         FA2(A1,A2,A3,A4,B1,B2,B3,B4) =
     F        (A1-A2) / ((B2-B1)*(B3-B2)) +
     F        (A3-A1) / ((B3-B1)*(B3-B2)) -
     F        (B1+B2+B3) * FA3(A1,A2,A3,A4,B1,B2,B3,B4)
         FA1(A1,A2,A3,A4,B1,B2,B3,B4) =
     F        (A1-A2) / (B1-B2) -
     F        (B1+B2) * FA2(A1,A2,A3,A4,B1,B2,B3,B4) -
     F        (B1*B1+B1*B2+B2*B2) * FA3(A1,A2,A3,A4,B1,B2,B3,B4)
         FA0(A1,A2,A3,A4,B1,B2,B3,B4) =
     F        A1 -
     F        B1 * (FA1(A1,A2,A3,A4,B1,B2,B3,B4) +
     F              B1 * (FA2(A1,A2,A3,A4,B1,B2,B3,B4) +
     F                    B1 * FA3(A1,A2,A3,A4,B1,B2,B3,B4)))
C ----------------------------------------------------------------------
C -- FCCCC0 GIVES THE VALUE OF THE FUNCTION AT POINT PX:              --
C -- FCCCC0(......,PX) = F(PX)                                        --
C ----------------------------------------------------------------------
        FCCCC0(A1,A2,A3,A4,B1,B2,B3,B4,PX) =
     F              FA0(A1,A2,A3,A4,B1,B2,B3,B4) +
     F              PX * (FA1(A1,A2,A3,A4,B1,B2,B3,B4) +
     F                    PX * (FA2(A1,A2,A3,A4,B1,B2,B3,B4) +
     F                          PX * FA3(A1,A2,A3,A4,B1,B2,B3,B4)))
C ----------------------------------------------------------------------
C -- FCCCC1 GIVES THE VALUE OF THE FIRST DERIVATIVE OF F(X) AT PX:    --
C -- FCCCC1(......,PX) = DF/DX (PX)                                   --
C ----------------------------------------------------------------------
        FCCCC1(A1,A2,A3,A4,B1,B2,B3,B4,PX) =
     F              FA1(A1,A2,A3,A4,B1,B2,B3,B4) +
     F              PX * (ZTWO * FA2(A1,A2,A3,A4,B1,B2,B3,B4) +
     F                    ZTHREE * PX * FA3(A1,A2,A3,A4,B1,B2,B3,B4))
C ----------------------------------------------------------------------
C -- FCCCC2 GIVES THE VALUE OF THE SECOND DERIVATIVE OF F(X) AT PX:   --
C -- FCCCC2(......,PX) = D2F/DX2 (PX)                                 --
C ----------------------------------------------------------------------
         FCCCC2(A1,A2,A3,A4,B1,B2,B3,B4,PX) =
     F             ZTWO * FA2(A1,A2,A3,A4,B1,B2,B3,B4) +
     F             zsix * FA3(A1,A2,A3,A4,B1,B2,B3,B4) * PX
C ----------------------------------------------------------------------
C -- FCCCC3 GIVES THE VALUE OF THE THIRD DERIVATIVE OF F(X) AT PX:     -
C -- FCCCC3(......,PX) = D3F/DX3 (PX)                                  -
C ----------------------------------------------------------------------
         FCCCC3(A1,A2,A3,A4,B1,B2,B3,B4,PX) =
     F                      zsix * FA3(A1,A2,A3,A4,B1,B2,B3,B4)
c.......................................................................
C
C-----------------------------------------------------------------------
C
C     0. INITIALIZATION
C
      iptr_ftauk = malloc_f(LENREAL*knin,infomalloc)
      if (infomalloc .ne. 0) return
      ITAUVAL = 0
      IF (PSIG(1) .NE. 0.) ITAUVAL = 1
      ZTAUEFF = abs(PSIG(1))
      ZXEXP0 = PXEXP0
      ZXEXPDL = PXEXPDL
      ZCOFEXP = 1.0
      IF (ZXEXP0.LT.PXIN(1) .OR. ZXEXP0.GT.PXIN(KNIN)) ZCOFEXP=0.0
C
      ISYM = 1
      IF (NBCLFT.GE.10 .OR. NBCRGT.GE.10 .OR. NBCLFT.EQ.2
     +  .OR. NBCRGT.EQ.2) ISYM = 0
C
      N = KNIN
      DO I=1,N
        PYINPP(I) = 0.0
      END DO
      DO I=1,MDAMAT
        DO J=1,N
          PAMAT(I,J) = 0.0
        END DO
      END DO
C
C     0.2 PRE-COMPUTE H_K AND 1./H_K, and zftauk
C
      DO K=1,N-1
        WHK(K)  = (PXIN(K+1) - PXIN(K))
      enddo
      DO K=1,N-1
        WHK(K)  = (PXIN(K+1) - PXIN(K))
        WOHK(K) = zone / WHK(K)
        ftauk(k)=psig(k)
c%OS        ftauk(k)=fun_ftauk(k)
      END DO
      WHK(N) = 0.0
      WOHK(N) = 0.0
      ftauk(n)=psig(N)
c%OS      ftauk(n)=fun_ftauk(N)
      if (PSIG(1).lt.0.) ftauk(1)=-10.*ftauk(1)
C
C     0.3 PREPARE BAND WIDTH
C
      IUP   = 2
      IDOWN = 2
      IF (ITAUVAL .EQ. 0) IUP   = 1
      IF (ITAUVAL .EQ. 0) IDOWN = 1
      IDIAG = IUP + 1
      IF (ISYM .EQ. 0) THEN
        IF (NBCLFT .GE. 10) IUP = IUP + 1
        IF (NBCRGT .GE. 10) IDOWN = IDOWN + 1
        IDIAG = IUP + IDOWN + 1
      ENDIF
      IF (MDAMAT .LT. IUP+1+2*(ISYM-1)*IDOWN) THEN
        PRINT *,' MDAMAT= ',MDAMAT,' < IUP+1+2*(ISYM-1)*IDOWN= ',
     +    IUP+1+2*(ISYM-1)*IDOWN
c%OS        STOP 'PAMAT'
        print *,'PAMAT'
        RETURN
      ENDIF
C
C     0.4 DETERMINE NEIGHBOURS: K-2, K-1, .., K+2
C     WHEN OUT OF BOUNDS, POINT TO INDEX N, AS WHK, WOHK(N)=0.0
C
      DO ISHIFT=-2,+2
        DO K=1,N
          KPM2(K,ISHIFT) = K + ISHIFT
        END DO
      END DO
C     OUT OF INTERVAL: SEND TO N
      KPM2(1,-2)   = N
      KPM2(1,-1)   = N
      KPM2(2,-2)   = N
      KPM2(N-1,+2) = N
      KPM2(N  ,+1) = N
      KPM2(N  ,+2) = N
C
C     1. CONSTRUCT MATRIX AND R.H.S
C     LAPACK SET-UP OF MATRIX:    A(I,J) -> PAMAT(I-J+IDIAG, J)
C
c.......................................................................
C     AS MATRIX SYMMETRIC, COMPUTE ONLY UPPER PART, THAT IS IF J.GE.I
C
      DO K=1,N-1
        IEFF = K + IDIAG
        IKP2 = KPM2(K,+2)
        IKP1 = KPM2(K,+1)
        IKM1 = KPM2(K,-1)
        IKM2 = KPM2(K,-2)
C     A(K,K)
        JK = K
        PAMAT(IEFF-JK,JK) = FAKK(FTAUK(IKP1),FTAUK(K),FTAUK(IKM1),WHK(K)
     +    ,WHK(IKM1),WOHK(K),WOHK(IKM1))
C
C     A(K,K+1)
        JKP1 = K + 1
        PAMAT(IEFF-JKP1,JKP1) = FAKKP1(FTAUK(IKP1),FTAUK(K),WHK(K),
     +    WOHK(IKP1),WOHK(K),WOHK(IKM1))
C     A(K,K-1)
c%OS        JKM1 = K - 1
c%OS        PAMAT(IEFF-JKM1,JKM1) = FAKKM1(FTAUK(K),FTAUK(IKM1),WHK(IKM1),
c%OS     +    WOHK(K),WOHK(IKM1),WOHK(IKM2))
C
        IF (ITAUVAL .EQ. 1) THEN
C     A(K,K+2)
          JKP2 = K + 2
          IF (JKP2 .LE. N)
     +      PAMAT(IEFF-JKP2,JKP2) = FAKKP2(FTAUK(IKP1),WOHK(IKP1),
     +      WOHK(K))
C     A(K,K-2)
c%OS          JKM2 = K - 2
c%OS          PAMAT(IEFF-JKM2,JKM2) = FAKKM2(FTAUK(IKM1),WOHK(IKM1),
c%OS     +      WOHK(IKM2))
        ENDIF
C     B(K)
        PYINPP(K) = FRHS(PYIN(IKP1),PYIN(K),PYIN(IKM1),WOHK(K),
     +    WOHK(IKM1))
C
      END DO
C
C     2. BOUNDARY CONDITIONS
C
C     2.1 IF NON-SYMMETRIC, COPY TOP PART TO BOTTOM BEFORE APPLYING
C     B.C.
C
      IF (ISYM .EQ. 0) THEN
        DO I=1,N
          IEFF = I + IDIAG
          DO J=I+1,MIN(I+MIN(IUP,IDOWN),N)
            JEFF = J + IDIAG
C     A(J,I) = A(I,J)
            PAMAT(JEFF-I,I) = PAMAT(IEFF-J,J)
          END DO
        END DO
      ENDIF
c%OS
c     debug, print matrix and rhs
c%OS      write(6,'(3a4,a)') 'i ','j1 ','j2 ',' i,j1  i,j1+1,..... i,j2'
c%OS      do i=1,n
c%OS        ieff = idiag + i
c%OS        j1 = i-idown
c%OS        j2 = i+iup
c%OSc%OS        j1 = max(i-idown,1)
c%OSc%OS        j2 = min(i+iup,n)
c%OS        write(6,'(3i4,1p10e13.4)') i,j1,j2,(pamat(ieff-j,j),j=j1,j2)
c%OS      end do
c%OS      write(6,'(a4,a12)') 'i','RHS'
c%OS      write(6,'(i4,1pe13.4)') (i,pyinpp(i),i=1,n)
c
c%OS
C
C     2.2 B.C. AT TWO LOCATIONS PXIN(IXBC(JBC)), JBC=1,2
C     IBCTYP(JBC) = 0, 1 OR 2 (TYPE OF B.C, SEE ABOVE).
C     SO FAR USES NODE CLOSEST TO XBCLFT/RGT FOR LOCATION
C     OF B.C., INSTEAD OF ACTUAL VALUE OF XBCLFT/RGT
C
      IXBC(1) = 1
      IXBC(2) = N
      IF (NBCLFT .GE. 10) THEN
        DO I=1,KNIN
          IF (PXIN(I) .GE. XBCLFT) GO TO 220
        END DO
 220    CONTINUE
        ZDELX = ABS(PXIN(I)-XBCLFT)
        IXBC(1) = I
        IF (I .GE. N) THEN
          IXBC(1) = N
          PRINT *,' WARNING: LEFT B.C. AT I=N: XBCLFT=',XBCLFT,
     +      '  PXIN(N)= ',PXIN(N)
        ELSE IF (ABS(PXIN(I-1)-XBCLFT).LE.ZDELX .AND. I.NE.1) THEN
          IXBC(1) = I-1
        ENDIF
      ENDIF
C
      IF (NBCRGT .GE. 10) THEN
        DO I=1,KNIN
          IF (PXIN(I) .GE. XBCRGT) GO TO 221
        END DO
 221    CONTINUE
        ZDELX = ABS(PXIN(I)-XBCRGT)
        IXBC(2) = I
        IF (I .LE. 1) THEN
          IXBC(2) = 1
          PRINT *,' WARNING: RIGHT B.C. AT I=1: XBCRGT=',XBCRGT,
     +      '  PXIN(1)= ',PXIN(1)
        ELSE IF (I .GT. N) THEN
          IXBC(2) = N
        ELSE IF (ABS(PXIN(I-1)-XBCRGT) .LE. ZDELX) THEN
          IXBC(2) = I-1
        ENDIF
      ENDIF
C
      ZYBC(1) = YBCLFT
      ZYBC(2) = YBCRGT
      IBCTYP(1) = MOD(NBCLFT,10)
      IBCTYP(2) = MOD(NBCRGT,10)
      IF (IXBC(1) .EQ. IXBC(2)) THEN
        PRINT *,' ERROR, B.C. AT SAME LOCATIONS: IXBC(1)=IXBC(2)= ',
     +    IXBC(1)
c%OS        STOP '1=2'
        RETURN
      ELSE IF (IXBC(1) .GT. IXBC(2)) THEN
        PRINT *,' WARNING, NEEDED TO SWITCH B.C. POINTS AS IXBC(1)= ',
     +    IXBC(1),' > IXBC(2)= ',IXBC(2)
        III = IXBC(1)
        IXBC(1) = IXBC(2)
        IXBC(2) = III
        ZYBC(1) = YBCRGT
        ZYBC(2) = YBCLFT
        IBCTYP(1) = MOD(NBCRGT,10)
        IBCTYP(2) = MOD(NBCLFT,10)
      ENDIF
C
C     2.3 MOVE EQUATIONS UP OR DOWN IF B.C. IS NOT AN END POINT
C
      IF (IXBC(1) .NE. 1) THEN
C
C     MOVE ROW EQ. K=2,..,IXBC(1) UP BY ONE
        IUPSOFAR = IUP - 1
        DO K=2,IXBC(1)
          DO J=MAX(1,K-IDOWN),MIN(N,K+IUPSOFAR)
            PAMAT(IDIAG+(K-1)-J,J) = PAMAT(IDIAG+K-J,J)
          END DO
          PYINPP(K-1) = PYINPP(K)
C     ZERO A((K-1),(K-1)-IDOWN)
          IF (K-1-IDOWN .GE. 1) PAMAT(IDIAG+IDOWN,K-1-IDOWN) = 0.0
        END DO
C     ZERO ROW IXBC(1) AND RHS
        K = IXBC(1)
        DO J=MAX(1,K-IDOWN),MIN(N,K+IUP)
          PAMAT(IDIAG+K-J,J) = 0.0
        END DO
        PYINPP(K) = 0.0
      ENDIF
C
      IF (IXBC(2) .NE. N) THEN
C     
C     MOVE EQ. K=IXBC(2),..,N-1 DOWN BY ONE
        IDWNSOFA = IDOWN - 1
        DO K=N-1,IXBC(2),-1
          DO J=MAX(1,K-IDWNSOFA),MIN(N,K+IUP)
            PAMAT(IDIAG+(K+1)-J,J) = PAMAT(IDIAG+K-J,J)
          END DO
          PYINPP(K+1) = PYINPP(K)
C     ZERO A((K+1),(K+1)+IUP)
          IF (K+1+IUP .LE. N) PAMAT(IDIAG-IUP,K+1+IUP) = 0.0
        END DO
C     ZERO ROW IXBC(2) AND RHS
        K = IXBC(2)
        DO J=MAX(1,K-IDOWN),MIN(N,K+IUP)
          PAMAT(IDIAG+K-J,J) = 0.0
        END DO
        PYINPP(K) = 0.0
      ENDIF
C
C     2.4 FOR ROW=IXBC(), MODIFY MATRIX AND RHS ACCORDING TO B.C. TYPE
C
      ZERO = 0.0
      DO JBC=1,2
        IK = IXBC(JBC)
        ZVALUE = ZYBC(JBC)
        IEFF = IK + IDIAG
        IKP2 = KPM2(IK,+2)
        IKP1 = KPM2(IK,+1)
        IKM1 = KPM2(IK,-1)
        IKM2 = KPM2(IK,-2)
        IF (IBCTYP(JBC) .EQ. 0) THEN
C
C     SYMMETRIZE => COL IK GOES TO RIGHT-HAND SIDE AND THEN ZEROED
C
          IF (ISYM .EQ. 1) THEN
            IDIAMIK = IDIAG - IK
            DO I=MAX(1,IK-IUP),IK-1
              PYINPP(I) = PYINPP(I) - ZVALUE * PAMAT(I+IDIAMIK,IK)
              PAMAT(I+IDIAMIK,IK) = 0.0
            END DO
            IDIAPIK = IDIAG + IK
            DO I=IK+1,MIN(N,IK+IUP)
              PYINPP(I) = PYINPP(I) - ZVALUE * PAMAT(IDIAPIK-I,I)
              PAMAT(IDIAPIK-I,I) = 0.0
            END DO
          ELSE
            IDIAMIK = IDIAG - IK
            DO I=MAX(1,IK-IUP),MIN(N,IK+IDOWN)
              PYINPP(I) = PYINPP(I) - ZVALUE * PAMAT(I+IDIAMIK,IK)
              PAMAT(I+IDIAMIK,IK) = 0.0
            END DO
C     ZERO ROW IK
            DO J=MAX(1,IK-IDOWN),MIN(N,IK+IUP)
              PAMAT(IEFF-J,J) = 0.0
            END DO
          ENDIF
C
C     REPLACE ROW IK BY EQUATION: G_K = ZVALUE
          PAMAT(IDIAG,IK) = 1.0
          PYINPP(IK) = ZVALUE
C
        ELSE IF (IBCTYP(JBC) .EQ. 1) THEN
C
C     1ST DERIVATIVE GIVEN
C
          ZYPEFF = ZVALUE
          IF (ZVALUE .GT. 1.E+31) THEN
C     FROM LGRANGIAN INTERPOLATION
            IKLFT = IK - 1
            IF (IK .EQ. 1) IKLFT = IK
            IF (IKLFT+3 .GT. N) IKLFT = N - 3
            ZYPEFF = FCCCC1(PYIN(IKLFT),PYIN(IKLFT+1),PYIN(IKLFT+2),
     +        PYIN(IKLFT+3),PXIN(IKLFT),PXIN(IKLFT+1),PXIN(IKLFT+2),
     +        PXIN(IKLFT+3),PXIN(IK))
          ELSE IF (ZVALUE .LT. -1.E+31) THEN
            IKLFT = IK
            IF (IK .EQ. N) IKLFT = IK - 1
            print *,PXIN(IKLFT+1)-PXIN(IKLFT)
            ZYPEFF = (PYIN(IKLFT+1)-PYIN(IKLFT))
     +        / (PXIN(IKLFT+1)-PXIN(IKLFT))
          ENDIF
          ZTOHKK1 = FTAUK(IK)*WOHK(IK)*WOHK(IKM1)
C     A(IK,IK)
          IF (IK .NE. N) PAMAT(IEFF-IK,IK) = FAKK(FTAUK(IKP1),FTAUK(IK),
     +      ZERO,WHK(IK),ZERO,WOHK(IK),ZERO) + ZTOHKK1
          IF (IK .EQ. N) PAMAT(IEFF-IK,IK) = FAKK(ZERO,FTAUK(IK),
     +      FTAUK(IKM1),ZERO,WHK(IKM1),ZERO,WOHK(IKM1)) + ZTOHKK1
C     A(IK,IK-1)
          JKM1 = IK - 1
          IF (ISYM.EQ.0 .AND. JKM1.GE.1) THEN
            IF (IK .NE. N) PAMAT(IEFF-JKM1,JKM1) = - ZTOHKK1
            IF (IK .EQ. N) PAMAT(IEFF-JKM1,JKM1) = FAKKM1(FTAUK(IK),
     +        FTAUK(IKM1),WHK(IKM1),ZERO,WOHK(IKM1),WOHK(IKM2))
          ENDIF
C     A(IK,IK+1)
          JKP1 = IK + 1
          IF (JKP1 .LE. N)
     +      PAMAT(IEFF-JKP1,JKP1) = FAKKP1(FTAUK(IKP1),
     +      FTAUK(IK),WHK(IK),WOHK(IKP1),WOHK(IK),ZERO)
C
          IF (ITAUVAL .EQ. 1) THEN
C     A(IK,IK+2)
            JKP2 = IK + 2
            IF (JKP2 .LE. N)
     +        PAMAT(IEFF-JKP2,JKP2) = FAKKP2(FTAUK(IKP1),WOHK(IKP1),
     +        WOHK(IK))
C     A(IK,IK-2)
            JKM2 = IK - 2
            IF (ISYM.EQ.0 .AND. JKM2.GE.1) THEN
              IF (IK .NE. N) PAMAT(IEFF-JKM2,JKM2) = 0.0
              IF (IK .EQ. N) PAMAT(IEFF-JKM2,JKM2) = FAKKM2(FTAUK(IKM1),
     +          WOHK(IKM1),WOHK(IKM2))
            ENDIF
          ENDIF
C     RHS
          ZSIGN = -1.
          IF (IK .EQ. N) ZSIGN = +1.
          IF (IK .NE. N) PYINPP(IK) = FRHS(PYIN(IKP1),PYIN(IK),ZERO,
     +      WOHK(IK),ZERO) - ZYPEFF
          IF (IK .EQ. N) PYINPP(IK) = FRHS(ZERO,PYIN(IK),PYIN(IKM1),
     +      ZERO,WOHK(IKM1)) + ZYPEFF
C
        ELSE IF (IBCTYP(JBC) .EQ. 2) THEN
C
C     FUNCTION IS GIVEN
C
C     A(IK,IK)
          PAMAT(IEFF-IK,IK) = - FTAUK(IK) * (WOHK(IK) + WOHK(IKM1))
C     A(IK,IK+1)
          JKP1 = IK + 1
          IF (JKP1 .LE. N)
     +      PAMAT(IEFF-JKP1,JKP1) = FTAUK(IK) * WOHK(IK)
C     A(IK,IK-1)
          JKM1 = IK - 1
          IF (ISYM.EQ.0 .AND. JKM1.GE.1)
     +      PAMAT(IEFF-JKM1,JKM1) = FTAUK(IK) * WOHK(IKM1)
C
          IF (ITAUVAL .EQ. 1) THEN
C     A(IK,IK+2)
            JKP2 = IK + 2
            IF (JKP2 .LE. N) PAMAT(IEFF-JKP2,JKP2) = 0.0
C     A(IK,IK-2)
            JKM2 = IK - 2
            IF (ISYM.EQ.0 .AND. JKM2.GE.1) PAMAT(IEFF-JKM2,JKM2) = 0.0
          ENDIF
C     RHS
          PYINPP(IK) = PYIN(IK) - ZVALUE
C
        ENDIF
C
      END DO
C
C     3. SOLVE SYSTEM
C
C     USE INTEGER WORK SPACE FROM KPM2(0) ARRAY FOR IPIVOT,
C     AS KPM2(K,0) NOT NEEDED NOR USED
C
      IDIMA = MDAMAT
      IDIMRHS = N
      IRHS = 1
c%OS
c     debug, print matrix and rhs
c%OS      write(6,'(3a4,a)') 'i ','j1 ','j2 ',' i,j1  i,j1+1,..... i,j2'
c%OS      do i=1,n
c%OS        ieff = idiag + i
c%OS        j1 = i-idown
c%OS        j2 = i+iup
c%OSc%OS        j1 = max(i-idown,1)
c%OSc%OS        j2 = min(i+iup,n)
c%OS        write(6,'(3i4,1p10e13.4)') i,j1,j2,(pamat(ieff-j,j),j=j1,j2)
c%OS      end do
c%OS      write(6,'(a4,a12)') 'i','RHS'
c%OS      write(6,'(i4,1pe13.4)') (i,pyinpp(i),i=1,n)
c
c%OS
      IF (ISYM .EQ. 1) THEN
        CALL SPBTRF('U',N,IUP,PAMAT,IDIMA,INFO)
      ELSE
        CALL SGBTRF(N,N,IDOWN,IUP,PAMAT,IDIMA,KPM2(1,0),INFO)
      ENDIF
      IF (INFO .EQ. 0) THEN
        IF (ISYM .EQ. 1) THEN
          CALL SPBTRS('U',N,IUP,IRHS,PAMAT,IDIMA,PYINPP,IDIMRHS,INFO2)
        ELSE
          CALL SGBTRS('N',N,IDOWN,IUP,IRHS,PAMAT,IDIMA,KPM2(1,0),PYINPP,
     +      IDIMRHS,INFO2)
        ENDIF
      ELSE
        PRINT *,' ERROR IN SP/GBTRF: INFO = ',INFO
c%OS        STOP 'INFO'
        RETURN
      ENDIF
C
C     4. COMPUTE NEW VALUES OF Y_K (NON-STANDARD CUBIC SPLINE ONLY)
C
      IF (ITAUVAL .EQ. 1) THEN
        DO K=1,N
          IKP1 = KPM2(K,+1)
          IKM1 = KPM2(K,-1)
          PYINNEW(K) = PYIN(K) - FTAUK(K) *
     +      ((PYINPP(IKP1)-PYINPP(K))*WOHK(K)
     +      - (PYINPP(K)-PYINPP(IKM1))*WOHK(IKM1))
        END DO
C
      ENDIF

      IF (INFO2 .LT. 0) THEN
        PRINT *,' ERROR IN SP/GBTRS: INFO2 = ',INFO2
c%OS        STOP 'INFO2'
        RETURN
      ENDIF
C
      call free(iptr_ftauk)
C
      RETURN
      END

C-----------------------------------------------------------------------
      SUBROUTINE SPLIBND(XA,YA,Y2A,N,X,Y,YP,YPP,KOPTIN)
      implicit none
      real*8 zsix, zthree, ztwo, zone
      PARAMETER(zsix=6., zthree=3., ztwo=2., zone=1.)
      real*8 XA, YA, Y2A, X, Y, YP, YPP
      integer N, KOPTIN
      DIMENSION XA(N),YA(N),Y2A(N)
C
C   ABS(KOPTIN):
C     KOPTIN = 0: STOP WITH ERROR MESSAGE IF OUT OF BOUND
C     KOPTIN = 1: LINEAR EXTRAPOLATION
C     KOPTIN = 2: USE QUADRATIC INTERPOLATION IF X OUT OF BOUND
C     KOPTIN = 3: USE CUBIC INTERPOLATION IF X OUT OF BOUND
C     KOPTIN = 21: USE QUADRATIC WITHIN ALFA*DELTA_X AND LINEAR FURTHER
C     KOPTIN = 31: USE CUBIC WITHIN ALFA*DELTA_X AND LINEAR    FURTHER
C     KOPTIN = 32: USE CUBIC WITHIN ALFA*DELTA_X AND QUADRATIC FURTHER
C
C   KOPTIN >=0 => INCONTDER = 1
C   KOPTIN < 0 => INCONTDER = 0
C     ICONTDER = 1: VALUE AND 1ST DER. CONTINUOUS AT END OF INTERVAL, THUS
C     .             USES CUBIC SPLINE OF LAST INTERVAL TO CONTINUE
C     ICONTDER = 0: ONLY VALUE CONTINUOUS AND USES VALUES AT LAST BUT ONE,
C     .             TWO, THREE POINTS TO EXTRAPOLATE (BETTER IF DER. AT EDGE
C     .             IS WILD)
C
C-----------------------------------------------------------------------
C
      real*8 ALFA
      PARAMETER(ALFA = 0.)
C
c.......................................................................
C*COMDECK CUCCCC
C ----------------------------------------------------------------------
C --     STATEMENT FUNCTION FOR CUBIC INTERPOLATION                   --
C --                         23.04.88            AR        CRPP       --
C --                                                                  --
C -- CUBIC INTERPOLATION OF A FUNCTION F(X)                           --
C -- THE EIGHT ARGUMENTS A1,A2,A3,A4,B1,B2,B3,B4 ARE DEFINED BY:      --
C -- F(B1) = A1 , F(B2) = A2 , F(B3) = A3 , F(B4) = A4                --
C ----------------------------------------------------------------------
C
      real*8 h, a, b, dellim, zxedg, zyedg, zypedg, zxdellim, zyxdl,
     +  zypxdl, zxedgm1, zyedgm1, zyxdlp, zypedgm1
      integer icontder, kopt, klo, khi, k, ikopt, k1, k2, k3, k4, klohi
      real*8 fc3, x1, f1, p1, x2,
     +  f2, p2, fc2, fc1, fc0, fqqq0, fqqq1, fqqq2,
     +  flinear, flinearp, fcccc0, fcccc1, fcccc2, fcccc3, fqdq0, fqdq1,
     +  fqdq2, fcdcd0, fcdcd1, fcdcd2, fcdcd3, fb1,
     +  fb2, fa2, fa3, fd2, fd1
      real*8 a1, a2, a3, a4, b1, b2, b3, b4, px
      real*8 fb0, fd0, fa0, fa1
         FA3(A1,A2,A3,A4,B1,B2,B3,B4) =
     F        (A1-A2) / ((B1-B2)*(B2-B4)*(B2-B3)) +
     F        (A1-A3) / ((B4-B3)*(B3-B1)*(B3-B2)) +
     F        (A1-A4) / ((B1-B4)*(B2-B4)*(B3-B4))
         FA2(A1,A2,A3,A4,B1,B2,B3,B4) =
     F        (A1-A2) / ((B2-B1)*(B3-B2)) +
     F        (A3-A1) / ((B3-B1)*(B3-B2)) -
     F        (B1+B2+B3) * FA3(A1,A2,A3,A4,B1,B2,B3,B4)
         FA1(A1,A2,A3,A4,B1,B2,B3,B4) =
     F        (A1-A2) / (B1-B2) -
     F        (B1+B2) * FA2(A1,A2,A3,A4,B1,B2,B3,B4) -
     F        (B1*B1+B1*B2+B2*B2) * FA3(A1,A2,A3,A4,B1,B2,B3,B4)
         FA0(A1,A2,A3,A4,B1,B2,B3,B4) =
     F        A1 -
     F        B1 * (FA1(A1,A2,A3,A4,B1,B2,B3,B4) +
     F              B1 * (FA2(A1,A2,A3,A4,B1,B2,B3,B4) +
     F                    B1 * FA3(A1,A2,A3,A4,B1,B2,B3,B4)))
C ----------------------------------------------------------------------
C -- FCCCC0 GIVES THE VALUE OF THE FUNCTION AT POINT PX:              --
C -- FCCCC0(......,PX) = F(PX)                                        --
C ----------------------------------------------------------------------
        FCCCC0(A1,A2,A3,A4,B1,B2,B3,B4,PX) =
     F              FA0(A1,A2,A3,A4,B1,B2,B3,B4) +
     F              PX * (FA1(A1,A2,A3,A4,B1,B2,B3,B4) +
     F                    PX * (FA2(A1,A2,A3,A4,B1,B2,B3,B4) +
     F                          PX * FA3(A1,A2,A3,A4,B1,B2,B3,B4)))
C ----------------------------------------------------------------------
C -- FCCCC1 GIVES THE VALUE OF THE FIRST DERIVATIVE OF F(X) AT PX:    --
C -- FCCCC1(......,PX) = DF/DX (PX)                                   --
C ----------------------------------------------------------------------
        FCCCC1(A1,A2,A3,A4,B1,B2,B3,B4,PX) =
     F              FA1(A1,A2,A3,A4,B1,B2,B3,B4) +
     F              PX * (ZTWO * FA2(A1,A2,A3,A4,B1,B2,B3,B4) +
     F                    ZTHREE * PX * FA3(A1,A2,A3,A4,B1,B2,B3,B4))
C ----------------------------------------------------------------------
C -- FCCCC2 GIVES THE VALUE OF THE SECOND DERIVATIVE OF F(X) AT PX:   --
C -- FCCCC2(......,PX) = D2F/DX2 (PX)                                 --
C ----------------------------------------------------------------------
         FCCCC2(A1,A2,A3,A4,B1,B2,B3,B4,PX) =
     F             ZTWO * FA2(A1,A2,A3,A4,B1,B2,B3,B4) +
     F             zsix * FA3(A1,A2,A3,A4,B1,B2,B3,B4) * PX
C ----------------------------------------------------------------------
C -- FCCCC3 GIVES THE VALUE OF THE THIRD DERIVATIVE OF F(X) AT PX:     -
C -- FCCCC3(......,PX) = D3F/DX3 (PX)                                  -
C ----------------------------------------------------------------------
         FCCCC3(A1,A2,A3,A4,B1,B2,B3,B4,PX) =
     F                      zsix * FA3(A1,A2,A3,A4,B1,B2,B3,B4)
C-----------------------------------------------------------------------
c.......................................................................
C*COMDECK CUCDCD
C ----------------------------------------------------------------------
C --     STATEMENT FUNCTION FOR CUBIC INTERPOLATION                   --
C --                         19.01.87            AR        CRPP       --
C --                                                                  --
C -- CUBIC INTERPOLATION OF A FUNCTION F(X)                           --
C -- THE SIX ARGUMENTS X1,F1,P1,X2,F2,P2 ARE DEFINED AS FOLLOWS:      --
C -- F(X1) = F1 , F(X2) = F2 , DF/DX(X1) = P1 , DF/DX(X2) = P2        --
C ----------------------------------------------------------------------
C
         FC3(X1,F1,P1,X2,F2,P2) =
     =      (ZTWO * (F2 - F1) / (X1 - X2) + (P1 + P2)) / 
     /      ((X1 - X2) * (X1 - X2))
         FC2(X1,F1,P1,X2,F2,P2) =
     =      (ZTHREE * (X1 + X2) * (F1 - F2) / (X1 - X2) - 
     *       P1 * (X1 + ZTWO * X2) - P2 * (X2 + ZTWO * X1)) /
     /      ((X1 - X2) * (X1 - X2))
         FC1(X1,F1,P1,X2,F2,P2) =
     =      (zsix * X1 * X2 * (F2 - F1) / (X1 - X2) + 
     *       X2 * P1 * (2 * X1 + X2) + X1 * P2 * (X1 + ZTWO * X2)) /
     /      ((X1 - X2) * (X1 - X2))
         FC0(X1,F1,P1,X2,F2,P2) =
     =      (F1 * X2**2 + F2 * X1**2 - X1 * X2 * (X2 * P1 + X1 * P2) +
     +       ZTWO * X1 * X2 * (F1 * X2 - F2 * X1) / (X1 - X2)) /
     /      ((X1 - X2) * (X1 - X2))
C ----------------------------------------------------------------------
C -- FCDCD0 GIVES THE VALUE OF THE FUNCTION AT POINT PX               --
C -- FCDCD0(......,PX) = F(PX)                                        --
C ----------------------------------------------------------------------
         FCDCD0(X1,F1,P1,X2,F2,P2,PX) =
     =              FC0(X1,F1,P1,X2,F2,P2) +
     +              PX * (FC1(X1,F1,P1,X2,F2,P2) +
     +                    PX * (FC2(X1,F1,P1,X2,F2,P2) +
     +                          PX * FC3(X1,F1,P1,X2,F2,P2)))
C ----------------------------------------------------------------------
C -- FCDCD1 GIVES THE VALUE OF THE FIRST DERIVATIVE OF F(X) AT PX:    --
C -- FCDCD1(......,PX) = DF/DX (PX)                                   --
C ----------------------------------------------------------------------
         FCDCD1(X1,F1,P1,X2,F2,P2,PX) =
     =              FC1(X1,F1,P1,X2,F2,P2) +
     +              PX * (ZTWO * FC2(X1,F1,P1,X2,F2,P2) +
     +                    ZTHREE * PX * FC3(X1,F1,P1,X2,F2,P2))
C ----------------------------------------------------------------------
C -- FCDCD2 GIVES THE VALUE OF THE SECOND DERIVATIVE OF F(X) AT PX:   --
C -- FCDCD2(......,PX) = D2F/DX2 (PX)                                 --
C ----------------------------------------------------------------------
         FCDCD2(X1,F1,P1,X2,F2,P2,PX) =
     =             ZTWO * FC2(X1,F1,P1,X2,F2,P2) +
     +             zsix * FC3(X1,F1,P1,X2,F2,P2) * PX
C ----------------------------------------------------------------------
C -- FCDCD3 GIVES THE VALUE OF THE THIRD DERIVATIVE OF F(X) AT PX:    --
C -- FCDCD3(......,PX) = D3F/DX3 (PX)                                 --
C ----------------------------------------------------------------------
         FCDCD3(X1,F1,P1,X2,F2,P2,PX) =
     =                      zsix * FC3(X1,F1,P1,X2,F2,P2)
C
C.......................................................................
C*COMDECK QUAQQQ
C ----------------------------------------------------------------------
C --     STATEMENT FUNCTION FOR QUADRATIC INTERPOLATION               --
C --                         19.01.87            AR        CRPP       --
C --                                                                  --
C -- QUADRATIC INTERPOLATION OF A FUNCTION F(X)                       --
C -- THE SIX PARAMETERS A1,A2,A3,B1,B2,B3 ARE DEFINED AS FOLLOWS:     --
C -- F(B1) = A1 , F(B2) = A2 , F(B3) = A3                             --
C ----------------------------------------------------------------------
C
         FB2(A1,A2,A3,B1,B2,B3) =
     F               ((A1-A2)/(B1-B2)-(A1-A3)/(B1-B3))/(B2-B3)
         FB1(A1,A2,A3,B1,B2,B3) = ((A1-A2)/(B1-B2))-
     F         FB2(A1,A2,A3,B1,B2,B3)*(B1+B2)
         FB0(A1,A2,A3,B1,B2,B3) = A1-FB1(A1,A2,A3,B1,B2,B3)*B1
     F         -FB2(A1,A2,A3,B1,B2,B3)*B1*B1
C ----------------------------------------------------------------------
C -- FQQQ0 GIVES THE VALUE OF THE FUNCTION AT THE POINT PX            --
C -- FQQQ0(......,PX) = F(PX)                                         --
C ----------------------------------------------------------------------
         FQQQ0(A1,A2,A3,B1,B2,B3,PX) = FB0(A1,A2,A3,B1,B2,B3) +
     F                                 PX * (FB1(A1,A2,A3,B1,B2,B3) +
     F                                 PX * FB2(A1,A2,A3,B1,B2,B3))
C ----------------------------------------------------------------------
C -- FQQQ1 GIVES THE VALUE OF THE FIRST DERIVATIVE OF F(X) AT PX      --
C -- FQQQ1(......,PX) = DF/DX (PX)                                    --
C ----------------------------------------------------------------------
         FQQQ1(A1,A2,A3,B1,B2,B3,PX) = FB1(A1,A2,A3,B1,B2,B3) +
     F     ZTWO * PX * FB2(A1,A2,A3,B1,B2,B3)
C ----------------------------------------------------------------------
C -- FQQQ2 GIVES THE VALUE OF THE SECOND DERIVATIVE OF F(X) AT PX     --
C -- FQQQ2(......,PX) = D2F/DX2 (PX)                                  --
C ----------------------------------------------------------------------
         FQQQ2(A1,A2,A3,B1,B2,B3) = ZTWO * FB2(A1,A2,A3,B1,B2,B3)
c.......................................................................
C*COMDECK QUAQDQ
C ----------------------------------------------------------------------
C --     STATEMENT FUNCTION FOR QUADRATIC INTERPOLATION               --
C --                         19.01.87            AR        CRPP       --
C --                                                                  --
C -- QUADRATIC INTERPOLATION OF A FUNCTION F(X)                       --
C -- THE FIVE PARAMETERS X1,F1,P1,X2,F2    ARE DEFINED AS FOLLOWS:    --
C -- F(X1) = F1 , DF/DX(X1) = P1 , F(X2) = F2                         --
C ----------------------------------------------------------------------
C
         FD2(X1,F1,P1,X2,F2) = ((F2-F1)/(X2-X1) - P1) / (X2-X1)
         FD1(X1,F1,P1,X2,F2) = P1 - ZTWO*X1*FD2(X1,F1,P1,X2,F2)
         FD0(X1,F1,P1,X2,F2) = F1 - X1*(X1*FD2(X1,F1,P1,X2,F2) +
     +                                     FD1(X1,F1,P1,X2,F2))
C ----------------------------------------------------------------------
C -- FQDQ0 GIVES THE VALUE OF THE FUNCTION AT POINT PX                --
C -- FQDQ0(......,PX) = F(PX)                                         --
C ----------------------------------------------------------------------
         FQDQ0(X1,F1,P1,X2,F2,PX) = FD0(X1,F1,P1,X2,F2) +
     F                              PX * (FD1(X1,F1,P1,X2,F2) +
     F                                    PX * FD2(X1,F1,P1,X2,F2))
C ----------------------------------------------------------------------
C -- FQDQ1 GIVES THE VALUE OF THE FIRST DERIVATIVE OF F(X) AT PX      --
C -- FQDQ1(......,PX) = DF/DX (PX)                                    --
C ----------------------------------------------------------------------
         FQDQ1(X1,F1,P1,X2,F2,PX) = FD1(X1,F1,P1,X2,F2) +
     F                              ZTWO* PX * FD2(X1,F1,P1,X2,F2)
C ----------------------------------------------------------------------
C -- FQDQ2 GIVES THE VALUE OF THE SECOND DERIVATIVE OF F(X) AT PX     --
C -- FQDQ2(......,PX) = D2F/DX2 (PX)                                  --
C ----------------------------------------------------------------------
         FQDQ2(X1,F1,P1,X2,F2) = ZTWO * FD2(X1,F1,P1,X2,F2)
C-----------------------------------------------------------------------
c.......................................................................
C     LINEAR
C
      FLINEAR(X1,F1,X2,F2,PX) = F2 + (PX-X2)/(X2-X1) * (F2-F1)
      FLINEARP(X1,F1,X2,F2) = (F2-F1) / (X2-X1)
C-----------------------------------------------------------------------
      ICONTDER = 1
      IF (KOPTIN .LT. 0) ICONTDER = 0
      KOPT=ABS(KOPTIN)
C
C     1. POINT INSIDE INTERVAL
C
      IF (X.LT.XA(1) .OR. X.GT.XA(N)) GO TO 200
C
      KLO=1
      KHI=N
1     IF (KHI-KLO.GT.1) THEN
        K=(KHI+KLO)/2
        IF(XA(K).GT.X)THEN
          KHI=K
        ELSE
          KLO=K
        ENDIF
      GOTO 1
      ENDIF
      H=XA(KHI)-XA(KLO)
c%OS      IF (H.EQ.0.) STOP 'BAD XA INPUT.'
      IF (H .EQ. 0.) THEN
        PRINT *,'BAD XA INPUT.'
        RETURN
      ENDIF
      A=(XA(KHI)-X)/H
      B=(X-XA(KLO))/H
      Y=A*YA(KLO)+B*YA(KHI)+
     *      ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/zsix
      YP=(YA(KHI)-YA(KLO))/H -
     - ((ZTHREE*A*A-zone)*Y2A(KLO)-(ZTHREE*B*B-ZONE)*Y2A(KHI) )*H/zsix
      YPP=A*Y2A(KLO)+B*Y2A(KHI)
C
      RETURN
C
C     2. POINT OUTSIDE INTERVAL
C
 200  CONTINUE
C
      IKOPT = KOPT
      IF (KOPT .EQ. 0) THEN
        PRINT *,' POINT X=',X,' IS OUTSIDE INTERVAL [',XA(1),','
     +    ,XA(N),']'
        IKOPT = 1
CC        STOP 'KOPT=0'
      ENDIF
C
      KLO = 1
      KHI = 2
      IF (X .GT. XA(N)) THEN
        KLO = N-1
        KHI = N
      ENDIF
      H=XA(KHI)-XA(KLO)
      DELLIM = ALFA * H
c%OS      IF (H.EQ.0.) STOP 'BAD XA INPUT.'
      IF (H .EQ. 0.) THEN
        PRINT *,'BAD XA INPUT.'
        RETURN
      ENDIF
C
C.......................................................................
C     2.1 LINEAR, IKOPT=1
C
      IF (IKOPT .EQ. 1) THEN
C
C     LINEAR EXTRAPOLATION
        IF (ICONTDER .EQ. 0) THEN
          Y = YA(KHI)
          YP = 0.0
          YPP = 0.0
c%OS          Y = FLINEAR(XA(KLO),YA(KLO),XA(KHI),YA(KHI),X)
c%OS          YP = FLINEARP(XA(KLO),YA(KLO),XA(KHI),YA(KHI))
c%OS          YPP = 0.0
        ELSE
C     COMPUTE VALUE AND 1ST DER. AT EDGE
          IF (KLO .EQ. 1) THEN
            ZXEDG = XA(KLO)
            ZYEDG = YA(KLO)
            ZYPEDG= (YA(KHI)-YA(KLO))/H-(ZTWO*Y2A(KLO)+Y2A(KHI))*H/zsix
          ELSE
            ZXEDG = XA(KHI)
            ZYEDG = YA(KHI)
            ZYPEDG= (YA(KHI)-YA(KLO))/H+(Y2A(KLO)+ZTWO*Y2A(KHI))*H/zsix
          ENDIF
          Y = FLINEAR(ZXEDG,ZYEDG,ZXEDG+ZONE,ZYEDG+ZYPEDG,X)
          YP = ZYPEDG
          YPP = 0.0
        ENDIF
C
C.......................................................................
C     2.2 LINEAR FAR END OF IKOPT=21 AND IKOPT=31
C
      ELSE IF (MAX(XA(KLO)-X,X-XA(KHI)).GT.DELLIM .AND.
     +    (IKOPT.EQ.21 .OR. IKOPT.EQ.31)) THEN
C
C     LINEAR EXTRAPOLATION OUTSIDE DELLIM
C     COMPUTE STARTING POINT AND DERIVATIVE FROM END OF QUADR. OR CUBIC
C     INTERPOLATION IN ALFA*DELTA_X INTERVAL
C
        IF (IKOPT .EQ. 21) THEN
C     QUADRATIC
          IF (ICONTDER .EQ. 0) THEN
            K1 = 1
            IF (KHI .EQ. N) K1 = N-2
            K2 = K1 + 1
            K3 = K1 + 2
            ZXDELLIM = XA(1) - DELLIM
            IF (KHI .EQ. N) ZXDELLIM = XA(N) + DELLIM
            ZYXDL = FQQQ0(YA(K1),YA(K2),YA(K3),XA(K1),XA(K2),XA(K3),
     +        ZXDELLIM)
            ZYPXDL= FQQQ1(YA(K1),YA(K2),YA(K3),XA(K1),XA(K2),XA(K3),
     +        ZXDELLIM)
          ELSE
C     COMPUTE VALUE AND 1ST DER. AT EDGE
            IF (KLO .EQ. 1) THEN
              ZXEDGM1 = XA(KHI)
              ZYEDGM1 = YA(KHI)
              ZXEDG = XA(KLO)
              ZYEDG = YA(KLO)
              ZYPEDG=(YA(KHI)-YA(KLO))/H-(ZTWO*Y2A(KLO)+Y2A(KHI))*H/zsix
            ELSE
              ZXEDGM1 = XA(KLO)
              ZYEDGM1 = YA(KLO)
              ZXEDG = XA(KHI)
              ZYEDG = YA(KHI)
              ZYPEDG=(YA(KHI)-YA(KLO))/H+(Y2A(KLO)+ZTWO*Y2A(KHI))*H/zsix
            ENDIF
            ZXDELLIM = XA(1) - DELLIM
            IF (KHI .EQ. N) ZXDELLIM = XA(N) + DELLIM
            ZYXDL = FQDQ0(ZXEDG,ZYEDG,ZYPEDG,ZXEDGM1,ZYEDGM1,ZXDELLIM)
            ZYPXDL= FQDQ1(ZXEDG,ZYEDG,ZYPEDG,ZXEDGM1,ZYEDGM1,ZXDELLIM)
          ENDIF
        ELSE IF (IKOPT .EQ. 31) THEN
C     CUBIC
          IF (ICONTDER .EQ. 0) THEN
            K1 = 1
            IF (KHI .EQ. N) K1 = N-3
            K2 = K1 + 1
            K3 = K1 + 2
            K4 = K1 + 3
            ZXDELLIM = XA(1) - DELLIM
            IF (KHI .EQ. N) ZXDELLIM = XA(N) + DELLIM
            IF (K1.LT.1 .OR. K4.GT.N) THEN
              IF (KHI .EQ. N) K1 = N-1
              K2 = K1 + 1
              ZYXDL=FLINEAR(XA(K1),YA(K1),XA(K2),YA(K2),ZXDELLIM)
              ZYXDLP=FLINEARP(XA(K1),YA(K1),XA(K2),YA(K2))
            ELSE
              ZYXDL = FCCCC0(YA(K1),YA(K2),YA(K3),YA(K4),XA(K1),XA(K2),
     +          XA(K3),XA(K4),ZXDELLIM)
              ZYPXDL= FCCCC1(YA(K1),YA(K2),YA(K3),YA(K4),XA(K1),XA(K2),
     +          XA(K3),XA(K4),ZXDELLIM)
            ENDIF
          ELSE
C     COMPUTE VALUE AND 1ST DER. AT EDGE AND EDGE-1
            IF (KLO .EQ. 1) THEN
              ZXEDGM1 = XA(KHI)
              ZYEDGM1 = YA(KHI)
              ZYPEDGM1 = (YA(KHI)-YA(KLO))/H + (Y2A(KLO)+ZTWO*Y2A(KHI))
     +          *H/zsix
              ZXEDG = XA(KLO)
              ZYEDG = YA(KLO)
              ZYPEDG=(YA(KHI)-YA(KLO))/H-(ZTWO*Y2A(KLO)+Y2A(KHI))*H/zsix
            ELSE
c%OS              IF (H.LE.1.E-02) THEN
c%OS                KLO=KLO-1
c%OS                H=XA(KHI)-XA(KLO)
c%OS              endif
              ZXEDGM1 = XA(KLO)
              ZYEDGM1 = YA(KLO)
              ZYPEDGM1 =(YA(KHI)-YA(KLO))/H-(ZTWO*Y2A(KLO)+Y2A(KHI))
     +          *H/zsix
              ZXEDG = XA(KHI)
              ZYEDG = YA(KHI)
              ZYPEDG=(YA(KHI)-YA(KLO))/H+(Y2A(KLO)+ZTWO*Y2A(KHI))*H/zsix
            ENDIF
            ZXDELLIM = XA(1) - DELLIM
            IF (KHI .EQ. N) ZXDELLIM = XA(N) + DELLIM
            ZYXDL = FCDCD0(ZXEDGM1,ZYEDGM1,ZYPEDGM1,ZXEDG,ZYEDG,ZYPEDG,
     +        ZXDELLIM)
            ZYPXDL= FCDCD1(ZXEDGM1,ZYEDGM1,ZYPEDGM1,ZXEDG,ZYEDG,ZYPEDG,
     +        ZXDELLIM)
          ENDIF
        ENDIF
C
        Y = FLINEAR(ZXDELLIM,ZYXDL,ZXDELLIM+ZONE,ZYXDL+ZYPXDL,X)
        YP = ZYPXDL
        YPP = 0.0
C
C.......................................................................
C     2.3 QUADRATIC, IKOPT=2 OR FIRST PART OF IKOPT=21
C
      ELSE IF (IKOPT.EQ.2 .OR.
     +    (MAX(XA(KLO)-X,X-XA(KHI)).LE.DELLIM .AND. IKOPT.EQ.21)) THEN
C
C     QUADRATIC
        IF (ICONTDER .EQ. 0) THEN
          K1 = 1
          IF (KHI .EQ. N) K1 = N-2
          K2 = K1 + 1
          K3 = K1 + 2
          Y =  FQQQ0(YA(K1),YA(K2),YA(K3),XA(K1),XA(K2),XA(K3),X)
          YP = FQQQ1(YA(K1),YA(K2),YA(K3),XA(K1),XA(K2),XA(K3),X)
          YPP =FQQQ2(YA(K1),YA(K2),YA(K3),XA(K1),XA(K2),XA(K3))
        ELSE
C     COMPUTE VALUE AND 1ST DER. AT EDGE
          IF (KLO .EQ. 1) THEN
            ZXEDGM1 = XA(KHI)
            ZYEDGM1 = YA(KHI)
            ZXEDG = XA(KLO)
            ZYEDG = YA(KLO)
            ZYPEDG= (YA(KHI)-YA(KLO))/H-(ZTWO*Y2A(KLO)+Y2A(KHI))*H/zsix
          ELSE
            ZXEDGM1 = XA(KLO)
            ZYEDGM1 = YA(KLO)
            ZXEDG = XA(KHI)
            ZYEDG = YA(KHI)
            ZYPEDG= (YA(KHI)-YA(KLO))/H+(Y2A(KLO)+ZTWO*Y2A(KHI))*H/zsix
          ENDIF
          Y  = FQDQ0(ZXEDG,ZYEDG,ZYPEDG,ZXEDGM1,ZYEDGM1,X)
          YP = FQDQ1(ZXEDG,ZYEDG,ZYPEDG,ZXEDGM1,ZYEDGM1,X)
          YPP= FQDQ2(ZXEDG,ZYEDG,ZYPEDG,ZXEDGM1,ZYEDGM1)
        ENDIF
C
C.......................................................................
C     2.4 QUADRATIC, FAR END PART OF IKOPT=32
C
      ELSE IF (MAX(XA(KLO)-X,X-XA(KHI)).GT.DELLIM .AND. IKOPT.EQ.32)THEN
C
C     QUADRATIC FROM X+ALFA*DELTA_X
C     COMPUTE STARTING POINT AND DERIVATIVE FROM END OF CUBIC
C     WARNING: MUST BE COMPATIBLE WITH ALFA=0
        IF (ICONTDER .EQ. 0) THEN
          K1 = 1
          IF (KHI .EQ. N) K1 = N-3
          K2 = K1 + 1
          K3 = K1 + 2
          K4 = K1 + 3
          ZXDELLIM = XA(1) - DELLIM
          IF (KHI .EQ. N) ZXDELLIM = XA(N) + DELLIM
          IF (K1.LT.1 .OR. K4.GT.N) THEN
            IF (KHI .EQ. N) K1 = N-1
            K2 = K1 + 1
            ZYXDL=FLINEAR(XA(K1),YA(K1),XA(K2),YA(K2),ZXDELLIM)
            ZYXDLP=FLINEARP(XA(K1),YA(K1),XA(K2),YA(K2))
          ELSE
            ZYXDL = FCCCC0(YA(K1),YA(K2),YA(K3),YA(K4),XA(K1),XA(K2),
     +        XA(K3),XA(K4),ZXDELLIM)
            ZYPXDL= FCCCC1(YA(K1),YA(K2),YA(K3),YA(K4),XA(K1),XA(K2),
     +        XA(K3),XA(K4),ZXDELLIM)
          ENDIF
        ELSE
C     COMPUTE VALUE AND 1ST DER. AT EDGE AND EDGE-1
          IF (KLO .EQ. 1) THEN
            ZXEDGM1 = XA(KHI)
            ZYEDGM1 = YA(KHI)
            ZYPEDGM1 = (YA(KHI)-YA(KLO))/H + (Y2A(KLO)+ZTWO*Y2A(KHI))
     +        *H/zsix
            ZXEDG = XA(KLO)
            ZYEDG = YA(KLO)
            ZYPEDG= (YA(KHI)-YA(KLO))/H-(ZTWO*Y2A(KLO)+Y2A(KHI))*H/zsix
          ELSE
            ZXEDGM1 = XA(KLO)
            ZYEDGM1 = YA(KLO)
            ZYPEDGM1 = (YA(KHI)-YA(KLO))/H - (ZTWO*Y2A(KLO)+Y2A(KHI))
     +        *H/zsix
            ZXEDG = XA(KHI)
            ZYEDG = YA(KHI)
            ZYPEDG= (YA(KHI)-YA(KLO))/H+(Y2A(KLO)+ZTWO*Y2A(KHI))*H/zsix
          ENDIF
          ZXDELLIM = XA(1) - DELLIM
          IF (KHI .EQ. N) ZXDELLIM = XA(N) + DELLIM
          ZYXDL = FCDCD0(ZXEDGM1,ZYEDGM1,ZYPEDGM1,ZXEDG,ZYEDG,ZYPEDG,
     +      ZXDELLIM)
          ZYPXDL= FCDCD1(ZXEDGM1,ZYEDGM1,ZYPEDGM1,ZXEDG,ZYEDG,ZYPEDG,
     +      ZXDELLIM)
        ENDIF
C     QUADRATIC FROM END OF INTERVAL, END+DELLIM AND 1ST DER. AT END+DELLIM
        KLOHI = 1
        IF (KHI .EQ. N) KLOHI = N
        IF (ABS(ALFA/H).GT.1E-04) THEN
          Y  = FQDQ0(ZXDELLIM,ZYXDL,ZYPXDL,XA(KLOHI),YA(KLOHI),X)
          YP = FQDQ1(ZXDELLIM,ZYXDL,ZYPXDL,XA(KLOHI),YA(KLOHI),X)
          YPP =FQDQ2(ZXDELLIM,ZYXDL,ZYPXDL,XA(KLOHI),YA(KLOHI))
        ELSE
          Y  = FQDQ0(ZXDELLIM,ZYXDL,ZYPXDL,XA(KLO),YA(KLO),X)
          YP = FQDQ1(ZXDELLIM,ZYXDL,ZYPXDL,XA(KLO),YA(KLO),X)
          YPP =FQDQ2(ZXDELLIM,ZYXDL,ZYPXDL,XA(KLO),YA(KLO))
        ENDIF
C
C.......................................................................
C     2.5 CUBIC, IKOPT=3 OR FIRST PART OF IKOPT=31 AND IKOPT=32
C
      ELSE IF (IKOPT.EQ.3 .OR.
     +    (MAX(XA(KLO)-X,X-XA(KHI)).LE.DELLIM .AND.
     +    (IKOPT.EQ.32 .OR. IKOPT.EQ.31)) ) THEN
C
C     CUBIC
        IF (ICONTDER .EQ. 0) THEN
          K1 = 1
          IF (KHI .EQ. N) K1 = N-3
          K2 = K1 + 1
          K3 = K1 + 2
          K4 = K1 + 3
          IF (K1.LT.1 .OR. K4.GT.N) THEN
            IF (KHI .EQ. N) K1 = N-1
            K2 = K1 + 1
            Y =FLINEAR(XA(K1),YA(K1),XA(K2),YA(K2),X)
            YP=FLINEARP(XA(K1),YA(K1),XA(K2),YA(K2))
            YPP=0.
          ELSE
            Y   =FCCCC0(YA(K1),YA(K2),YA(K3),YA(K4),XA(K1),XA(K2),XA(K3)
     +        ,XA(K4),X)
            YP  =FCCCC1(YA(K1),YA(K2),YA(K3),YA(K4),XA(K1),XA(K2),XA(K3)
     +        ,XA(K4),X)
            YPP =FCCCC2(YA(K1),YA(K2),YA(K3),YA(K4),XA(K1),XA(K2),XA(K3)
     +        ,XA(K4),X)
          ENDIF
        ELSE
C     COMPUTE VALUE AND 1ST DER. AT EDGE AND EDGE-1
          IF (KLO .EQ. 1) THEN
            ZXEDGM1 = XA(KHI)
            ZYEDGM1 = YA(KHI)
            ZYPEDGM1=(YA(KHI)-YA(KLO))/H+(Y2A(KLO)+ZTWO*Y2A(KHI))*H/zsix
            ZXEDG = XA(KLO)
            ZYEDG = YA(KLO)
            ZYPEDG= (YA(KHI)-YA(KLO))/H-(ZTWO*Y2A(KLO)+Y2A(KHI))*H/zsix
          ELSE
c%OS            IF (H.LE.1.E-02) THEN
c%OS              KLO=KLO-1
c%OS              H=XA(KHI)-XA(KLO)
c%OS            endif
            ZXEDGM1 = XA(KLO)
            ZYEDGM1 = YA(KLO)
            ZYPEDGM1=(YA(KHI)-YA(KLO))/H-(ZTWO*Y2A(KLO)+Y2A(KHI))*H/zsix
            ZXEDG = XA(KHI)
            ZYEDG = YA(KHI)
            ZYPEDG= (YA(KHI)-YA(KLO))/H+(Y2A(KLO)+ZTWO*Y2A(KHI))*H/zsix
          ENDIF
          Y  = FCDCD0(ZXEDGM1,ZYEDGM1,ZYPEDGM1,ZXEDG,ZYEDG,ZYPEDG,X)
          YP = FCDCD1(ZXEDGM1,ZYEDGM1,ZYPEDGM1,ZXEDG,ZYEDG,ZYPEDG,X)
          YPP= FCDCD2(ZXEDGM1,ZYEDGM1,ZYPEDGM1,ZXEDG,ZYEDG,ZYPEDG,X)
        ENDIF
      ENDIF
C
      RETURN
      END
c
c
      SUBROUTINE SGBTRF( M, N, KL, KU, AB, LDAB, IPIV, INFO )
*
*  -- LAPACK routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      implicit none
      INTEGER            INFO, KL, KU, LDAB, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      REAL*8               AB( LDAB, * )
*     ..
*
*  Purpose
*  =======
*
*  SGBTRF computes an LU factorization of a real m-by-n band matrix A
*  using partial pivoting with row interchanges.
*
*  This is the blocked version of the algorithm, calling Level 3 BLAS.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  KL      (input) INTEGER
*          The number of subdiagonals within the band of A.  KL >= 0.
*
*  KU      (input) INTEGER
*          The number of superdiagonals within the band of A.  KU >= 0.
*
*  AB      (input/output) REAL array, dimension (LDAB,N)
*          On entry, the matrix A in band storage, in rows KL+1 to
*          2*KL+KU+1; rows 1 to KL of the array need not be set.
*          The j-th column of A is stored in the j-th column of the
*          array AB as follows:
*          AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
*
*          On exit, details of the factorization: U is stored as an
*          upper triangular band matrix with KL+KU superdiagonals in
*          rows 1 to KL+KU+1, and the multipliers used during the
*          factorization are stored in rows KL+KU+2 to 2*KL+KU+1.
*          See below for further details.
*
*  LDAB    (input) INTEGER
*          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
*
*  IPIV    (output) INTEGER array, dimension (min(M,N))
*          The pivot indices; for 1 <= i <= min(M,N), row i of the
*          matrix was interchanged with row IPIV(i).
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*          > 0: if INFO = +i, U(i,i) is exactly zero. The factorization
*               has been completed, but the factor U is exactly
*               singular, and division by zero will occur if it is used
*               to solve a system of equations.
*
*  Further Details
*  ===============
*
*  The band storage scheme is illustrated by the following example, when
*  M = N = 6, KL = 2, KU = 1:
*
*  On entry:                       On exit:
*
*      *    *    *    +    +    +       *    *    *   u14  u25  u36
*      *    *    +    +    +    +       *    *   u13  u24  u35  u46
*      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
*     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
*     a21  a32  a43  a54  a65   *      m21  m32  m43  m54  m65   *
*     a31  a42  a53  a64   *    *      m31  m42  m53  m64   *    *
*
*  Array elements marked * are not used by the routine; elements marked
*  + need not be set on entry, but are required by the routine to store
*  elements of U because of fill-in resulting from the row interchanges.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL*8               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
      INTEGER            NBMAX, LDWORK
      PARAMETER          ( NBMAX = 64, LDWORK = NBMAX+1 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, I2, I3, II, IP, J, J2, J3, JB, JJ, JM, JP,
     $                   JU, K2, KM, KV, NB, NW
      REAL*8               TEMP
*     ..
*     .. Local Arrays ..
      REAL*8               WORK13( LDWORK, NBMAX ),
     $                   WORK31( LDWORK, NBMAX )
*     ..
*     .. External Functions ..
      INTEGER            ILAENV, ISAMAX
      EXTERNAL           ILAENV, ISAMAX
*     ..
*     .. External Subroutines ..
      EXTERNAL           SCOPY, SGBTF2, SGEMM, SGER, SLASWP, SSCAL,
     $                   SSWAP, STRSM, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     KV is the number of superdiagonals in the factor U, allowing for
*     fill-in
*
      KV = KU + KL
*
*     Test the input parameters.
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KL.LT.0 ) THEN
         INFO = -3
      ELSE IF( KU.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDAB.LT.KL+KV+1 ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGBTRF', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
*
*     Determine the block size for this environment
*
      NB = ILAENV( 1, 'SGBTRF', ' ', M, N, KL, KU )
*
*     The block size must not exceed the limit set by the size of the
*     local arrays WORK13 and WORK31.
*
      NB = MIN( NB, NBMAX )
*
      IF( NB.LE.1 .OR. NB.GT.KL ) THEN
*
*        Use unblocked code
*
         CALL SGBTF2( M, N, KL, KU, AB, LDAB, IPIV, INFO )
      ELSE
*
*        Use blocked code
*
*        Zero the superdiagonal elements of the work array WORK13
*
         DO 20 J = 1, NB
            DO 10 I = 1, J - 1
               WORK13( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
*
*        Zero the subdiagonal elements of the work array WORK31
*
         DO 40 J = 1, NB
            DO 30 I = J + 1, NB
               WORK31( I, J ) = ZERO
   30       CONTINUE
   40    CONTINUE
*
*        Gaussian elimination with partial pivoting
*
*        Set fill-in elements in columns KU+2 to KV to zero
*
         DO 60 J = KU + 2, MIN( KV, N )
            DO 50 I = KV - J + 2, KL
               AB( I, J ) = ZERO
   50       CONTINUE
   60    CONTINUE
*
*        JU is the index of the last column affected by the current
*        stage of the factorization
*
         JU = 1
*
         DO 180 J = 1, MIN( M, N ), NB
            JB = MIN( NB, MIN( M, N )-J+1 )
*
*           The active part of the matrix is partitioned
*
*              A11   A12   A13
*              A21   A22   A23
*              A31   A32   A33
*
*           Here A11, A21 and A31 denote the current block of JB columns
*           which is about to be factorized. The number of rows in the
*           partitioning are JB, I2, I3 respectively, and the numbers
*           of columns are JB, J2, J3. The superdiagonal elements of A13
*           and the subdiagonal elements of A31 lie outside the band.
*
            I2 = MIN( KL-JB, M-J-JB+1 )
            I3 = MIN( JB, M-J-KL+1 )
*
*           J2 and J3 are computed after JU has been updated.
*
*           Factorize the current block of JB columns
*
            DO 80 JJ = J, J + JB - 1
*
*              Set fill-in elements in column JJ+KV to zero
*
               IF( JJ+KV.LE.N ) THEN
                  DO 70 I = 1, KL
                     AB( I, JJ+KV ) = ZERO
   70             CONTINUE
               END IF
*
*              Find pivot and test for singularity. KM is the number of
*              subdiagonal elements in the current column.
*
               KM = MIN( KL, M-JJ )
               JP = ISAMAX( KM+1, AB( KV+1, JJ ), 1 )
               IPIV( JJ ) = JP + JJ - J
               IF( AB( KV+JP, JJ ).NE.ZERO ) THEN
                  JU = MAX( JU, MIN( JJ+KU+JP-1, N ) )
                  IF( JP.NE.1 ) THEN
*
*                    Apply interchange to columns J to J+JB-1
*
                     IF( JP+JJ-1.LT.J+KL ) THEN
*
                        CALL SSWAP( JB, AB( KV+1+JJ-J, J ), LDAB-1,
     $                              AB( KV+JP+JJ-J, J ), LDAB-1 )
                     ELSE
*
*                       The interchange affects columns J to JJ-1 of A31
*                       which are stored in the work array WORK31
*
                        CALL SSWAP( JJ-J, AB( KV+1+JJ-J, J ), LDAB-1,
     $                              WORK31( JP+JJ-J-KL, 1 ), LDWORK )
                        CALL SSWAP( J+JB-JJ, AB( KV+1, JJ ), LDAB-1,
     $                              AB( KV+JP, JJ ), LDAB-1 )
                     END IF
                  END IF
*
*                 Compute multipliers
*
                  CALL SSCAL( KM, ONE / AB( KV+1, JJ ), AB( KV+2, JJ ),
     $                        1 )
*
*                 Update trailing submatrix within the band and within
*                 the current block. JM is the index of the last column
*                 which needs to be updated.
*
                  JM = MIN( JU, J+JB-1 )
                  IF( JM.GT.JJ )
     $               CALL SGER( KM, JM-JJ, -ONE, AB( KV+2, JJ ), 1,
     $                          AB( KV, JJ+1 ), LDAB-1,
     $                          AB( KV+1, JJ+1 ), LDAB-1 )
               ELSE
*
*                 If pivot is zero, set INFO to the index of the pivot
*                 unless a zero pivot has already been found.
*
                  IF( INFO.EQ.0 )
     $               INFO = JJ
               END IF
*
*              Copy current column of A31 into the work array WORK31
*
               NW = MIN( JJ-J+1, I3 )
               IF( NW.GT.0 )
     $            CALL SCOPY( NW, AB( KV+KL+1-JJ+J, JJ ), 1,
     $                        WORK31( 1, JJ-J+1 ), 1 )
   80       CONTINUE
            IF( J+JB.LE.N ) THEN
*
*              Apply the row interchanges to the other blocks.
*
               J2 = MIN( JU-J+1, KV ) - JB
               J3 = MAX( 0, JU-J-KV+1 )
*
*              Use SLASWP to apply the row interchanges to A12, A22, and
*              A32.
*
               CALL SLASWP( J2, AB( KV+1-JB, J+JB ), LDAB-1, 1, JB,
     $                      IPIV( J ), 1 )
*
*              Adjust the pivot indices.
*
               DO 90 I = J, J + JB - 1
                  IPIV( I ) = IPIV( I ) + J - 1
   90          CONTINUE
*
*              Apply the row interchanges to A13, A23, and A33
*              columnwise.
*
               K2 = J - 1 + JB + J2
               DO 110 I = 1, J3
                  JJ = K2 + I
                  DO 100 II = J + I - 1, J + JB - 1
                     IP = IPIV( II )
                     IF( IP.NE.II ) THEN
                        TEMP = AB( KV+1+II-JJ, JJ )
                        AB( KV+1+II-JJ, JJ ) = AB( KV+1+IP-JJ, JJ )
                        AB( KV+1+IP-JJ, JJ ) = TEMP
                     END IF
  100             CONTINUE
  110          CONTINUE
*
*              Update the relevant part of the trailing submatrix
*
               IF( J2.GT.0 ) THEN
*
*                 Update A12
*
                  CALL STRSM( 'Left', 'Lower', 'No transpose', 'Unit',
     $                        JB, J2, ONE, AB( KV+1, J ), LDAB-1,
     $                        AB( KV+1-JB, J+JB ), LDAB-1 )
*
                  IF( I2.GT.0 ) THEN
*
*                    Update A22
*
                     CALL SGEMM( 'No transpose', 'No transpose', I2, J2,
     $                           JB, -ONE, AB( KV+1+JB, J ), LDAB-1,
     $                           AB( KV+1-JB, J+JB ), LDAB-1, ONE,
     $                           AB( KV+1, J+JB ), LDAB-1 )
                  END IF
*
                  IF( I3.GT.0 ) THEN
*
*                    Update A32
*
                     CALL SGEMM( 'No transpose', 'No transpose', I3, J2,
     $                           JB, -ONE, WORK31, LDWORK,
     $                           AB( KV+1-JB, J+JB ), LDAB-1, ONE,
     $                           AB( KV+KL+1-JB, J+JB ), LDAB-1 )
                  END IF
               END IF
*
               IF( J3.GT.0 ) THEN
*
*                 Copy the lower triangle of A13 into the work array
*                 WORK13
*
                  DO 130 JJ = 1, J3
                     DO 120 II = JJ, JB
                        WORK13( II, JJ ) = AB( II-JJ+1, JJ+J+KV-1 )
  120                CONTINUE
  130             CONTINUE
*
*                 Update A13 in the work array
*
                  CALL STRSM( 'Left', 'Lower', 'No transpose', 'Unit',
     $                        JB, J3, ONE, AB( KV+1, J ), LDAB-1,
     $                        WORK13, LDWORK )
*
                  IF( I2.GT.0 ) THEN
*
*                    Update A23
*
                     CALL SGEMM( 'No transpose', 'No transpose', I2, J3,
     $                           JB, -ONE, AB( KV+1+JB, J ), LDAB-1,
     $                           WORK13, LDWORK, ONE, AB( 1+JB, J+KV ),
     $                           LDAB-1 )
                  END IF
*
                  IF( I3.GT.0 ) THEN
*
*                    Update A33
*
                     CALL SGEMM( 'No transpose', 'No transpose', I3, J3,
     $                           JB, -ONE, WORK31, LDWORK, WORK13,
     $                           LDWORK, ONE, AB( 1+KL, J+KV ), LDAB-1 )
                  END IF
*
*                 Copy the lower triangle of A13 back into place
*
                  DO 150 JJ = 1, J3
                     DO 140 II = JJ, JB
                        AB( II-JJ+1, JJ+J+KV-1 ) = WORK13( II, JJ )
  140                CONTINUE
  150             CONTINUE
               END IF
            ELSE
*
*              Adjust the pivot indices.
*
               DO 160 I = J, J + JB - 1
                  IPIV( I ) = IPIV( I ) + J - 1
  160          CONTINUE
            END IF
*
*           Partially undo the interchanges in the current block to
*           restore the upper triangular form of A31 and copy the upper
*           triangle of A31 back into place
*
            DO 170 JJ = J + JB - 1, J, -1
               JP = IPIV( JJ ) - JJ + 1
               IF( JP.NE.1 ) THEN
*
*                 Apply interchange to columns J to JJ-1
*
                  IF( JP+JJ-1.LT.J+KL ) THEN
*
*                    The interchange does not affect A31
*
                     CALL SSWAP( JJ-J, AB( KV+1+JJ-J, J ), LDAB-1,
     $                           AB( KV+JP+JJ-J, J ), LDAB-1 )
                  ELSE
*
*                    The interchange does affect A31
*
                     CALL SSWAP( JJ-J, AB( KV+1+JJ-J, J ), LDAB-1,
     $                           WORK31( JP+JJ-J-KL, 1 ), LDWORK )
                  END IF
               END IF
*
*              Copy the current column of A31 back into place
*
               NW = MIN( I3, JJ-J+1 )
               IF( NW.GT.0 )
     $            CALL SCOPY( NW, WORK31( 1, JJ-J+1 ), 1,
     $                        AB( KV+KL+1-JJ+J, JJ ), 1 )
  170       CONTINUE
  180    CONTINUE
      END IF
*
      RETURN
*
*     End of SGBTRF
*
      END

      SUBROUTINE SGBTRS( TRANS, N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB,
     $                   INFO )
*
*  -- LAPACK routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993 
*
*     .. Scalar Arguments ..
      implicit none
      CHARACTER          TRANS
      INTEGER            INFO, KL, KU, LDAB, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      REAL*8               AB( LDAB, * ), B( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  SGBTRS solves a system of linear equations
*     A * X = B  or  A' * X = B
*  with a general band matrix A using the LU factorization computed
*  by SGBTRF.
*
*  Arguments
*  =========
*
*  TRANS   (input) CHARACTER*1
*          Specifies the form of the system of equations.
*          = 'N':  A * X = B  (No transpose)
*          = 'T':  A'* X = B  (Transpose)
*          = 'C':  A'* X = B  (Conjugate transpose = Transpose)
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  KL      (input) INTEGER
*          The number of subdiagonals within the band of A.  KL >= 0.
*
*  KU      (input) INTEGER
*          The number of superdiagonals within the band of A.  KU >= 0.
*
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the matrix B.  NRHS >= 0.
*
*  AB      (input) REAL array, dimension (LDAB,N)
*          Details of the LU factorization of the band matrix A, as
*          computed by SGBTRF.  U is stored as an upper triangular band
*          matrix with KL+KU superdiagonals in rows 1 to KL+KU+1, and
*          the multipliers used during the factorization are stored in
*          rows KL+KU+2 to 2*KL+KU+1.
*
*  LDAB    (input) INTEGER
*          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
*
*  IPIV    (input) INTEGER array, dimension (N)
*          The pivot indices; for 1 <= i <= N, row i of the matrix was
*          interchanged with row IPIV(i).
*
*  B       (input/output) REAL array, dimension (LDB,NRHS)
*          On entry, the right hand side matrix B.
*          On exit, the solution matrix X.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      REAL*8               ONE
      PARAMETER          ( ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LNOTI, NOTRAN
      INTEGER            I, J, KD, L, LM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           SGEMV, SGER, SSWAP, STBSV, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      NOTRAN = LSAME( TRANS, 'N' )
      IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT.
     $    LSAME( TRANS, 'C' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KL.LT.0 ) THEN
         INFO = -3
      ELSE IF( KU.LT.0 ) THEN
         INFO = -4
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDAB.LT.( 2*KL+KU+1 ) ) THEN
         INFO = -7
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGBTRS', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. NRHS.EQ.0 )
     $   RETURN
*
      KD = KU + KL + 1
      LNOTI = KL.GT.0
*
      IF( NOTRAN ) THEN
*
*        Solve  A*X = B.
*
*        Solve L*X = B, overwriting B with X.
*
*        L is represented as a product of permutations and unit lower
*        triangular matrices L = P(1) * L(1) * ... * P(n-1) * L(n-1),
*        where each transformation L(i) is a rank-one modification of
*        the identity matrix.
*
         IF( LNOTI ) THEN
            DO 10 J = 1, N - 1
               LM = MIN( KL, N-J )
               L = IPIV( J )
               IF( L.NE.J )
     $            CALL SSWAP( NRHS, B( L, 1 ), LDB, B( J, 1 ), LDB )
               CALL SGER( LM, NRHS, -ONE, AB( KD+1, J ), 1, B( J, 1 ),
     $                    LDB, B( J+1, 1 ), LDB )
   10       CONTINUE
         END IF
*
         DO 20 I = 1, NRHS
*
*           Solve U*X = B, overwriting B with X.
*
            CALL STBSV( 'Upper', 'No transpose', 'Non-unit', N, KL+KU,
     $                  AB, LDAB, B( 1, I ), 1 )
   20    CONTINUE
*
      ELSE
*
*        Solve A'*X = B.
*
         DO 30 I = 1, NRHS
*
*           Solve U'*X = B, overwriting B with X.
*
            CALL STBSV( 'Upper', 'Transpose', 'Non-unit', N, KL+KU, AB,
     $                  LDAB, B( 1, I ), 1 )
   30    CONTINUE
*
*        Solve L'*X = B, overwriting B with X.
*
         IF( LNOTI ) THEN
            DO 40 J = N - 1, 1, -1
               LM = MIN( KL, N-J )
               CALL SGEMV( 'Transpose', LM, NRHS, -ONE, B( J+1, 1 ),
     $                     LDB, AB( KD+1, J ), 1, ONE, B( J, 1 ), LDB )
               L = IPIV( J )
               IF( L.NE.J )
     $            CALL SSWAP( NRHS, B( L, 1 ), LDB, B( J, 1 ), LDB )
   40       CONTINUE
         END IF
      END IF
      RETURN
*
*     End of SGBTRS
*
      END

      SUBROUTINE SGBTF2( M, N, KL, KU, AB, LDAB, IPIV, INFO )
*
*  -- LAPACK routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      implicit none
      INTEGER            INFO, KL, KU, LDAB, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      REAL*8               AB( LDAB, * )
*     ..
*
*  Purpose
*  =======
*
*  SGBTF2 computes an LU factorization of a real m-by-n band matrix A
*  using partial pivoting with row interchanges.
*
*  This is the unblocked version of the algorithm, calling Level 2 BLAS.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  KL      (input) INTEGER
*          The number of subdiagonals within the band of A.  KL >= 0.
*
*  KU      (input) INTEGER
*          The number of superdiagonals within the band of A.  KU >= 0.
*
*  AB      (input/output) REAL array, dimension (LDAB,N)
*          On entry, the matrix A in band storage, in rows KL+1 to
*          2*KL+KU+1; rows 1 to KL of the array need not be set.
*          The j-th column of A is stored in the j-th column of the
*          array AB as follows:
*          AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
*
*          On exit, details of the factorization: U is stored as an
*          upper triangular band matrix with KL+KU superdiagonals in
*          rows 1 to KL+KU+1, and the multipliers used during the
*          factorization are stored in rows KL+KU+2 to 2*KL+KU+1.
*          See below for further details.
*
*  LDAB    (input) INTEGER
*          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
*
*  IPIV    (output) INTEGER array, dimension (min(M,N))
*          The pivot indices; for 1 <= i <= min(M,N), row i of the
*          matrix was interchanged with row IPIV(i).
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*          > 0: if INFO = +i, U(i,i) is exactly zero. The factorization
*               has been completed, but the factor U is exactly
*               singular, and division by zero will occur if it is used
*               to solve a system of equations.
*
*  Further Details
*  ===============
*
*  The band storage scheme is illustrated by the following example, when
*  M = N = 6, KL = 2, KU = 1:
*
*  On entry:                       On exit:
*
*      *    *    *    +    +    +       *    *    *   u14  u25  u36
*      *    *    +    +    +    +       *    *   u13  u24  u35  u46
*      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
*     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
*     a21  a32  a43  a54  a65   *      m21  m32  m43  m54  m65   *
*     a31  a42  a53  a64   *    *      m31  m42  m53  m64   *    *
*
*  Array elements marked * are not used by the routine; elements marked
*  + need not be set on entry, but are required by the routine to store
*  elements of U, because of fill-in resulting from the row
*  interchanges.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL*8               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J, JP, JU, KM, KV
*     ..
*     .. External Functions ..
      INTEGER            ISAMAX
      EXTERNAL           ISAMAX
*     ..
*     .. External Subroutines ..
      EXTERNAL           SGER, SSCAL, SSWAP, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     KV is the number of superdiagonals in the factor U, allowing for
*     fill-in.
*
      KV = KU + KL
*
*     Test the input parameters.
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KL.LT.0 ) THEN
         INFO = -3
      ELSE IF( KU.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDAB.LT.KL+KV+1 ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGBTF2', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
*
*     Gaussian elimination with partial pivoting
*
*     Set fill-in elements in columns KU+2 to KV to zero.
*
      DO 20 J = KU + 2, MIN( KV, N )
         DO 10 I = KV - J + 2, KL
            AB( I, J ) = ZERO
   10    CONTINUE
   20 CONTINUE
*
*     JU is the index of the last column affected by the current stage
*     of the factorization.
*
      JU = 1
*
      DO 40 J = 1, MIN( M, N )
*
*        Set fill-in elements in column J+KV to zero.
*
         IF( J+KV.LE.N ) THEN
            DO 30 I = 1, KL
               AB( I, J+KV ) = ZERO
   30       CONTINUE
         END IF
*
*        Find pivot and test for singularity. KM is the number of
*        subdiagonal elements in the current column.
*
         KM = MIN( KL, M-J )
         JP = ISAMAX( KM+1, AB( KV+1, J ), 1 )
         IPIV( J ) = JP + J - 1
         IF( AB( KV+JP, J ).NE.ZERO ) THEN
            JU = MAX( JU, MIN( J+KU+JP-1, N ) )
*
*           Apply interchange to columns J to JU.
*
            IF( JP.NE.1 )
     $         CALL SSWAP( JU-J+1, AB( KV+JP, J ), LDAB-1,
     $                     AB( KV+1, J ), LDAB-1 )
*
            IF( KM.GT.0 ) THEN
*
*              Compute multipliers.
*
               CALL SSCAL( KM, ONE / AB( KV+1, J ), AB( KV+2, J ), 1 )
*
*              Update trailing submatrix within the band.
*
               IF( JU.GT.J )
     $            CALL SGER( KM, JU-J, -ONE, AB( KV+2, J ), 1,
     $                       AB( KV, J+1 ), LDAB-1, AB( KV+1, J+1 ),
     $                       LDAB-1 )
            END IF
         ELSE
*
*           If pivot is zero, set INFO to the index of the pivot
*           unless a zero pivot has already been found.
*
            IF( INFO.EQ.0 )
     $         INFO = J
         END IF
   40 CONTINUE
      RETURN
*
*     End of SGBTF2
*
      END

      SUBROUTINE SLASWP( N, A, LDA, K1, K2, IPIV, INCX )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      implicit none
      INTEGER            INCX, K1, K2, LDA, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      REAL*8               A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  SLASWP performs a series of row interchanges on the matrix A.
*  One row interchange is initiated for each of rows K1 through K2 of A.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.
*
*  A       (input/output) REAL array, dimension (LDA,N)
*          On entry, the matrix of column dimension N to which the row
*          interchanges will be applied.
*          On exit, the permuted matrix.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.
*
*  K1      (input) INTEGER
*          The first element of IPIV for which a row interchange will
*          be done.
*
*  K2      (input) INTEGER
*          The last element of IPIV for which a row interchange will
*          be done.
*
*  IPIV    (input) INTEGER array, dimension (M*abs(INCX))
*          The vector of pivot indices.  Only the elements in positions
*          K1 through K2 of IPIV are accessed.
*          IPIV(K) = L implies rows K and L are to be interchanged.
*
*  INCX    (input) INTEGER
*          The increment between successive values of IPIV.  If IPIV
*          is negative, the pivots are applied in reverse order.
*
* =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, IP, IX
*     ..
*     .. External Subroutines ..
      EXTERNAL           SSWAP
*     ..
*     .. Executable Statements ..
*
*     Interchange row I with row IPIV(I) for each of rows K1 through K2.
*
      IF( INCX.EQ.0 )
     $   RETURN
      IF( INCX.GT.0 ) THEN
         IX = K1
      ELSE
         IX = 1 + ( 1-K2 )*INCX
      END IF
      IF( INCX.EQ.1 ) THEN
         DO 10 I = K1, K2
            IP = IPIV( I )
            IF( IP.NE.I )
     $         CALL SSWAP( N, A( I, 1 ), LDA, A( IP, 1 ), LDA )
   10    CONTINUE
      ELSE IF( INCX.GT.1 ) THEN
         DO 20 I = K1, K2
            IP = IPIV( IX )
            IF( IP.NE.I )
     $         CALL SSWAP( N, A( I, 1 ), LDA, A( IP, 1 ), LDA )
            IX = IX + INCX
   20    CONTINUE
      ELSE IF( INCX.LT.0 ) THEN
         DO 30 I = K2, K1, -1
            IP = IPIV( IX )
            IF( IP.NE.I )
     $         CALL SSWAP( N, A( I, 1 ), LDA, A( IP, 1 ), LDA )
            IX = IX + INCX
   30    CONTINUE
      END IF
*
      RETURN
*
*     End of SLASWP
*
      END
      INTEGER          FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3,
     $                 N4 )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      implicit none
      CHARACTER*( * )    NAME, OPTS
      INTEGER            ISPEC, N1, N2, N3, N4
*     ..
*
*  Purpose
*  =======
*
*  ILAENV is called from the LAPACK routines to choose problem-dependent
*  parameters for the local environment.  See ISPEC for a description of
*  the parameters.
*
*  This version provides a set of parameters which should give good,
*  but not optimal, performance on many of the currently available
*  computers.  Users are encouraged to modify this subroutine to set
*  the tuning parameters for their particular machine using the option
*  and problem size information in the arguments.
*
*  This routine will not function correctly if it is converted to all
*  lower case.  Converting it to all upper case is allowed.
*
*  Arguments
*  =========
*
*  ISPEC   (input) INTEGER
*          Specifies the parameter to be returned as the value of
*          ILAENV.
*          = 1: the optimal blocksize; if this value is 1, an unblocked
*               algorithm will give the best performance.
*          = 2: the minimum block size for which the block routine
*               should be used; if the usable block size is less than
*               this value, an unblocked routine should be used.
*          = 3: the crossover point (in a block routine, for N less
*               than this value, an unblocked routine should be used)
*          = 4: the number of shifts, used in the nonsymmetric
*               eigenvalue routines
*          = 5: the minimum column dimension for blocking to be used;
*               rectangular blocks must have dimension at least k by m,
*               where k is given by ILAENV(2,...) and m by ILAENV(5,...)
*          = 6: the crossover point for the SVD (when reducing an m by n
*               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
*               this value, a QR factorization is used first to reduce
*               the matrix to a triangular form.)
*          = 7: the number of processors
*          = 8: the crossover point for the multishift QR and QZ methods
*               for nonsymmetric eigenvalue problems.
*
*  NAME    (input) CHARACTER*(*)
*          The name of the calling subroutine, in either upper case or
*          lower case.
*
*  OPTS    (input) CHARACTER*(*)
*          The character options to the subroutine NAME, concatenated
*          into a single character string.  For example, UPLO = 'U',
*          TRANS = 'T', and DIAG = 'N' for a triangular routine would
*          be specified as OPTS = 'UTN'.
*
*  N1      (input) INTEGER
*  N2      (input) INTEGER
*  N3      (input) INTEGER
*  N4      (input) INTEGER
*          Problem dimensions for the subroutine NAME; these may not all
*          be required.
*
* (ILAENV) (output) INTEGER
*          >= 0: the value of the parameter specified by ISPEC
*          < 0:  if ILAENV = -k, the k-th argument had an illegal value.
*
*  Further Details
*  ===============
*
*  The following conventions have been used when calling ILAENV from the
*  LAPACK routines:
*  1)  OPTS is a concatenation of all of the character options to
*      subroutine NAME, in the same order that they appear in the
*      argument list for NAME, even if they are not used in determining
*      the value of the parameter specified by ISPEC.
*  2)  The problem dimensions N1, N2, N3, N4 are specified in the order
*      that they appear in the argument list for NAME.  N1 is used
*      first, N2 second, and so on, and unused problem dimensions are
*      passed a value of -1.
*  3)  The parameter value returned by ILAENV is checked for validity in
*      the calling subroutine.  For example, ILAENV is used to retrieve
*      the optimal blocksize for STRTRI as follows:
*
*      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
*      IF( NB.LE.1 ) NB = MAX( 1, N )
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            CNAME, SNAME
      CHARACTER*1        C1
      CHARACTER*2        C2, C4
      CHARACTER*3        C3
      CHARACTER*6        SUBNAM
      INTEGER            I, IC, IZ, NB, NBMIN, NX
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CHAR, ICHAR, INT, MIN, REAL
*     ..
*     .. Executable Statements ..
*
      GO TO ( 100, 100, 100, 400, 500, 600, 700, 800 ) ISPEC
*
*     Invalid value for ISPEC
*
      ILAENV = -1
      RETURN
*
  100 CONTINUE
*
*     Convert NAME to upper case if the first character is lower case.
*
      ILAENV = 1
      SUBNAM = NAME
      IC = ICHAR( SUBNAM( 1:1 ) )
      IZ = ICHAR( 'Z' )
      IF( IZ.EQ.90 .OR. IZ.EQ.122 ) THEN
*
*        ASCII character set
*
         IF( IC.GE.97 .AND. IC.LE.122 ) THEN
            SUBNAM( 1:1 ) = CHAR( IC-32 )
            DO 10 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( IC.GE.97 .AND. IC.LE.122 )
     $            SUBNAM( I:I ) = CHAR( IC-32 )
   10       CONTINUE
         END IF
*
      ELSE IF( IZ.EQ.233 .OR. IZ.EQ.169 ) THEN
*
*        EBCDIC character set
*
         IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.
     $       ( IC.GE.145 .AND. IC.LE.153 ) .OR.
     $       ( IC.GE.162 .AND. IC.LE.169 ) ) THEN
            SUBNAM( 1:1 ) = CHAR( IC+64 )
            DO 20 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.
     $             ( IC.GE.145 .AND. IC.LE.153 ) .OR.
     $             ( IC.GE.162 .AND. IC.LE.169 ) )
     $            SUBNAM( I:I ) = CHAR( IC+64 )
   20       CONTINUE
         END IF
*
      ELSE IF( IZ.EQ.218 .OR. IZ.EQ.250 ) THEN
*
*        Prime machines:  ASCII+128
*
         IF( IC.GE.225 .AND. IC.LE.250 ) THEN
            SUBNAM( 1:1 ) = CHAR( IC-32 )
            DO 30 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( IC.GE.225 .AND. IC.LE.250 )
     $            SUBNAM( I:I ) = CHAR( IC-32 )
   30       CONTINUE
         END IF
      END IF
*
      C1 = SUBNAM( 1:1 )
      SNAME = C1.EQ.'S' .OR. C1.EQ.'D'
      CNAME = C1.EQ.'C' .OR. C1.EQ.'Z'
      IF( .NOT.( CNAME .OR. SNAME ) )
     $   RETURN
      C2 = SUBNAM( 2:3 )
      C3 = SUBNAM( 4:6 )
      C4 = C3( 2:3 )
*
      GO TO ( 110, 200, 300 ) ISPEC
*
  110 CONTINUE
*
*     ISPEC = 1:  block size
*
*     In these examples, separate code is provided for setting NB for
*     real and complex.  We assume that NB will take the same value in
*     single or double precision.
*
      NB = 1
*
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.
     $            C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'PO' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NB = 1
         ELSE IF( SNAME .AND. C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            NB = 64
         ELSE IF( C3.EQ.'TRD' ) THEN
            NB = 1
         ELSE IF( C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         END IF
      ELSE IF( C2.EQ.'GB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'PB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'TR' ) THEN
         IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'LA' ) THEN
         IF( C3.EQ.'UUM' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'ST' ) THEN
         IF( C3.EQ.'EBZ' ) THEN
            NB = 1
         END IF
      END IF
      ILAENV = NB
      RETURN
*
  200 CONTINUE
*
*     ISPEC = 2:  minimum block size
*
      NBMIN = 2
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.
     $       C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 8
            ELSE
               NBMIN = 8
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         END IF
      END IF
      ILAENV = NBMIN
      RETURN
*
  300 CONTINUE
*
*     ISPEC = 3:  crossover point
*
      NX = 0
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.
     $       C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NX = 1
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NX = 1
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NX = 128
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NX = 128
            END IF
         END IF
      END IF
      ILAENV = NX
      RETURN
*
  400 CONTINUE
*
*     ISPEC = 4:  number of shifts (used by xHSEQR)
*
      ILAENV = 6
      RETURN
*
  500 CONTINUE
*
*     ISPEC = 5:  minimum column dimension (not used)
*
      ILAENV = 2
      RETURN
*
  600 CONTINUE 
*
*     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)
*
      ILAENV = INT( REAL( MIN( N1, N2 ))*1.6 )
      RETURN
*
  700 CONTINUE
*
*     ISPEC = 7:  number of processors (not used)
*
      ILAENV = 1
      RETURN
*
  800 CONTINUE
*
*     ISPEC = 8:  crossover point for multishift (used by xHSEQR)
*
      ILAENV = 50
      RETURN
*
*     End of ILAENV
*
      END
      SUBROUTINE XERBLA( SRNAME, INFO )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      implicit none
      CHARACTER*6        SRNAME
      INTEGER            INFO
*     ..
*
*  Purpose
*  =======
*
*  XERBLA  is an error handler for the LAPACK routines.
*  It is called by an LAPACK routine if an input parameter has an
*  invalid value.  A message is printed and execution stops.
*
*  Installers may consider modifying the STOP statement in order to
*  call system-specific exception-handling facilities.
*
*  Arguments
*  =========
*
*  SRNAME  (input) CHARACTER*6
*          The name of the routine which called XERBLA.
*
*  INFO    (input) INTEGER
*          The position of the invalid parameter in the parameter list
*          of the calling routine.
*
* =====================================================================
*
*     .. Executable Statements ..
*
      WRITE( *, FMT = 9999 )SRNAME, INFO
*
      STOP
*
 9999 FORMAT( ' ** On entry to ', A6, ' parameter number ', I2, ' had ',
     $      'an illegal value' )
*
*     End of XERBLA
*
      END
      LOGICAL          FUNCTION LSAME( CA, CB )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      implicit none
      CHARACTER          CA, CB
*     ..
*
*  Purpose
*  =======
*
*  LSAME returns .TRUE. if CA is the same letter as CB regardless of
*  case.
*
*  Arguments
*  =========
*
*  CA      (input) CHARACTER*1
*  CB      (input) CHARACTER*1
*          CA and CB specify the single characters to be compared.
*
* =====================================================================
*
*     .. Intrinsic Functions ..
      INTRINSIC          ICHAR
*     ..
*     .. Local Scalars ..
      INTEGER            INTA, INTB, ZCODE
*     ..
*     .. Executable Statements ..
*
*     Test if the characters are equal
*
      LSAME = CA.EQ.CB
      IF( LSAME )
     $   RETURN
*
*     Now test for equivalence if both characters are alphabetic.
*
      ZCODE = ICHAR( 'Z' )
*
*     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
*     machines, on which ICHAR returns a value with bit 8 set.
*     ICHAR('A') on Prime machines returns 193 which is the same as
*     ICHAR('A') on an EBCDIC machine.
*
      INTA = ICHAR( CA )
      INTB = ICHAR( CB )
*
      IF( ZCODE.EQ.90 .OR. ZCODE.EQ.122 ) THEN
*
*        ASCII is assumed - ZCODE is the ASCII code of either lower or
*        upper case 'Z'.
*
         IF( INTA.GE.97 .AND. INTA.LE.122 ) INTA = INTA - 32
         IF( INTB.GE.97 .AND. INTB.LE.122 ) INTB = INTB - 32
*
      ELSE IF( ZCODE.EQ.233 .OR. ZCODE.EQ.169 ) THEN
*
*        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
*        upper case 'Z'.
*
         IF( INTA.GE.129 .AND. INTA.LE.137 .OR.
     $       INTA.GE.145 .AND. INTA.LE.153 .OR.
     $       INTA.GE.162 .AND. INTA.LE.169 ) INTA = INTA + 64
         IF( INTB.GE.129 .AND. INTB.LE.137 .OR.
     $       INTB.GE.145 .AND. INTB.LE.153 .OR.
     $       INTB.GE.162 .AND. INTB.LE.169 ) INTB = INTB + 64
*
      ELSE IF( ZCODE.EQ.218 .OR. ZCODE.EQ.250 ) THEN
*
*        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
*        plus 128 of either lower or upper case 'Z'.
*
         IF( INTA.GE.225 .AND. INTA.LE.250 ) INTA = INTA - 32
         IF( INTB.GE.225 .AND. INTB.LE.250 ) INTB = INTB - 32
      END IF
      LSAME = INTA.EQ.INTB
*
*     RETURN
*
*     End of LSAME
*
      END

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
c
c     all_sgbtrf_s.blas.f
c
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

      integer function isamax(n,sx,incx)
c
c     finds the index of element having max. absolute value.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      implicit none
      real*8 sx(*),smax
      integer i,incx,ix,n
c
      isamax = 0
      if( n.lt.1 .or. incx.le.0 ) return
      isamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      smax = abs(sx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(abs(sx(ix)).le.smax) go to 5
         isamax = i
         smax = abs(sx(ix))
    5    ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 smax = abs(sx(1))
      do 30 i = 2,n
         if(abs(sx(i)).le.smax) go to 30
         isamax = i
         smax = abs(sx(i))
   30 continue
      return
      end
      subroutine scopy(n,sx,incx,sy,incy)
c
c     copies a vector, x, to a vector, y.
c     uses unrolled loops for increments equal to 1.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      implicit none
      real*8 sx(*),sy(*)
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        sy(iy) = sx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,7)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        sy(i) = sx(i)
   30 continue
      if( n .lt. 7 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,7
        sy(i) = sx(i)
        sy(i + 1) = sx(i + 1)
        sy(i + 2) = sx(i + 2)
        sy(i + 3) = sx(i + 3)
        sy(i + 4) = sx(i + 4)
        sy(i + 5) = sx(i + 5)
        sy(i + 6) = sx(i + 6)
   50 continue
      return
      end
      SUBROUTINE SGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB,
     $                   BETA, C, LDC )
*     .. Scalar Arguments ..
      implicit none
      CHARACTER*1        TRANSA, TRANSB
      INTEGER            M, N, K, LDA, LDB, LDC
      REAL*8               ALPHA, BETA
*     .. Array Arguments ..
      REAL*8               A( LDA, * ), B( LDB, * ), C( LDC, * )
*     ..
*
*  Purpose
*  =======
*
*  SGEMM  performs one of the matrix-matrix operations
*
*     C := alpha*op( A )*op( B ) + beta*C,
*
*  where  op( X ) is one of
*
*     op( X ) = X   or   op( X ) = X',
*
*  alpha and beta are scalars, and A, B and C are matrices, with op( A )
*  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
*
*  Parameters
*  ==========
*
*  TRANSA - CHARACTER*1.
*           On entry, TRANSA specifies the form of op( A ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSA = 'N' or 'n',  op( A ) = A.
*
*              TRANSA = 'T' or 't',  op( A ) = A'.
*
*              TRANSA = 'C' or 'c',  op( A ) = A'.
*
*           Unchanged on exit.
*
*  TRANSB - CHARACTER*1.
*           On entry, TRANSB specifies the form of op( B ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSB = 'N' or 'n',  op( B ) = B.
*
*              TRANSB = 'T' or 't',  op( B ) = B'.
*
*              TRANSB = 'C' or 'c',  op( B ) = B'.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry,  M  specifies  the number  of rows  of the  matrix
*           op( A )  and of the  matrix  C.  M  must  be at least  zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry,  N  specifies the number  of columns of the matrix
*           op( B ) and the number of columns of the matrix C. N must be
*           at least zero.
*           Unchanged on exit.
*
*  K      - INTEGER.
*           On entry,  K  specifies  the number of columns of the matrix
*           op( A ) and the number of rows of the matrix op( B ). K must
*           be at least  zero.
*           Unchanged on exit.
*
*  ALPHA  - REAL            .
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - REAL             array of DIMENSION ( LDA, ka ), where ka is
*           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
*           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
*           part of the array  A  must contain the matrix  A,  otherwise
*           the leading  k by m  part of the array  A  must contain  the
*           matrix A.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
*           LDA must be at least  max( 1, m ), otherwise  LDA must be at
*           least  max( 1, k ).
*           Unchanged on exit.
*
*  B      - REAL             array of DIMENSION ( LDB, kb ), where kb is
*           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
*           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
*           part of the array  B  must contain the matrix  B,  otherwise
*           the leading  n by k  part of the array  B  must contain  the
*           matrix B.
*           Unchanged on exit.
*
*  LDB    - INTEGER.
*           On entry, LDB specifies the first dimension of B as declared
*           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
*           LDB must be at least  max( 1, k ), otherwise  LDB must be at
*           least  max( 1, n ).
*           Unchanged on exit.
*
*  BETA   - REAL            .
*           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
*           supplied as zero then C need not be set on input.
*           Unchanged on exit.
*
*  C      - REAL             array of DIMENSION ( LDC, n ).
*           Before entry, the leading  m by n  part of the array  C must
*           contain the matrix  C,  except when  beta  is zero, in which
*           case C need not be set on entry.
*           On exit, the array  C  is overwritten by the  m by n  matrix
*           ( alpha*op( A )*op( B ) + beta*C ).
*
*  LDC    - INTEGER.
*           On entry, LDC specifies the first dimension of C as declared
*           in  the  calling  (sub)  program.   LDC  must  be  at  least
*           max( 1, m ).
*           Unchanged on exit.
*
*
*  Level 3 Blas routine.
*
*  -- Written on 8-February-1989.
*     Jack Dongarra, Argonne National Laboratory.
*     Iain Duff, AERE Harwell.
*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*     Sven Hammarling, Numerical Algorithms Group Ltd.
*
*
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     .. Local Scalars ..
      LOGICAL            NOTA, NOTB
      INTEGER            I, INFO, J, L, NCOLA, NROWA, NROWB
      REAL*8               TEMP
*     .. Parameters ..
      REAL*8               ONE         , ZERO
      PARAMETER        ( ONE = 1.0E+0, ZERO = 0.0E+0 )
*     ..
*     .. Executable Statements ..
*
*     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
*     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
*     and  columns of  A  and the  number of  rows  of  B  respectively.
*
      NOTA  = LSAME( TRANSA, 'N' )
      NOTB  = LSAME( TRANSB, 'N' )
      IF( NOTA )THEN
         NROWA = M
         NCOLA = K
      ELSE
         NROWA = K
         NCOLA = M
      END IF
      IF( NOTB )THEN
         NROWB = K
      ELSE
         NROWB = N
      END IF
*
*     Test the input parameters.
*
      INFO = 0
      IF(      ( .NOT.NOTA                 ).AND.
     $         ( .NOT.LSAME( TRANSA, 'C' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'T' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.NOTB                 ).AND.
     $         ( .NOT.LSAME( TRANSB, 'C' ) ).AND.
     $         ( .NOT.LSAME( TRANSB, 'T' ) )      )THEN
         INFO = 2
      ELSE IF( M  .LT.0               )THEN
         INFO = 3
      ELSE IF( N  .LT.0               )THEN
         INFO = 4
      ELSE IF( K  .LT.0               )THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 8
      ELSE IF( LDB.LT.MAX( 1, NROWB ) )THEN
         INFO = 10
      ELSE IF( LDC.LT.MAX( 1, M     ) )THEN
         INFO = 13
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'SGEMM ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.
     $    ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
*
*     And if  alpha.eq.zero.
*
      IF( ALPHA.EQ.ZERO )THEN
         IF( BETA.EQ.ZERO )THEN
            DO 20, J = 1, N
               DO 10, I = 1, M
                  C( I, J ) = ZERO
   10          CONTINUE
   20       CONTINUE
         ELSE
            DO 40, J = 1, N
               DO 30, I = 1, M
                  C( I, J ) = BETA*C( I, J )
   30          CONTINUE
   40       CONTINUE
         END IF
         RETURN
      END IF
*
*     Start the operations.
*
      IF( NOTB )THEN
         IF( NOTA )THEN
*
*           Form  C := alpha*A*B + beta*C.
*
            DO 90, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 50, I = 1, M
                     C( I, J ) = ZERO
   50             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 60, I = 1, M
                     C( I, J ) = BETA*C( I, J )
   60             CONTINUE
               END IF
               DO 80, L = 1, K
                  IF( B( L, J ).NE.ZERO )THEN
                     TEMP = ALPHA*B( L, J )
                     DO 70, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
   70                CONTINUE
                  END IF
   80          CONTINUE
   90       CONTINUE
         ELSE
*
*           Form  C := alpha*A'*B + beta*C
*
            DO 120, J = 1, N
               DO 110, I = 1, M
                  TEMP = ZERO
                  DO 100, L = 1, K
                     TEMP = TEMP + A( L, I )*B( L, J )
  100             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  110          CONTINUE
  120       CONTINUE
         END IF
      ELSE
         IF( NOTA )THEN
*
*           Form  C := alpha*A*B' + beta*C
*
            DO 170, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 130, I = 1, M
                     C( I, J ) = ZERO
  130             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 140, I = 1, M
                     C( I, J ) = BETA*C( I, J )
  140             CONTINUE
               END IF
               DO 160, L = 1, K
                  IF( B( J, L ).NE.ZERO )THEN
                     TEMP = ALPHA*B( J, L )
                     DO 150, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
  150                CONTINUE
                  END IF
  160          CONTINUE
  170       CONTINUE
         ELSE
*
*           Form  C := alpha*A'*B' + beta*C
*
            DO 200, J = 1, N
               DO 190, I = 1, M
                  TEMP = ZERO
                  DO 180, L = 1, K
                     TEMP = TEMP + A( L, I )*B( J, L )
  180             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  190          CONTINUE
  200       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of SGEMM .
*
      END
      SUBROUTINE SGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
*     .. Scalar Arguments ..
      implicit none
      REAL*8               ALPHA, BETA
      INTEGER            INCX, INCY, LDA, M, N
      CHARACTER*1        TRANS
*     .. Array Arguments ..
      REAL*8               A( LDA, * ), X( * ), Y( * )
*     ..
*
*  Purpose
*  =======
*
*  SGEMV  performs one of the matrix-vector operations
*
*     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
*
*  where alpha and beta are scalars, x and y are vectors and A is an
*  m by n matrix.
*
*  Parameters
*  ==========
*
*  TRANS  - CHARACTER*1.
*           On entry, TRANS specifies the operation to be performed as
*           follows:
*
*              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
*
*              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
*
*              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of the matrix A.
*           M must be at least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - REAL            .
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - REAL             array of DIMENSION ( LDA, n ).
*           Before entry, the leading m by n part of the array A must
*           contain the matrix of coefficients.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           max( 1, m ).
*           Unchanged on exit.
*
*  X      - REAL             array of DIMENSION at least
*           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
*           and at least
*           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
*           Before entry, the incremented array X must contain the
*           vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*  BETA   - REAL            .
*           On entry, BETA specifies the scalar beta. When BETA is
*           supplied as zero then Y need not be set on input.
*           Unchanged on exit.
*
*  Y      - REAL             array of DIMENSION at least
*           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
*           and at least
*           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
*           Before entry with BETA non-zero, the incremented array Y
*           must contain the vector y. On exit, Y is overwritten by the
*           updated vector y.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y. INCY must not be zero.
*           Unchanged on exit.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*
*     .. Parameters ..
      REAL*8               ONE         , ZERO
      PARAMETER        ( ONE = 1.0E+0, ZERO = 0.0E+0 )
*     .. Local Scalars ..
      REAL*8               TEMP
      INTEGER            I, INFO, IX, IY, J, JX, JY, KX, KY, LENX, LENY
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF     ( .NOT.LSAME( TRANS, 'N' ).AND.
     $         .NOT.LSAME( TRANS, 'T' ).AND.
     $         .NOT.LSAME( TRANS, 'C' )      )THEN
         INFO = 1
      ELSE IF( M.LT.0 )THEN
         INFO = 2
      ELSE IF( N.LT.0 )THEN
         INFO = 3
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 6
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 8
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'SGEMV ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.
     $    ( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
*
*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
*     up the start points in  X  and  Y.
*
      IF( LSAME( TRANS, 'N' ) )THEN
         LENX = N
         LENY = M
      ELSE
         LENX = M
         LENY = N
      END IF
      IF( INCX.GT.0 )THEN
         KX = 1
      ELSE
         KX = 1 - ( LENX - 1 )*INCX
      END IF
      IF( INCY.GT.0 )THEN
         KY = 1
      ELSE
         KY = 1 - ( LENY - 1 )*INCY
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
*     First form  y := beta*y.
*
      IF( BETA.NE.ONE )THEN
         IF( INCY.EQ.1 )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 10, I = 1, LENY
                  Y( I ) = ZERO
   10          CONTINUE
            ELSE
               DO 20, I = 1, LENY
                  Y( I ) = BETA*Y( I )
   20          CONTINUE
            END IF
         ELSE
            IY = KY
            IF( BETA.EQ.ZERO )THEN
               DO 30, I = 1, LENY
                  Y( IY ) = ZERO
                  IY      = IY   + INCY
   30          CONTINUE
            ELSE
               DO 40, I = 1, LENY
                  Y( IY ) = BETA*Y( IY )
                  IY      = IY           + INCY
   40          CONTINUE
            END IF
         END IF
      END IF
      IF( ALPHA.EQ.ZERO )
     $   RETURN
      IF( LSAME( TRANS, 'N' ) )THEN
*
*        Form  y := alpha*A*x + y.
*
         JX = KX
         IF( INCY.EQ.1 )THEN
            DO 60, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  DO 50, I = 1, M
                     Y( I ) = Y( I ) + TEMP*A( I, J )
   50             CONTINUE
               END IF
               JX = JX + INCX
   60       CONTINUE
         ELSE
            DO 80, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  IY   = KY
                  DO 70, I = 1, M
                     Y( IY ) = Y( IY ) + TEMP*A( I, J )
                     IY      = IY      + INCY
   70             CONTINUE
               END IF
               JX = JX + INCX
   80       CONTINUE
         END IF
      ELSE
*
*        Form  y := alpha*A'*x + y.
*
         JY = KY
         IF( INCX.EQ.1 )THEN
            DO 100, J = 1, N
               TEMP = ZERO
               DO 90, I = 1, M
                  TEMP = TEMP + A( I, J )*X( I )
   90          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  100       CONTINUE
         ELSE
            DO 120, J = 1, N
               TEMP = ZERO
               IX   = KX
               DO 110, I = 1, M
                  TEMP = TEMP + A( I, J )*X( IX )
                  IX   = IX   + INCX
  110          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  120       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of SGEMV .
*
      END
      SUBROUTINE SGER  ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
*     .. Scalar Arguments ..
      implicit none
      REAL*8               ALPHA
      INTEGER            INCX, INCY, LDA, M, N
*     .. Array Arguments ..
      REAL*8               A( LDA, * ), X( * ), Y( * )
*     ..
*
*  Purpose
*  =======
*
*  SGER   performs the rank 1 operation
*
*     A := alpha*x*y' + A,
*
*  where alpha is a scalar, x is an m element vector, y is an n element
*  vector and A is an m by n matrix.
*
*  Parameters
*  ==========
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of the matrix A.
*           M must be at least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - REAL            .
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  X      - REAL             array of dimension at least
*           ( 1 + ( m - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the m
*           element vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*  Y      - REAL             array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCY ) ).
*           Before entry, the incremented array Y must contain the n
*           element vector y.
*           Unchanged on exit.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y. INCY must not be zero.
*           Unchanged on exit.
*
*  A      - REAL             array of DIMENSION ( LDA, n ).
*           Before entry, the leading m by n part of the array A must
*           contain the matrix of coefficients. On exit, A is
*           overwritten by the updated matrix.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           max( 1, m ).
*           Unchanged on exit.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*
*     .. Parameters ..
      REAL*8               ZERO
      PARAMETER        ( ZERO = 0.0E+0 )
*     .. Local Scalars ..
      REAL*8               TEMP
      INTEGER            I, INFO, IX, J, JY, KX
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF     ( M.LT.0 )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 5
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 7
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'SGER  ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.( ALPHA.EQ.ZERO ) )
     $   RETURN
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
      IF( INCY.GT.0 )THEN
         JY = 1
      ELSE
         JY = 1 - ( N - 1 )*INCY
      END IF
      IF( INCX.EQ.1 )THEN
         DO 20, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               DO 10, I = 1, M
                  A( I, J ) = A( I, J ) + X( I )*TEMP
   10          CONTINUE
            END IF
            JY = JY + INCY
   20    CONTINUE
      ELSE
         IF( INCX.GT.0 )THEN
            KX = 1
         ELSE
            KX = 1 - ( M - 1 )*INCX
         END IF
         DO 40, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               IX   = KX
               DO 30, I = 1, M
                  A( I, J ) = A( I, J ) + X( IX )*TEMP
                  IX        = IX        + INCX
   30          CONTINUE
            END IF
            JY = JY + INCY
   40    CONTINUE
      END IF
*
      RETURN
*
*     End of SGER  .
*
      END
      subroutine sscal(n,sa,sx,incx)
c
c     scales a vector by a constant.
c     uses unrolled loops for increment equal to 1.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      implicit none
      real*8 sa,sx(*)
      integer i,incx,m,mp1,n,nincx
c
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
        sx(i) = sa*sx(i)
   10 continue
      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        sx(i) = sa*sx(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        sx(i) = sa*sx(i)
        sx(i + 1) = sa*sx(i + 1)
        sx(i + 2) = sa*sx(i + 2)
        sx(i + 3) = sa*sx(i + 3)
        sx(i + 4) = sa*sx(i + 4)
   50 continue
      return
      end
      subroutine sswap (n,sx,incx,sy,incy)
c
c     interchanges two vectors.
c     uses unrolled loops for increments equal to 1.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      implicit none
      real*8 sx(*),sy(*),stemp
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c       code for unequal increments or equal increments not equal
c         to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        stemp = sx(ix)
        sx(ix) = sy(iy)
        sy(iy) = stemp
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c       code for both increments equal to 1
c
c
c       clean-up loop
c
   20 m = mod(n,3)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        stemp = sx(i)
        sx(i) = sy(i)
        sy(i) = stemp
   30 continue
      if( n .lt. 3 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,3
        stemp = sx(i)
        sx(i) = sy(i)
        sy(i) = stemp
        stemp = sx(i + 1)
        sx(i + 1) = sy(i + 1)
        sy(i + 1) = stemp
        stemp = sx(i + 2)
        sx(i + 2) = sy(i + 2)
        sy(i + 2) = stemp
   50 continue
      return
      end
      SUBROUTINE STBSV ( UPLO, TRANS, DIAG, N, K, A, LDA, X, INCX )
*     .. Scalar Arguments ..
      implicit none
      INTEGER            INCX, K, LDA, N
      CHARACTER*1        DIAG, TRANS, UPLO
*     .. Array Arguments ..
      REAL*8               A( LDA, * ), X( * )
*     ..
*
*  Purpose
*  =======
*
*  STBSV  solves one of the systems of equations
*
*     A*x = b,   or   A'*x = b,
*
*  where b and x are n element vectors and A is an n by n unit, or
*  non-unit, upper or lower triangular band matrix, with ( k + 1 )
*  diagonals.
*
*  No test for singularity or near-singularity is included in this
*  routine. Such tests must be performed before calling this routine.
*
*  Parameters
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the matrix is an upper or
*           lower triangular matrix as follows:
*
*              UPLO = 'U' or 'u'   A is an upper triangular matrix.
*
*              UPLO = 'L' or 'l'   A is a lower triangular matrix.
*
*           Unchanged on exit.
*
*  TRANS  - CHARACTER*1.
*           On entry, TRANS specifies the equations to be solved as
*           follows:
*
*              TRANS = 'N' or 'n'   A*x = b.
*
*              TRANS = 'T' or 't'   A'*x = b.
*
*              TRANS = 'C' or 'c'   A'*x = b.
*
*           Unchanged on exit.
*
*  DIAG   - CHARACTER*1.
*           On entry, DIAG specifies whether or not A is unit
*           triangular as follows:
*
*              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
*
*              DIAG = 'N' or 'n'   A is not assumed to be unit
*                                  triangular.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the order of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  K      - INTEGER.
*           On entry with UPLO = 'U' or 'u', K specifies the number of
*           super-diagonals of the matrix A.
*           On entry with UPLO = 'L' or 'l', K specifies the number of
*           sub-diagonals of the matrix A.
*           K must satisfy  0 .le. K.
*           Unchanged on exit.
*
*  A      - REAL             array of DIMENSION ( LDA, n ).
*           Before entry with UPLO = 'U' or 'u', the leading ( k + 1 )
*           by n part of the array A must contain the upper triangular
*           band part of the matrix of coefficients, supplied column by
*           column, with the leading diagonal of the matrix in row
*           ( k + 1 ) of the array, the first super-diagonal starting at
*           position 2 in row k, and so on. The top left k by k triangle
*           of the array A is not referenced.
*           The following program segment will transfer an upper
*           triangular band matrix from conventional full matrix storage
*           to band storage:
*
*                 DO 20, J = 1, N
*                    M = K + 1 - J
*                    DO 10, I = MAX( 1, J - K ), J
*                       A( M + I, J ) = matrix( I, J )
*              10    CONTINUE
*              20 CONTINUE
*
*           Before entry with UPLO = 'L' or 'l', the leading ( k + 1 )
*           by n part of the array A must contain the lower triangular
*           band part of the matrix of coefficients, supplied column by
*           column, with the leading diagonal of the matrix in row 1 of
*           the array, the first sub-diagonal starting at position 1 in
*           row 2, and so on. The bottom right k by k triangle of the
*           array A is not referenced.
*           The following program segment will transfer a lower
*           triangular band matrix from conventional full matrix storage
*           to band storage:
*
*                 DO 20, J = 1, N
*                    M = 1 - J
*                    DO 10, I = J, MIN( N, J + K )
*                       A( M + I, J ) = matrix( I, J )
*              10    CONTINUE
*              20 CONTINUE
*
*           Note that when DIAG = 'U' or 'u' the elements of the array A
*           corresponding to the diagonal elements of the matrix are not
*           referenced, but are assumed to be unity.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           ( k + 1 ).
*           Unchanged on exit.
*
*  X      - REAL             array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the n
*           element right-hand side vector b. On exit, X is overwritten
*           with the solution vector x.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*
*     .. Parameters ..
      REAL*8               ZERO
      PARAMETER        ( ZERO = 0.0E+0 )
*     .. Local Scalars ..
      REAL*8               TEMP
      INTEGER            I, INFO, IX, J, JX, KPLUS1, KX, L
      LOGICAL            NOUNIT
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF     ( .NOT.LSAME( UPLO , 'U' ).AND.
     $         .NOT.LSAME( UPLO , 'L' )      )THEN
         INFO = 1
      ELSE IF( .NOT.LSAME( TRANS, 'N' ).AND.
     $         .NOT.LSAME( TRANS, 'T' ).AND.
     $         .NOT.LSAME( TRANS, 'C' )      )THEN
         INFO = 2
      ELSE IF( .NOT.LSAME( DIAG , 'U' ).AND.
     $         .NOT.LSAME( DIAG , 'N' )      )THEN
         INFO = 3
      ELSE IF( N.LT.0 )THEN
         INFO = 4
      ELSE IF( K.LT.0 )THEN
         INFO = 5
      ELSE IF( LDA.LT.( K + 1 ) )THEN
         INFO = 7
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'STBSV ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( N.EQ.0 )
     $   RETURN
*
      NOUNIT = LSAME( DIAG, 'N' )
*
*     Set up the start point in X if the increment is not unity. This
*     will be  ( N - 1 )*INCX  too small for descending loops.
*
      IF( INCX.LE.0 )THEN
         KX = 1 - ( N - 1 )*INCX
      ELSE IF( INCX.NE.1 )THEN
         KX = 1
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed by sequentially with one pass through A.
*
      IF( LSAME( TRANS, 'N' ) )THEN
*
*        Form  x := inv( A )*x.
*
         IF( LSAME( UPLO, 'U' ) )THEN
            KPLUS1 = K + 1
            IF( INCX.EQ.1 )THEN
               DO 20, J = N, 1, -1
                  IF( X( J ).NE.ZERO )THEN
                     L = KPLUS1 - J
                     IF( NOUNIT )
     $                  X( J ) = X( J )/A( KPLUS1, J )
                     TEMP = X( J )
                     DO 10, I = J - 1, MAX( 1, J - K ), -1
                        X( I ) = X( I ) - TEMP*A( L + I, J )
   10                CONTINUE
                  END IF
   20          CONTINUE
            ELSE
               KX = KX + ( N - 1 )*INCX
               JX = KX
               DO 40, J = N, 1, -1
                  KX = KX - INCX
                  IF( X( JX ).NE.ZERO )THEN
                     IX = KX
                     L  = KPLUS1 - J
                     IF( NOUNIT )
     $                  X( JX ) = X( JX )/A( KPLUS1, J )
                     TEMP = X( JX )
                     DO 30, I = J - 1, MAX( 1, J - K ), -1
                        X( IX ) = X( IX ) - TEMP*A( L + I, J )
                        IX      = IX      - INCX
   30                CONTINUE
                  END IF
                  JX = JX - INCX
   40          CONTINUE
            END IF
         ELSE
            IF( INCX.EQ.1 )THEN
               DO 60, J = 1, N
                  IF( X( J ).NE.ZERO )THEN
                     L = 1 - J
                     IF( NOUNIT )
     $                  X( J ) = X( J )/A( 1, J )
                     TEMP = X( J )
                     DO 50, I = J + 1, MIN( N, J + K )
                        X( I ) = X( I ) - TEMP*A( L + I, J )
   50                CONTINUE
                  END IF
   60          CONTINUE
            ELSE
               JX = KX
               DO 80, J = 1, N
                  KX = KX + INCX
                  IF( X( JX ).NE.ZERO )THEN
                     IX = KX
                     L  = 1  - J
                     IF( NOUNIT )
     $                  X( JX ) = X( JX )/A( 1, J )
                     TEMP = X( JX )
                     DO 70, I = J + 1, MIN( N, J + K )
                        X( IX ) = X( IX ) - TEMP*A( L + I, J )
                        IX      = IX      + INCX
   70                CONTINUE
                  END IF
                  JX = JX + INCX
   80          CONTINUE
            END IF
         END IF
      ELSE
*
*        Form  x := inv( A')*x.
*
         IF( LSAME( UPLO, 'U' ) )THEN
            KPLUS1 = K + 1
            IF( INCX.EQ.1 )THEN
               DO 100, J = 1, N
                  TEMP = X( J )
                  L    = KPLUS1 - J
                  DO 90, I = MAX( 1, J - K ), J - 1
                     TEMP = TEMP - A( L + I, J )*X( I )
   90             CONTINUE
                  IF( NOUNIT )
     $               TEMP = TEMP/A( KPLUS1, J )
                  X( J ) = TEMP
  100          CONTINUE
            ELSE
               JX = KX
               DO 120, J = 1, N
                  TEMP = X( JX )
                  IX   = KX
                  L    = KPLUS1  - J
                  DO 110, I = MAX( 1, J - K ), J - 1
                     TEMP = TEMP - A( L + I, J )*X( IX )
                     IX   = IX   + INCX
  110             CONTINUE
                  IF( NOUNIT )
     $               TEMP = TEMP/A( KPLUS1, J )
                  X( JX ) = TEMP
                  JX      = JX   + INCX
                  IF( J.GT.K )
     $               KX = KX + INCX
  120          CONTINUE
            END IF
         ELSE
            IF( INCX.EQ.1 )THEN
               DO 140, J = N, 1, -1
                  TEMP = X( J )
                  L    = 1      - J
                  DO 130, I = MIN( N, J + K ), J + 1, -1
                     TEMP = TEMP - A( L + I, J )*X( I )
  130             CONTINUE
                  IF( NOUNIT )
     $               TEMP = TEMP/A( 1, J )
                  X( J ) = TEMP
  140          CONTINUE
            ELSE
               KX = KX + ( N - 1 )*INCX
               JX = KX
               DO 160, J = N, 1, -1
                  TEMP = X( JX )
                  IX   = KX
                  L    = 1       - J
                  DO 150, I = MIN( N, J + K ), J + 1, -1
                     TEMP = TEMP - A( L + I, J )*X( IX )
                     IX   = IX   - INCX
  150             CONTINUE
                  IF( NOUNIT )
     $               TEMP = TEMP/A( 1, J )
                  X( JX ) = TEMP
                  JX      = JX   - INCX
                  IF( ( N - J ).GE.K )
     $               KX = KX - INCX
  160          CONTINUE
            END IF
         END IF
      END IF
*
      RETURN
*
*     End of STBSV .
*
      END
      SUBROUTINE STRSM ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA,
     $                   B, LDB )
*     .. Scalar Arguments ..
      implicit none
      CHARACTER*1        SIDE, UPLO, TRANSA, DIAG
      INTEGER            M, N, LDA, LDB
      REAL*8               ALPHA
*     .. Array Arguments ..
      REAL*8               A( LDA, * ), B( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  STRSM  solves one of the matrix equations
*
*     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
*
*  where alpha is a scalar, X and B are m by n matrices, A is a unit, or
*  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
*
*     op( A ) = A   or   op( A ) = A'.
*
*  The matrix X is overwritten on B.
*
*  Parameters
*  ==========
*
*  SIDE   - CHARACTER*1.
*           On entry, SIDE specifies whether op( A ) appears on the left
*           or right of X as follows:
*
*              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
*
*              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
*
*           Unchanged on exit.
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the matrix A is an upper or
*           lower triangular matrix as follows:
*
*              UPLO = 'U' or 'u'   A is an upper triangular matrix.
*
*              UPLO = 'L' or 'l'   A is a lower triangular matrix.
*
*           Unchanged on exit.
*
*  TRANSA - CHARACTER*1.
*           On entry, TRANSA specifies the form of op( A ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSA = 'N' or 'n'   op( A ) = A.
*
*              TRANSA = 'T' or 't'   op( A ) = A'.
*
*              TRANSA = 'C' or 'c'   op( A ) = A'.
*
*           Unchanged on exit.
*
*  DIAG   - CHARACTER*1.
*           On entry, DIAG specifies whether or not A is unit triangular
*           as follows:
*
*              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
*
*              DIAG = 'N' or 'n'   A is not assumed to be unit
*                                  triangular.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of B. M must be at
*           least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of B.  N must be
*           at least zero.
*           Unchanged on exit.
*
*  ALPHA  - REAL            .
*           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
*           zero then  A is not referenced and  B need not be set before
*           entry.
*           Unchanged on exit.
*
*  A      - REAL             array of DIMENSION ( LDA, k ), where k is m
*           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
*           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
*           upper triangular part of the array  A must contain the upper
*           triangular matrix  and the strictly lower triangular part of
*           A is not referenced.
*           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
*           lower triangular part of the array  A must contain the lower
*           triangular matrix  and the strictly upper triangular part of
*           A is not referenced.
*           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
*           A  are not referenced either,  but are assumed to be  unity.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
*           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
*           then LDA must be at least max( 1, n ).
*           Unchanged on exit.
*
*  B      - REAL             array of DIMENSION ( LDB, n ).
*           Before entry,  the leading  m by n part of the array  B must
*           contain  the  right-hand  side  matrix  B,  and  on exit  is
*           overwritten by the solution matrix  X.
*
*  LDB    - INTEGER.
*           On entry, LDB specifies the first dimension of B as declared
*           in  the  calling  (sub)  program.   LDB  must  be  at  least
*           max( 1, m ).
*           Unchanged on exit.
*
*
*  Level 3 Blas routine.
*
*
*  -- Written on 8-February-1989.
*     Jack Dongarra, Argonne National Laboratory.
*     Iain Duff, AERE Harwell.
*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*     Sven Hammarling, Numerical Algorithms Group Ltd.
*
*
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     .. Local Scalars ..
      LOGICAL            LSIDE, NOUNIT, UPPER
      INTEGER            I, INFO, J, K, NROWA
      REAL*8               TEMP
*     .. Parameters ..
      REAL*8               ONE         , ZERO
      PARAMETER        ( ONE = 1.0E+0, ZERO = 0.0E+0 )
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      LSIDE  = LSAME( SIDE  , 'L' )
      IF( LSIDE )THEN
         NROWA = M
      ELSE
         NROWA = N
      END IF
      NOUNIT = LSAME( DIAG  , 'N' )
      UPPER  = LSAME( UPLO  , 'U' )
*
      INFO   = 0
      IF(      ( .NOT.LSIDE                ).AND.
     $         ( .NOT.LSAME( SIDE  , 'R' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.UPPER                ).AND.
     $         ( .NOT.LSAME( UPLO  , 'L' ) )      )THEN
         INFO = 2
      ELSE IF( ( .NOT.LSAME( TRANSA, 'N' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'T' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'C' ) )      )THEN
         INFO = 3
      ELSE IF( ( .NOT.LSAME( DIAG  , 'U' ) ).AND.
     $         ( .NOT.LSAME( DIAG  , 'N' ) )      )THEN
         INFO = 4
      ELSE IF( M  .LT.0               )THEN
         INFO = 5
      ELSE IF( N  .LT.0               )THEN
         INFO = 6
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 9
      ELSE IF( LDB.LT.MAX( 1, M     ) )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'STRSM ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( N.EQ.0 )
     $   RETURN
*
*     And when  alpha.eq.zero.
*
      IF( ALPHA.EQ.ZERO )THEN
         DO 20, J = 1, N
            DO 10, I = 1, M
               B( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
         RETURN
      END IF
*
*     Start the operations.
*
      IF( LSIDE )THEN
         IF( LSAME( TRANSA, 'N' ) )THEN
*
*           Form  B := alpha*inv( A )*B.
*
            IF( UPPER )THEN
               DO 60, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 30, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
   30                CONTINUE
                  END IF
                  DO 50, K = M, 1, -1
                     IF( B( K, J ).NE.ZERO )THEN
                        IF( NOUNIT )
     $                     B( K, J ) = B( K, J )/A( K, K )
                        DO 40, I = 1, K - 1
                           B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
   40                   CONTINUE
                     END IF
   50             CONTINUE
   60          CONTINUE
            ELSE
               DO 100, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 70, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
   70                CONTINUE
                  END IF
                  DO 90 K = 1, M
                     IF( B( K, J ).NE.ZERO )THEN
                        IF( NOUNIT )
     $                     B( K, J ) = B( K, J )/A( K, K )
                        DO 80, I = K + 1, M
                           B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
   80                   CONTINUE
                     END IF
   90             CONTINUE
  100          CONTINUE
            END IF
         ELSE
*
*           Form  B := alpha*inv( A' )*B.
*
            IF( UPPER )THEN
               DO 130, J = 1, N
                  DO 120, I = 1, M
                     TEMP = ALPHA*B( I, J )
                     DO 110, K = 1, I - 1
                        TEMP = TEMP - A( K, I )*B( K, J )
  110                CONTINUE
                     IF( NOUNIT )
     $                  TEMP = TEMP/A( I, I )
                     B( I, J ) = TEMP
  120             CONTINUE
  130          CONTINUE
            ELSE
               DO 160, J = 1, N
                  DO 150, I = M, 1, -1
                     TEMP = ALPHA*B( I, J )
                     DO 140, K = I + 1, M
                        TEMP = TEMP - A( K, I )*B( K, J )
  140                CONTINUE
                     IF( NOUNIT )
     $                  TEMP = TEMP/A( I, I )
                     B( I, J ) = TEMP
  150             CONTINUE
  160          CONTINUE
            END IF
         END IF
      ELSE
         IF( LSAME( TRANSA, 'N' ) )THEN
*
*           Form  B := alpha*B*inv( A ).
*
            IF( UPPER )THEN
               DO 210, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 170, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
  170                CONTINUE
                  END IF
                  DO 190, K = 1, J - 1
                     IF( A( K, J ).NE.ZERO )THEN
                        DO 180, I = 1, M
                           B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
  180                   CONTINUE
                     END IF
  190             CONTINUE
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( J, J )
                     DO 200, I = 1, M
                        B( I, J ) = TEMP*B( I, J )
  200                CONTINUE
                  END IF
  210          CONTINUE
            ELSE
               DO 260, J = N, 1, -1
                  IF( ALPHA.NE.ONE )THEN
                     DO 220, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
  220                CONTINUE
                  END IF
                  DO 240, K = J + 1, N
                     IF( A( K, J ).NE.ZERO )THEN
                        DO 230, I = 1, M
                           B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
  230                   CONTINUE
                     END IF
  240             CONTINUE
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( J, J )
                     DO 250, I = 1, M
                       B( I, J ) = TEMP*B( I, J )
  250                CONTINUE
                  END IF
  260          CONTINUE
            END IF
         ELSE
*
*           Form  B := alpha*B*inv( A' ).
*
            IF( UPPER )THEN
               DO 310, K = N, 1, -1
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( K, K )
                     DO 270, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  270                CONTINUE
                  END IF
                  DO 290, J = 1, K - 1
                     IF( A( J, K ).NE.ZERO )THEN
                        TEMP = A( J, K )
                        DO 280, I = 1, M
                           B( I, J ) = B( I, J ) - TEMP*B( I, K )
  280                   CONTINUE
                     END IF
  290             CONTINUE
                  IF( ALPHA.NE.ONE )THEN
                     DO 300, I = 1, M
                        B( I, K ) = ALPHA*B( I, K )
  300                CONTINUE
                  END IF
  310          CONTINUE
            ELSE
               DO 360, K = 1, N
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( K, K )
                     DO 320, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  320                CONTINUE
                  END IF
                  DO 340, J = K + 1, N
                     IF( A( J, K ).NE.ZERO )THEN
                        TEMP = A( J, K )
                        DO 330, I = 1, M
                           B( I, J ) = B( I, J ) - TEMP*B( I, K )
  330                   CONTINUE
                     END IF
  340             CONTINUE
                  IF( ALPHA.NE.ONE )THEN
                     DO 350, I = 1, M
                        B( I, K ) = ALPHA*B( I, K )
  350                CONTINUE
                  END IF
  360          CONTINUE
            END IF
         END IF
      END IF
*
      RETURN
*
*     End of STRSM .
*
      END


      SUBROUTINE SPBTRF( UPLO, N, KD, AB, LDAB, INFO )
*
*  -- LAPACK routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993 
*
*     .. Scalar Arguments ..
      implicit none
      CHARACTER          UPLO
      INTEGER            INFO, KD, LDAB, N
*     ..
*     .. Array Arguments ..
      REAL*8               AB( LDAB, * )
*     ..
*
*  Purpose
*  =======
*
*  SPBTRF computes the Cholesky factorization of a real symmetric
*  positive definite band matrix A.
*
*  The factorization has the form
*     A = U**T * U,  if UPLO = 'U', or
*     A = L  * L**T,  if UPLO = 'L',
*  where U is an upper triangular matrix and L is lower triangular.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangle of A is stored;
*          = 'L':  Lower triangle of A is stored.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  KD      (input) INTEGER
*          The number of superdiagonals of the matrix A if UPLO = 'U',
*          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.
*
*  AB      (input/output) REAL array, dimension (LDAB,N)
*          On entry, the upper or lower triangle of the symmetric band
*          matrix A, stored in the first KD+1 rows of the array.  The
*          j-th column of A is stored in the j-th column of the array AB
*          as follows:
*          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
*          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
*
*          On exit, if INFO = 0, the triangular factor U or L from the
*          Cholesky factorization A = U**T*U or A = L*L**T of the band
*          matrix A, in the same storage format as A.
*
*  LDAB    (input) INTEGER
*          The leading dimension of the array AB.  LDAB >= KD+1.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, the leading minor of order i is not
*                positive definite, and the factorization could not be
*                completed.
*
*  Further Details
*  ===============
*
*  The band storage scheme is illustrated by the following example, when
*  N = 6, KD = 2, and UPLO = 'U':
*
*  On entry:                       On exit:
*
*      *    *   a13  a24  a35  a46      *    *   u13  u24  u35  u46
*      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
*     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
*
*  Similarly, if UPLO = 'L' the format of A is as follows:
*
*  On entry:                       On exit:
*
*     a11  a22  a33  a44  a55  a66     l11  l22  l33  l44  l55  l66
*     a21  a32  a43  a54  a65   *      l21  l32  l43  l54  l65   *
*     a31  a42  a53  a64   *    *      l31  l42  l53  l64   *    *
*
*  Array elements marked * are not used by the routine.
*
*  Contributed by
*  Peter Mayes and Giuseppe Radicati, IBM ECSEC, Rome, March 23, 1989
*
*  =====================================================================
*
*     .. Parameters ..
      REAL*8               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
      INTEGER            NBMAX, LDWORK
      PARAMETER          ( NBMAX = 32, LDWORK = NBMAX+1 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, I2, I3, IB, II, J, JJ, NB
*     ..
*     .. Local Arrays ..
      REAL*8               WORK( LDWORK, NBMAX )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           SGEMM, SPBTF2, SPOTF2, SSYRK, STRSM, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF( ( .NOT.LSAME( UPLO, 'U' ) ) .AND.
     $    ( .NOT.LSAME( UPLO, 'L' ) ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KD.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDAB.LT.KD+1 ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SPBTRF', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Determine the block size for this environment
*
      NB = ILAENV( 1, 'SPBTRF', UPLO, N, KD, -1, -1 )
*
*     The block size must not exceed the semi-bandwidth KD, and must not
*     exceed the limit set by the size of the local array WORK.
*
      NB = MIN( NB, NBMAX )
*
      IF( NB.LE.1 .OR. NB.GT.KD ) THEN
*
*        Use unblocked code
*
         CALL SPBTF2( UPLO, N, KD, AB, LDAB, INFO )
      ELSE
*
*        Use blocked code
*
         IF( LSAME( UPLO, 'U' ) ) THEN
*
*           Compute the Cholesky factorization of a symmetric band
*           matrix, given the upper triangle of the matrix in band
*           storage.
*
*           Zero the upper triangle of the work array.
*
            DO 20 J = 1, NB
               DO 10 I = 1, J - 1
                  WORK( I, J ) = ZERO
   10          CONTINUE
   20       CONTINUE
*
*           Process the band matrix one diagonal block at a time.
*
            DO 70 I = 1, N, NB
               IB = MIN( NB, N-I+1 )
*
*              Factorize the diagonal block
*
               CALL SPOTF2( UPLO, IB, AB( KD+1, I ), LDAB-1, II )
               IF( II.NE.0 ) THEN
                  INFO = I + II - 1
                  GO TO 150
               END IF
               IF( I+IB.LE.N ) THEN
*
*                 Update the relevant part of the trailing submatrix.
*                 If A11 denotes the diagonal block which has just been
*                 factorized, then we need to update the remaining
*                 blocks in the diagram:
*
*                    A11   A12   A13
*                          A22   A23
*                                A33
*
*                 The numbers of rows and columns in the partitioning
*                 are IB, I2, I3 respectively. The blocks A12, A22 and
*                 A23 are empty if IB = KD. The upper triangle of A13
*                 lies outside the band.
*
                  I2 = MIN( KD-IB, N-I-IB+1 )
                  I3 = MIN( IB, N-I-KD+1 )
*
                  IF( I2.GT.0 ) THEN
*
*                    Update A12
*
                     CALL STRSM( 'Left', 'Upper', 'Transpose',
     $                           'Non-unit', IB, I2, ONE, AB( KD+1, I ),
     $                           LDAB-1, AB( KD+1-IB, I+IB ), LDAB-1 )
*
*                    Update A22
*
                     CALL SSYRK( 'Upper', 'Transpose', I2, IB, -ONE,
     $                           AB( KD+1-IB, I+IB ), LDAB-1, ONE,
     $                           AB( KD+1, I+IB ), LDAB-1 )
                  END IF
*
                  IF( I3.GT.0 ) THEN
*
*                    Copy the lower triangle of A13 into the work array.
*
                     DO 40 JJ = 1, I3
                        DO 30 II = JJ, IB
                           WORK( II, JJ ) = AB( II-JJ+1, JJ+I+KD-1 )
   30                   CONTINUE
   40                CONTINUE
*
*                    Update A13 (in the work array).
*
                     CALL STRSM( 'Left', 'Upper', 'Transpose',
     $                           'Non-unit', IB, I3, ONE, AB( KD+1, I ),
     $                           LDAB-1, WORK, LDWORK )
*
*                    Update A23
*
                     IF( I2.GT.0 )
     $                  CALL SGEMM( 'Transpose', 'No Transpose', I2, I3,
     $                              IB, -ONE, AB( KD+1-IB, I+IB ),
     $                              LDAB-1, WORK, LDWORK, ONE,
     $                              AB( 1+IB, I+KD ), LDAB-1 )
*
*                    Update A33
*
                     CALL SSYRK( 'Upper', 'Transpose', I3, IB, -ONE,
     $                           WORK, LDWORK, ONE, AB( KD+1, I+KD ),
     $                           LDAB-1 )
*
*                    Copy the lower triangle of A13 back into place.
*
                     DO 60 JJ = 1, I3
                        DO 50 II = JJ, IB
                           AB( II-JJ+1, JJ+I+KD-1 ) = WORK( II, JJ )
   50                   CONTINUE
   60                CONTINUE
                  END IF
               END IF
   70       CONTINUE
         ELSE
*
*           Compute the Cholesky factorization of a symmetric band
*           matrix, given the lower triangle of the matrix in band
*           storage.
*
*           Zero the lower triangle of the work array.
*
            DO 90 J = 1, NB
               DO 80 I = J + 1, NB
                  WORK( I, J ) = ZERO
   80          CONTINUE
   90       CONTINUE
*
*           Process the band matrix one diagonal block at a time.
*
            DO 140 I = 1, N, NB
               IB = MIN( NB, N-I+1 )
*
*              Factorize the diagonal block
*
               CALL SPOTF2( UPLO, IB, AB( 1, I ), LDAB-1, II )
               IF( II.NE.0 ) THEN
                  INFO = I + II - 1
                  GO TO 150
               END IF
               IF( I+IB.LE.N ) THEN
*
*                 Update the relevant part of the trailing submatrix.
*                 If A11 denotes the diagonal block which has just been
*                 factorized, then we need to update the remaining
*                 blocks in the diagram:
*
*                    A11
*                    A21   A22
*                    A31   A32   A33
*
*                 The numbers of rows and columns in the partitioning
*                 are IB, I2, I3 respectively. The blocks A21, A22 and
*                 A32 are empty if IB = KD. The lower triangle of A31
*                 lies outside the band.
*
                  I2 = MIN( KD-IB, N-I-IB+1 )
                  I3 = MIN( IB, N-I-KD+1 )
*
                  IF( I2.GT.0 ) THEN
*
*                    Update A21
*
                     CALL STRSM( 'Right', 'Lower', 'Transpose',
     $                           'Non-unit', I2, IB, ONE, AB( 1, I ),
     $                           LDAB-1, AB( 1+IB, I ), LDAB-1 )
*
*                    Update A22
*
                     CALL SSYRK( 'Lower', 'No Transpose', I2, IB, -ONE,
     $                           AB( 1+IB, I ), LDAB-1, ONE,
     $                           AB( 1, I+IB ), LDAB-1 )
                  END IF
*
                  IF( I3.GT.0 ) THEN
*
*                    Copy the upper triangle of A31 into the work array.
*
                     DO 110 JJ = 1, IB
                        DO 100 II = 1, MIN( JJ, I3 )
                           WORK( II, JJ ) = AB( KD+1-JJ+II, JJ+I-1 )
  100                   CONTINUE
  110                CONTINUE
*
*                    Update A31 (in the work array).
*
                     CALL STRSM( 'Right', 'Lower', 'Transpose',
     $                           'Non-unit', I3, IB, ONE, AB( 1, I ),
     $                           LDAB-1, WORK, LDWORK )
*
*                    Update A32
*
                     IF( I2.GT.0 )
     $                  CALL SGEMM( 'No transpose', 'Transpose', I3, I2,
     $                              IB, -ONE, WORK, LDWORK,
     $                              AB( 1+IB, I ), LDAB-1, ONE,
     $                              AB( 1+KD-IB, I+IB ), LDAB-1 )
*
*                    Update A33
*
                     CALL SSYRK( 'Lower', 'No Transpose', I3, IB, -ONE,
     $                           WORK, LDWORK, ONE, AB( 1, I+KD ),
     $                           LDAB-1 )
*
*                    Copy the upper triangle of A31 back into place.
*
                     DO 130 JJ = 1, IB
                        DO 120 II = 1, MIN( JJ, I3 )
                           AB( KD+1-JJ+II, JJ+I-1 ) = WORK( II, JJ )
  120                   CONTINUE
  130                CONTINUE
                  END IF
               END IF
  140       CONTINUE
         END IF
      END IF
      RETURN
*
  150 CONTINUE
      RETURN
*
*     End of SPBTRF
*
      END

      SUBROUTINE SPBTRS( UPLO, N, KD, NRHS, AB, LDAB, B, LDB, INFO )
*
*  -- LAPACK routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      implicit none
      CHARACTER          UPLO
      INTEGER            INFO, KD, LDAB, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      REAL*8               AB( LDAB, * ), B( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  SPBTRS solves a system of linear equations A*X = B with a symmetric
*  positive definite band matrix A using the Cholesky factorization
*  A = U**T*U or A = L*L**T computed by SPBTRF.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangular factor stored in AB;
*          = 'L':  Lower triangular factor stored in AB.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  KD      (input) INTEGER
*          The number of superdiagonals of the matrix A if UPLO = 'U',
*          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.
*
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the matrix B.  NRHS >= 0.
*
*  AB      (input) REAL array, dimension (LDAB,N)
*          The triangular factor U or L from the Cholesky factorization
*          A = U**T*U or A = L*L**T of the band matrix A, stored in the
*          first KD+1 rows of the array.  The j-th column of U or L is
*          stored in the j-th column of the array AB as follows:
*          if UPLO ='U', AB(kd+1+i-j,j) = U(i,j) for max(1,j-kd)<=i<=j;
*          if UPLO ='L', AB(1+i-j,j)    = L(i,j) for j<=i<=min(n,j+kd).
*
*  LDAB    (input) INTEGER
*          The leading dimension of the array AB.  LDAB >= KD+1.
*
*  B       (input/output) REAL array, dimension (LDB,NRHS)
*          On entry, the right hand side matrix B.
*          On exit, the solution matrix X.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            J
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           STBSV, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KD.LT.0 ) THEN
         INFO = -3
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDAB.LT.KD+1 ) THEN
         INFO = -6
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SPBTRS', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. NRHS.EQ.0 )
     $   RETURN
*
      IF( UPPER ) THEN
*
*        Solve A*X = B where A = U'*U.
*
         DO 10 J = 1, NRHS
*
*           Solve U'*X = B, overwriting B with X.
*
            CALL STBSV( 'Upper', 'Transpose', 'Non-unit', N, KD, AB,
     $                  LDAB, B( 1, J ), 1 )
*
*           Solve U*X = B, overwriting B with X.
*
            CALL STBSV( 'Upper', 'No transpose', 'Non-unit', N, KD, AB,
     $                  LDAB, B( 1, J ), 1 )
   10    CONTINUE
      ELSE
*
*        Solve A*X = B where A = L*L'.
*
         DO 20 J = 1, NRHS
*
*           Solve L*X = B, overwriting B with X.
*
            CALL STBSV( 'Lower', 'No transpose', 'Non-unit', N, KD, AB,
     $                  LDAB, B( 1, J ), 1 )
*
*           Solve L'*X = B, overwriting B with X.
*
            CALL STBSV( 'Lower', 'Transpose', 'Non-unit', N, KD, AB,
     $                  LDAB, B( 1, J ), 1 )
   20    CONTINUE
      END IF
*
      RETURN
*
*     End of SPBTRS
*
      END

      SUBROUTINE SPBTF2( UPLO, N, KD, AB, LDAB, INFO )
*
*  -- LAPACK routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      implicit none
      CHARACTER          UPLO
      INTEGER            INFO, KD, LDAB, N
*     ..
*     .. Array Arguments ..
      REAL*8               AB( LDAB, * )
*     ..
*
*  Purpose
*  =======
*
*  SPBTF2 computes the Cholesky factorization of a real symmetric
*  positive definite band matrix A.
*
*  The factorization has the form
*     A = U' * U ,  if UPLO = 'U', or
*     A = L  * L',  if UPLO = 'L',
*  where U is an upper triangular matrix, U' is the transpose of U, and
*  L is lower triangular.
*
*  This is the unblocked version of the algorithm, calling Level 2 BLAS.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          Specifies whether the upper or lower triangular part of the
*          symmetric matrix A is stored:
*          = 'U':  Upper triangular
*          = 'L':  Lower triangular
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  KD      (input) INTEGER
*          The number of super-diagonals of the matrix A if UPLO = 'U',
*          or the number of sub-diagonals if UPLO = 'L'.  KD >= 0.
*
*  AB      (input/output) REAL array, dimension (LDAB,N)
*          On entry, the upper or lower triangle of the symmetric band
*          matrix A, stored in the first KD+1 rows of the array.  The
*          j-th column of A is stored in the j-th column of the array AB
*          as follows:
*          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
*          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
*
*          On exit, if INFO = 0, the triangular factor U or L from the
*          Cholesky factorization A = U'*U or A = L*L' of the band
*          matrix A, in the same storage format as A.
*
*  LDAB    (input) INTEGER
*          The leading dimension of the array AB.  LDAB >= KD+1.
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -k, the k-th argument had an illegal value
*          > 0: if INFO = k, the leading minor of order k is not
*               positive definite, and the factorization could not be
*               completed.
*
*  Further Details
*  ===============
*
*  The band storage scheme is illustrated by the following example, when
*  N = 6, KD = 2, and UPLO = 'U':
*
*  On entry:                       On exit:
*
*      *    *   a13  a24  a35  a46      *    *   u13  u24  u35  u46
*      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
*     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
*
*  Similarly, if UPLO = 'L' the format of A is as follows:
*
*  On entry:                       On exit:
*
*     a11  a22  a33  a44  a55  a66     l11  l22  l33  l44  l55  l66
*     a21  a32  a43  a54  a65   *      l21  l32  l43  l54  l65   *
*     a31  a42  a53  a64   *    *      l31  l42  l53  l64   *    *
*
*  Array elements marked * are not used by the routine.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL*8               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            J, KLD, KN
      REAL*8               AJJ
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           SSCAL, SSYR, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KD.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDAB.LT.KD+1 ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SPBTF2', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      KLD = MAX( 1, LDAB-1 )
*
      IF( UPPER ) THEN
*
*        Compute the Cholesky factorization A = U'*U.
*
         DO 10 J = 1, N
*
*           Compute U(J,J) and test for non-positive-definiteness.
*
            AJJ = AB( KD+1, J )
            IF( AJJ.LE.ZERO )
     $         GO TO 30
            AJJ = SQRT( AJJ )
            AB( KD+1, J ) = AJJ
*
*           Compute elements J+1:J+KN of row J and update the
*           trailing submatrix within the band.
*
            KN = MIN( KD, N-J )
            IF( KN.GT.0 ) THEN
               CALL SSCAL( KN, ONE / AJJ, AB( KD, J+1 ), KLD )
               CALL SSYR( 'Upper', KN, -ONE, AB( KD, J+1 ), KLD,
     $                    AB( KD+1, J+1 ), KLD )
            END IF
   10    CONTINUE
      ELSE
*
*        Compute the Cholesky factorization A = L*L'.
*
         DO 20 J = 1, N
*
*           Compute L(J,J) and test for non-positive-definiteness.
*
            AJJ = AB( 1, J )
            IF( AJJ.LE.ZERO )
     $         GO TO 30
            AJJ = SQRT( AJJ )
            AB( 1, J ) = AJJ
*
*           Compute elements J+1:J+KN of column J and update the
*           trailing submatrix within the band.
*
            KN = MIN( KD, N-J )
            IF( KN.GT.0 ) THEN
               CALL SSCAL( KN, ONE / AJJ, AB( 2, J ), 1 )
               CALL SSYR( 'Lower', KN, -ONE, AB( 2, J ), 1,
     $                    AB( 1, J+1 ), KLD )
            END IF
   20    CONTINUE
      END IF
      RETURN
*
   30 CONTINUE
      INFO = J
      RETURN
*
*     End of SPBTF2
*
      END
      SUBROUTINE SPOTF2( UPLO, N, A, LDA, INFO )
*
*  -- LAPACK routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      implicit none
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      REAL*8               A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  SPOTF2 computes the Cholesky factorization of a real symmetric
*  positive definite matrix A.
*
*  The factorization has the form
*     A = U' * U ,  if UPLO = 'U', or
*     A = L  * L',  if UPLO = 'L',
*  where U is an upper triangular matrix and L is lower triangular.
*
*  This is the unblocked version of the algorithm, calling Level 2 BLAS.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          Specifies whether the upper or lower triangular part of the
*          symmetric matrix A is stored.
*          = 'U':  Upper triangular
*          = 'L':  Lower triangular
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) REAL array, dimension (LDA,N)
*          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
*          n by n upper triangular part of A contains the upper
*          triangular part of the matrix A, and the strictly lower
*          triangular part of A is not referenced.  If UPLO = 'L', the
*          leading n by n lower triangular part of A contains the lower
*          triangular part of the matrix A, and the strictly upper
*          triangular part of A is not referenced.
*
*          On exit, if INFO = 0, the factor U or L from the Cholesky
*          factorization A = U'*U  or A = L*L'.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -k, the k-th argument had an illegal value
*          > 0: if INFO = k, the leading minor of order k is not
*               positive definite, and the factorization could not be
*               completed.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL*8               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            J
      REAL*8               AJJ
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      REAL*8               SDOT
      EXTERNAL           LSAME, SDOT
*     ..
*     .. External Subroutines ..
      EXTERNAL           SGEMV, SSCAL, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SPOTF2', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      IF( UPPER ) THEN
*
*        Compute the Cholesky factorization A = U'*U.
*
         DO 10 J = 1, N
*
*           Compute U(J,J) and test for non-positive-definiteness.
*
            AJJ = A( J, J ) - SDOT( J-1, A( 1, J ), 1, A( 1, J ), 1 )
            IF( AJJ.LE.ZERO ) THEN
               A( J, J ) = AJJ
               GO TO 30
            END IF
            AJJ = SQRT( AJJ )
            A( J, J ) = AJJ
*
*           Compute elements J+1:N of row J.
*
            IF( J.LT.N ) THEN
               CALL SGEMV( 'Transpose', J-1, N-J, -ONE, A( 1, J+1 ),
     $                     LDA, A( 1, J ), 1, ONE, A( J, J+1 ), LDA )
               CALL SSCAL( N-J, ONE / AJJ, A( J, J+1 ), LDA )
            END IF
   10    CONTINUE
      ELSE
*
*        Compute the Cholesky factorization A = L*L'.
*
         DO 20 J = 1, N
*
*           Compute L(J,J) and test for non-positive-definiteness.
*
            AJJ = A( J, J ) - SDOT( J-1, A( J, 1 ), LDA, A( J, 1 ),
     $            LDA )
            IF( AJJ.LE.ZERO ) THEN
               A( J, J ) = AJJ
               GO TO 30
            END IF
            AJJ = SQRT( AJJ )
            A( J, J ) = AJJ
*
*           Compute elements J+1:N of column J.
*
            IF( J.LT.N ) THEN
               CALL SGEMV( 'No transpose', N-J, J-1, -ONE, A( J+1, 1 ),
     $                     LDA, A( J, 1 ), LDA, ONE, A( J+1, J ), 1 )
               CALL SSCAL( N-J, ONE / AJJ, A( J+1, J ), 1 )
            END IF
   20    CONTINUE
      END IF
      GO TO 40
*
   30 CONTINUE
      INFO = J
*
   40 CONTINUE
      RETURN
*
*     End of SPOTF2
*
      END

      real*8 function sdot(n,sx,incx,sy,incy)
c
c     forms the dot product of two vectors.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      implicit none
      real*8 sx(*),sy(*),stemp
      integer i,incx,incy,ix,iy,m,mp1,n
c
      stemp = 0.0e0
      sdot = 0.0e0
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        stemp = stemp + sx(ix)*sy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      sdot = stemp
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        stemp = stemp + sx(i)*sy(i)
   30 continue
      if( n .lt. 5 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        stemp = stemp + sx(i)*sy(i) + sx(i + 1)*sy(i + 1) +
     *   sx(i + 2)*sy(i + 2) + sx(i + 3)*sy(i + 3) + sx(i + 4)*sy(i + 4)
   50 continue
   60 sdot = stemp
      return
      end

      SUBROUTINE SSYR  ( UPLO, N, ALPHA, X, INCX, A, LDA )
*     .. Scalar Arguments ..
      implicit none
      REAL*8               ALPHA
      INTEGER            INCX, LDA, N
      CHARACTER*1        UPLO
*     .. Array Arguments ..
      REAL*8               A( LDA, * ), X( * )
*     ..
*
*  Purpose
*  =======
*
*  SSYR   performs the symmetric rank 1 operation
*
*     A := alpha*x*x' + A,
*
*  where alpha is a real scalar, x is an n element vector and A is an
*  n by n symmetric matrix.
*
*  Parameters
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the upper or lower
*           triangular part of the array A is to be referenced as
*           follows:
*
*              UPLO = 'U' or 'u'   Only the upper triangular part of A
*                                  is to be referenced.
*
*              UPLO = 'L' or 'l'   Only the lower triangular part of A
*                                  is to be referenced.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the order of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - REAL            .
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  X      - REAL             array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the n
*           element vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*  A      - REAL             array of DIMENSION ( LDA, n ).
*           Before entry with  UPLO = 'U' or 'u', the leading n by n
*           upper triangular part of the array A must contain the upper
*           triangular part of the symmetric matrix and the strictly
*           lower triangular part of A is not referenced. On exit, the
*           upper triangular part of the array A is overwritten by the
*           upper triangular part of the updated matrix.
*           Before entry with UPLO = 'L' or 'l', the leading n by n
*           lower triangular part of the array A must contain the lower
*           triangular part of the symmetric matrix and the strictly
*           upper triangular part of A is not referenced. On exit, the
*           lower triangular part of the array A is overwritten by the
*           lower triangular part of the updated matrix.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           max( 1, n ).
*           Unchanged on exit.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*
*     .. Parameters ..
      REAL*8               ZERO
      PARAMETER        ( ZERO = 0.0E+0 )
*     .. Local Scalars ..
      REAL*8               TEMP
      INTEGER            I, INFO, IX, J, JX, KX
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF     ( .NOT.LSAME( UPLO, 'U' ).AND.
     $         .NOT.LSAME( UPLO, 'L' )      )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, N ) )THEN
         INFO = 7
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'SSYR  ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( N.EQ.0 ).OR.( ALPHA.EQ.ZERO ) )
     $   RETURN
*
*     Set the start point in X if the increment is not unity.
*
      IF( INCX.LE.0 )THEN
         KX = 1 - ( N - 1 )*INCX
      ELSE IF( INCX.NE.1 )THEN
         KX = 1
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through the triangular part
*     of A.
*
      IF( LSAME( UPLO, 'U' ) )THEN
*
*        Form  A  when A is stored in upper triangle.
*
         IF( INCX.EQ.1 )THEN
            DO 20, J = 1, N
               IF( X( J ).NE.ZERO )THEN
                  TEMP = ALPHA*X( J )
                  DO 10, I = 1, J
                     A( I, J ) = A( I, J ) + X( I )*TEMP
   10             CONTINUE
               END IF
   20       CONTINUE
         ELSE
            JX = KX
            DO 40, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  IX   = KX
                  DO 30, I = 1, J
                     A( I, J ) = A( I, J ) + X( IX )*TEMP
                     IX        = IX        + INCX
   30             CONTINUE
               END IF
               JX = JX + INCX
   40       CONTINUE
         END IF
      ELSE
*
*        Form  A  when A is stored in lower triangle.
*
         IF( INCX.EQ.1 )THEN
            DO 60, J = 1, N
               IF( X( J ).NE.ZERO )THEN
                  TEMP = ALPHA*X( J )
                  DO 50, I = J, N
                     A( I, J ) = A( I, J ) + X( I )*TEMP
   50             CONTINUE
               END IF
   60       CONTINUE
         ELSE
            JX = KX
            DO 80, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  IX   = JX
                  DO 70, I = J, N
                     A( I, J ) = A( I, J ) + X( IX )*TEMP
                     IX        = IX        + INCX
   70             CONTINUE
               END IF
               JX = JX + INCX
   80       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of SSYR  .
*
      END

      SUBROUTINE SSYRK ( UPLO, TRANS, N, K, ALPHA, A, LDA,
     $                   BETA, C, LDC )
*     .. Scalar Arguments ..
      implicit none
      CHARACTER*1        UPLO, TRANS
      INTEGER            N, K, LDA, LDC
      REAL*8               ALPHA, BETA
*     .. Array Arguments ..
      REAL*8               A( LDA, * ), C( LDC, * )
*     ..
*
*  Purpose
*  =======
*
*  SSYRK  performs one of the symmetric rank k operations
*
*     C := alpha*A*A' + beta*C,
*
*  or
*
*     C := alpha*A'*A + beta*C,
*
*  where  alpha and beta  are scalars, C is an  n by n  symmetric matrix
*  and  A  is an  n by k  matrix in the first case and a  k by n  matrix
*  in the second case.
*
*  Parameters
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On  entry,   UPLO  specifies  whether  the  upper  or  lower
*           triangular  part  of the  array  C  is to be  referenced  as
*           follows:
*
*              UPLO = 'U' or 'u'   Only the  upper triangular part of  C
*                                  is to be referenced.
*
*              UPLO = 'L' or 'l'   Only the  lower triangular part of  C
*                                  is to be referenced.
*
*           Unchanged on exit.
*
*  TRANS  - CHARACTER*1.
*           On entry,  TRANS  specifies the operation to be performed as
*           follows:
*
*              TRANS = 'N' or 'n'   C := alpha*A*A' + beta*C.
*
*              TRANS = 'T' or 't'   C := alpha*A'*A + beta*C.
*
*              TRANS = 'C' or 'c'   C := alpha*A'*A + beta*C.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry,  N specifies the order of the matrix C.  N must be
*           at least zero.
*           Unchanged on exit.
*
*  K      - INTEGER.
*           On entry with  TRANS = 'N' or 'n',  K  specifies  the number
*           of  columns   of  the   matrix   A,   and  on   entry   with
*           TRANS = 'T' or 't' or 'C' or 'c',  K  specifies  the  number
*           of rows of the matrix  A.  K must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - REAL            .
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - REAL             array of DIMENSION ( LDA, ka ), where ka is
*           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
*           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
*           part of the array  A  must contain the matrix  A,  otherwise
*           the leading  k by n  part of the array  A  must contain  the
*           matrix A.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
*           then  LDA must be at least  max( 1, n ), otherwise  LDA must
*           be at least  max( 1, k ).
*           Unchanged on exit.
*
*  BETA   - REAL            .
*           On entry, BETA specifies the scalar beta.
*           Unchanged on exit.
*
*  C      - REAL             array of DIMENSION ( LDC, n ).
*           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n
*           upper triangular part of the array C must contain the upper
*           triangular part  of the  symmetric matrix  and the strictly
*           lower triangular part of C is not referenced.  On exit, the
*           upper triangular part of the array  C is overwritten by the
*           upper triangular part of the updated matrix.
*           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n
*           lower triangular part of the array C must contain the lower
*           triangular part  of the  symmetric matrix  and the strictly
*           upper triangular part of C is not referenced.  On exit, the
*           lower triangular part of the array  C is overwritten by the
*           lower triangular part of the updated matrix.
*
*  LDC    - INTEGER.
*           On entry, LDC specifies the first dimension of C as declared
*           in  the  calling  (sub)  program.   LDC  must  be  at  least
*           max( 1, n ).
*           Unchanged on exit.
*
*
*  Level 3 Blas routine.
*
*  -- Written on 8-February-1989.
*     Jack Dongarra, Argonne National Laboratory.
*     Iain Duff, AERE Harwell.
*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*     Sven Hammarling, Numerical Algorithms Group Ltd.
*
*
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I, INFO, J, L, NROWA
      REAL*8               TEMP
*     .. Parameters ..
      REAL*8               ONE ,         ZERO
      PARAMETER        ( ONE = 1.0E+0, ZERO = 0.0E+0 )
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      IF( LSAME( TRANS, 'N' ) )THEN
         NROWA = N
      ELSE
         NROWA = K
      END IF
      UPPER = LSAME( UPLO, 'U' )
*
      INFO = 0
      IF(      ( .NOT.UPPER               ).AND.
     $         ( .NOT.LSAME( UPLO , 'L' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.LSAME( TRANS, 'N' ) ).AND.
     $         ( .NOT.LSAME( TRANS, 'T' ) ).AND.
     $         ( .NOT.LSAME( TRANS, 'C' ) )      )THEN
         INFO = 2
      ELSE IF( N  .LT.0               )THEN
         INFO = 3
      ELSE IF( K  .LT.0               )THEN
         INFO = 4
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 7
      ELSE IF( LDC.LT.MAX( 1, N     ) )THEN
         INFO = 10
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'SSYRK ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( N.EQ.0 ).OR.
     $    ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
*
*     And when  alpha.eq.zero.
*
      IF( ALPHA.EQ.ZERO )THEN
         IF( UPPER )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 20, J = 1, N
                  DO 10, I = 1, J
                     C( I, J ) = ZERO
   10             CONTINUE
   20          CONTINUE
            ELSE
               DO 40, J = 1, N
                  DO 30, I = 1, J
                     C( I, J ) = BETA*C( I, J )
   30             CONTINUE
   40          CONTINUE
            END IF
         ELSE
            IF( BETA.EQ.ZERO )THEN
               DO 60, J = 1, N
                  DO 50, I = J, N
                     C( I, J ) = ZERO
   50             CONTINUE
   60          CONTINUE
            ELSE
               DO 80, J = 1, N
                  DO 70, I = J, N
                     C( I, J ) = BETA*C( I, J )
   70             CONTINUE
   80          CONTINUE
            END IF
         END IF
         RETURN
      END IF
*
*     Start the operations.
*
      IF( LSAME( TRANS, 'N' ) )THEN
*
*        Form  C := alpha*A*A' + beta*C.
*
         IF( UPPER )THEN
            DO 130, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 90, I = 1, J
                     C( I, J ) = ZERO
   90             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 100, I = 1, J
                     C( I, J ) = BETA*C( I, J )
  100             CONTINUE
               END IF
               DO 120, L = 1, K
                  IF( A( J, L ).NE.ZERO )THEN
                     TEMP = ALPHA*A( J, L )
                     DO 110, I = 1, J
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
  110                CONTINUE
                  END IF
  120          CONTINUE
  130       CONTINUE
         ELSE
            DO 180, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 140, I = J, N
                     C( I, J ) = ZERO
  140             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 150, I = J, N
                     C( I, J ) = BETA*C( I, J )
  150             CONTINUE
               END IF
               DO 170, L = 1, K
                  IF( A( J, L ).NE.ZERO )THEN
                     TEMP      = ALPHA*A( J, L )
                     DO 160, I = J, N
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
  160                CONTINUE
                  END IF
  170          CONTINUE
  180       CONTINUE
         END IF
      ELSE
*
*        Form  C := alpha*A'*A + beta*C.
*
         IF( UPPER )THEN
            DO 210, J = 1, N
               DO 200, I = 1, J
                  TEMP = ZERO
                  DO 190, L = 1, K
                     TEMP = TEMP + A( L, I )*A( L, J )
  190             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  200          CONTINUE
  210       CONTINUE
         ELSE
            DO 240, J = 1, N
               DO 230, I = J, N
                  TEMP = ZERO
                  DO 220, L = 1, K
                     TEMP = TEMP + A( L, I )*A( L, J )
  220             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  230          CONTINUE
  240       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of SSYRK .
*
      END
c
c.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
c
      subroutine copy8tor(px8,px,kn)
c
c   copy real*8 array to real array
c
      implicit none
      integer kn, i
      real*8 px8(kn)
      real*8   px (kn)
c
      do i=1,kn
        px(i) = px8(i)
      end do
c
      return
      end
c     
      subroutine copyrto8(px,px8,kn)
c
c   copy real*8 array to real array
c
      implicit none
      integer kn, i
      real*8 px8(kn)
      real*8 px (kn)
c
      do i=1,kn
        px8(i) = px(i)
      end do
c
      return
      end
      integer*4 function malloc_f(nbytes,infomalloc)
      implicit none
      integer nbytes, malloc, nbyt_tot, maxmem,infomalloc
      parameter(maxmem = 500e6)
      common /memoryuse/ nbyt_tot
c%OS      save nbyt_tot
      if (nbyt_tot + nbytes .gt. maxmem) then
        print *,' Extra ',nbytes,
     +    ' bytes not allocated as total memory = ',nbyt_tot + nbytes,
     +    ' would be larger than allowed maxmem = ',maxmem
        infomalloc=1
      else
        nbyt_tot = nbyt_tot + nbytes
c%OS        print *,' malloc bytes so far: ',nbyt_tot/1e6,' Mb'
C   SUN
c%OS        malloc_f = malloc(nbytes)
C   IBM
        malloc_f = malloc(%val(nbytes))
c
        infomalloc=0
      endif
c%OS      print *,malloc_f
      return
      end
