    program cbspltest
      USE prec_rkind
      USE interpos
      implicit none
      integer, parameter :: nin=1000, nout=3000
      real(rkind) :: xin(nin), yin(nin), yinspl(nin), yinsplpp(nin)
      real(rkind) :: xout(nout), yout(nout), youtp(nout), youtpp(nout)
      integer :: i, ioptder, iextrapol, nbc(2), iflag, icount, icount2, irate, icmax
      real(rkind) :: taus, ybc(6), sigma(nin)
      REAL(RKIND), PARAMETER :: TWOPI=6.283185307179586476925286766559005768394_rkind
!
      do i=1,nin
         xin(i)=real(i-1,rkind)/real(nin-1,rkind)
         yin(i)=1._rkind + cos(xin(i)**2 * TWOPI)
      end do
      do i=1,nout
         xout(i)=real(i-1,rkind)/real(nout-1,rkind)
      end do
      !
      ! prepare parameters for cbsplgen0
      ioptder = 2
      iextrapol = 1
      taus = 1.0e-05_rkind
      sigma = taus
      nbc = 0
      ybc = 0.0_rkind
      ! more complicated option
!      nbc(1)=2
!      ybc(1)=1._rkind
      ! call spline with tension
      call system_clock(icount, irate, icmax)
      write(7,*) ' icount, irate, icmax= ',icount, irate, icmax
      call cbsplgen0(xin,yin,nin,xout,yout,youtp,youtpp,nout, &
           & ioptder,iextrapol,sigma,nbc,ybc,iflag)
      call system_clock(icount2, irate, icmax)
      write(7,*) ' icount2, irate, icmax, icount2-icount= ',icount2, irate, icmax, icount2-icount
      write(7,*) ' dtime = ',real(icount2-icount,rkind) / real(irate,rkind)
!
      print *,'% iflag = ',iflag
      print *,' in=[ '
      write(*,'(1p2e15.6)') (xin(i), yin(i), i=1,nin)
      print *,'];'
      print *,' out=[ '
      write(*,'(1p4e15.6)') (xout(i), yout(i), youtp(i), youtpp(i), i=1,nout)
      print *,'];'

      ! test using splibnd
      ! first calculated spline function and 2nd der: yinspl and yinsplpp on xin
      call system_clock(icount)
      call cbsplgen0(xin,yin,nin,xin,yinspl,youtp,yinsplpp,nin, &
           & ioptder,iextrapol,sigma,nbc,ybc,iflag)
      call system_clock(icount2)
      write(7,*) ' icount2-icount, dtime = ',icount2-icount,real(icount2-icount,rkind) / real(irate,rkind)
      ! Then can calculate spline on any new x value
      call system_clock(icount)
      do i=1,nout
         CALL SPLIBND(xin,yinspl,yinsplpp,nin,xout(i),yout(i), youtp(i), youtpp(i), &
     &      iextrapol)
      end do
      call system_clock(icount2)
      write(7,*) ' icount2-icount, dtime = ',icount2-icount,real(icount2-icount,rkind) / real(irate,rkind)
!
      print *,' outsplibnd=[ '
      write(*,'(1p4e15.6)') (xout(i), yout(i), youtp(i), youtpp(i), i=1,nout)
      print *,'];'
!
    end program cbspltest
