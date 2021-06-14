MODULE interpos_module
  !
  ! Examples of calls to this generic routine:
  !
  ! captured by interpos_def: (not clear can have optional youtp and then nout (?), if yes include kopt as 1st option
  ! call interpos(xin,yin,nin,xout,yout,nout)
  ! call interpos(xin,yin,nin,xout,yout,nout,taus)
  ! call interpos(xin,yin,nin,xout,yout,youtp,nout[,taus])
  ! call interpos(xin,yin,nin,xout,yout,youtp,youtpp,nout[,taus])
  ! call interpos(xin,yin,nin,xout,yout,youtp,youtpp,youtint,nout[,taus])
  ! below can be included if optional in middle possible, but then probalem if keywords given?
  ! call interpos(xin,yin,nin,xout,yout,youtp,youtpp,youtint,nout[,taus,nbc,ybc,sigma])
  !
  ! captured by interpos_def2: (not clear can have optional youtp and then nout (?)
  ! call interpos(xin,yin,yinnew,yinpp,nin[,taus])  to just compute arrays for later multiple spline interpolations
  ! call interpos(xin,yin,yinnew,yinpp,nin,xout,yout,nout[,taus])
  ! call interpos(xin,yin,yinnew,yinpp,nin,xout,yout,youtp,nout[,taus])
  ! call interpos(xin,yin,yinnew,yinpp,nin,xout,yout,youtp,youtpp,nout[,taus])
  ! call interpos(xin,yin,yinnew,yinpp,nin,xout,yout,youtp,youtpp,youtint,nout[,taus])
  !
  ! captured by interpos_def3: yinpp already computed, so just compute interpolation on xout
  ! call interpos(xin,yinnew,yinpp,nin,xout,yout,nout[,taus])
  ! call interpos(xin,yinnew,yinpp,nin,xout,yout[,youtp,youtpp,youtint],nout[,taus])
  !
  ! captured by interpos_full: (not clear can have optional youtp and then nout (?)
  ! call interpos(xin,yin,yinnew,yinpp,nin,xout,yout,youtp,youtpp,youtint,nout[,taus],nbc,ybc,sigma)
  !
  ! To test if gives keyword can you check with "present(youtpp)" if needed?
  !
  implicit none
  public :: interpos
  interface interpos
     module procedure interpos_def_xoutscal, interpos_defp_xoutscal, interpos_defint, interpos_defp, &
          & interpos_sfull, interpos_vfull, &
          & interpos_ssimple, interpos_vsimple
  end interface

contains

!!$  SUBROUTINE interpos_def(XIN,YIN,nin,xout,yout,youtp,youtpp,youtint,nout,taus)
!!$
!!$  !
!!$  SUBROUTINE interpos_def2(XIN,YIN,YNEW,YINPP,nin,xout,yout,nout,taus)
!!$
!!$  SUBROUTINE interpos_def2p(XIN,YIN,YNEW,YINPP,nin,xout,yout,youtp,nout,taus)
!!$  !
!!$  SUBROUTINE interpos_def2pp(XIN,YIN,YNEW,YINPP,nin,xout,yout,youtp,youtpp,nout,taus)
!!$  !
!!$  SUBROUTINE interpos_def2ppint(XIN,YIN,YNEW,YINPP,nin,xout,yout,youtp,youtpp,youtint,nout,taus)
!!$  !
!!$    SUBROUTINE interpos_sfull(XIN,YIN,YNEW,YINPP,nin,xout,yout,youtp,youtpp,youtint, &
!!$         &  nout,ioptder,iextrapo,taus,nbc,ybc,IFLAG)
!!$  !
!!$  SUBROUTINE interpos_sfull(XIN,YIN,YNEW,YINPP,nin,xout,yout,youtp,youtpp,youtint, &
!!$  &  nout,ioptder,iextrapo,taus,nbc,ybc,IFLAG)
!!$  !
!!$  SUBROUTINE interpos_sfull(XIN,YIN,YNEW,YINPP,nin,xout,yout,youtp,youtpp,youtint, &
!!$  &  nout,ioptder,iextrapo,taus,nbc,ybc,IFLAG)
!!$  !
!!$  SUBROUTINE interpos_sfull(XIN,YIN,YNEW,YINPP,nin,xout,yout,youtp,youtpp,youtint, &
!!$  &  nout,ioptder,iextrapo,taus,nbc,ybc,IFLAG)
!!$  !
!!$  SUBROUTINE interpos_sfull(XIN,YIN,YNEW,YINPP,nin,xout,yout,youtp,youtpp,youtint, &
!!$  &  nout,ioptder,iextrapo,taus,nbc,ybc,IFLAG)
!!$  !
!!$  SUBROUTINE interpos_def(XIN,YIN,nin,xout,yout,nout,taus)

  SUBROUTINE interpos_defp_xoutscal(XIN,YIN,nin,xout,yout,youtp,nout,taus)
    USE prec_rkind
    implicit none
    integer :: nin,nout
    REAL(RKIND) :: xin(nin), yin(nin), xout
    REAL(RKIND) :: yout, youtp
    REAL(RKIND), optional ::  taus
    !
    REAL(RKIND) ::  ynew(NIN), yinpp(NIN), ybc(6),psig(nin)
    REAL(RKIND) ::  youtpp, youtint
    INTEGER ioptder, NBCLFT,NBCRGT, iextrapo, nbc(2)
    REAL(RKIND) :: XBCLFT,XBCRGT, YBCLFT,YBCRGT
    REAL(RKIND) :: PXEXP0,PXEXPDL,zz
    !
    integer iflag
    !
    ! default values
    !
    nbc=0
    ybc=0._rkind
    ioptder=1
    iextrapo=32
    !
!!$    if (present(youtp)) then
!!$      print *,' youtp is present'
!!$    else
!!$      print *,' youtp is not present'
!!$    end if
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
    PXEXP0 = ybc(5)
    PXEXPDL= ybc(6)
    CALL CBSPLGEN(XIN,YIN,YNEW,YINPP,Nin,XOUT,YOUT,YOUTP,YOUTPP,YOUTINT, &
         &    nout,ioptder,PSIG,NBCLFT, &
         &    NBCRGT,XBCLFT,XBCRGT,YBCLFT,YBCRGT,iextrapo,PXEXP0,PXEXPDL &
         &    ,IFLAG)
    !

  end SUBROUTINE interpos_defp_xoutscal

  SUBROUTINE interpos_def_xoutscal(XIN,YIN,nin,xout,yout,nout,taus)
    USE prec_rkind
    implicit none
    integer :: nin,nout
    REAL(RKIND) :: xin(nin), yin(nin), xout
    REAL(RKIND) :: yout
    REAL(RKIND), optional ::  taus
    !
    REAL(RKIND) ::  ynew(NIN), yinpp(NIN), ybc(6),psig(nin)
    REAL(RKIND) ::  youtp, youtpp, youtint
    INTEGER ioptder, NBCLFT,NBCRGT, iextrapo, nbc(2)
    REAL(RKIND) :: XBCLFT,XBCRGT, YBCLFT,YBCRGT
    REAL(RKIND) :: PXEXP0,PXEXPDL,zz
    !
    integer iflag
    !
    ! default values
    !
    nbc=0
    ybc=0._rkind
    ioptder=0
    iextrapo=32
    !
!!$    if (present(youtp)) then
!!$      print *,' youtp is present'
!!$    else
!!$      print *,' youtp is not present'
!!$    end if
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
    PXEXP0 = ybc(5)
    PXEXPDL= ybc(6)
    CALL CBSPLGEN(XIN,YIN,YNEW,YINPP,Nin,XOUT,YOUT,YOUTP,YOUTPP,YOUTINT, &
         &    nout,ioptder,PSIG,NBCLFT, &
         &    NBCRGT,XBCLFT,XBCRGT,YBCLFT,YBCRGT,iextrapo,PXEXP0,PXEXPDL &
         &    ,IFLAG)
    !

  end SUBROUTINE interpos_def_xoutscal

  SUBROUTINE interpos_def(XIN,YIN,nin,xout,yout,nout,taus)
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
    !
    integer iflag
    !
    ! default values
    !
    nbc=0
    ybc=0._rkind
    ioptder=0
    iextrapo=32
    !
!!$    if (present(youtp)) then
!!$      print *,' youtp is present'
!!$    else
!!$      print *,' youtp is not present'
!!$    end if
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
    PXEXP0 = ybc(5)
    PXEXPDL= ybc(6)
    CALL CBSPLGEN(XIN,YIN,YNEW,YINPP,Nin,XOUT,YOUT,YOUTP,YOUTPP,YOUTINT, &
         &    nout,ioptder,PSIG,NBCLFT, &
         &    NBCRGT,XBCLFT,XBCRGT,YBCLFT,YBCRGT,iextrapo,PXEXP0,PXEXPDL &
         &    ,IFLAG)
    !

  end SUBROUTINE interpos_def

  SUBROUTINE interpos_defp(XIN,YIN,nin,xout,yout,youtp,nout,taus)
    USE prec_rkind
    implicit none
    integer :: nin,nout
    REAL(RKIND) :: xin(nin), yin(nin), xout(nout)
    REAL(RKIND) :: yout(nout), youtp(nout)
    REAL(RKIND), optional ::  taus
    !
    REAL(RKIND) ::  ynew(NIN), yinpp(NIN), ybc(6),psig(nin)
    REAL(RKIND) ::  youtpp(nout), youtint(nout)
    INTEGER ioptder, NBCLFT,NBCRGT, iextrapo, nbc(2)
    REAL(RKIND) :: XBCLFT,XBCRGT, YBCLFT,YBCRGT
    REAL(RKIND) :: PXEXP0,PXEXPDL,zz
    !
    integer iflag
    !
    ! default values
    !
    nbc=0
    ybc=0._rkind
    ioptder=1
    iextrapo=32
    !
!!$    if (present(youtp)) then
!!$      print *,' youtp is present'
!!$    else
!!$      print *,' youtp is not present'
!!$    end if
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
    PXEXP0 = ybc(5)
    PXEXPDL= ybc(6)
    CALL CBSPLGEN(XIN,YIN,YNEW,YINPP,Nin,XOUT,YOUT,YOUTP,YOUTPP,YOUTINT, &
         &    nout,ioptder,PSIG,NBCLFT, &
         &    NBCRGT,XBCLFT,XBCRGT,YBCLFT,YBCRGT,iextrapo,PXEXP0,PXEXPDL &
         &    ,IFLAG)
    !

  end SUBROUTINE interpos_defp

  SUBROUTINE interpos_defint(XIN,YIN,nin,xout,yout,nout,taus,youtint)
    USE prec_rkind
    implicit none
    integer :: nin,nout
    REAL(RKIND) :: xin(nin), yin(nin), xout(nout)
    REAL(RKIND) :: yout(nout)
    REAL(RKIND),optional :: youtint(nout)
    REAL(RKIND), optional ::  taus
    !
    REAL(RKIND) :: youtp2(nout), youtpp2(nout), youtint2(nout)
    REAL(RKIND) ::  ynew(NIN), yinpp(NIN), ybc(6),psig(nin)
    REAL(RKIND) ::  youtpp(nout)
    INTEGER ioptder, NBCLFT,NBCRGT, iextrapo, nbc(2)
    REAL(RKIND) :: XBCLFT,XBCRGT, YBCLFT,YBCRGT
    REAL(RKIND) :: PXEXP0,PXEXPDL,zz
    !
    integer iflag
    !
    ! default values
    !
    nbc=0
    ybc=0._rkind
    ioptder=3
    iextrapo=32
    !
!!$    if (present(youtp)) then
!!$      print *,' youtp is present'
!!$    else
!!$      print *,' youtp is not present'
!!$    end if
    if ((present(taus)) .and. (taus .ge. 0._rkind)) then
      print *,'% in defint: taus= ',taus
      psig=taus
    else
      zz=minval(xin(2:nin)-xin(1:nin-1))
      psig=zz**3
      print *,'% in defint: taus defined as = ',psig(1)
    end if
    NBCLFT = nbc(1)
    NBCRGT = nbc(2)
    YBCLFT = ybc(1)
    YBCRGT = ybc(2)
    if (NBCLFT .ge. 10) XBCLFT = ybc(3)
    if (NBCRGT .ge. 10) XBCRGT = ybc(4)
    !
    PXEXP0 = ybc(5)
    PXEXPDL= ybc(6)
    CALL CBSPLGEN(XIN,YIN,YNEW,YINPP,Nin,XOUT,YOUT,YOUTP2,YOUTPP2,YOUTINT2, &
         &    nout,ioptder,PSIG,NBCLFT, &
         &    NBCRGT,XBCLFT,XBCRGT,YBCLFT,YBCRGT,iextrapo,PXEXP0,PXEXPDL &
         &    ,IFLAG)
    !
    if (present(youtint)) then
       youtint=youtint2
       print *,' youtint is present'
    else
       print *,' youtint is not present'
    end if

  end SUBROUTINE interpos_defint

  SUBROUTINE interpos_sfull(XIN,YIN,YNEW,YINPP,nin,xout,yout,youtp,youtpp,youtint, &
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
  !
  integer iflag
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
  PXEXP0 = ybc(5)
  PXEXPDL= ybc(6)
  CALL CBSPLGEN(XIN,YIN,YNEW,YINPP,Nin,XOUT,YOUT,YOUTP,YOUTPP,YOUTINT, &
    &    nout,ioptder,PSIG,NBCLFT, &
    &    NBCRGT,XBCLFT,XBCRGT,YBCLFT,YBCRGT,iextrapo,PXEXP0,PXEXPDL &
    &    ,IFLAG)
  !
END SUBROUTINE interpos_sfull

  SUBROUTINE interpos_vfull(XIN,YIN,YNEW,YINPP,nin,xout,yout,youtp,youtpp,youtint, &
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
  !
  integer iflag
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
  PXEXP0 = ybc(5)
  PXEXPDL= ybc(6)
  CALL CBSPLGEN(XIN,YIN,YNEW,YINPP,Nin,XOUT,YOUT,YOUTP,YOUTPP,YOUTINT, &
    &    nout,ioptder,PSIG,NBCLFT, &
    &    NBCRGT,XBCLFT,XBCRGT,YBCLFT,YBCRGT,iextrapo,PXEXP0,PXEXPDL &
    &    ,IFLAG)
  !
END SUBROUTINE interpos_vfull

  SUBROUTINE interpos_ssimple(XIN,YIN,nin,xout,yout,youtp,youtpp,youtint, &
  &  nout,ioptder,iextrapo,taus,nbc,ybc)
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
  !
  integer iflag
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
  PXEXP0 = ybc(5)
  PXEXPDL= ybc(6)
  CALL CBSPLGEN(XIN,YIN,YNEW,YINPP,Nin,XOUT,YOUT,YOUTP,YOUTPP,YOUTINT, &
    &    nout,ioptder,PSIG,NBCLFT, &
    &    NBCRGT,XBCLFT,XBCRGT,YBCLFT,YBCRGT,iextrapo,PXEXP0,PXEXPDL &
    &    ,IFLAG)
  !
END SUBROUTINE interpos_ssimple

  SUBROUTINE interpos_vsimple(XIN,YIN,nin,xout,yout,youtp,youtpp,youtint, &
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
  PXEXP0 = ybc(5)
  PXEXPDL= ybc(6)
  CALL CBSPLGEN(XIN,YIN,YNEW,YINPP,Nin,XOUT,YOUT,YOUTP,YOUTPP,YOUTINT, &
    &    nout,ioptder,PSIG,NBCLFT, &
    &    NBCRGT,XBCLFT,XBCRGT,YBCLFT,YBCRGT,iextrapo,PXEXP0,PXEXPDL &
    &    ,IFLAG)
  !
END SUBROUTINE interpos_vsimple

end MODULE interpos_module
