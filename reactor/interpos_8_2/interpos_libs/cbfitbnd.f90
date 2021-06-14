!-----------------------------------------------------------------------
SUBROUTINE CBFITBND(PXIN,PYIN,PYINNEW,KNIN,PYINPP,PSIG, &
  &  NBCLFT,NBCRGT,XBCLFT,XBCRGT,YBCLFT, &
  &  YBCRGT,PXEXP0,PXEXPDL)
  !
  !     NON-PERIODIC B.C
  !
  !     PREPARE SECOND DERIVATIVE OF CUBIC SPLINE INTERPOLATION AND NEW
  !     VALUES OF Y AT NODES YNEW FITTED SUCH THAT CHI**2 + TAUS*F''**2
  !     IS MINIMIZED ACCORDING TO HIRSHMAN ET AL, PHYS. PLASMAS 1 (1994) 2280.
  !     SIG = TAU*SIGMA_K/min(SIGMA_K) OF PAPER.
  !
  !     SETTING TAUS=0., ONE FINDS THE USUAL CUBIC SPLINE INT. WITH CHI**2=0
  !     TAUS LARGE => FIT CLOSER TO STRAIGHT LINE (SECOND DERIV.=0)
  !
  !     IF LAPACK ROUTINES NOT AVAILABLE, USE spgbtrf_s.f
  !     (LAPACK SOURCE COPIED FROM NETLIB.ORG)
  !
  !     IF TAUS=0, PYINNEW NOT USED => PYINNEW(1) OR PYINNEW=PYIN IS OK
  !
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
  !     IF END POINTS ARE USED FOR THE B.C. AND TYPE 0 OR 1, THEN USE SYMMETRY
  !     OF MATRIX
  !
  !     IF XBCLFT OR XBCRGT ARE USED, IMPOSE B.C. ON NODE CLOSEST TO XBCLFT OR XBCRGT
  !
  !     TENSION TAUS(K) IS GIVEN WITH AN EXPONENTIAL FORM TO BE ABLE TO LOCALIZE
  !     IT:
  !     .     TAU_K = PTAUS * EXP( -COF * ((X-X0)/DX)**2)
  !
  !     WHERE X0 = PXEXP0 AND DX = PXEXPDL, AND:
  !     COF = 1. IF PXEXP0 IN [PXIN(1),PXIN(KNIN)], 0 OTHERWISE
  !     THUS SETTING PXEXP0 OUTSIDE DOMAIN GIVES A CST TAU_K VALUE
  !
  !     TOTAL DIMENSION OF PAMAT. THE REQUIRED SPACE DEPENDS IF
  !     .         PTAUS IS ZERO OR NOT AND ON THE B.C. (SYMMETRIC OR NOT)
  !     THUS, MDMATOT CAN VARY BETWEEN 2*KNIN AND 10*KNIN (MAX. VALUE NEEDED)
  !
  !-----------------------------------------------------------------------
  !
  USE prec_rkind
  implicit none
  ! arguments
  integer :: KNIN, NBCLFT, NBCRGT
  real(rkind)  :: PXIN(KNIN), PYIN(KNIN), PSIG(KNIN), &
    & xbclft, xbcrgt, ybclft, ybcrgt, pxexp0, pxexpdl
  real(rkind) :: PYINNEW(KNIN), PYINPP(KNIN)
  !
  integer :: idamat, kpm2, ixbc, ibctyp, iik, &
    &  itauval, isym, n, i, j, k, iup, idown, idiag, ishift, ieff, ikp2 &
    &  ,ikp1, ikm1, ikm2, jk, jkp1, jkp2, jeff, iii, iupsofar, idwnsofa &
    &  , jbc, ik, idiamik, idiapik, iklft, irhs, info, &
    &  info2, jkm1, jkm2
  DIMENSION :: KPM2(KNIN,-2:+2), IXBC(2), IBCTYP(2)
  !OS      pointer(iptr_ftauk,ftauk)
  REAL(RKIND), DIMENSION(KNIN) :: ftauk
  !
  REAL(RKIND) :: WHK(KNIN), WOHK(KNIN), ZYBC(2)
  REAL(RKIND) :: xtkm1, xohkm1, xohkm2, xtk, xhkm1, &
    &  xohk, xtkp1, xhk, xohkp1, xykp1, xyk, &
    &  xykm1, ztaueff, zcofexp, zxexp0, zxexpdl, a1, a2, a3, a4, b1, &
    &  b2, b3, b4, px, fakk, fakkp1, fakkp2, frhs, zdelx, zero, zvalue, &
    &  zypeff, ztohkk1, zsign, fa2, fa3, fa0, fa1, &
    &  fakkm1, fakkm2, fun_ftauk, fcccc0, fcccc1, fcccc2, fcccc3
  integer nel
  REAL(RKIND), ALLOCATABLE :: PAMAT(:,:)
  !
  !
  !     FUNCTIONS FOR MATRIX COEFFICIENTS
  !
  REAL(RKIND) :: zsix, zthree, ztwo, zone
  PARAMETER(zsix=6._RKIND, zthree=3._RKIND, ztwo=2._RKIND, zone=1._RKIND)
  FAKKM2(XTKM1,XOHKM1,XOHKM2) = XTKM1*XOHKM1*XOHKM2
  FAKKM1(XTK,XTKM1,XHKM1,XOHK,XOHKM1,XOHKM2) = XHKM1/zsix &
    &  - XOHKM1*(XTK*XOHK+(XTK+XTKM1)*XOHKM1 + XTKM1*XOHKM2)
  FAKK(XTKP1,XTK,XTKM1,XHK,XHKM1,XOHK,XOHKM1) = (XHK+XHKM1)/ZTHREE &
    &  + XOHK*XOHK*(XTKP1+XTK) &
    &  + XOHKM1*(ZTWO*XTK*XOHK+(XTK+XTKM1)*XOHKM1)
  FAKKP1(XTKP1,XTK,XHK,XOHKP1,XOHK,XOHKM1) = XHK/zsix &
    &  - XOHK*(XTKP1*XOHKP1+(XTK+XTKP1)*XOHK + XTK*XOHKM1)
  FAKKP2(XTKP1,XOHKP1,XOHK) = XTKP1*XOHKP1*XOHK
  !
  FRHS(XYKP1,XYK,XYKM1,XOHK,XOHKM1) = (XYKP1-XYK)*XOHK &
    &  - (XYK-XYKM1)*XOHKM1
  !
  !     WHEN ONE WANTS AN ARRAY FOR TAU*SIGMA_K**2, THEN ONE SHOULD REPLACE
  !     THE FUNCTION FTAU BY AN ARRAY
  !
  fun_FTAUK(IIK)= ZTAUEFF*EXP(-ZCOFEXP*(PXIN(IIK)-ZXEXP0)**2 &
    &  /ZXEXPDL**2)
  !%OS      FTAUK(IIK) = ZTAUEFF
  !
  !.......................................................................
  !*COMDECK CUCCCC
  ! ----------------------------------------------------------------------
  ! --     STATEMENT FUNCTION FOR CUBIC INTERPOLATION                   --
  ! --                         23.04.88            AR        CRPP       --
  ! --                                                                  --
  ! -- CUBIC INTERPOLATION OF A FUNCTION F(X)                           --
  ! -- THE EIGHT ARGUMENTS A1,A2,A3,A4,B1,B2,B3,B4 ARE DEFINED BY:      --
  ! -- F(B1) = A1 , F(B2) = A2 , F(B3) = A3 , F(B4) = A4                --
  ! ----------------------------------------------------------------------
  !
  FA3(A1,A2,A3,A4,B1,B2,B3,B4) = &
    &        (A1-A2) / ((B1-B2)*(B2-B4)*(B2-B3)) + &
    &        (A1-A3) / ((B4-B3)*(B3-B1)*(B3-B2)) + &
    &        (A1-A4) / ((B1-B4)*(B2-B4)*(B3-B4))
  FA2(A1,A2,A3,A4,B1,B2,B3,B4) = &
    &        (A1-A2) / ((B2-B1)*(B3-B2)) + &
    &        (A3-A1) / ((B3-B1)*(B3-B2)) - &
    &        (B1+B2+B3) * FA3(A1,A2,A3,A4,B1,B2,B3,B4)
  FA1(A1,A2,A3,A4,B1,B2,B3,B4) = &
    &        (A1-A2) / (B1-B2) - &
    &        (B1+B2) * FA2(A1,A2,A3,A4,B1,B2,B3,B4) - &
    &        (B1*B1+B1*B2+B2*B2) * FA3(A1,A2,A3,A4,B1,B2,B3,B4)
  FA0(A1,A2,A3,A4,B1,B2,B3,B4) = &
    &        A1 - &
    &        B1 * (FA1(A1,A2,A3,A4,B1,B2,B3,B4) + &
    &              B1 * (FA2(A1,A2,A3,A4,B1,B2,B3,B4) + &
    &                    B1 * FA3(A1,A2,A3,A4,B1,B2,B3,B4)))
  ! ----------------------------------------------------------------------
  ! -- FCCCC0 GIVES THE VALUE OF THE FUNCTION AT POINT PX:              --
  ! -- FCCCC0(......,PX) = F(PX)                                        --
  ! ----------------------------------------------------------------------
  FCCCC0(A1,A2,A3,A4,B1,B2,B3,B4,PX) = &
    &              FA0(A1,A2,A3,A4,B1,B2,B3,B4) + &
    &              PX * (FA1(A1,A2,A3,A4,B1,B2,B3,B4) + &
    &                    PX * (FA2(A1,A2,A3,A4,B1,B2,B3,B4) + &
    &                          PX * FA3(A1,A2,A3,A4,B1,B2,B3,B4)))
  ! ----------------------------------------------------------------------
  ! -- FCCCC1 GIVES THE VALUE OF THE FIRST DERIVATIVE OF F(X) AT PX:    --
  ! -- FCCCC1(......,PX) = DF/DX (PX)                                   --
  ! ----------------------------------------------------------------------
  FCCCC1(A1,A2,A3,A4,B1,B2,B3,B4,PX) = &
    &              FA1(A1,A2,A3,A4,B1,B2,B3,B4) + &
    &              PX * (ZTWO * FA2(A1,A2,A3,A4,B1,B2,B3,B4) + &
    &                    ZTHREE * PX * FA3(A1,A2,A3,A4,B1,B2,B3,B4))
  ! ----------------------------------------------------------------------
  ! -- FCCCC2 GIVES THE VALUE OF THE SECOND DERIVATIVE OF F(X) AT PX:   --
  ! -- FCCCC2(......,PX) = D2F/DX2 (PX)                                 --
  ! ----------------------------------------------------------------------
  FCCCC2(A1,A2,A3,A4,B1,B2,B3,B4,PX) = &
    &             ZTWO * FA2(A1,A2,A3,A4,B1,B2,B3,B4) + &
    &             zsix * FA3(A1,A2,A3,A4,B1,B2,B3,B4) * PX
  ! ----------------------------------------------------------------------
  ! -- FCCCC3 GIVES THE VALUE OF THE THIRD DERIVATIVE OF F(X) AT PX:     -
  ! -- FCCCC3(......,PX) = D3F/DX3 (PX)                                  -
  ! ----------------------------------------------------------------------
  FCCCC3(A1,A2,A3,A4,B1,B2,B3,B4,PX) = &
    &                      zsix * FA3(A1,A2,A3,A4,B1,B2,B3,B4)
  !.......................................................................
  !
  !-----------------------------------------------------------------------
  !
  !     0. INITIALIZATION
  !
  ITAUVAL = 1
  IF (PSIG(1) .EQ. 0._RKIND) ITAUVAL = 0
  ZTAUEFF = abs(PSIG(1))
  ZXEXP0 = PXEXP0
  ZXEXPDL = PXEXPDL
  ZCOFEXP = 1.0_RKIND
  IF (ZXEXP0.LT.PXIN(1) .OR. ZXEXP0.GT.PXIN(KNIN)) ZCOFEXP=0.0_RKIND
  !
  ISYM = 1
  IF (NBCLFT.GE.10 .OR. NBCRGT.GE.10 .OR. NBCLFT.EQ.2 &
    &  .OR. NBCRGT.EQ.2) ISYM = 0
  !
  !
  !     0.3 PREPARE BAND WIDTH
  !
  !     PREPARE MATRIX DIMENSION.
  !     MINIMUM REQUIRED:
  !     IUP = 1 = IDOWN; IF PTAUS.NE.0 => IUP = IDOWN = 2
  !     IF SYMMETRIC, USE ONLY UPPER BAND AND DIAGONAL =>IDAMAT=IUP+1
  !     IF ASYMMETRIC => IDAMAT = 2*IDOWN + IUP + 1
  !     IF B.C. NOT AT END OF INTERVAL => IUP = IUP + 1, AND/OR IDOWN=IDOWN+1
  !
  !     => ALTOGETHER, MINIMUM VALUE: IDOWN=1, IUP=1, SYM. =>IDAMAT_MAX = 2
  !     => ALTOGETHER, MAXIMUM VALUE: IDOWN=3, IUP=3 =>IDAMAT_MAX = 10
  !
  !
  IF (ITAUVAL .EQ. 0) THEN
    IUP   = 1
    IDOWN = 1
    IDAMAT = 2
  ELSE
    IUP   = 2
    IDOWN = 2
    IDAMAT = 3
  END IF
  IDIAG = IUP + 1
  IF (ISYM .EQ. 0) THEN
    IDAMAT = 3*(IDAMAT-1)+1
    IF (NBCLFT .GE. 10) THEN
      IUP = IUP + 1
      IDAMAT = IDAMAT + 1
    END IF
    IF (NBCRGT .GE. 10) THEN
      IDOWN = IDOWN + 1
      IDAMAT = IDAMAT + 2
    END IF
    IDIAG = IUP + IDOWN + 1
  ENDIF
  !
  ALLOCATE(PAMAT(IDAMAT,KNIN))
  !
  PYINPP = 0._RKIND
  PAMAT = 0._RKIND
  !
  !     0.2 PRE-COMPUTE H_K AND 1./H_K, and zftauk
  !
  DO K=1,KNIN-1
    WHK(K)  = (PXIN(K+1) - PXIN(K))
  enddo
  DO K=1,KNIN-1
    WHK(K)  = (PXIN(K+1) - PXIN(K))
    WOHK(K) = zone / WHK(K)
    ftauk(k)=psig(k)
    !%OS        ftauk(k)=fun_ftauk(k)
  END DO
  WHK(KNIN) = 0.0_RKIND
  WOHK(KNIN) = 0.0_RKIND
  ftauk(knin)=psig(KNIN)
  !%OS      ftauk(knin)=fun_ftauk(KNIN)
  if (PSIG(1).lt.0._RKIND) ftauk(1)=-10._RKIND*ftauk(1)
  !
  !     0.4 DETERMINE NEIGHBOURS: K-2, K-1, .., K+2
  !     WHEN OUT OF BOUNDS, POINT TO INDEX KNIN, AS WHK, WOHK(KNIN)=0.0
  !
  DO ISHIFT=-2,+2
    DO K=1,KNIN
      KPM2(K,ISHIFT) = K + ISHIFT
    END DO
  END DO
  !     OUT OF INTERVAL: SEND TO KNIN
  KPM2(1,-2)   = KNIN
  KPM2(1,-1)   = KNIN
  KPM2(2,-2)   = KNIN
  KPM2(KNIN-1,+2) = KNIN
  KPM2(KNIN  ,+1) = KNIN
  KPM2(KNIN  ,+2) = KNIN
  !
  !     1. CONSTRUCT MATRIX AND R.H.S
  !     LAPACK SET-UP OF MATRIX:    A(I,J) -> PAMAT(I-J+IDIAG, J)
  !
  !.......................................................................
  !     AS MATRIX SYMMETRIC, COMPUTE ONLY UPPER PART, THAT IS IF J.GE.I
  !
  DO K=1,KNIN-1
    IEFF = K + IDIAG
    IKP2 = KPM2(K,+2)
    IKP1 = KPM2(K,+1)
    IKM1 = KPM2(K,-1)
    IKM2 = KPM2(K,-2)
    !     A(K,K)
    JK = K
    PAMAT(IEFF-JK,JK) = FAKK(FTAUK(IKP1),FTAUK(K),FTAUK(IKM1),WHK(K) &
      &    ,WHK(IKM1),WOHK(K),WOHK(IKM1))
    !
    !     A(K,K+1)
    JKP1 = K + 1
    PAMAT(IEFF-JKP1,JKP1) = FAKKP1(FTAUK(IKP1),FTAUK(K),WHK(K), &
      &    WOHK(IKP1),WOHK(K),WOHK(IKM1))
    !     A(K,K-1)
    !%OS        JKM1 = K - 1
    !%OS        PAMAT(IEFF-JKM1,JKM1) = FAKKM1(FTAUK(K),FTAUK(IKM1),WHK(IKM1), &
    !%OS     &    WOHK(K),WOHK(IKM1),WOHK(IKM2))
    !
    IF (ITAUVAL .EQ. 1) THEN
      !     A(K,K+2)
      JKP2 = K + 2
      IF (JKP2 .LE. KNIN) &
        &      PAMAT(IEFF-JKP2,JKP2) = FAKKP2(FTAUK(IKP1),WOHK(IKP1), &
        &      WOHK(K))
      !     A(K,K-2)
      !%OS          JKM2 = K - 2
      !%OS          PAMAT(IEFF-JKM2,JKM2) = FAKKM2(FTAUK(IKM1),WOHK(IKM1), &
      !%OS     &      WOHK(IKM2))
    ENDIF
    !     B(K)
    PYINPP(K) = FRHS(PYIN(IKP1),PYIN(K),PYIN(IKM1),WOHK(K), &
      &    WOHK(IKM1))
    !
  END DO
  !
  !     2. BOUNDARY CONDITIONS
  !
  !     2.1 IF NON-SYMMETRIC, COPY TOP PART TO BOTTOM BEFORE APPLYING
  !     B.C.
  !
  IF (ISYM .EQ. 0) THEN
    DO I=1,KNIN
      IEFF = I + IDIAG
      DO J=I+1,MIN(I+MIN(IUP,IDOWN),KNIN)
        JEFF = J + IDIAG
        !     A(J,I) = A(I,J)
        PAMAT(JEFF-I,I) = PAMAT(IEFF-J,J)
      END DO
    END DO
  ENDIF
  !%OS
  !     debug, print matrix and rhs
  !%OS      write(6,'(3a4,a)') 'i ','j1 ','j2 ',' i,j1  i,j1+1,..... i,j2'
  !%OS      do i=1,knin
  !%OS        ieff = idiag + i
  !%OS        j1 = i-idown
  !%OS        j2 = i+iup
  !%OSc%OS        j1 = max(i-idown,1)
  !%OSc%OS        j2 = min(i+iup,knin)
  !%OS        write(6,'(3i4,1p10e13.4)') i,j1,j2,(pamat(ieff-j,j),j=j1,j2)
  !%OS      end do
  !%OS      write(6,'(a4,a12)') 'i','RHS'
  !%OS      write(6,'(i4,1pe13.4)') (i,pyinpp(i),i=1,knin)
  !
  !%OS
  !
  !     2.2 B.C. AT TWO LOCATIONS PXIN(IXBC(JBC)), JBC=1,2
  !     IBCTYP(JBC) = 0, 1 OR 2 (TYPE OF B.C, SEE ABOVE).
  !     SO FAR USES NODE CLOSEST TO XBCLFT/RGT FOR LOCATION
  !     OF B.C., INSTEAD OF ACTUAL VALUE OF XBCLFT/RGT
  !
  IXBC(1) = 1
  IXBC(2) = KNIN
  IF (NBCLFT .GE. 10) THEN
    DO I=1,KNIN
      IF (PXIN(I) .GE. XBCLFT) GO TO 220
    END DO
220 CONTINUE
    ZDELX = ABS(PXIN(I)-XBCLFT)
    IXBC(1) = I
    IF (I .GE. KNIN) THEN
      IXBC(1) = KNIN
      PRINT *,' WARNING: LEFT B.C. AT I=KNIN: XBCLFT=',XBCLFT, &
        &      '  PXIN(KNIN)= ',PXIN(KNIN)
    ELSE IF (ABS(PXIN(I-1)-XBCLFT).LE.ZDELX .AND. I.NE.1) THEN
      IXBC(1) = I-1
    ENDIF
  ENDIF
  !
  IF (NBCRGT .GE. 10) THEN
    DO I=1,KNIN
      IF (PXIN(I) .GE. XBCRGT) GO TO 221
    END DO
221 CONTINUE
    ZDELX = ABS(PXIN(I)-XBCRGT)
    IXBC(2) = I
    IF (I .LE. 1) THEN
      IXBC(2) = 1
      PRINT *,' WARNING: RIGHT B.C. AT I=1: XBCRGT=',XBCRGT, &
        &      '  PXIN(1)= ',PXIN(1)
    ELSE IF (I .GT. KNIN) THEN
      IXBC(2) = KNIN
    ELSE IF (ABS(PXIN(I-1)-XBCRGT) .LE. ZDELX) THEN
      IXBC(2) = I-1
    ENDIF
  ENDIF
  !
  ZYBC(1) = YBCLFT
  ZYBC(2) = YBCRGT
  nel=10
  IBCTYP(1) = MOD(NBCLFT,nel)
  IBCTYP(2) = MOD(NBCRGT,nel)
  IF (IXBC(1) .EQ. IXBC(2)) THEN
    PRINT *,' ERROR, B.C. AT SAME LOCATIONS: IXBC(1)=IXBC(2)= ', &
      &    IXBC(1)
    !%OS        STOP '1=2'
    RETURN
  ELSE IF (IXBC(1) .GT. IXBC(2)) THEN
    PRINT *,' WARNING, NEEDED TO SWITCH B.C. POINTS AS IXBC(1)= ', &
      &    IXBC(1),' > IXBC(2)= ',IXBC(2)
    III = IXBC(1)
    IXBC(1) = IXBC(2)
    IXBC(2) = III
    ZYBC(1) = YBCRGT
    ZYBC(2) = YBCLFT
    nel=10
    IBCTYP(1) = MOD(NBCRGT,nel)
    IBCTYP(2) = MOD(NBCLFT,nel)
  ENDIF
  !
  !     2.3 MOVE EQUATIONS UP OR DOWN IF B.C. IS NOT AN END POINT
  !
  IF (IXBC(1) .NE. 1) THEN
    !
    !     MOVE ROW EQ. K=2,..,IXBC(1) UP BY ONE
    IUPSOFAR = IUP - 1
    DO K=2,IXBC(1)
      DO J=MAX(1,K-IDOWN),MIN(KNIN,K+IUPSOFAR)
        PAMAT(IDIAG+(K-1)-J,J) = PAMAT(IDIAG+K-J,J)
      END DO
      PYINPP(K-1) = PYINPP(K)
      !     ZERO A((K-1),(K-1)-IDOWN)
      IF (K-1-IDOWN .GE. 1) PAMAT(IDIAG+IDOWN,K-1-IDOWN) = 0.0_RKIND
    END DO
    !     ZERO ROW IXBC(1) AND RHS
    K = IXBC(1)
    DO J=MAX(1,K-IDOWN),MIN(KNIN,K+IUP)
      PAMAT(IDIAG+K-J,J) = 0.0_RKIND
    END DO
    PYINPP(K) = 0.0_RKIND
  ENDIF
  !
  IF (IXBC(2) .NE. KNIN) THEN
    !     
    !     MOVE EQ. K=IXBC(2),..,KNIN-1 DOWN BY ONE
    IDWNSOFA = IDOWN - 1
    DO K=KNIN-1,IXBC(2),-1
      DO J=MAX(1,K-IDWNSOFA),MIN(KNIN,K+IUP)
        PAMAT(IDIAG+(K+1)-J,J) = PAMAT(IDIAG+K-J,J)
      END DO
      PYINPP(K+1) = PYINPP(K)
      !     ZERO A((K+1),(K+1)+IUP)
      IF (K+1+IUP .LE. KNIN) PAMAT(IDIAG-IUP,K+1+IUP) = 0.0_RKIND
    END DO
    !     ZERO ROW IXBC(2) AND RHS
    K = IXBC(2)
    DO J=MAX(1,K-IDOWN),MIN(KNIN,K+IUP)
      PAMAT(IDIAG+K-J,J) = 0.0_RKIND
    END DO
    PYINPP(K) = 0.0_RKIND
  ENDIF
  !
  !     2.4 FOR ROW=IXBC(), MODIFY MATRIX AND RHS ACCORDING TO B.C. TYPE
  !
  ZERO = 0.0_RKIND
  DO JBC=1,2
    IK = IXBC(JBC)
    ZVALUE = ZYBC(JBC)
    IEFF = IK + IDIAG
    IKP2 = KPM2(IK,+2)
    IKP1 = KPM2(IK,+1)
    IKM1 = KPM2(IK,-1)
    IKM2 = KPM2(IK,-2)
    IF (IBCTYP(JBC) .EQ. 0) THEN
      !
      !     SYMMETRIZE => COL IK GOES TO RIGHT-HAND SIDE AND THEN ZEROED
      !
      IF (ISYM .EQ. 1) THEN
        IDIAMIK = IDIAG - IK
        DO I=MAX(1,IK-IUP),IK-1
          PYINPP(I) = PYINPP(I) - ZVALUE * PAMAT(I+IDIAMIK,IK)
          PAMAT(I+IDIAMIK,IK) = 0.0_RKIND
        END DO
        IDIAPIK = IDIAG + IK
        DO I=IK+1,MIN(KNIN,IK+IUP)
          PYINPP(I) = PYINPP(I) - ZVALUE * PAMAT(IDIAPIK-I,I)
          PAMAT(IDIAPIK-I,I) = 0.0_RKIND
        END DO
      ELSE
        IDIAMIK = IDIAG - IK
        DO I=MAX(1,IK-IUP),MIN(KNIN,IK+IDOWN)
          PYINPP(I) = PYINPP(I) - ZVALUE * PAMAT(I+IDIAMIK,IK)
          PAMAT(I+IDIAMIK,IK) = 0.0_RKIND
        END DO
        !     ZERO ROW IK
        DO J=MAX(1,IK-IDOWN),MIN(KNIN,IK+IUP)
          PAMAT(IEFF-J,J) = 0.0_RKIND
        END DO
      ENDIF
      !
      !     REPLACE ROW IK BY EQUATION: G_K = ZVALUE
      PAMAT(IDIAG,IK) = 1.0_RKIND
      PYINPP(IK) = ZVALUE
      !
    ELSE IF (IBCTYP(JBC) .EQ. 1) THEN
      !
      !     1ST DERIVATIVE GIVEN
      !
      ZYPEFF = ZVALUE
      IF (ZVALUE .GT. 1.E+31_RKIND) THEN
        !     FROM LGRANGIAN INTERPOLATION
        IKLFT = IK - 1
        IF (IK .EQ. 1) IKLFT = IK
        IF (IKLFT+3 .GT. KNIN) IKLFT = KNIN - 3
        ZYPEFF = FCCCC1(PYIN(IKLFT),PYIN(IKLFT+1),PYIN(IKLFT+2), &
          &        PYIN(IKLFT+3),PXIN(IKLFT),PXIN(IKLFT+1),PXIN(IKLFT+2), &
          &        PXIN(IKLFT+3),PXIN(IK))
      ELSE IF (ZVALUE .LT. -1.E+31_RKIND) THEN
        IKLFT = IK
        IF (IK .EQ. KNIN) IKLFT = IK - 1
        ZYPEFF = (PYIN(IKLFT+1)-PYIN(IKLFT)) &
          &        / (PXIN(IKLFT+1)-PXIN(IKLFT))
      ENDIF
      ZTOHKK1 = FTAUK(IK)*WOHK(IK)*WOHK(IKM1)
      !     A(IK,IK)
      IF (IK .NE. KNIN) PAMAT(IEFF-IK,IK) = FAKK(FTAUK(IKP1),FTAUK(IK), &
        &      ZERO,WHK(IK),ZERO,WOHK(IK),ZERO) + ZTOHKK1
      IF (IK .EQ. KNIN) PAMAT(IEFF-IK,IK) = FAKK(ZERO,FTAUK(IK), &
        &      FTAUK(IKM1),ZERO,WHK(IKM1),ZERO,WOHK(IKM1)) + ZTOHKK1
      !     A(IK,IK-1)
      JKM1 = IK - 1
      IF (ISYM.EQ.0 .AND. JKM1.GE.1) THEN
        IF (IK .NE. KNIN) PAMAT(IEFF-JKM1,JKM1) = - ZTOHKK1
        IF (IK .EQ. KNIN) PAMAT(IEFF-JKM1,JKM1) = FAKKM1(FTAUK(IK), &
          &        FTAUK(IKM1),WHK(IKM1),ZERO,WOHK(IKM1),WOHK(IKM2))
      ENDIF
      !     A(IK,IK+1)
      JKP1 = IK + 1
      IF (JKP1 .LE. KNIN) &
        &      PAMAT(IEFF-JKP1,JKP1) = FAKKP1(FTAUK(IKP1), &
        &      FTAUK(IK),WHK(IK),WOHK(IKP1),WOHK(IK),ZERO)
      !
      IF (ITAUVAL .EQ. 1) THEN
        !     A(IK,IK+2)
        JKP2 = IK + 2
        IF (JKP2 .LE. KNIN) &
          &        PAMAT(IEFF-JKP2,JKP2) = FAKKP2(FTAUK(IKP1),WOHK(IKP1), &
          &        WOHK(IK))
        !     A(IK,IK-2)
        JKM2 = IK - 2
        IF (ISYM.EQ.0 .AND. JKM2.GE.1) THEN
          IF (IK .NE. KNIN) PAMAT(IEFF-JKM2,JKM2) = 0.0_RKIND
          IF (IK .EQ. KNIN) PAMAT(IEFF-JKM2,JKM2) = FAKKM2(FTAUK(IKM1), &
            &          WOHK(IKM1),WOHK(IKM2))
        ENDIF
      ENDIF
      !     RHS
      ZSIGN = -1._RKIND
      IF (IK .EQ. KNIN) ZSIGN = +1._RKIND
      IF (IK .NE. KNIN) PYINPP(IK) = FRHS(PYIN(IKP1),PYIN(IK),ZERO, &
        &      WOHK(IK),ZERO) - ZYPEFF
      IF (IK .EQ. KNIN) PYINPP(IK) = FRHS(ZERO,PYIN(IK),PYIN(IKM1), &
        &      ZERO,WOHK(IKM1)) + ZYPEFF
      !
    ELSE IF (IBCTYP(JBC) .EQ. 2) THEN
      !
      !     FUNCTION IS GIVEN
      !
      !     A(IK,IK)
      PAMAT(IEFF-IK,IK) = - FTAUK(IK) * (WOHK(IK) + WOHK(IKM1))
      !     A(IK,IK+1)
      JKP1 = IK + 1
      IF (JKP1 .LE. KNIN) &
        &      PAMAT(IEFF-JKP1,JKP1) = FTAUK(IK) * WOHK(IK)
      !     A(IK,IK-1)
      JKM1 = IK - 1
      IF (ISYM.EQ.0 .AND. JKM1.GE.1) &
        &      PAMAT(IEFF-JKM1,JKM1) = FTAUK(IK) * WOHK(IKM1)
      !
      IF (ITAUVAL .EQ. 1) THEN
        !     A(IK,IK+2)
        JKP2 = IK + 2
        IF (JKP2 .LE. KNIN) PAMAT(IEFF-JKP2,JKP2) = 0.0_RKIND
        !     A(IK,IK-2)
        JKM2 = IK - 2
        IF (ISYM.EQ.0 .AND. JKM2.GE.1) PAMAT(IEFF-JKM2,JKM2) = 0.0_RKIND
      ENDIF
      !     RHS
      PYINPP(IK) = PYIN(IK) - ZVALUE
      !
    ENDIF
    !
  END DO
  !
  !     3. SOLVE SYSTEM
  !
  !     USE INTEGER WORK SPACE FROM KPM2(0) ARRAY FOR IPIVOT,
  !     AS KPM2(K,0) NOT NEEDED NOR USED
  !
  !%OS
  !     debug, print matrix and rhs
  !%OS      write(6,'(3a4,a)') 'i ','j1 ','j2 ',' i,j1  i,j1+1,..... i,j2'
  !%OS      do i=1,knin
  !%OS        ieff = idiag + i
  !%OS        j1 = i-idown
  !%OS        j2 = i+iup
  !%OSc%OS        j1 = max(i-idown,1)
  !%OSc%OS        j2 = min(i+iup,knin)
  !%OS        write(6,'(3i4,1p10e13.4)') i,j1,j2,(pamat(ieff-j,j),j=j1,j2)
  !%OS      end do
  !%OS      write(6,'(a4,a12)') 'i','RHS'
  !%OS      write(6,'(i4,1pe13.4)') (i,pyinpp(i),i=1,knin)
  !
  !%OS
  !%OS      print *,'isym= ',isym
  IF (ISYM .EQ. 1) THEN
    !   Single precision (using routines in spgbtrf_s.f compiled with -r8 for example)
    !   Double precision (using the lapack libraries from the compiler)
    CALL DPBTRF('U',KNIN,IUP,PAMAT,IDAMAT,INFO)
  ELSE
    CALL DGBTRF(KNIN,KNIN,IDOWN,IUP,PAMAT,IDAMAT,KPM2(1,0),INFO)
  ENDIF
  IRHS = 1
  IF (INFO .EQ. 0) THEN
    IF (ISYM .EQ. 1) THEN
      CALL DPBTRS('U',KNIN,IUP,IRHS,PAMAT,IDAMAT,PYINPP,KNIN,INFO2)
    ELSE
      CALL DGBTRS('N',KNIN,IDOWN,IUP,IRHS,PAMAT,IDAMAT,KPM2(1,0),PYINPP, &
        &      KNIN,INFO2)
    ENDIF
  ELSE
    PRINT *,' ERROR IN SP/GBTRF: INFO = ',INFO
    !%OS        STOP 'INFO'
    RETURN
  ENDIF
  !
  !     4. COMPUTE NEW VALUES OF Y_K (NON-STANDARD CUBIC SPLINE ONLY)
  !
  IF (ITAUVAL .EQ. 1) THEN
    DO K=1,KNIN
      IKP1 = KPM2(K,+1)
      IKM1 = KPM2(K,-1)
      PYINNEW(K) = PYIN(K) - FTAUK(K) * &
        &      ((PYINPP(IKP1)-PYINPP(K))*WOHK(K) &
        &      - (PYINPP(K)-PYINPP(IKM1))*WOHK(IKM1))
    END DO
    !
  ENDIF

  IF (INFO2 .LT. 0) THEN
    PRINT *,' ERROR IN SP/GBTRS: INFO2 = ',INFO2
    !%OS        STOP 'INFO2'
    RETURN
  ENDIF
  !
  DEALLOCATE(PAMAT)
  !
  RETURN
END SUBROUTINE CBFITBND
