!......................................................................
!     
SUBROUTINE INTQUADRATIC(PXIN,PYIN,KNIN,PXOUT,PYOUT,PYOUTP,PYOUTPP,PYOUTINT,KNOUT,KOPTDER,KOPTXPOL,NBC)
  !
  !   COMPUTE LINEAR INTERPOLATION OF (PXIN,PYIN) ON (PXOUT,PYOUT).
  !
  !   KOPTDER = 0: COMPUTE ONLY PYOUT (PYOUTP, PYOUTPP NOT USED)
  !   KOPTDER = 1: COMPUTE ALSO 1ST DER. IN PYOUTP (PYOUTPP NOT USED)
  !   KOPTDER = 2: AS 1 AND 2ND DER. IN PYOUTPP
  !   KOPTDER = 3: AS 2 AND INTEGRAL FROM (PXIN(1) IN PYOUTPP
  !
  !   KOPTXPOL = 0: SEND MESSAGE IF NEED TO EXTRAPOLATE
  !   KOPTXPOL = 1: LINEAR EXTRAPOLATION WITH CONTINUOUS DERIVATIVE
  !   KOPTXPOL = -1: LINEAR EXTRAPOLATION WITH LAST THREE Y VALUES
  !   KOPTXPOL = 2: QUADRATIC EXTRAPOLATION WITH CONTINUOUS DERIVATIVE
  !   KOPTXPOL = -2: QUADRATIC EXTRAPOLATION WITH LAST THREE Y VALUES
  !   KOPTXPOL = 21: QUADRATIC EXTRAPOLATION WITHIN ONE DELTA_X, THEN LINEAR
  !   KOPTXPOL = -21: SAME AS 21 BUT USES VALUES INSTEAD OF CLOSEST DERIVATIVE FOR EXTRAPOLATION
  !   KOPTXPOL = 10: Y=Y_EDGE FOR EXTRAPOLATION => CONSTANT EXTRAPOLATION
  !   KOPTXPOL = -10: Y=0 FOR EXTRAPOLATION
  !
  ! SOLVES FOR QUADRATICS WITH CONTINUOUS DERIVATIVE AND 
  !   NBC = 0:  2ND DERIVATIVE =0 AT LEFT (DEFAULT)
  !   NBC = 1:  2ND DERIVATIVE =0 AT RIGHT
  !
  USE PREC_RKIND
  IMPLICIT NONE
  REAL(RKIND) :: ZSIX, ZTHREE, ZTWO, ZONE
  PARAMETER(ZSIX=6._RKIND, ZTHREE=3._RKIND, ZTWO=2._RKIND, ZONE=1._RKIND)
  REAL(RKIND) :: ALFA
  PARAMETER(ALFA = 1._RKIND)
  ! arguments
  INTEGER :: KNIN, KNOUT, KOPTDER, KOPTXPOL
  INTEGER, OPTIONAL ::  NBC
  REAL(RKIND) :: PXIN(KNIN), PYIN(KNIN), PXOUT(KNOUT)
  REAL(RKIND) :: PYOUT(KNOUT), PYOUTP(KNOUT), PYOUTPP(KNOUT), PYOUTINT(KNOUT)
  !
  REAL(RKIND) :: BCOEF(KNIN-1), ZYINP(KNIN), ZYINT_XIN(KNIN), ZYIN_XPOL, ZYINT, &
       & ZXLFTDEL, ZXRGTDEL, ZX1, ZXN, ZY1, ZYN, ZY1P, ZX2, ZXNN, ZYNP, ZY2, ZYNN
  INTEGER I,J,IOPTXPOL, INBC, ICONTDER, K, KLO, KHI
  !
  !
  ! VARIABLES RELATED TO FUNCTIONS:
  REAL(RKIND) :: FQQQ0, FQQQ1, FQQQ2, &
    &  FLINEAR, FLINEARP,  &
    &  FQQQM1, FQDQM1, FLINEARM1, FLINXYP, FLINXYPM1, FPARABOLP
  REAL(RKIND) :: A1, A2, A3, A4, B1, B2, B3, B4, PX
  REAL(RKIND) :: FB0, FB1, FB2, FA0, FA1, FD2, FD1, FD0, FQDQ0, FQDQ1, FQDQ2
  REAL(RKIND) :: X1, F1, P1, X2, F2
  REAL(RKIND) :: X3, F3
  !
  !.......................................................................
  !*COMDECK QUAQQQ
  ! ----------------------------------------------------------------------
  ! --     STATEMENT FUNCTION FOR QUADRATIC INTERPOLATION               --
  ! --                         19.01.87            AR        CRPP       --
  ! --                                                                  --
  ! -- QUADRATIC INTERPOLATION OF A FUNCTION F(X)                       --
  ! -- THE SIX PARAMETERS A1,A2,A3,B1,B2,B3 ARE DEFINED AS FOLLOWS:     --
  ! -- F(B1) = A1 , F(B2) = A2 , F(B3) = A3                             --
  ! ----------------------------------------------------------------------
  !
  FB2(A1,A2,A3,B1,B2,B3) = &
    &               ((A1-A2)/(B1-B2)-(A1-A3)/(B1-B3))/(B2-B3)
  FB1(A1,A2,A3,B1,B2,B3) = ((A1-A2)/(B1-B2))- &
    &         FB2(A1,A2,A3,B1,B2,B3)*(B1+B2)
  FB0(A1,A2,A3,B1,B2,B3) = A1-FB1(A1,A2,A3,B1,B2,B3)*B1 &
    &         -FB2(A1,A2,A3,B1,B2,B3)*B1*B1
  ! ----------------------------------------------------------------------
  ! -- FQQQ0 GIVES THE VALUE OF THE FUNCTION AT THE POINT PX            --
  ! -- FQQQ0(......,PX) = F(PX)                                         --
  ! ----------------------------------------------------------------------
  FQQQ0(A1,A2,A3,B1,B2,B3,PX) = FB0(A1,A2,A3,B1,B2,B3) + &
    &                                 PX * (FB1(A1,A2,A3,B1,B2,B3) + &
    &                                 PX * FB2(A1,A2,A3,B1,B2,B3))
  ! ----------------------------------------------------------------------
  ! -- FQQQ1 GIVES THE VALUE OF THE FIRST DERIVATIVE OF F(X) AT PX      --
  ! -- FQQQ1(......,PX) = DF/DX (PX)                                    --
  ! ----------------------------------------------------------------------
  FQQQ1(A1,A2,A3,B1,B2,B3,PX) = FB1(A1,A2,A3,B1,B2,B3) + &
    &     ZTWO * PX * FB2(A1,A2,A3,B1,B2,B3)
  ! ----------------------------------------------------------------------
  ! -- FQQQ2 GIVES THE VALUE OF THE SECOND DERIVATIVE OF F(X) AT PX     --
  ! -- FQQQ2(......,PX) = D2F/DX2 (PX)                                  --
  ! ----------------------------------------------------------------------
  FQQQ2(A1,A2,A3,B1,B2,B3) = ZTWO * FB2(A1,A2,A3,B1,B2,B3)
  ! ----------------------------------------------------------------------
  ! -- FQQQM1 GIVES THE VALUE OF THE INTEGRAL OF F(X) FROM B1 TO PX:
  !
  FQQQM1(A1,A2,A3,B1,B2,B3,PX) = &
    & (PX-B1)*(FB0(A1,A2,A3,B1,B2,B3) + &
    &  0.5_RKIND*(PX+B1)*FB1(A1,A2,A3,B1,B2,B3) + &
    &  FB2(A1,A2,A3,B1,B2,B3)/3._RKIND*(PX*(PX+B1)+B1*B1))
  !.......................................................................
  !*COMDECK QUAQDQ
  ! ----------------------------------------------------------------------
  ! --     STATEMENT FUNCTION FOR QUADRATIC INTERPOLATION               --
  ! --                         19.01.87            AR        CRPP       --
  ! --                                                                  --
  ! -- QUADRATIC INTERPOLATION OF A FUNCTION F(X)                       --
  ! -- THE FIVE PARAMETERS X1,F1,P1,X2,F2    ARE DEFINED AS FOLLOWS:    --
  ! -- F(X1) = F1 , DF/DX(X1) = P1 , F(X2) = F2                         --
  ! ----------------------------------------------------------------------
  !
  FD2(X1,F1,P1,X2,F2) = ((F2-F1)/(X2-X1) - P1) / (X2-X1)
  FD1(X1,F1,P1,X2,F2) = P1 - ZTWO*X1*FD2(X1,F1,P1,X2,F2)
  FD0(X1,F1,P1,X2,F2) = F1 - X1*(X1*FD2(X1,F1,P1,X2,F2) + &
    &                                     FD1(X1,F1,P1,X2,F2))
  ! ----------------------------------------------------------------------
  ! -- FQDQ0 GIVES THE VALUE OF THE FUNCTION AT POINT PX                --
  ! -- FQDQ0(......,PX) = F(PX)                                         --
  ! ----------------------------------------------------------------------
  FQDQ0(X1,F1,P1,X2,F2,PX) = FD0(X1,F1,P1,X2,F2) + &
    &                              PX * (FD1(X1,F1,P1,X2,F2) + &
    &                                    PX * FD2(X1,F1,P1,X2,F2))
  ! ----------------------------------------------------------------------
  ! -- FQDQ1 GIVES THE VALUE OF THE FIRST DERIVATIVE OF F(X) AT PX      --
  ! -- FQDQ1(......,PX) = DF/DX (PX)                                    --
  ! ----------------------------------------------------------------------
  FQDQ1(X1,F1,P1,X2,F2,PX) = FD1(X1,F1,P1,X2,F2) + &
    &                              ZTWO* PX * FD2(X1,F1,P1,X2,F2)
  ! ----------------------------------------------------------------------
  ! -- FQDQ2 GIVES THE VALUE OF THE SECOND DERIVATIVE OF F(X) AT PX     --
  ! -- FQDQ2(......,PX) = D2F/DX2 (PX)                                  --
  ! ----------------------------------------------------------------------
  FQDQ2(X1,F1,P1,X2,F2) = ZTWO * FD2(X1,F1,P1,X2,F2)
  ! ----------------------------------------------------------------------
  ! -- FQDQM1 GIVES THE VALUE OF THE INTEGRAL OF F(X) FROM X1 TO PX:
  !
  FQDQM1(X1,F1,P1,X2,F2,PX) = &
    & (PX-X1)*(FD0(X1,F1,P1,X2,F2) + &
    &  0.5_RKIND*(PX+X1)*FD1(X1,F1,P1,X2,F2) + &
    &  FD2(X1,F1,P1,X2,F2)/3._RKIND*(PX*(PX+X1)+X1*X1))
  !-----------------------------------------------------------------------
  FPARABOLP(X1,X2,X3,F1,F2,F3,PX) = &
    &  ((PX-X1)+(PX-X2))*F3/((X3-X1)*(X3-X2))+ &
    &  ((PX-X1)+(PX-X3))*F2/((X2-X1)*(X2-X3))+ &
    &  ((PX-X2)+(PX-X3))*F1/((X1-X2)*(X1-X3))
  !.......................................................................
  !     LINEAR
  !
  FLINEAR(X1,F1,X2,F2,PX) = F1 + (PX-X1)/(X2-X1) * (F2-F1)
  FLINEARP(X1,F1,X2,F2) = (F2-F1) / (X2-X1)
  FLINEARM1(X1,F1,X2,F2,PX) = (PX-X1)*(F1+0.5_RKIND*(PX-X1)*(F2-F1)/(X2-X1))
  FLINXYP(X1,F1,P1,PX) = P1*(PX-X1) + F1
  FLINXYPM1(X1,F1,P1,PX) = (PX-X1)*(F1+0.5_RKIND*P1*(PX-X1))
  !
  !-----------------------------------------------------------------------
  ! 0. DEFAULTS
  !
  ICONTDER = 1
  IF (KOPTXPOL .LT. 0) ICONTDER = 0
  IOPTXPOL=ABS(KOPTXPOL)
  INBC = 0
  IF (PRESENT(NBC)) INBC = NBC
  !
  ! DETERMINE DERIVATIVE ZYINP(I) SUCH THAT QUADRATICS ARE THEN DEFINED BY XI, YI, YPI, XI+1, YI+1 WITHIN [XI,XI+1]
  ! VARIOUS OPTIONS ^CHOSEN THROUGH NBC
  !
!!$  ACOEF(1:KNIN-1) = PYIN(1:KNIN-1)
!!$  CCOEF(1:KNIN-1) = PYIN(2:KNIN)
  IF (INBC .EQ. 1) THEN
    ! FIND COEFFICIENTS: Y(X)=AI*(XI+1-X)^2 + BI*(XI+1-X)*(X-XI) + CI*(X-XI)^2 FOR [XI,XI+1], I=1,N-1
    ! AI = YI ; CI=YI+1
    ! BI/HI + BI+1/HI+1 = 2 YI+1 (1/HI+1/HI+1) FOR CONTINUITY OF DERIVATIVES
    ! NEED ONE BOUNDARY CONDITION
    ! d2y/dx2=0 at knin
    BCOEF(KNIN-1) = PYIN(KNIN-1) + PYIN(KNIN)
    DO I=KNIN-2,1,-1
      BCOEF(I) = (PXIN(I+1)-PXIN(I))*(-BCOEF(I+1)/(PXIN(I+2)-PXIN(I+1)) + & 
        & 2._RKIND*PYIN(I+1)*(1._RKIND/(PXIN(I+1)-PXIN(I))+1._RKIND/(PXIN(I+2)-PXIN(I+1))))
    END DO
    ! d2y/dx2=0 at 1
    !equiv    BCOEF(1) = PYIN(1) + PYIN(2)
    !equiv    DO I=2,KNIN-1
    !equiv      BCOEF(I) = (PXIN(I+1)-PXIN(I))*(-BCOEF(I-1)/(PXIN(I)-PXIN(I-1)) + & 
    !equiv        & 2._RKIND*PYIN(I)*(1._RKIND/(PXIN(I)-PXIN(I-1))+1._RKIND/(PXIN(I+1)-PXIN(I))))
    !equiv    END DO
    DO I=1,KNIN-1
      ZYINP(I) = (BCOEF(I)-2._RKIND*PYIN(I))/(PXIN(I+1)-PXIN(I))
    END DO
    ZYINP(KNIN) = (-BCOEF(KNIN-1)+2._RKIND*PYIN(KNIN))/(PXIN(KNIN)-PXIN(KNIN-1))
  ELSE
    ! USE 3 POINTS TO DEFINE QUADRATICS EXCEPT FOR LAST INTERVALS
    DO I=1,KNIN-2
      ZYINP(I) = FPARABOLP(PXIN(I),PXIN(I+1),PXIN(I+2),PYIN(I),PYIN(I+1),PYIN(I+2),PXIN(I))
    END DO
    ZYINP(KNIN-1) = FPARABOLP(PXIN(KNIN-2),PXIN(KNIN-1),PXIN(KNIN),PYIN(KNIN-2),PYIN(KNIN-1),PYIN(KNIN),PXIN(KNIN-1))
    ZYINP(KNIN) = FPARABOLP(PXIN(KNIN-2),PXIN(KNIN-1),PXIN(KNIN),PYIN(KNIN-2),PYIN(KNIN-1),PYIN(KNIN),PXIN(KNIN))
  END IF
  !write(*,'(1p2e14.5)') (pxin(i),zyinp(i),i=1,knin)
  !
  ! COMPUTE INT. UP TO EACH INPUT INTERVAL
  IF (KOPTDER .GE. 3) THEN
    ZYINT_XIN(1) = 0._RKIND
    DO I=1,KNIN-1
      ZYINT_XIN(I+1) = ZYINT_XIN(I) + FQDQM1(PXIN(I),PYIN(I),ZYINP(I),PXIN(I+1),PYIN(I+1),PXIN(I+1))
    END DO
  END IF
  !
  ! LOOP OVER PXOUT POINTS WHICH CAN BE IN RANDOM ORDER
  DO 100 J=1,KNOUT
    IF ((PXOUT(J) .LT. PXIN(1)) .OR. (PXOUT(J) .GT. PXIN(KNIN))) GO TO 200
    !
    ! 1.1 POINTS INSIDE INTERVAL [XIN(1),XIN(KNIN)]
    ! FIND PXIN INTERVAL BY BI-SECTION
    KLO=1
    KHI=KNIN
10  CONTINUE
    IF (KHI-KLO.GT.1) THEN
      K=(KHI+KLO)/2
      IF(PXIN(K) .GT. PXOUT(J))THEN
        KHI=K
      ELSE
        KLO=K
      ENDIF
      GOTO 10
    ENDIF
    SELECT CASE (KOPTDER)
    CASE (0)
      PYOUT(J)=FQDQ0(PXIN(KLO),PYIN(KLO),ZYINP(KLO),PXIN(KHI),PYIN(KHI),PXOUT(J))
    CASE (1)
      PYOUT(J)=FQDQ0(PXIN(KLO),PYIN(KLO),ZYINP(KLO),PXIN(KHI),PYIN(KHI),PXOUT(J))
      PYOUTP(J)=FQDQ1(PXIN(KLO),PYIN(KLO),ZYINP(KLO),PXIN(KHI),PYIN(KHI),PXOUT(J))
    CASE (2)
      PYOUT(J)=FQDQ0(PXIN(KLO),PYIN(KLO),ZYINP(KLO),PXIN(KHI),PYIN(KHI),PXOUT(J))
      PYOUTP(J)=FQDQ1(PXIN(KLO),PYIN(KLO),ZYINP(KLO),PXIN(KHI),PYIN(KHI),PXOUT(J))
      PYOUTPP(J)=FQDQ2(PXIN(KLO),PYIN(KLO),ZYINP(KLO),PXIN(KHI),PYIN(KHI))
    CASE (3)
      PYOUT(J)=FQDQ0(PXIN(KLO),PYIN(KLO),ZYINP(KLO),PXIN(KHI),PYIN(KHI),PXOUT(J))
      PYOUTP(J)=FQDQ1(PXIN(KLO),PYIN(KLO),ZYINP(KLO),PXIN(KHI),PYIN(KHI),PXOUT(J))
      PYOUTPP(J)=FQDQ2(PXIN(KLO),PYIN(KLO),ZYINP(KLO),PXIN(KHI),PYIN(KHI))
      PYOUTINT(J) = ZYINT_XIN(KLO) + FQDQM1(PXIN(KLO),PYIN(KLO),ZYINP(KLO),PXIN(KHI),PYIN(KHI),PXOUT(J))
    END SELECT
    !
    GO TO 100
    !
    !     2 POINT OUTSIDE INTERVAL
    !
200 CONTINUE
    !
    !     2.1 IF KOPTXPOL=0, PRINT WARNING AND RETURN OR STOP
    !
    IF (IOPTXPOL .EQ. 0) THEN
      PRINT *,' PXOUT(',J,')=',PXOUT(J),' OUTSIDE INTERVAL [',PXIN(1),',',PXIN(KNIN),']'
      RETURN
      !        STOP 'IOPTXPOL=0'
    ENDIF
    !
    !     2.2 COMPUTE VALUES FOR POINTS ON THE LEFT OF PXIN(1)
    !           EXTRAPOLATION DEPENDS ON VALUE OF KOPTXPOL
    !
    IF (PXOUT(J) .LT. PXIN(1)) THEN
      ! ZY1 REFERS TO FIRST KNOWN POINTS AND ZY2 TO SECOND KNOWN POINT, 
      ! TYPICALLY AT PXIN(1) AND PXIN(2) OR PXIN(1)-ALFA*H AND PXIN(1)
      ZXLFTDEL = PXIN(1) - ALFA*(PXIN(2) - PXIN(1))
      !
      !   2.2.1 SPECIAL PART [XIN(1)-ALFA*H,XIN(1)] IF IOPTXPOL>20
      IF ((PXOUT(J) .GE. ZXLFTDEL) .AND. (IOPTXPOL.EQ.21)) THEN
        ! QUADRATIC PART
        IF (ICONTDER .EQ. 1) THEN
          PYOUT(J) = FQDQ0(PXIN(1),PYIN(1),ZYINP(1),PXIN(2),PYIN(2),PXOUT(J))
          IF (KOPTDER .GE. 1) PYOUTP(J) = FQDQ1(PXIN(1),PYIN(1),ZYINP(1),PXIN(2),PYIN(2),PXOUT(J))
          IF (KOPTDER .GE. 2) PYOUTPP(J) = FQDQ2(PXIN(1),PYIN(1),ZYINP(1),PXIN(2),PYIN(2))
          ! INTEGRATES DIRECTLY FROM PXIN(1)
          IF (KOPTDER .GE. 3) PYOUTINT(J) = FQDQM1(PXIN(1),PYIN(1),ZYINP(1),PXIN(2),PYIN(2),PXOUT(J))
        ELSE
          ! KOPTXPOL = -21
          PYOUT(J) = FQQQ0(PYIN(1),PYIN(2),PYIN(3),PXIN(1),PXIN(2),PXIN(3),PXOUT(J))
          IF (KOPTDER .GE. 1) PYOUTP(J) = FQQQ1(PYIN(1),PYIN(2),PYIN(3),PXIN(1),PXIN(2),PXIN(3),PXOUT(J))
          IF (KOPTDER .GE. 2) PYOUTPP(J) = FQQQ2(PYIN(1),PYIN(2),PYIN(3),PXIN(1),PXIN(2),PXIN(3))
          IF (KOPTDER .GE. 3) PYOUTINT(J) = FQQQM1(PYIN(1),PYIN(2),PYIN(3),PXIN(1),PXIN(2),PXIN(3),PXOUT(J))
        END IF
        !
      ELSE
        !
        ! 2.2.2 EXTRAPOLATION FAR LEFT: X<XIN(1)-ALFA*H OR X<XIN(1) IF NO ALFA PART CONSIDERED
        !
        IF ((IOPTXPOL .EQ. 1) .OR. (IOPTXPOL .EQ. 21)) THEN
          ! LINEAR EXTRAPOLATION
          SELECT CASE (KOPTXPOL)
          CASE (1)
            ZX1 = PXIN(1)
            ZY1 = PYIN(1)
            ZY1P = ZYINP(1)
            ZYINT = 0._RKIND
          CASE(-1)
            ZX1 = PXIN(1)
            ZY1 = PYIN(1)
            ZY1P = FLINEARP(PXIN(1),PYIN(1),PXIN(2),PYIN(2))
            ZYINT = 0._RKIND
          CASE (21)
            ZX1 = ZXLFTDEL
            ZY1 = FQDQ0(PXIN(1),PYIN(1),ZYINP(1),PXIN(2),PYIN(2),ZXLFTDEL)
            ZY1P = FQDQ1(PXIN(1),PYIN(1),ZYINP(1),PXIN(2),PYIN(2),ZXLFTDEL)
            ZYINT = FQDQM1(PXIN(1),PYIN(1),ZYINP(1),PXIN(2),PYIN(2),ZXLFTDEL)
          CASE (-21)
            ZX1 = ZXLFTDEL
            ZY1 = FQQQ0(PYIN(1),PYIN(2),PYIN(3),PXIN(1),PXIN(2),PXIN(3),ZXLFTDEL)
            ZY1P = FQQQ1(PYIN(1),PYIN(2),PYIN(3),PXIN(1),PXIN(2),PXIN(3),ZXLFTDEL)
            ZYINT = FQQQM1(PYIN(1),PYIN(2),PYIN(3),PXIN(1),PXIN(2),PXIN(3),ZXLFTDEL)
          END SELECT
          PYOUT(J) = FLINXYP(ZX1,ZY1,ZY1P,PXOUT(J))
          IF (KOPTDER .GE. 1) PYOUTP(J) = ZY1P
          IF (KOPTDER .GE. 2) PYOUTPP(J) = 0._RKIND
          IF (KOPTDER .GE. 3) PYOUTINT(J) = ZYINT + FLINXYPM1(ZX1,ZY1,ZY1P,PXOUT(J))
        ELSE IF (IOPTXPOL .EQ. 10) THEN
          ! CONSTANT OUTSIDE PXIN
          IF (ICONTDER .EQ. 1) THEN
            ! KOPTXPOL = +10
            ZYIN_XPOL = PYIN(1)
          ELSE
            ! KOPTXPOL = -10
            ZYIN_XPOL = 0._RKIND
          END IF
          PYOUT(J) = ZYIN_XPOL
          IF (KOPTDER .GE. 1) PYOUTP(J) = 0._RKIND
          IF (KOPTDER .GE. 2) PYOUTPP(J) = 0._RKIND
          ! INTEGRATES FROM PXIN(1) DIRECTLY
          IF (KOPTDER .GE. 3) PYOUTINT(J) = (PXOUT(J)-PXIN(1))*ZYIN_XPOL
          !
        ELSE
          ! QUADRATIC EXTRAPOLATION
          ZX1 = PXIN(1)
          ZY1 = PYIN(1)
          ZX2 = PXIN(2)
          ZY2 = PYIN(2)
          ZYINT = 0._RKIND
          IF (ICONTDER .EQ. 1) THEN
            ! KOPTXPOL = 2
            ZY1P = ZYINP(1)
          ELSE
            ! KOPTXPOL = -2
            ZY1P = FQQQ1(PYIN(1),PYIN(2),PYIN(3),PXIN(1),PXIN(2),PXIN(3),PXIN(1))
          ENDIF
          PYOUT(J) = FQDQ0(ZX1,ZY1,ZY1P,ZX2,ZY2,PXOUT(J))
          IF (KOPTDER .GE. 1) PYOUTP(J) = FQDQ1(ZX1,ZY1,ZY1P,ZX2,ZY2,PXOUT(J))
          IF (KOPTDER .GE. 2) PYOUTPP(J) = FQDQ2(ZX1,ZY1,ZY1P,ZX2,ZY2)
          IF (KOPTDER .GE. 3) PYOUTINT(J) = ZYINT + FQDQM1(ZX1,ZY1,ZY1P,ZX2,ZY2,PXOUT(J))
        END IF
      END IF
      !
    ELSE ! PXOUT(J) .GT. PXIN(KNIN)
      !
      !     2.3 COMPUTE VALUES FOR POINTS ON THE RIGHT OF PXIN(1)
      !         EXTRAPOLATION DEPENDS ON VALUE OF KOPTXPOL
      !
      ! ZYN REFERS TO FIRST KNOWN POINTS AND ZYNN TO SECOND KNOWN POINT, 
      ! TYPICALLY AT PXIN(KNIN) AND PXIN(KNIN-1) OR PXIN(KNIN)+ALFA*H AND PXIN(KNIN)
      ! YP(PXIN(KNIN)) FROM SPLINE
      ZXRGTDEL = PXIN(KNIN) + ALFA * (PXIN(KNIN) - PXIN(KNIN-1))
      !
      !   2.3.1 SPECIAL PART ]XIN(KNIN),XIN(KNIN)+ALFA*H] IF IOPTXPOL>20
      IF ((PXOUT(J) .LE. ZXRGTDEL) .AND. (IOPTXPOL.EQ.21)) THEN
        ZYINT = ZYINT_XIN(KNIN)
        ! QUADRATIC
        IF (ICONTDER .EQ. 1) THEN
          ! KOPTXPOL = +21
          PYOUT(J) = FQDQ0(PXIN(KNIN),PYIN(KNIN),ZYINP(KNIN),PXIN(KNIN-1),PYIN(KNIN-1),PXOUT(J))
          IF (KOPTDER .GE. 1) PYOUTP(J) = FQDQ1(PXIN(KNIN),PYIN(KNIN),ZYINP(KNIN),PXIN(KNIN-1),PYIN(KNIN-1),PXOUT(J))
          IF (KOPTDER .GE. 2) PYOUTPP(J) = FQDQ2(PXIN(KNIN),PYIN(KNIN),ZYINP(KNIN),PXIN(KNIN-1),PYIN(KNIN-1))
          ! INTEGRATES FROM PXIN(KNIN) SO ADD INTEGRAL UP TO PXIN(KNIN)
          IF (KOPTDER .GE. 3) PYOUTINT(J) = ZYINT + FQDQM1(PXIN(KNIN),PYIN(KNIN),ZYINP(KNIN),PXIN(KNIN-1),PYIN(KNIN-1),PXOUT(J))
        ELSE
          ! KOPTXPOL = -21
          PYOUT(J) = FQQQ0(PYIN(KNIN),PYIN(KNIN-1),PYIN(KNIN-2),PXIN(KNIN),PXIN(KNIN-1),PXIN(KNIN-2),PXOUT(J))
          IF (KOPTDER .GE. 1) PYOUTP(J) = &
               & FQQQ1(PYIN(KNIN),PYIN(KNIN-1),PYIN(KNIN-2),PXIN(KNIN),PXIN(KNIN-1),PXIN(KNIN-2),PXOUT(J))
          IF (KOPTDER .GE. 2) PYOUTPP(J) = &
               & FQQQ2(PYIN(KNIN),PYIN(KNIN-1),PYIN(KNIN-2),PXIN(KNIN),PXIN(KNIN-1),PXIN(KNIN-2))
          IF (KOPTDER .GE. 3) PYOUTINT(J) = ZYINT + &
               & FQQQM1(PYIN(KNIN),PYIN(KNIN-1),PYIN(KNIN-2),PXIN(KNIN),PXIN(KNIN-1),PXIN(KNIN-2),PXOUT(J))
        END IF
        !
      ELSE
        !
        ! 2.3.2 EXTRAPOLATION FAR RIGHT: X>XIN(KNIN)+ALFA*H OR X>XIN(KNIN) IF NO ALFA PART CONSIDERED
        !
        IF ((IOPTXPOL .EQ. 1) .OR. (IOPTXPOL .EQ. 21)) THEN
          ! LINEAR EXTRAPOLATION
          SELECT CASE (KOPTXPOL)
          CASE (1)
            ZXN = PXIN(KNIN)
            ZYN = PYIN(KNIN)
            ZYNP = ZYINP(KNIN)
            ZYINT = ZYINT_XIN(KNIN)
          CASE(-1)
            ZXN = PXIN(KNIN)
            ZYN = PYIN(KNIN)
            ZYNP = FLINEARP(PXIN(KNIN),PYIN(KNIN),PXIN(KNIN-1),PYIN(KNIN-1))
            ZYINT = ZYINT_XIN(KNIN)
          CASE (21)
            ZYINT = ZYINT_XIN(KNIN) + FQDQM1(PXIN(KNIN),PYIN(KNIN),ZYINP(KNIN),PXIN(KNIN-1),PYIN(KNIN-1),ZXRGTDEL)
            ZXN = ZXRGTDEL
            ZYN = FQDQ0(PXIN(KNIN),PYIN(KNIN),ZYINP(KNIN),PXIN(KNIN-1),PYIN(KNIN-1),ZXRGTDEL)
            ZYNP = FQDQ1(PXIN(KNIN),PYIN(KNIN),ZYINP(KNIN),PXIN(KNIN-1),PYIN(KNIN-1),ZXRGTDEL)
          CASE (-21)
            ZYINT = ZYINT_XIN(KNIN) + &
              & FQQQM1(PYIN(KNIN),PYIN(KNIN-1),PYIN(KNIN-2),PXIN(KNIN),PXIN(KNIN-1),PXIN(KNIN-2),ZXRGTDEL)
            ZXN = ZXRGTDEL
            ZYN = FQQQ0(PYIN(KNIN),PYIN(KNIN-1),PYIN(KNIN-2),PXIN(KNIN),PXIN(KNIN-1),PXIN(KNIN-2),ZXRGTDEL)
            ZYNP = FQQQ1(PYIN(KNIN),PYIN(KNIN-1),PYIN(KNIN-2),PXIN(KNIN),PXIN(KNIN-1),PXIN(KNIN-2),ZXRGTDEL)
          END SELECT
          PYOUT(J) = FLINXYP(ZXN,ZYN,ZYNP,PXOUT(J))
          IF (KOPTDER .GE. 1) PYOUTP(J) = ZYINP(KNIN)
          IF (KOPTDER .GE. 2) PYOUTPP(J) = 0._RKIND
          IF (KOPTDER .GE. 3) PYOUTINT(J) = ZYINT + FLINXYPM1(ZXN,ZYN,ZYNP,PXOUT(J))
        ELSE IF (IOPTXPOL .EQ. 10) THEN
          ! CONSTANT OUTSIDE PXIN
          IF (ICONTDER .EQ. 1) THEN
            ! KOPTXPOL = +10
            ZYIN_XPOL = PYIN(KNIN)
          ELSE
            ! KOPTXPOL = -10
            ZYIN_XPOL = 0._RKIND
          END IF
          PYOUT(J) = ZYIN_XPOL
          IF (KOPTDER .GE. 1) PYOUTP(J) = 0._RKIND
          IF (KOPTDER .GE. 2) PYOUTPP(J) = 0._RKIND
          IF (KOPTDER .GE. 3) PYOUTINT(J) = ZYINT_XIN(KNIN) + (PXOUT(J)-PXIN(KNIN))*ZYIN_XPOL
          !
        ELSE
          ! QUADRATIC EXTRAPOLATION
          ZXN = PXIN(KNIN)
          ZYN = PYIN(KNIN)
          ZXNN = PXIN(KNIN-1)
          ZYNN = PYIN(KNIN-1)
          ZYINT = ZYINT_XIN(KNIN)
          IF (KOPTXPOL .EQ. -2) ZYINP(KNIN) = &
               & FQQQ1(PYIN(KNIN),PYIN(KNIN-1),PYIN(KNIN-2),PXIN(KNIN),PXIN(KNIN-1),PXIN(KNIN-2),PXIN(KNIN))
          PYOUT(J) = FQDQ0(ZXN,ZYN,ZYINP(KNIN),ZXNN,ZYNN,PXOUT(J))
          IF (KOPTDER .GE. 1) PYOUTP(J) = FQDQ1(ZXN,ZYN,ZYINP(KNIN),ZXNN,ZYNN,PXOUT(J))
          IF (KOPTDER .GE. 2) PYOUTPP(J) = FQDQ2(ZXN,ZYN,ZYINP(KNIN),ZXNN,ZYNN)
          IF (KOPTDER .GE. 3) PYOUTINT(J) = ZYINT + FQDQM1(ZXN,ZYN,ZYINP(KNIN),ZXNN,ZYNN,PXOUT(J))
        END IF
      END IF
    END IF
    !
100 END DO
  RETURN 
END SUBROUTINE INTQUADRATIC
