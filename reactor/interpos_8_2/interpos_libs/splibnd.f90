!-----------------------------------------------------------------------
      SUBROUTINE SPLIBND(XA,YA,Y2A,N,X,Y,YP,YPP,KOPTIN)
      USE prec_rkind
      implicit none
      REAL(RKIND) :: zsix, zthree, ztwo, zone
      PARAMETER(zsix=6._RKIND, zthree=3._RKIND, ztwo=2._RKIND, zone=1._RKIND)
      REAL(RKIND) :: XA, YA, Y2A, X, Y, YP, YPP
      integer N, KOPTIN
      DIMENSION :: XA(N),YA(N),Y2A(N)
!
!   ABS(KOPTIN):
!     KOPTIN = 0: STOP WITH ERROR MESSAGE IF OUT OF BOUND
!     KOPTIN = 1: LINEAR EXTRAPOLATION
!     KOPTIN = 2: USE QUADRATIC INTERPOLATION IF X OUT OF BOUND
!     KOPTIN = 3: USE CUBIC INTERPOLATION IF X OUT OF BOUND
!     KOPTIN = 21: USE QUADRATIC WITHIN ALFA*DELTA_X AND LINEAR FURTHER
!     KOPTIN = 31: USE CUBIC WITHIN ALFA*DELTA_X AND LINEAR    FURTHER
!     KOPTIN = 32: USE CUBIC WITHIN ALFA*DELTA_X AND QUADRATIC FURTHER
!
!   KOPTIN >=0 => INCONTDER = 1
!   KOPTIN < 0 => INCONTDER = 0
!     ICONTDER = 1: VALUE AND 1ST DER. CONTINUOUS AT END OF INTERVAL, THUS
!     .             USES CUBIC SPLINE OF LAST INTERVAL TO CONTINUE
!     ICONTDER = 0: ONLY VALUE CONTINUOUS AND USES VALUES AT LAST BUT ONE,
!     .             TWO, THREE POINTS TO EXTRAPOLATE (BETTER IF DER. AT EDGE
!     .             IS WILD)
!
!-----------------------------------------------------------------------
!
      REAL(RKIND) :: ALFA
      PARAMETER(ALFA = 1._RKIND)
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
      REAL(RKIND) :: h, a, b, dellim, zxedg, zyedg, zypedg, zxdellim, zyxdl, &
     &  zypxdl, zxedgm1, zyedgm1, zyxdlp, zypedgm1
      integer icontder, kopt, klo, khi, k, ikopt, k1, k2, k3, k4, klohi
      REAL(RKIND) :: fc3, x1, f1, p1, x2, &
     &  f2, p2, fc2, fc1, fc0, fqqq0, fqqq1, fqqq2, &
     &  flinear, flinearp, fcccc0, fcccc1, fcccc2, fcccc3, fqdq0, fqdq1, &
     &  fqdq2, fcdcd0, fcdcd1, fcdcd2, fcdcd3, fb1, &
     &  fb2, fa2, fa3, fd2, fd1
      REAL(RKIND) :: a1, a2, a3, a4, b1, b2, b3, b4, px
      REAL(RKIND) :: fb0, fd0, fa0, fa1
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
!-----------------------------------------------------------------------
!.......................................................................
!*COMDECK CUCDCD
! ----------------------------------------------------------------------
! --     STATEMENT FUNCTION FOR CUBIC INTERPOLATION                   --
! --                         19.01.87            AR        CRPP       --
! --                                                                  --
! -- CUBIC INTERPOLATION OF A FUNCTION F(X)                           --
! -- THE SIX ARGUMENTS X1,F1,P1,X2,F2,P2 ARE DEFINED AS FOLLOWS:      --
! -- F(X1) = F1 , F(X2) = F2 , DF/DX(X1) = P1 , DF/DX(X2) = P2        --
! ----------------------------------------------------------------------
!
         FC3(X1,F1,P1,X2,F2,P2) = &
     &      (ZTWO * (F2 - F1) / (X1 - X2) + (P1 + P2)) / &
     &      ((X1 - X2) * (X1 - X2))
         FC2(X1,F1,P1,X2,F2,P2) = &
     &      (ZTHREE * (X1 + X2) * (F1 - F2) / (X1 - X2) - &
     &       P1 * (X1 + ZTWO * X2) - P2 * (X2 + ZTWO * X1)) / &
     &      ((X1 - X2) * (X1 - X2))
         FC1(X1,F1,P1,X2,F2,P2) = &
     &      (zsix * X1 * X2 * (F2 - F1) / (X1 - X2) + &
     &       X2 * P1 * (2 * X1 + X2) + X1 * P2 * (X1 + ZTWO * X2)) / &
     &      ((X1 - X2) * (X1 - X2))
         FC0(X1,F1,P1,X2,F2,P2) = &
     &      (F1 * X2**2 + F2 * X1**2 - X1 * X2 * (X2 * P1 + X1 * P2) + &
     &       ZTWO * X1 * X2 * (F1 * X2 - F2 * X1) / (X1 - X2)) / &
     &      ((X1 - X2) * (X1 - X2))
! ----------------------------------------------------------------------
! -- FCDCD0 GIVES THE VALUE OF THE FUNCTION AT POINT PX               --
! -- FCDCD0(......,PX) = F(PX)                                        --
! ----------------------------------------------------------------------
         FCDCD0(X1,F1,P1,X2,F2,P2,PX) = &
     &              FC0(X1,F1,P1,X2,F2,P2) + &
     &              PX * (FC1(X1,F1,P1,X2,F2,P2) + &
     &                    PX * (FC2(X1,F1,P1,X2,F2,P2) + &
     &                          PX * FC3(X1,F1,P1,X2,F2,P2)))
! ----------------------------------------------------------------------
! -- FCDCD1 GIVES THE VALUE OF THE FIRST DERIVATIVE OF F(X) AT PX:    --
! -- FCDCD1(......,PX) = DF/DX (PX)                                   --
! ----------------------------------------------------------------------
         FCDCD1(X1,F1,P1,X2,F2,P2,PX) = &
     &              FC1(X1,F1,P1,X2,F2,P2) + &
     &              PX * (ZTWO * FC2(X1,F1,P1,X2,F2,P2) + &
     &                    ZTHREE * PX * FC3(X1,F1,P1,X2,F2,P2))
! ----------------------------------------------------------------------
! -- FCDCD2 GIVES THE VALUE OF THE SECOND DERIVATIVE OF F(X) AT PX:   --
! -- FCDCD2(......,PX) = D2F/DX2 (PX)                                 --
! ----------------------------------------------------------------------
         FCDCD2(X1,F1,P1,X2,F2,P2,PX) = &
     &             ZTWO * FC2(X1,F1,P1,X2,F2,P2) + &
     &             zsix * FC3(X1,F1,P1,X2,F2,P2) * PX
! ----------------------------------------------------------------------
! -- FCDCD3 GIVES THE VALUE OF THE THIRD DERIVATIVE OF F(X) AT PX:    --
! -- FCDCD3(......,PX) = D3F/DX3 (PX)                                 --
! ----------------------------------------------------------------------
         FCDCD3(X1,F1,P1,X2,F2,P2,PX) = &
     &                      zsix * FC3(X1,F1,P1,X2,F2,P2)
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
!-----------------------------------------------------------------------
!.......................................................................
!     LINEAR
!
      FLINEAR(X1,F1,X2,F2,PX) = F2 + (PX-X2)/(X2-X1) * (F2-F1)
      FLINEARP(X1,F1,X2,F2) = (F2-F1) / (X2-X1)
!-----------------------------------------------------------------------
      ICONTDER = 1
      IF (KOPTIN .LT. 0) ICONTDER = 0
      KOPT=ABS(KOPTIN)
!
!     1. POINT INSIDE INTERVAL
!
      IF (X.LT.XA(1) .OR. X.GT.XA(N)) GO TO 200
!
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
!%OS      IF (H.EQ.0._RKIND) STOP 'BAD XA INPUT.'
      IF (H .EQ. 0._RKIND) THEN
        PRINT *,'BAD XA INPUT.'
        RETURN
      ENDIF
      A=(XA(KHI)-X)/H
      B=(X-XA(KLO))/H
      Y=A*YA(KLO)+B*YA(KHI)+ &
     &      ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/zsix
      YP=(YA(KHI)-YA(KLO))/H - &
     & ((ZTHREE*A*A-zone)*Y2A(KLO)-(ZTHREE*B*B-ZONE)*Y2A(KHI) )*H/zsix
      YPP=A*Y2A(KLO)+B*Y2A(KHI)
!
      RETURN
!
!     2. POINT OUTSIDE INTERVAL
!
 200  CONTINUE
!
      IKOPT = KOPT
      IF (KOPT .EQ. 0) THEN
        PRINT *,' POINT X=',X,' IS OUTSIDE INTERVAL [',XA(1),',' &
     &    ,XA(N),']'
        IKOPT = 1
!C        STOP 'KOPT=0'
      ENDIF
!
      KLO = 1
      KHI = 2
      IF (X .GT. XA(N)) THEN
        KLO = N-1
        KHI = N
      ENDIF
      H=XA(KHI)-XA(KLO)
      DELLIM = ALFA * H
!%OS      IF (H.EQ.0._RKIND) STOP 'BAD XA INPUT.'
      IF (H .EQ. 0._RKIND) THEN
        PRINT *,'BAD XA INPUT.'
        RETURN
      ENDIF
!
!.......................................................................
!     2.1 LINEAR, IKOPT=1
!
      IF (IKOPT .EQ. 1) THEN
!
!     LINEAR EXTRAPOLATION
        IF (ICONTDER .EQ. 0) THEN
          Y = YA(KHI)
          YP = 0.0_RKIND
          YPP = 0.0_RKIND
!%OS          Y = FLINEAR(XA(KLO),YA(KLO),XA(KHI),YA(KHI),X)
!%OS          YP = FLINEARP(XA(KLO),YA(KLO),XA(KHI),YA(KHI))
!%OS          YPP = 0.0_RKIND
        ELSE
!     COMPUTE VALUE AND 1ST DER. AT EDGE
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
          YPP = 0.0_RKIND
        ENDIF
!
!.......................................................................
!     2.2 LINEAR FAR END OF IKOPT=21 AND IKOPT=31
!
      ELSE IF (MAX(XA(KLO)-X,X-XA(KHI)).GT.DELLIM .AND. &
     &    (IKOPT.EQ.21 .OR. IKOPT.EQ.31)) THEN
!
!     LINEAR EXTRAPOLATION OUTSIDE DELLIM
!     COMPUTE STARTING POINT AND DERIVATIVE FROM END OF QUADR. OR CUBIC
!     INTERPOLATION IN ALFA*DELTA_X INTERVAL
!
        IF (IKOPT .EQ. 21) THEN
!     QUADRATIC
          IF (ICONTDER .EQ. 0) THEN
            K1 = 1
            IF (KHI .EQ. N) K1 = N-2
            K2 = K1 + 1
            K3 = K1 + 2
            ZXDELLIM = XA(1) - DELLIM
            IF (KHI .EQ. N) ZXDELLIM = XA(N) + DELLIM
            ZYXDL = FQQQ0(YA(K1),YA(K2),YA(K3),XA(K1),XA(K2),XA(K3), &
     &        ZXDELLIM)
            ZYPXDL= FQQQ1(YA(K1),YA(K2),YA(K3),XA(K1),XA(K2),XA(K3), &
     &        ZXDELLIM)
          ELSE
!     COMPUTE VALUE AND 1ST DER. AT EDGE
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
!     CUBIC
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
              ZYXDL = FCCCC0(YA(K1),YA(K2),YA(K3),YA(K4),XA(K1),XA(K2), &
     &          XA(K3),XA(K4),ZXDELLIM)
              ZYPXDL= FCCCC1(YA(K1),YA(K2),YA(K3),YA(K4),XA(K1),XA(K2), &
     &          XA(K3),XA(K4),ZXDELLIM)
            ENDIF
          ELSE
!     COMPUTE VALUE AND 1ST DER. AT EDGE AND EDGE-1
            IF (KLO .EQ. 1) THEN
              ZXEDGM1 = XA(KHI)
              ZYEDGM1 = YA(KHI)
              ZYPEDGM1 = (YA(KHI)-YA(KLO))/H + (Y2A(KLO)+ZTWO*Y2A(KHI)) &
     &          *H/zsix
              ZXEDG = XA(KLO)
              ZYEDG = YA(KLO)
              ZYPEDG=(YA(KHI)-YA(KLO))/H-(ZTWO*Y2A(KLO)+Y2A(KHI))*H/zsix
            ELSE
!%OS              IF (H.LE.1.E-02_RKIND) THEN
!%OS                KLO=KLO-1
!%OS                H=XA(KHI)-XA(KLO)
!%OS              endif
              ZXEDGM1 = XA(KLO)
              ZYEDGM1 = YA(KLO)
              ZYPEDGM1 =(YA(KHI)-YA(KLO))/H-(ZTWO*Y2A(KLO)+Y2A(KHI)) &
     &          *H/zsix
              ZXEDG = XA(KHI)
              ZYEDG = YA(KHI)
              ZYPEDG=(YA(KHI)-YA(KLO))/H+(Y2A(KLO)+ZTWO*Y2A(KHI))*H/zsix
            ENDIF
            ZXDELLIM = XA(1) - DELLIM
            IF (KHI .EQ. N) ZXDELLIM = XA(N) + DELLIM
            ZYXDL = FCDCD0(ZXEDGM1,ZYEDGM1,ZYPEDGM1,ZXEDG,ZYEDG,ZYPEDG, &
     &        ZXDELLIM)
            ZYPXDL= FCDCD1(ZXEDGM1,ZYEDGM1,ZYPEDGM1,ZXEDG,ZYEDG,ZYPEDG, &
     &        ZXDELLIM)
          ENDIF
        ENDIF
!
        Y = FLINEAR(ZXDELLIM,ZYXDL,ZXDELLIM+ZONE,ZYXDL+ZYPXDL,X)
        YP = ZYPXDL
        YPP = 0.0_RKIND
!
!.......................................................................
!     2.3 QUADRATIC, IKOPT=2 OR FIRST PART OF IKOPT=21
!
      ELSE IF (IKOPT.EQ.2 .OR. &
     &    (MAX(XA(KLO)-X,X-XA(KHI)).LE.DELLIM .AND. IKOPT.EQ.21)) THEN
!
!     QUADRATIC
        IF (ICONTDER .EQ. 0) THEN
          K1 = 1
          IF (KHI .EQ. N) K1 = N-2
          K2 = K1 + 1
          K3 = K1 + 2
          Y =  FQQQ0(YA(K1),YA(K2),YA(K3),XA(K1),XA(K2),XA(K3),X)
          YP = FQQQ1(YA(K1),YA(K2),YA(K3),XA(K1),XA(K2),XA(K3),X)
          YPP =FQQQ2(YA(K1),YA(K2),YA(K3),XA(K1),XA(K2),XA(K3))
        ELSE
!     COMPUTE VALUE AND 1ST DER. AT EDGE
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
!
!.......................................................................
!     2.4 QUADRATIC, FAR END PART OF IKOPT=32
!
      ELSE IF (MAX(XA(KLO)-X,X-XA(KHI)).GT.DELLIM .AND. IKOPT.EQ.32)THEN
!
!     QUADRATIC FROM X+ALFA*DELTA_X
!     COMPUTE STARTING POINT AND DERIVATIVE FROM END OF CUBIC
!     WARNING: MUST BE COMPATIBLE WITH ALFA=0
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
            ZYXDL = FCCCC0(YA(K1),YA(K2),YA(K3),YA(K4),XA(K1),XA(K2), &
     &        XA(K3),XA(K4),ZXDELLIM)
            ZYPXDL= FCCCC1(YA(K1),YA(K2),YA(K3),YA(K4),XA(K1),XA(K2), &
     &        XA(K3),XA(K4),ZXDELLIM)
          ENDIF
        ELSE
!     COMPUTE VALUE AND 1ST DER. AT EDGE AND EDGE-1
          IF (KLO .EQ. 1) THEN
            ZXEDGM1 = XA(KHI)
            ZYEDGM1 = YA(KHI)
            ZYPEDGM1 = (YA(KHI)-YA(KLO))/H + (Y2A(KLO)+ZTWO*Y2A(KHI)) &
     &        *H/zsix
            ZXEDG = XA(KLO)
            ZYEDG = YA(KLO)
            ZYPEDG= (YA(KHI)-YA(KLO))/H-(ZTWO*Y2A(KLO)+Y2A(KHI))*H/zsix
          ELSE
            ZXEDGM1 = XA(KLO)
            ZYEDGM1 = YA(KLO)
            ZYPEDGM1 = (YA(KHI)-YA(KLO))/H - (ZTWO*Y2A(KLO)+Y2A(KHI)) &
     &        *H/zsix
            ZXEDG = XA(KHI)
            ZYEDG = YA(KHI)
            ZYPEDG= (YA(KHI)-YA(KLO))/H+(Y2A(KLO)+ZTWO*Y2A(KHI))*H/zsix
          ENDIF
          ZXDELLIM = XA(1) - DELLIM
          IF (KHI .EQ. N) ZXDELLIM = XA(N) + DELLIM
          ZYXDL = FCDCD0(ZXEDGM1,ZYEDGM1,ZYPEDGM1,ZXEDG,ZYEDG,ZYPEDG, &
     &      ZXDELLIM)
          ZYPXDL= FCDCD1(ZXEDGM1,ZYEDGM1,ZYPEDGM1,ZXEDG,ZYEDG,ZYPEDG, &
     &      ZXDELLIM)
        ENDIF
!     QUADRATIC FROM END OF INTERVAL, END+DELLIM AND 1ST DER. AT END+DELLIM
        KLOHI = 1
        IF (KHI .EQ. N) KLOHI = N
        IF (ABS(ALFA/H).GT.1E-04_RKIND) THEN
          Y  = FQDQ0(ZXDELLIM,ZYXDL,ZYPXDL,XA(KLOHI),YA(KLOHI),X)
          YP = FQDQ1(ZXDELLIM,ZYXDL,ZYPXDL,XA(KLOHI),YA(KLOHI),X)
          YPP =FQDQ2(ZXDELLIM,ZYXDL,ZYPXDL,XA(KLOHI),YA(KLOHI))
        ELSE
          Y  = FQDQ0(ZXDELLIM,ZYXDL,ZYPXDL,XA(KLO),YA(KLO),X)
          YP = FQDQ1(ZXDELLIM,ZYXDL,ZYPXDL,XA(KLO),YA(KLO),X)
          YPP =FQDQ2(ZXDELLIM,ZYXDL,ZYPXDL,XA(KLO),YA(KLO))
        ENDIF
!
!.......................................................................
!     2.5 CUBIC, IKOPT=3 OR FIRST PART OF IKOPT=31 AND IKOPT=32
!
      ELSE IF (IKOPT.EQ.3 .OR. &
     &    (MAX(XA(KLO)-X,X-XA(KHI)).LE.DELLIM .AND. &
     &    (IKOPT.EQ.32 .OR. IKOPT.EQ.31)) ) THEN
!
!     CUBIC
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
            YPP=0._RKIND
          ELSE
            Y   =FCCCC0(YA(K1),YA(K2),YA(K3),YA(K4),XA(K1),XA(K2),XA(K3) &
     &        ,XA(K4),X)
            YP  =FCCCC1(YA(K1),YA(K2),YA(K3),YA(K4),XA(K1),XA(K2),XA(K3) &
     &        ,XA(K4),X)
            YPP =FCCCC2(YA(K1),YA(K2),YA(K3),YA(K4),XA(K1),XA(K2),XA(K3) &
     &        ,XA(K4),X)
          ENDIF
        ELSE
!     COMPUTE VALUE AND 1ST DER. AT EDGE AND EDGE-1
          IF (KLO .EQ. 1) THEN
            ZXEDGM1 = XA(KHI)
            ZYEDGM1 = YA(KHI)
            ZYPEDGM1=(YA(KHI)-YA(KLO))/H+(Y2A(KLO)+ZTWO*Y2A(KHI))*H/zsix
            ZXEDG = XA(KLO)
            ZYEDG = YA(KLO)
            ZYPEDG= (YA(KHI)-YA(KLO))/H-(ZTWO*Y2A(KLO)+Y2A(KHI))*H/zsix
          ELSE
!%OS            IF (H.LE.1.E-02_RKIND) THEN
!%OS              KLO=KLO-1
!%OS              H=XA(KHI)-XA(KLO)
!%OS            endif
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
!
      RETURN
      END SUBROUTINE SPLIBND
