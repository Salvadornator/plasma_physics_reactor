!-----------------------------------------------------------------------
SUBROUTINE SPLIBNDA(PXIN,PYIN,PYINPP,KNIN,PXOUT,PY,PYP,PYPP,PYINT,KNOUT,KOPTXPOL,KOPTDER)
  USE prec_rkind
  implicit none
  REAL(RKIND) :: zone, ztwo, zthree, zfour, zsix
  PARAMETER(ZONE=1._RKIND, ZTWO=2._RKIND, ZTHREE=3._RKIND, ZFOUR=4._RKIND, ZSIX=6._RKIND)
  REAL(RKIND) :: PXIN(KNIN),PYIN(KNIN),PYINPP(KNIN)
  REAL(RKIND) :: PXOUT(KNOUT), PY(KNOUT), PYP(KNOUT), PYPP(KNOUT), PYINT(KNOUT)
  integer KNIN, KNOUT, KOPTXPOL, KOPTDER
  !
  !     KOPTDER = 0: COMPUTE ONLY FUNCTION PY POINTS PXOUT
  !     KOPTDER = 1: COMPUTE ALSO 1ST DERIVATIVE PYP
  !     KOPTDER = 2: COMPUTE ALSO 2ND DERIVATIVE PYPP
  !     KOPTDER = 3: COMPUTE ALSO INTEGRAL OF Y FROM XIN(1) TO XOUT(J) IN PYINT(J)
  !
  !   ABS(KOPTXPOL):
  !     KOPTXPOL = 0: STOP WITH ERROR MESSAGE IF OUT OF BOUND
  !     KOPTXPOL = 1: LINEAR EXTRAPOLATION: +1 FROM YEDGE,YPEDGE ; -1: FROM Y(LAST TWO POINTS)
  !     KOPTXPOL = 10: CONSTANT EXTRAPOLATION: +10: Y=YEDGE, -10: Y=0. OUTSIDE PXIN INTERVAL
  !     KOPTXPOL = 2: USE QUADRATIC EXTRAPOLATION IF X OUT OF BOUND
  !     KOPTXPOL = 3: USE CUBIC EXTRAPOLATION IF X OUT OF BOUND
  !     KOPTXPOL = 21: USE QUADRATIC WITHIN ALFA*DELTA_X AND LINEAR FURTHER
  !     KOPTXPOL = 31: USE CUBIC WITHIN ALFA*DELTA_X AND LINEAR    FURTHER
  !     KOPTXPOL = 32: USE CUBIC WITHIN ALFA*DELTA_X AND QUADRATIC FURTHER
  !
  !     KOPTXPOL > 0: VALUE AND 1ST DER. CONTINUOUS AT END OF INTERVAL, THUS
  !     .             USES CUBIC SPLINE OF LAST INTERVAL TO CONTINUE
  !     KOPTXPOL < 0: ONLY Y VALUE CONTINUOUS AND USES VALUES AT LAST BUT ONE,
  !     .             TWO, THREE POINTS TO EXTRAPOLATE (BETTER IF DER. AT EDGE
  !     .             IS WILD)
  !
  !-----------------------------------------------------------------------
  ! LOCAL VARIABLES:
  REAL(RKIND) :: ALFA
  PARAMETER(ALFA = 1._RKIND)
  REAL(RKIND) :: H, H2, A, B,  &
    & ZX1, ZY1, ZY1P, zx2, zy2,  &
    & ZXN, ZYN, ZYNP, ZXNN, ZYNN, &
    & ZXLFTDEL, ZYLFTDEL, ZYPLFTDEL, ZYINTLFTDL, &
    & ZXRGTDEL, ZYRGTDEL, ZYPRGTDEL, ZYINTRGTDL, &
    & ZYINT, ZYINT_XIN(KNIN), ZYIN_XPOL
  INTEGER ICONTDER
  INTEGER :: J1ST_XIN(KNIN),JLAST_XIN(KNIN), JLEFT(KNOUT), JRIGHT(KNOUT), &
    & JLEFT_DEL(KNOUT), JRIGHT_DEL(KNOUT), IOPTXPOL, I, J, K, KLO, KHI
  !
  ! VARIABLES RELATED TO FUNCTIONS:
  REAL(RKIND) :: FC3, X1, F1, P1, X2, PP1, PP2, HH, HPX1, &
       &  F2, P2, FC2, FC1, FC0, FQQQ0, FQQQ1, FQQQ2, &
       &  FLINEAR, FLINEARP, FCCCC0, FCCCC1, FCCCC2, FCCCC3, FQDQ0, FQDQ1, &
       &  FQDQ2, FCDCD0, FCDCD1, FCDCD2, FCDCD3, FB1, &
       &  FB2, FA2, FA3, FD2, FD1, &
       &  FCCCCM1, FCDCDM1, FQQQM1, FQDQM1, FLINEARM1, FLINXYP, FLINXYPM1, &
       &  FC0D, FC1C, FC2B, FC3A, &
       &  F2D3H, F2D2H, F2D1H, F2D0H,FC2DC2DH0,  FC2DC2DH1, FC2DC2DH2, FC2DC2DH3, FC2DC2DHM1, &
       &  FD0H, FD1H, FD2H, FD3H, FCDCDH0, FCDCDH1, FCDCDH2, FCDCDH3, FCDCDHM1, &
       &  FCH0, FCH1, FCH2, FCH3, FCCCCH0, FCCCCH1, X3, X4, F3, F4, H21, H31, H41, &
       &  FCCCCH2, FCCCCH3, FCCCCHM1, &
       &  FQH0, FQH1, FQH2, FQH3, FQQQH0, FQQQH1, FQQQH2, FQQQHM1, &
       &  FQDH0, FQDH1, FQDH2, FQDQH0, FQDQH1, FQDQH2, FQDQHM1, &
       &  FQ2DH0, FQ2DH1, FQ2DH2, FQ2DQH0, FQ2DQH1, FQ2DQH2, FQ2DQHM1
  REAL(RKIND) :: A1, A2, A3, A4, B1, B2, B3, B4, PX
  REAL(RKIND) :: FB0, FD0, FA0, FA1
  !
  ! REQUIRES VARIABLES ZTWO, ZTHREE, ZFOUR AND ZSIX
  ! (DEFINED AS 2._RKIND ETC TYPICALLY BUT IS USER DEPENDENT FOR PRECISION)
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
    &             ZSIX * FA3(A1,A2,A3,A4,B1,B2,B3,B4) * PX
  ! ----------------------------------------------------------------------
  ! -- FCCCC3 GIVES THE VALUE OF THE THIRD DERIVATIVE OF F(X) AT PX:     -
  ! -- FCCCC3(......,PX) = D3F/DX3 (PX)                                  -
  ! ----------------------------------------------------------------------
  FCCCC3(A1,A2,A3,A4,B1,B2,B3,B4,PX) = &
    &                      ZSIX * FA3(A1,A2,A3,A4,B1,B2,B3,B4)
  ! ----------------------------------------------------------------------
  ! -- FCCCCM1 GIVES THE VALUE OF THE INTEGRAL OF F(X) FROM B1 TO PX:     -
  FCCCCM1(A1,A2,A3,A4,B1,B2,B3,B4,PX) = &
    &  (PX-B1)*(FA0(A1,A2,A3,A4,B1,B2,B3,B4) + &
    &  (PX+B1)*FA1(A1,A2,A3,A4,B1,B2,B3,B4)/ZTWO + &
    &  FA2(A1,A2,A3,A4,B1,B2,B3,B4)/ZTHREE*(PX*(PX+B1)+B1*B1) + &
    &  (PX+B1)*(PX*PX+B1*B1)*FA3(A1,A2,A3,A4,B1,B2,B3,B4)/ZFOUR)
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
    &      (P1 + P2 - ZTWO*FC2(X1,F1,P1,X2,F2,P2)*(X1+X2) - &
    & ZTHREE*FC3(X1,F1,P1,X2,F2,P2)*(X1*X1+X2*X2))/ZTWO
  FC0(X1,F1,P1,X2,F2,P2) = &
    &      (F1 + F2 - FC1(X1,F1,P1,X2,F2,P2)*(X1+X2) - &
    & FC2(X1,F1,P1,X2,F2,P2)*(X1*X1+X2*X2) - &
    & FC3(X1,F1,P1,X2,F2,P2)*(X1**3+X2**3))/ZTWO
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
    &             ZSIX * FC3(X1,F1,P1,X2,F2,P2) * PX
  ! ----------------------------------------------------------------------
  ! -- FCDCD3 GIVES THE VALUE OF THE THIRD DERIVATIVE OF F(X) AT PX:    --
  ! -- FCDCD3(......,PX) = D3F/DX3 (PX)                                 --
  ! ----------------------------------------------------------------------
  FCDCD3(X1,F1,P1,X2,F2,P2,PX) = &
    &                      ZSIX * FC3(X1,F1,P1,X2,F2,P2)
  ! ----------------------------------------------------------------------
  ! -- FCDCDM1 GIVES THE VALUE OF THE INTEGRAL OF F(X) FROM X1 TO PX:
  FCDCDM1(X1,F1,P1,X2,F2,P2,PX) = &
    &  (PX-X1)*(FC0(X1,F1,P1,X2,F2,P2) + &
    &    (PX+X1) * FC1(X1,F1,P1,X2,F2,P2)/ZTWO + &
    &    (PX*(PX+X1)+X1*X1) * FC2(X1,F1,P1,X2,F2,P2)/ZTHREE + &
    &    (PX+X1)*(PX*PX+X1*X1) * FC3(X1,F1,P1,X2,F2,P2)/ZFOUR)
  !-----------------------------------------------------------------------
  !.......................................................................
  !*COMDECK CUC2DC2DH
  ! ----------------------------------------------------------------------
  ! --     STATEMENT FUNCTION FOR CUBIC INTERPOLATION                   --
  ! --                         10.10.2014  OS CRPP                      --
  ! --                                                                  --
  ! -- CUBIC INTERPOLATION OF A FUNCTION F(X-X(i))                           --
  ! -- THE SEVEN ARGUMENTS X1,F1,PP1,X2,F2,PP2,HH ARE DEFINED AS FOLLOWS: 
  ! -- F(X1) = F1 , F(X2) = F2 , D2F/DX2(X1) = PP1 , DF2/D2X(X2) = PP2, HH = X2-X1
  ! -- WRITTEN FROM F(X)=A (X-XI)^3 + B (X-XI)^2 + C (X-XI)^2 + D
  ! -- (BETTER WHEN X(I+1)-X(I) SMALL)
  ! ----------------------------------------------------------------------
  ! ADD "H" FOR "DELTA X" INSTEAD OF F(X)
  F2D3H(X1,F1,PP1,X2,F2,PP2,HH) = (PP2-PP1) / HH / ZSIX
  F2D2H(X1,F1,PP1,X2,F2,PP2,HH) = PP1 / ZTWO
  F2D1H(X1,F1,PP1,X2,F2,PP2,HH) = (F2-F1)/HH - HH/ZSIX*(ZTWO*PP1+PP2)
  F2D0H(X1,F1,PP1,X2,F2,PP2,HH) = F1
  ! ----------------------------------------------------------------------
  ! -- FC2DC2DH0 GIVES THE VALUE OF THE FUNCTION AT POINT PX with HPX1=PX-X1 --
  ! -- FC2DC2DH0(......,HPX1) = F(PX)                                        --
  ! ----------------------------------------------------------------------
  FC2DC2DH0(X1,F1,PP1,X2,F2,PP2,HH,HPX1) = &
       &           F2D0H(X1,F1,PP1,X2,F2,PP2,HH) + &
       &             HPX1 * (F2D1H(X1,F1,PP1,X2,F2,PP2,HH) + &
       &               HPX1 * (F2D2H(X1,F1,PP1,X2,F2,PP2,HH) + &
       &                 HPX1 * F2D3H(X1,F1,PP1,X2,F2,PP2,HH)))
  ! ----------------------------------------------------------------------
  ! -- FC2DC2DH1 GIVES THE VALUE OF THE 1st derivative AT POINT PX with HPX1=PX-X1
  ! -- FC2DC2DH1(......,HPX1) = dF(PX)/dx , similar for 2nd, 3RD derivative
  ! ----------------------------------------------------------------------
  FC2DC2DH1(X1,F1,PP1,X2,F2,PP2,HH,HPX1) = &
       &           (F2D1H(X1,F1,PP1,X2,F2,PP2,HH) + &
       &              HPX1 * (ZTWO*F2D2H(X1,F1,PP1,X2,F2,PP2,HH) + &
       &                HPX1 * ZTHREE*F2D3H(X1,F1,PP1,X2,F2,PP2,HH)))
  FC2DC2DH2(X1,F1,PP1,X2,F2,PP2,HH,HPX1) = &
       &           ZTWO * F2D2H(X1,F1,PP1,X2,F2,PP2,HH) + &
       &             ZSIX * HPX1 * F2D3H(X1,F1,PP1,X2,F2,PP2,HH)
  FC2DC2DH3(X1,F1,PP1,X2,F2,PP2,HH,HPX1) = &
       &             ZSIX * F2D3H(X1,F1,PP1,X2,F2,PP2,HH)
  ! Integral from X1 to PX
  FC2DC2DHM1(X1,F1,PP1,X2,F2,PP2,HH,HPX1) = &
       &           HPX1 * (F2D0H(X1,F1,PP1,X2,F2,PP2,HH) + &
       &             HPX1 * (F2D1H(X1,F1,PP1,X2,F2,PP2,HH)/ZTWO + &
       &               HPX1 * (F2D2H(X1,F1,PP1,X2,F2,PP2,HH)/ZTHREE + &
       &                 HPX1 * F2D3H(X1,F1,PP1,X2,F2,PP2,HH)/ZFOUR)))
  !-----------------------------------------------------------------------
  !.......................................................................
  !*COMDECK CUCDCDH
  ! ----------------------------------------------------------------------
  ! --     STATEMENT FUNCTION FOR CUBIC INTERPOLATION                   --
  ! --                         10.10.2014  OS CRPP                      --
  ! --                                                                  --
  ! -- CUBIC INTERPOLATION OF A FUNCTION F(X-X(i))                           --
  ! -- THE ARGUMENTS X1,F1,P1,X2,F2,P2,HH ARE DEFINED AS FOLLOWS: 
  ! -- F(X1) = F1 , F(X2) = F2 , DF/DX(X1) = P1 , DF/DX(X2) = P2, HH = X2-X1
  ! -- WRITTEN FROM F(X)=A (X-XI)^3 + B (X-XI)^2 + C (X-XI)^2 + D
  ! -- (BETTER WHEN X(I+1)-X(I) SMALL)
  ! ----------------------------------------------------------------------
  FD0H(X1,F1,P1,X2,F2,P2,HH) = F1
  FD1H(X1,F1,P1,X2,F2,P2,HH) = P1
  FD2H(X1,F1,P1,X2,F2,P2,HH) = (ZTHREE*(F2-F1) - HH*(P2+ZTWO*P1))/HH**2
  ! FD3H(X1,F1,P1,X2,F2,P2,HH) = (ZTWO*(F1-F2) + HH*(P1+P2))/HH**3
  FD3H(X1,F1,P1,X2,F2,P2,HH) = (P2-P1-ZTWO*FD2H(X1,F1,P1,X2,F2,P2,HH)*HH)/(HH*HH*ZTHREE)
  ! ----------------------------------------------------------------------
  ! -- FCDCDH0 GIVES THE VALUE OF THE FUNCTION AT POINT PX with HPX1=PX-X1 --
  ! -- FCDCDH0(......,HPX1) = F(PX)                                        --
  ! ----------------------------------------------------------------------
  FCDCDH0(X1,F1,P1,X2,F2,P2,HH,HPX1) = &
       &           FD0H(X1,F1,P1,X2,F2,P2,HH) + &
       &             HPX1 * (FD1H(X1,F1,P1,X2,F2,P2,HH) + &
       &               HPX1 * (FD2H(X1,F1,P1,X2,F2,P2,HH) + &
       &                 HPX1 * FD3H(X1,F1,P1,X2,F2,P2,HH)))
  ! ----------------------------------------------------------------------
  ! -- FCDCDH1 GIVES THE VALUE OF THE 1st derivative AT POINT PX with HPX1=PX-X1
  ! -- FCDCDH1(......,HPX1) = dF(PX)/dx , similar for 2nd, 3rd derivatives
  ! ----------------------------------------------------------------------
  FCDCDH1(X1,F1,P1,X2,F2,P2,HH,HPX1) = &
       &           (FD1H(X1,F1,P1,X2,F2,P2,HH) + &
       &              HPX1 * (ZTWO * FD2H(X1,F1,P1,X2,F2,P2,HH) + &
       &                HPX1 * ZTHREE * FD3H(X1,F1,P1,X2,F2,P2,HH)))
  FCDCDH2(X1,F1,P1,X2,F2,P2,HH,HPX1) = &
       &           ZTWO * FD2H(X1,F1,P1,X2,F2,P2,HH) + &
       &             ZSIX * HPX1 * FD3H(X1,F1,P1,X2,F2,P2,HH)
  FCDCDH3(X1,F1,P1,X2,F2,P2,HH,HPX1) = &
       &             ZSIX * FD3H(X1,F1,P1,X2,F2,P2,HH)
  ! Integral from X1 to PX
  FCDCDHM1(X1,F1,P1,X2,F2,P2,HH,HPX1) = &
       &           HPX1 * (FD0H(X1,F1,P1,X2,F2,P2,HH) + &
       &             HPX1 * (FD1H(X1,F1,P1,X2,F2,P2,HH)/ZTWO + &
       &               HPX1 * (FD2H(X1,F1,P1,X2,F2,P2,HH)/ZTHREE + &
       &                 HPX1 * FD3H(X1,F1,P1,X2,F2,P2,HH)/ZFOUR)))
  !-----------------------------------------------------------------------
  !.......................................................................
  !*COMDECK CUCCCCH
  ! ----------------------------------------------------------------------
  ! --     STATEMENT FUNCTION FOR CUBIC INTERPOLATION                   --
  ! --                         10.10.2014  OS CRPP                      --
  ! --                                                                  --
  ! -- CUBIC INTERPOLATION OF A FUNCTION F(PX)                           --
  ! -- THE ARGUMENTS X1,F1,X2,F2,X3,F3,X4,F4,H21,H31,H41 ARE DEFINED AS FOLLOWS: 
  ! -- F(X1) = F1 , F(X2) = F2 , F(X3) = F3 , F(X4) = F4, 
  ! -- H21 = X2-X1, H31 = X3-X1, H41 = X4-X1, HPX1=PX-X1
  ! -- Written from f(x)=a (x-x1)^3 + b (x-x1)^2 + c (x-x1)^2 + d
  ! -- (better when x(i+1)-x(i) small)
  ! ----------------------------------------------------------------------
  !
  FCH0(X1,F1,X2,F2,X3,F3,X4,F4,H21,H31,H41) = F1
  FCH1(X1,F1,X2,F2,X3,F3,X4,F4,H21,H31,H41) = &
       & (H41**2*(H41-H21)*(H31**3*(F2-F1)-H21**3*(F3-F1))-H31**2*(H31-H21)*(H41**3*(F2-F1)-H21**3*(F4-F1))) / &
       & (H41*(H41-H21)*(H31**2-H21**2)-H31*(H41**2-H21**2)*(H31-H21))/H21/H31/H41
!!$  FCH2(X1,F1,X2,F2,X3,F3,X4,F4,H21,H31,H41) = (H31**3*(F2-F1)-H21**3*(F3-F1)-(H31**3-H21**2*H31)* &
!!$       & FCH1(X1,F1,X2,F2,X3,F3,X4,F4,H21,H31,H41)*H21) / H21**2 / (H31-H21) / H31**2
  FCH2(X1,F1,X2,F2,X3,F3,X4,F4,H21,H31,H41) = &
       & (H31**3/H21*(F2-F1)-H21**2*(F3-F1) + H41**3/H21*(F2-F1)-H21**2*(F4-F1) &
       & - (H31*(H31**2-H21**2)+H41*(H41**2-H21**2))*FCH1(X1,F1,X2,F2,X3,F3,X4,F4,H21,H31,H41)) &
       & / H21 / (H31**2*(H31-H21)+H41**2*(H41-H21))
  FCH3(X1,F1,X2,F2,X3,F3,X4,F4,H21,H31,H41) = &
       & (F2+F3+F4-ZTHREE*F1-(H21**2+H31**2+H41**2)*FCH2(X1,F1,X2,F2,X3,F3,X4,F4,H21,H31,H41) &
       & - (H21+H31+H41)*FCH1(X1,F1,X2,F2,X3,F3,X4,F4,H21,H31,H41))/(H21**3+H31**3+H41**3)
  ! ----------------------------------------------------------------------
  ! -- FCCCCH0 GIVES THE VALUE OF THE FUNCTION AT POINT PX with HPX1=PX-X1 --
  ! -- FCCCCH0(......,HPX1) = F(PX)                                        --
  ! ----------------------------------------------------------------------
  FCCCCH0(X1,F1,X2,F2,X3,F3,X4,F4,H21,H31,H41,HPX1) = &
       &       FCH0(X1,F1,X2,F2,X3,F3,X4,F4,H21,H31,H41) + &
       &         HPX1 * (FCH1(X1,F1,X2,F2,X3,F3,X4,F4,H21,H31,H41) + &
       &           HPX1 * (FCH2(X1,F1,X2,F2,X3,F3,X4,F4,H21,H31,H41) + &
       &             HPX1 * FCH3(X1,F1,X2,F2,X3,F3,X4,F4,H21,H31,H41)))
  ! -- FCCCCH1(......,HPX1) = dF(PX)/dx
  FCCCCH1(X1,F1,X2,F2,X3,F3,X4,F4,H21,H31,H41,HPX1) = &
       &       FCH1(X1,F1,X2,F2,X3,F3,X4,F4,H21,H31,H41) + &
       &         HPX1 * (ZTWO*FCH2(X1,F1,X2,F2,X3,F3,X4,F4,H21,H31,H41) + &
       &           HPX1 * ZTHREE*FCH3(X1,F1,X2,F2,X3,F3,X4,F4,H21,H31,H41))
  ! -- FCCCCH2(......,HPX1) = d2F(PX)/dx^2
  FCCCCH2(X1,F1,X2,F2,X3,F3,X4,F4,H21,H31,H41,HPX1) = &
       &         ZTWO * FCH2(X1,F1,X2,F2,X3,F3,X4,F4,H21,H31,H41) + &
       &           ZSIX * HPX1 * FCH3(X1,F1,X2,F2,X3,F3,X4,F4,H21,H31,H41)
  ! -- FCCCCH3(......,HPX1) = d3F(PX)/dx^3
  FCCCCH3(X1,F1,X2,F2,X3,F3,X4,F4,H21,H31,H41,HPX1) = &
       &           ZSIX * FCH3(X1,F1,X2,F2,X3,F3,X4,F4,H21,H31,H41)
  ! -- FCCCCHM1(......,HPX1) = Integral of F(X) from X1 to PX
  FCCCCHM1(X1,F1,X2,F2,X3,F3,X4,F4,H21,H31,H41,HPX1) = &
       &       HPX1 * (FCH0(X1,F1,X2,F2,X3,F3,X4,F4,H21,H31,H41) + &
       &         HPX1 * (FCH1(X1,F1,X2,F2,X3,F3,X4,F4,H21,H31,H41)/ZTWO + &
       &           HPX1 * (FCH2(X1,F1,X2,F2,X3,F3,X4,F4,H21,H31,H41)/ZTHREE + &
       &             HPX1 * FCH3(X1,F1,X2,F2,X3,F3,X4,F4,H21,H31,H41)/ZFOUR)))
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
    &  (PX+B1)*FB1(A1,A2,A3,B1,B2,B3)/ZTWO + &
    &  FB2(A1,A2,A3,B1,B2,B3)/ZTHREE*(PX*(PX+B1)+B1*B1))
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
    &  (PX+X1)*FD1(X1,F1,P1,X2,F2)/ZTWO + &
    &  FD2(X1,F1,P1,X2,F2)/ZTHREE*(PX*(PX+X1)+X1*X1))
  !
  !.......................................................................
  !*COMDECK QUAQQQH
  ! ----------------------------------------------------------------------
  ! --     STATEMENT FUNCTION FOR QUADRATIC INTERPOLATION               --
  ! --                         19.01.87            AR        CRPP       --
  ! --                                                                  --
  ! -- QUADRATIC INTERPOLATION OF A FUNCTION F(X)                       --
  ! -- THE PARAMETERS X1,F1,X2,F2,X3,F3,X4,F4,H21,H31 ARE DEFINED AS FOLLOWS:     --
  ! -- F(X1) = F1 , F(X2) = F2 , F(X3) = F3 
  ! -- H21 = X2-X1, H31 = X3-X1, HPX1=PX-X1
  ! -- Written from f(x)=a (x-x1)^2 + b (x-x1) + c 
  ! ----------------------------------------------------------------------
  !
  FQH0(F1,F2,F3,H21,H31) = F1
  FQH1(F1,F2,F3,H21,H31) = &
       & ((F2-F1)*H31*H31 - (F3-F1)*H21*H21) / (H21*H31*(H31-H21))
  FQH2(F1,F2,F3,H21,H31) = &
       & (F3+F2-2*F1 - FQH1(F1,F2,F3,H21,H31)*(H21+H31)) &
       & / (H21*H21+H31*H31)
  ! -- FQQQH0(......,HPX1) = F(HPX1)
  FQQQH0(F1,F2,F3,H21,H31,HPX1) = &
       & FQH0(F1,F2,F3,H21,H31) + &
       &   HPX1 * (FQH1(F1,F2,F3,H21,H31) + &
       &     HPX1 * FQH2(F1,F2,F3,H21,H31))
  ! -- FQQQH1(......,HPX1) = DF/DX (HPX1)
  FQQQH1(F1,F2,F3,H21,H31,HPX1) = &
       &   FQH1(F1,F2,F3,H21,H31) + &
       &     ZTWO * HPX1 * FQH2(F1,F2,F3,H21,H31)
  ! -- FQQQH2(......,HPX1) = D2F/DX2 (HPX1)
  FQQQH2(F1,F2,F3,H21,H31,HPX1) = &
       & ZTWO * FQH2(F1,F2,F3,H21,H31)
  ! -- FQQQHM1 GIVES THE VALUE OF THE INTEGRAL OF F(X) FROM X1 TO PX:
  FQQQHM1(F1,F2,F3,H21,H31,HPX1) = &
       & HPX1 * (FQH0(F1,F2,F3,H21,H31) + &
       &   HPX1 * (FQH1(F1,F2,F3,H21,H31)/ZTWO + &
       &     HPX1 * FQH2(F1,F2,F3,H21,H31)/ZTHREE))
  !
  !.......................................................................
  !*COMDECK QUAQDQH
  ! ----------------------------------------------------------------------
  ! --     STATEMENT FUNCTION FOR QUADRATIC INTERPOLATION               --
  ! --                         19.01.87            AR        CRPP       --
  ! --                                                                  --
  ! -- QUADRATIC INTERPOLATION OF A FUNCTION F(X)                       --
  ! -- THE PARAMETERS X1,F1,X2,F2,X3,F3,X4,F4,H21,H31 ARE DEFINED AS FOLLOWS:     --
  ! -- F(X1) = F1 , DF(X1)/DX = P1 , F(X2) = F2 
  ! -- H21 = X2-X1, HPX1=PX-X1
  ! -- Written from f(x)=a (x-x1)^2 + b (x-x1) + c 
  ! ----------------------------------------------------------------------
  !
  FQDH0(F1,P1,F2,H21) = F1
  FQDH1(F1,P1,F2,H21) = P1
  FQDH2(F1,P1,F2,H21) = (F2 - P1*H21 - F1)/(H21*H21)
  ! -- FQDQH0(......,HPX1) = F(HPX1)
  FQDQH0(F1,P1,F2,H21,HPX1) = &
       & FQDH0(F1,P1,F2,H21) + &
       &   HPX1 * (FQDH1(F1,P1,F2,H21) + &
       &     HPX1 * FQDH2(F1,P1,F2,H21))
  ! -- FQQQH1(......,HPX1) = DF/DX (HPX1)
  FQDQH1(F1,P1,F2,H21,HPX1) = &
       &   FQDH1(F1,P1,F2,H21) + &
       &     ZTWO * HPX1 * FQDH2(F1,P1,F2,H21)
  ! -- FQDQH2(......,HPX1) = D2F/DX2 (HPX1)
  FQDQH2(F1,P1,F2,H21,HPX1) = ZTWO * FQDH2(F1,P1,F2,H21)
  ! -- FQQQHM1 GIVES THE VALUE OF THE INTEGRAL OF F(X) FROM X1 TO PX:
  FQDQHM1(F1,P1,F2,H21,HPX1) = &
       & HPX1 * (FQDH0(F1,P1,F2,H21) + &
       &   HPX1 * (FQDH1(F1,P1,F2,H21)/ZTWO + &
       &     HPX1 * FQDH2(F1,P1,F2,H21)/ZTHREE))
  !
  !.......................................................................
  !*COMDECK QUAQ2DQH
  ! ----------------------------------------------------------------------
  ! --     STATEMENT FUNCTION FOR QUADRATIC INTERPOLATION               --
  ! --                         19.01.87            AR        CRPP       --
  ! --                                                                  --
  ! -- QUADRATIC INTERPOLATION OF A FUNCTION F(X)                       --
  ! -- THE PARAMETERS X1,F1,X2,F2,X3,F3,X4,F4,H21,H31 ARE DEFINED AS FOLLOWS:     --
  ! -- F(X1) = F1 , D2F(X1)/DX^2 = PP1 , F(X2) = F2 
  ! -- H21 = X2-X1, HPX1=PX-X1
  ! -- Written from f(x)=a (x-x1)^2 + b (x-x1) + c 
  ! ----------------------------------------------------------------------
  !
  FQ2DH0(F1,PP1,F2,H21) = F1
  FQ2DH1(F1,PP1,F2,H21) = (F2-F1)/H21 - PP1*H21/ZTWO
  FQ2DH2(F1,PP1,F2,H21) = PP1 / ZTWO
  ! -- FQ2DQH0(......,HPX1) = F(HPX1)
  FQ2DQH0(F1,PP1,F2,H21,HPX1) = &
       & FQ2DH0(F1,PP1,F2,H21) + &
       &   HPX1 * (FQ2DH1(F1,PP1,F2,H21) + &
       &     HPX1 * FQ2DH2(F1,PP1,F2,H21))
  ! -- FQQQH1(......,HPX1) = DF/DX (HPX1)
  FQ2DQH1(F1,PP1,F2,H21,HPX1) = &
       &   FQ2DH1(F1,PP1,F2,H21) + &
       &     ZTWO * HPX1 * FQ2DH2(F1,PP1,F2,H21)
  ! -- FQ2DQH2(......,HPX1) = D2F/DX2 (HPX1)
  FQ2DQH2(F1,PP1,F2,H21,HPX1) = ZTWO * FQ2DH2(F1,PP1,F2,H21)
  ! -- FQQQHM1 GIVES THE VALUE OF THE INTEGRAL OF F(X) FROM X1 TO PX:
  FQ2DQHM1(F1,PP1,F2,H21,HPX1) = &
       & HPX1 * (FQ2DH0(F1,PP1,F2,H21) + &
       &   HPX1 * (FQ2DH1(F1,PP1,F2,H21)/ZTWO + &
       &     HPX1 * FQ2DH2(F1,PP1,F2,H21)/ZTHREE))
  !-----------------------------------------------------------------------
  !.......................................................................
  !     LINEAR
  !
  FLINEAR(X1,F1,X2,F2,PX) = F1 + (PX-X1)/(X2-X1) * (F2-F1)
  FLINEARP(X1,F1,X2,F2) = (F2-F1) / (X2-X1)
  FLINEARM1(X1,F1,X2,F2,PX) = (PX-X1)*(F1+(PX-X1)*(F2-F1)/(X2-X1)/ZTWO)
  FLINXYP(X1,F1,P1,PX) = P1*(PX-X1) + F1
  FLINXYPM1(X1,F1,P1,PX) = (PX-X1)*(F1+P1*(PX-X1)/ZTWO)

  !-----------------------------------------------------------------------
  ! 0. DEFAULTS
  !
  ICONTDER = 1
  IF (KOPTXPOL .LT. 0) ICONTDER = 0
  IOPTXPOL=ABS(KOPTXPOL)
  !
  ! COMPUTE INT. UP TO EACH INPUT INTERVAL
  IF (KOPTDER .GE. 3) THEN
    ZYINT_XIN(1) = 0._RKIND
    DO I=2,KNIN
      H = PXIN(I) - PXIN(I-1)
      ZYINT_XIN(I) = ZYINT_XIN(I-1) + H/2._RKIND*(PYIN(I) + PYIN(I-1) &
        & - H*H/12._RKIND*(PYINPP(I)+PYINPP(I-1)))
    END DO
  END IF
  !
  ! LOOP OVER XOUT POINTS WHICH CAN BE IN RANDOM ORDER
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
    H=PXIN(KHI)-PXIN(KLO)
    H2 = H * H
    SELECT CASE (KOPTDER)
    CASE (0)
      A=(PXIN(KHI)-PXOUT(J))/H
      B=(PXOUT(J)-PXIN(KLO))/H
      PY(J)=A*PYIN(KLO)+B*PYIN(KHI)+ ((A**3-A)*PYINPP(KLO)+(B**3-B)*PYINPP(KHI))*H2/6._RKIND
    CASE (1)
      A=(PXIN(KHI)-PXOUT(J))/H
      B=(PXOUT(J)-PXIN(KLO))/H
      PY(J)=A*PYIN(KLO)+B*PYIN(KHI)+ ((A**3-A)*PYINPP(KLO)+(B**3-B)*PYINPP(KHI))*H2/6._RKIND
      PYP(J)=(PYIN(KHI)-PYIN(KLO))/H - &
        & ((3._RKIND*A*A-1._RKIND)*PYINPP(KLO)-(3._RKIND*B*B-1._RKIND)*PYINPP(KHI) )*H/6._RKIND
    CASE (2)
      A=(PXIN(KHI)-PXOUT(J))/H
      B=(PXOUT(J)-PXIN(KLO))/H
      PY(J)=A*PYIN(KLO)+B*PYIN(KHI)+ ((A**3-A)*PYINPP(KLO)+(B**3-B)*PYINPP(KHI))*H2/6._RKIND
      PYP(J)=(PYIN(KHI)-PYIN(KLO))/H - &
        & ((3._RKIND*A*A-1._RKIND)*PYINPP(KLO)-(3._RKIND*B*B-1._RKIND)*PYINPP(KHI) )*H/6._RKIND
      PYPP(J)=A*PYINPP(KLO)+B*PYINPP(KHI)
    CASE (3)
      A=(PXIN(KHI)-PXOUT(J))/H
      B=ZONE - A
      A2 = A*A
      B2 = B*B
      PY(J)=A*PYIN(KLO)+B*PYIN(KHI)+ ((A**3-A)*PYINPP(KLO)+(B**3-B)*PYINPP(KHI))*H2/6._RKIND
      PYP(J)=(PYIN(KHI)-PYIN(KLO))/H - &
        & ((3._RKIND*A*A-1._RKIND)*PYINPP(KLO)-(3._RKIND*B*B-1._RKIND)*PYINPP(KHI) )*H/6._RKIND
      PYPP(J)=A*PYINPP(KLO)+B*PYINPP(KHI)
      PYINT(J) = ZYINT_XIN(KLO) + H/2._RKIND*(PYIN(KLO)*(1._RKIND-A2) + PYIN(KHI)*B2 &
        &  - (PYINPP(KLO)*(1._RKIND-A2)**2 + PYINPP(KHI)*B2*(2._RKIND-B2)) * H2/12._RKIND)
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
      PRINT *,' POINT PXOUT(',J,')=',PXOUT(J),' OUTSIDE INTERVAL [',PXIN(1),',',PXIN(KNIN),']'
      RETURN
      !        STOP 'IOPTXPOL=0'
    ENDIF
    !
    !     2.2 COMPUTE VALUES FOR POINTS ON THE LEFT OF PXIN(1)
    !           EXTRAPOLATION DEPENDS ON VALUE OF KOPTXPOL
    !
    IF (PXOUT(J) .LT. PXIN(1)) THEN
      H = PXIN(2) - PXIN(1)
      ! YP(PXIN(1)) FROM SPLINE
      ! ZY1 REFERS TO FIRST KNOWN POINTS AND ZY2 TO SECOND KNOWN POINT, 
      ! TYPICALLY AT PXIN(1) AND PXIN(2) OR PXIN(1)-ALFA*H AND PXIN(1)
      ZY1P = (PYIN(2)-PYIN(1))/H - (2._RKIND*PYINPP(1)+PYINPP(2))*H/6._RKIND
      ZXLFTDEL = PXIN(1) - ALFA*H
      !
      !   2.2.1 SPECIAL PART [XIN(1)-ALFA*H,XIN(1)] IF IOPTXPOL>20
      IF ((PXOUT(J) .GE. ZXLFTDEL) .AND. (IOPTXPOL.GE.21)) THEN
        IF (IOPTXPOL .EQ. 21) THEN
          ! QUADRATIC
          IF (ICONTDER .EQ. 1) THEN
            ! KOPTXPOL = +21
            PY(J) = FQDQ0(PXIN(1),PYIN(1),ZY1P,PXIN(2),PYIN(2),PXOUT(J))
            IF (KOPTDER .GE. 1) PYP(J) = FQDQ1(PXIN(1),PYIN(1),ZY1P,PXIN(2),PYIN(2),PXOUT(J))
            IF (KOPTDER .GE. 2) PYPP(J) = FQDQ2(PXIN(1),PYIN(1),ZY1P,PXIN(2),PYIN(2))
            ! INTEGRATES DIRECTLY FROM PXIN(1)
            IF (KOPTDER .GE. 3) PYINT(J) = FQDQM1(PXIN(1),PYIN(1),ZY1P,PXIN(2),PYIN(2),PXOUT(J))
          ELSE
            ! KOPTXPOL = -21
            PY(J) = FQQQ0(PYIN(1),PYIN(2),PYIN(3),PXIN(1),PXIN(2),PXIN(3),PXOUT(J))
            IF (KOPTDER .GE. 1) PYP(J) = FQQQ1(PYIN(1),PYIN(2),PYIN(3),PXIN(1),PXIN(2),PXIN(3),PXOUT(J))
            IF (KOPTDER .GE. 2) PYPP(J) = FQQQ2(PYIN(1),PYIN(2),PYIN(3),PXIN(1),PXIN(2),PXIN(3))
            IF (KOPTDER .GE. 3) PYINT(J) = FQQQM1(PYIN(1),PYIN(2),PYIN(3),PXIN(1),PXIN(2),PXIN(3),PXOUT(J))
          END IF
        ELSEIF (IOPTXPOL .GE. 31) THEN
          ! CUBIC PART OF 31, 32
          IF (ICONTDER .EQ. 1) THEN
            ! KOPTXPOL = +31 OR +32
            PY(J) = FC2DC2DH0(PXIN(1),PYIN(1),PYINPP(1),PXIN(2),PYIN(2),PYINPP(2),H,PXOUT(J)-PXIN(1))
            IF (KOPTDER .GE. 1) PYP(J) = FC2DC2DH1(PXIN(1),PYIN(1),PYINPP(1),PXIN(2),PYIN(2),PYINPP(2), &
                 & H,PXOUT(J)-PXIN(1))
            IF (KOPTDER .GE. 2) PYPP(J) = FC2DC2DH2(PXIN(1),PYIN(1),PYINPP(1),PXIN(2),PYIN(2),PYINPP(2), &
                 & H,PXOUT(J)-PXIN(1))
            IF (KOPTDER .GE. 3) PYINT(J) = FC2DC2DHM1(PXIN(1),PYIN(1),PYINPP(2),PXIN(2),PYIN(2),PYINPP(2), &
                 & H,PXOUT(J)-PXIN(1))
          ELSE
            ! KOPTXPOL = -31 OR -32
            PY(J) = FCCCC0(PYIN(1),PYIN(2),PYIN(3),PYIN(4),PXIN(1),PXIN(2),PXIN(3),PXIN(4),PXOUT(J))
            IF (KOPTDER .GE. 1) PYP(J) = FCCCC1(PYIN(1),PYIN(2),PYIN(3),PYIN(4),PXIN(1),PXIN(2),PXIN(3),PXIN(4),PXOUT(J))
            IF (KOPTDER .GE. 2) PYPP(J) = FCCCC2(PYIN(1),PYIN(2),PYIN(3),PYIN(4),PXIN(1),PXIN(2),PXIN(3),PXIN(4),PXOUT(J))
            IF (KOPTDER .GE. 3) PYINT(J) = FCCCCM1(PYIN(1),PYIN(2),PYIN(3),PYIN(4),PXIN(1),PXIN(2),PXIN(3),PXIN(4),PXOUT(J))
          END IF
        ELSE
          PRINT *,'OPTION  KOPTXPOL= ',KOPTXPOL,' NOT YET DEFINED'
          PRINT *,'KOPTXPOL 1'
          RETURN ! AVOID STOPS
        END IF
        !
      ELSE
        !
        ! 2.2.2 EXTRAPOLATION FAR LEFT: X<XIN(1)-ALFA*H OR X<XIN(1) IF NO ALFA PART CONSIDERED
        !
        IF ((IOPTXPOL .EQ. 1) .OR. (IOPTXPOL .EQ. 21) .OR. (IOPTXPOL .EQ. 31)) THEN
          ! LINEAR EXTRAPOLATION
          SELECT CASE (KOPTXPOL)
          CASE (1)
            ZX1 = PXIN(1)
            ZY1 = PYIN(1)
            ZY1P = ZY1P
            ZYINT = 0._RKIND
          CASE(-1)
            ZX1 = PXIN(1)
            ZY1 = PYIN(1)
            ZY1P = FLINEARP(PXIN(1),PYIN(1),PXIN(2),PYIN(2))
            ZYINT = 0._RKIND
          CASE (21)
            ZX1 = ZXLFTDEL
            ZY1 = FQDQ0(PXIN(1),PYIN(1),ZY1P,PXIN(2),PYIN(2),ZXLFTDEL)
            ZYINT = FQDQM1(PXIN(1),PYIN(1),ZY1P,PXIN(2),PYIN(2),ZXLFTDEL)
            ZY1P = FQDQ1(PXIN(1),PYIN(1),ZY1P,PXIN(2),PYIN(2),ZXLFTDEL)
          CASE (-21)
            ZX1 = ZXLFTDEL
            ZY1 = FQQQ0(PYIN(1),PYIN(2),PYIN(3),PXIN(1),PXIN(2),PXIN(3),ZXLFTDEL)
            ZY1P = FQQQ1(PYIN(1),PYIN(2),PYIN(3),PXIN(1),PXIN(2),PXIN(3),ZXLFTDEL)
            ZYINT = FQQQM1(PYIN(1),PYIN(2),PYIN(3),PXIN(1),PXIN(2),PXIN(3),ZXLFTDEL)
          CASE (31)
            ZYINT = FC2DC2DHM1(PXIN(1),PYIN(1),PYINPP(1),PXIN(2),PYIN(2),PYINPP(2),H,ZXLFTDEL-PXIN(1))
            ZX1 = ZXLFTDEL
            ZY1 = FC2DC2DH0(PXIN(1),PYIN(1),PYINPP(1),PXIN(2),PYIN(2),PYINPP(2),H,ZXLFTDEL-PXIN(1))
            ZY1P = FC2DC2DH1(PXIN(1),PYIN(1),PYINPP(1),PXIN(2),PYIN(2),PYINPP(2),H,ZXLFTDEL-PXIN(1))
          CASE (-31)
            ZX1 = ZXLFTDEL
            ZY1 = FCCCC0(PYIN(1),PYIN(2),PYIN(3),PYIN(4),PXIN(1),PXIN(2),PXIN(3),PXIN(4),ZXLFTDEL)
            ZY1P = FCCCC1(PYIN(1),PYIN(2),PYIN(3),PYIN(4),PXIN(1),PXIN(2),PXIN(3),PXIN(4),ZXLFTDEL)
            ZYINT = FCCCCM1(PYIN(1),PYIN(2),PYIN(3),PYIN(4),PXIN(1),PXIN(2),PXIN(3),PXIN(4),ZXLFTDEL)
          END SELECT
          PY(J) = FLINXYP(ZX1,ZY1,ZY1P,PXOUT(J))
          IF (KOPTDER .GE. 1) PYP(J) = ZY1P
          IF (KOPTDER .GE. 2) PYPP(J) = 0._RKIND
          IF (KOPTDER .GE. 3) PYINT(J) = ZYINT + FLINXYPM1(ZX1,ZY1,ZY1P,PXOUT(J))
        ELSE IF (IOPTXPOL .EQ. 10) THEN
          ! CONSTANT OUTSIDE PXIN
          IF (ICONTDER .EQ. 1) THEN
            ! KOPTXPOL = +10
            ZYIN_XPOL = PYIN(1)
          ELSE
            ! KOPTXPOL = -10
            ZYIN_XPOL = 0._RKIND
          END IF
          PY(J) = ZYIN_XPOL
          IF (KOPTDER .GE. 1) PYP(J) = 0._RKIND
          IF (KOPTDER .GE. 2) PYPP(J) = 0._RKIND
          ! INTEGRATES FROM PXIN(1) DIRECTLY
          IF (KOPTDER .GE. 3) PYINT(J) = (PXOUT(J)-PXIN(1))*ZYIN_XPOL
          !
        ELSE IF ((IOPTXPOL .EQ. 2) .OR. (IOPTXPOL .EQ. 32)) THEN
          ! QUADRATIC EXTRAPOLATION
          SELECT CASE (KOPTXPOL)
          CASE (2)
            ZX1 = PXIN(1)
            ZY1 = PYIN(1)
            ZX2 = PXIN(2)
            ZY2 = PYIN(2)
            ZYINT = 0._RKIND
          CASE (-2)
            ZX1 = PXIN(1)
            ZY1 = PYIN(1)
            ZY1P = FQQQ1(PYIN(1),PYIN(2),PYIN(3),PXIN(1),PXIN(2),PXIN(3),PXIN(1))
            ZX2 = PXIN(2)
            ZY2 = PYIN(2)
            ZYINT = 0._RKIND
          CASE (32)
            ZYINT = FC2DC2DHM1(PXIN(1),PYIN(1),PYINPP(1),PXIN(2),PYIN(2),PYINPP(2),H,ZXLFTDEL-PXIN(1))
            ZX1 = ZXLFTDEL
            ZY1 = FC2DC2DH0(PXIN(1),PYIN(1),PYINPP(1),PXIN(2),PYIN(2),PYINPP(2),H,ZXLFTDEL-PXIN(1))
            ZY1P = FC2DC2DH1(PXIN(1),PYIN(1),PYINPP(1),PXIN(2),PYIN(2),PYINPP(2),H,ZXLFTDEL-PXIN(1))
            ZX2 = PXIN(1)
            ZY2 = PYIN(1)
          CASE (-32)
            ZX1 = ZXLFTDEL
            ZY1 = FCCCC0(PYIN(1),PYIN(2),PYIN(3),PYIN(4),PXIN(1),PXIN(2),PXIN(3),PXIN(4),ZXLFTDEL)
            ZY1P = FCCCC1(PYIN(1),PYIN(2),PYIN(3),PYIN(4),PXIN(1),PXIN(2),PXIN(3),PXIN(4),ZXLFTDEL)
            ZX2 = PXIN(1)
            ZY2 = PYIN(1)
            ZYINT = FCCCCM1(PYIN(1),PYIN(2),PYIN(3),PYIN(4),PXIN(1),PXIN(2),PXIN(3),PXIN(4),ZXLFTDEL)
          END SELECT
          PY(J) = FQDQ0(ZX1,ZY1,ZY1P,ZX2,ZY2,PXOUT(J))
          IF (KOPTDER .GE. 1) PYP(J) = FQDQ1(ZX1,ZY1,ZY1P,ZX2,ZY2,PXOUT(J))
          IF (KOPTDER .GE. 2) PYPP(J) = FQDQ2(ZX1,ZY1,ZY1P,ZX2,ZY2)
          IF (KOPTDER .GE. 3) PYINT(J) = ZYINT + FQDQM1(ZX1,ZY1,ZY1P,ZX2,ZY2,PXOUT(J))
        ELSE IF (IOPTXPOL .EQ. 3) THEN
          ! CUBIC EXTRAPOLATION
          SELECT CASE (KOPTXPOL)
          CASE (3)
            PY(J) = FC2DC2DH0(PXIN(1),PYIN(1),PYINPP(1),PXIN(2),PYIN(2),PYINPP(2),H,PXOUT(J)-PXIN(1))
            IF (KOPTDER .GE. 1) PYP(J) = &
                 & FC2DC2DH1(PXIN(1),PYIN(1),PYINPP(1),PXIN(2),PYIN(2),PYINPP(2),H,PXOUT(J)-PXIN(1));
            IF (KOPTDER .GE. 2) PYPP(J) = &
                 & FC2DC2DH2(PXIN(1),PYIN(1),PYINPP(1),PXIN(2),PYIN(2),PYINPP(2),H,PXOUT(J)-PXIN(1))
            IF (KOPTDER .GE. 3) PYINT(J) = &
                 & FC2DC2DHM1(PXIN(1),PYIN(1),PYINPP(1),PXIN(2),PYIN(2),PYINPP(2),H,PXOUT(J)-PXIN(1))
          CASE (-3)
            PY(J) = FCCCC0(PYIN(1),PYIN(2),PYIN(3),PYIN(4),PXIN(1),PXIN(2),PXIN(3),PXIN(4),PXOUT(J))
            IF (KOPTDER .GE. 1) PYP(J) = FCCCC1(PYIN(1),PYIN(2),PYIN(3),PYIN(4),PXIN(1),PXIN(2),PXIN(3),PXIN(4),PXOUT(J))
            IF (KOPTDER .GE. 2) PYPP(J) = FCCCC2(PYIN(1),PYIN(2),PYIN(3),PYIN(4),PXIN(1),PXIN(2),PXIN(3),PXIN(4),PXOUT(J))
            IF (KOPTDER .GE. 3) PYINT(J) = FCCCCM1(PYIN(1),PYIN(2),PYIN(3),PYIN(4),PXIN(1),PXIN(2),PXIN(3),PXIN(4),PXOUT(J))
          END SELECT
        ELSE
          PRINT *,'OPTION  IOPTXPOL= ',IOPTXPOL,' NOT YET DEFINED'
          STOP 'KOPTXPOL 2'
        END IF
      END IF
      !
    ELSE ! PXOUT(J) .GT. PXIN(KNIN)
      !
      !     2.3 COMPUTE VALUES FOR POINTS ON THE RIGHT OF PXIN(1)
      !         EXTRAPOLATION DEPENDS ON VALUE OF KOPTXPOL
      !
      H = PXIN(KNIN) - PXIN(KNIN-1)
      ! ZYN REFERS TO FIRST KNOWN POINTS AND ZYNN TO SECOND KNOWN POINT, 
      ! TYPICALLY AT PXIN(KNIN) AND PXIN(KNIN-1) OR PXIN(KNIN)+ALFA*H AND PXIN(KNIN)
      ! YP(PXIN(KNIN)) FROM SPLINE
      ZYNP = (PYIN(KNIN)-PYIN(KNIN-1))/H + (PYINPP(KNIN-1)+2._RKIND*PYINPP(KNIN))*H/6._RKIND
      ZXRGTDEL = PXIN(KNIN) + ALFA*H
      !
      !   2.3.1 SPECIAL PART ]XIN(KNIN),XIN(KNIN)+ALFA*H] IF IOPTXPOL>20
      IF ((PXOUT(J) .LE. ZXRGTDEL) .AND. (IOPTXPOL.GE.21)) THEN
        ZYINT = ZYINT_XIN(KNIN)
        IF (IOPTXPOL .EQ. 21) THEN
          ! QUADRATIC
          IF (ICONTDER .EQ. 1) THEN
            ! KOPTXPOL = +21
            PY(J) = FQDQ0(PXIN(KNIN),PYIN(KNIN),ZYNP,PXIN(KNIN-1),PYIN(KNIN-1),PXOUT(J))
            IF (KOPTDER .GE. 1) PYP(J) = FQDQ1(PXIN(KNIN),PYIN(KNIN),ZYNP,PXIN(KNIN-1),PYIN(KNIN-1),PXOUT(J))
            IF (KOPTDER .GE. 2) PYPP(J) = FQDQ2(PXIN(KNIN),PYIN(KNIN),ZYNP,PXIN(KNIN-1),PYIN(KNIN-1))
            ! INTEGRATES FROM PXIN(KNIN) SO ADD INTEGRAL UP TO PXIN(KNIN)
            IF (KOPTDER .GE. 3) PYINT(J) = ZYINT + FQDQM1(PXIN(KNIN),PYIN(KNIN),ZYNP,PXIN(KNIN-1),PYIN(KNIN-1),PXOUT(J))
          ELSE
            ! KOPTXPOL = -21
            PY(J) = FQQQ0(PYIN(KNIN),PYIN(KNIN-1),PYIN(KNIN-2),PXIN(KNIN),PXIN(KNIN-1),PXIN(KNIN-2),PXOUT(J))
            IF (KOPTDER .GE. 1) PYP(J) = &
              & FQQQ1(PYIN(KNIN),PYIN(KNIN-1),PYIN(KNIN-2),PXIN(KNIN),PXIN(KNIN-1),PXIN(KNIN-2),PXOUT(J))
            IF (KOPTDER .GE. 2) PYPP(J) = &
              & FQQQ2(PYIN(KNIN),PYIN(KNIN-1),PYIN(KNIN-2),PXIN(KNIN),PXIN(KNIN-1),PXIN(KNIN-2))
            IF (KOPTDER .GE. 3) PYINT(J) = ZYINT + &
              & FQQQM1(PYIN(KNIN),PYIN(KNIN-1),PYIN(KNIN-2),PXIN(KNIN),PXIN(KNIN-1),PXIN(KNIN-2),PXOUT(J))
          END IF
        ELSEIF (IOPTXPOL .GE. 31) THEN
          ! CUBIC PART OF 31, 32
          IF (ICONTDER .EQ. 1) THEN
            ! KOPTXPOL = +31 OR +32
            PY(J) = FC2DC2DH0(PXIN(KNIN),PYIN(KNIN),PYINPP(KNIN), &
                 &      PXIN(KNIN-1),PYIN(KNIN-1),PYINPP(KNIN-1),-H,PXOUT(J)-PXIN(KNIN))
            IF (KOPTDER .GE. 1) PYP(J) = FC2DC2DH1(PXIN(KNIN),PYIN(KNIN),PYINPP(KNIN), &
                 & PXIN(KNIN-1),PYIN(KNIN-1),PYINPP(KNIN-1),-H,PXOUT(J)-PXIN(KNIN))
            IF (KOPTDER .GE. 2) PYPP(J) = FC2DC2DH2(PXIN(KNIN),PYIN(KNIN),PYINPP(KNIN), &
                 & PXIN(KNIN-1),PYIN(KNIN-1),PYINPP(KNIN-1),-H,PXOUT(J)-PXIN(KNIN))
            IF (KOPTDER .GE. 3) PYINT(J) = ZYINT + FC2DC2DHM1(PXIN(KNIN),PYIN(KNIN),PYINPP(KNIN), &
                 & PXIN(KNIN-1),PYIN(KNIN-1),PYINPP(KNIN-1),-H,PXOUT(J)-PXIN(KNIN))
          ELSE
            ! KOPTXPOL = -31 OR -32
            PY(J) = FCCCC0(PYIN(KNIN),PYIN(KNIN-1),PYIN(KNIN-2),PYIN(KNIN-3), &
              & PXIN(KNIN),PXIN(KNIN-1),PXIN(KNIN-2),PXIN(KNIN-3),PXOUT(J))
            IF (KOPTDER .GE. 1) PYP(J) = FCCCC1(PYIN(KNIN),PYIN(KNIN-1),PYIN(KNIN-2),PYIN(KNIN-3), &
              & PXIN(KNIN),PXIN(KNIN-1),PXIN(KNIN-2),PXIN(KNIN-3),PXOUT(J))
            IF (KOPTDER .GE. 2) PYPP(J) = FCCCC2(PYIN(KNIN),PYIN(KNIN-1),PYIN(KNIN-2),PYIN(KNIN-3), &
              & PXIN(KNIN),PXIN(KNIN-1),PXIN(KNIN-2),PXIN(KNIN-3),PXOUT(J))
            IF (KOPTDER .GE. 3) PYINT(J) = ZYINT + FCCCCM1(PYIN(KNIN),PYIN(KNIN-1),PYIN(KNIN-2),PYIN(KNIN-3), &
              & PXIN(KNIN),PXIN(KNIN-1),PXIN(KNIN-2),PXIN(KNIN-3),PXOUT(J))
          END IF
        ELSE
          PRINT *,'OPTION  KOPTXPOL= ',KOPTXPOL,' NOT YET DEFINED'
          STOP 'KOPTXPOL 3'
        END IF
        !
      ELSE
        !
        ! 2.3.2 EXTRAPOLATION FAR RIGHT: X>XIN(KNIN)+ALFA*H OR X>XIN(KNIN) IF NO ALFA PART CONSIDERED
        !
        IF ((IOPTXPOL .EQ. 1) .OR. (IOPTXPOL .EQ. 21) .OR. (IOPTXPOL .EQ. 31)) THEN
          ! LINEAR EXTRAPOLATION
          SELECT CASE (KOPTXPOL)
          CASE (1)
            ZXN = PXIN(KNIN)
            ZYN = PYIN(KNIN)
            ZYINT = ZYINT_XIN(KNIN)
          CASE(-1)
            ZXN = PXIN(KNIN)
            ZYN = PYIN(KNIN)
            ZYNP = FLINEARP(PXIN(KNIN),PYIN(KNIN),PXIN(KNIN-1),PYIN(KNIN-1))
            ZYINT = ZYINT_XIN(KNIN)
          CASE (21)
            ZYINT = ZYINT_XIN(KNIN) + FQDQM1(PXIN(KNIN),PYIN(KNIN),ZYNP,PXIN(KNIN-1),PYIN(KNIN-1),ZXRGTDEL)
            ZXN = ZXRGTDEL
            ZYN = FQDQ0(PXIN(KNIN),PYIN(KNIN),ZYNP,PXIN(KNIN-1),PYIN(KNIN-1),ZXRGTDEL)
            ZYNP = FQDQ1(PXIN(KNIN),PYIN(KNIN),ZYNP,PXIN(KNIN-1),PYIN(KNIN-1),ZXRGTDEL)
          CASE (-21)
            ZXN = ZXRGTDEL
            ZYN = FQQQ0(PYIN(KNIN),PYIN(KNIN-1),PYIN(KNIN-2),PXIN(KNIN),PXIN(KNIN-1),PXIN(KNIN-2),ZXRGTDEL)
            ZYNP = FQQQ1(PYIN(KNIN),PYIN(KNIN-1),PYIN(KNIN-2),PXIN(KNIN),PXIN(KNIN-1),PXIN(KNIN-2),ZXRGTDEL)
            ZYINT = ZYINT_XIN(KNIN) + &
              & FQQQM1(PYIN(KNIN),PYIN(KNIN-1),PYIN(KNIN-2),PXIN(KNIN),PXIN(KNIN-1),PXIN(KNIN-2),ZXRGTDEL)
          CASE (31)
            ZXN = ZXRGTDEL
            ZYINT = ZYINT_XIN(KNIN) + FC2DC2DHM1(PXIN(KNIN),PYIN(KNIN),PYINPP(KNIN), &
                 & PXIN(KNIN-1),PYIN(KNIN-1),PYINPP(KNIN-1),-H,ZXRGTDEL-PXIN(KNIN))
            ZYN = FC2DC2DH0(PXIN(KNIN),PYIN(KNIN),PYINPP(KNIN), &
                 & PXIN(KNIN-1),PYIN(KNIN-1),PYINPP(KNIN-1),-H,ZXRGTDEL-PXIN(KNIN))
            ZYNP = FC2DC2DH1(PXIN(KNIN),PYIN(KNIN),PYINPP(KNIN), &
                 & PXIN(KNIN-1),PYIN(KNIN-1),PYINPP(KNIN-1),-H,ZXRGTDEL-PXIN(KNIN))
          CASE (-31)
            ZXN = ZXRGTDEL
            ZYN = FCCCC0(PYIN(KNIN),PYIN(KNIN-1),PYIN(KNIN-2),PYIN(KNIN-3), &
              & PXIN(KNIN),PXIN(KNIN-1),PXIN(KNIN-2),PXIN(KNIN-3),ZXRGTDEL)
            ZYNP = FCCCC1(PYIN(KNIN),PYIN(KNIN-1),PYIN(KNIN-2),PYIN(KNIN-3), &
              & PXIN(KNIN),PXIN(KNIN-1),PXIN(KNIN-2),PXIN(KNIN-3),ZXRGTDEL)
            ZYINT = ZYINT_XIN(KNIN) + FCCCCM1(PYIN(KNIN),PYIN(KNIN-1),PYIN(KNIN-2),PYIN(KNIN-3), &
              & PXIN(KNIN),PXIN(KNIN-1),PXIN(KNIN-2),PXIN(KNIN-3),ZXRGTDEL)
          END SELECT
          PY(J) = FLINXYP(ZXN,ZYN,ZYNP,PXOUT(J))
          IF (KOPTDER .GE. 1) PYP(J) = ZYNP
          IF (KOPTDER .GE. 2) PYPP(J) = 0._RKIND
          IF (KOPTDER .GE. 3) PYINT(J) = ZYINT + FLINXYPM1(ZXN,ZYN,ZYNP,PXOUT(J))
        ELSE IF (IOPTXPOL .EQ. 10) THEN
          ! CONSTANT OUTSIDE PXIN
          IF (ICONTDER .EQ. 1) THEN
            ! KOPTXPOL = +10
            ZYIN_XPOL = PYIN(KNIN)
          ELSE
            ! KOPTXPOL = -10
            ZYIN_XPOL = 0._RKIND
          END IF
          PY(J) = ZYIN_XPOL
          IF (KOPTDER .GE. 1) PYP(J) = 0._RKIND
          IF (KOPTDER .GE. 2) PYPP(J) = 0._RKIND
          IF (KOPTDER .GE. 3) PYINT(J) = ZYINT_XIN(KNIN) + (PXOUT(J)-PXIN(KNIN))*ZYIN_XPOL
          !
        ELSE IF ((IOPTXPOL .EQ. 2) .OR. (IOPTXPOL .EQ. 32)) THEN
          ! QUADRATIC EXTRAPOLATION
          SELECT CASE (KOPTXPOL)
          CASE (2)
            ZXN = PXIN(KNIN)
            ZYN = PYIN(KNIN)
            ZXNN = PXIN(KNIN-1)
            ZYNN = PYIN(KNIN-1)
            ZYINT = ZYINT_XIN(KNIN)
          CASE (-2)
            ZXN = PXIN(KNIN)
            ZYN = PYIN(KNIN)
            ZYNP = FQQQ1(PYIN(KNIN),PYIN(KNIN-1),PYIN(KNIN-2),PXIN(KNIN),PXIN(KNIN-1),PXIN(KNIN-2),PXIN(KNIN))
            ZXNN = PXIN(KNIN-1)
            ZYNN = PYIN(KNIN-1)
            ZYINT = ZYINT_XIN(KNIN)
          CASE (32)
            ZXN = ZXRGTDEL
            ZYINT = ZYINT_XIN(KNIN) + FC2DC2DHM1(PXIN(KNIN),PYIN(KNIN),PYINPP(KNIN), &
                 & PXIN(KNIN-1),PYIN(KNIN-1),PYINPP(KNIN-1),-H,ZXRGTDEL-PXIN(KNIN))
            ZXNN = PXIN(KNIN)
            ZYNN = PYIN(KNIN)
            ZYN = FC2DC2DH0(PXIN(KNIN),PYIN(KNIN),PYINPP(KNIN), &
                 & PXIN(KNIN-1),PYIN(KNIN-1),PYINPP(KNIN-1),-H,ZXRGTDEL-PXIN(KNIN))
            ZYNP = FC2DC2DH1(PXIN(KNIN),PYIN(KNIN),PYINPP(KNIN), &
                 & PXIN(KNIN-1),PYIN(KNIN-1),PYINPP(KNIN-1),-H,ZXRGTDEL-PXIN(KNIN))
          CASE (-32)
            ZXN = ZXRGTDEL
            ZYN = FCCCC0(PYIN(KNIN),PYIN(KNIN-1),PYIN(KNIN-2),PYIN(KNIN-3), &
              & PXIN(KNIN),PXIN(KNIN-1),PXIN(KNIN-2),PXIN(KNIN-3),ZXRGTDEL)
            ZYNP = FCCCC1(PYIN(KNIN),PYIN(KNIN-1),PYIN(KNIN-2),PYIN(KNIN-3), &
              & PXIN(KNIN),PXIN(KNIN-1),PXIN(KNIN-2),PXIN(KNIN-3),ZXRGTDEL)
            ZXNN = PXIN(KNIN)
            ZYNN = PYIN(KNIN)
            ZYINT = ZYINT_XIN(KNIN) + FCCCCM1(PYIN(KNIN),PYIN(KNIN-1),PYIN(KNIN-2),PYIN(KNIN-3), &
              & PXIN(KNIN),PXIN(KNIN-1),PXIN(KNIN-2),PXIN(KNIN-3),ZXRGTDEL)
          END SELECT
          PY(J) = FQDQ0(ZXN,ZYN,ZYNP,ZXNN,ZYNN,PXOUT(J))
          IF (KOPTDER .GE. 1) PYP(J) = FQDQ1(ZXN,ZYN,ZYNP,ZXNN,ZYNN,PXOUT(J))
          IF (KOPTDER .GE. 2) PYPP(J) = FQDQ2(ZXN,ZYN,ZYNP,ZXNN,ZYNN)
          IF (KOPTDER .GE. 3) PYINT(J) = ZYINT + FQDQM1(ZXN,ZYN,ZYNP,ZXNN,ZYNN,PXOUT(J))
        ELSE IF (IOPTXPOL .EQ. 3) THEN
          ! CUBIC EXTRAPOLATION
          SELECT CASE (KOPTXPOL)
          CASE (3)
            PY(J) = FC2DC2DH0(PXIN(KNIN),PYIN(KNIN),PYINPP(KNIN), &
                 & PXIN(KNIN-1),PYIN(KNIN-1),PYINPP(KNIN-1),-H,PXOUT(J)-PXIN(KNIN))
            IF (KOPTDER .GE. 1) PYP(J) = FC2DC2DH1(PXIN(KNIN),PYIN(KNIN),PYINPP(KNIN), &
                 & PXIN(KNIN-1),PYIN(KNIN-1),PYINPP(KNIN-1),-H,PXOUT(J)-PXIN(KNIN))
            IF (KOPTDER .GE. 2) PYPP(J) = FC2DC2DH2(PXIN(KNIN),PYIN(KNIN),PYINPP(KNIN), &
                 & PXIN(KNIN-1),PYIN(KNIN-1),PYINPP(KNIN-1),-H,PXOUT(J)-PXIN(KNIN))
            IF (KOPTDER .GE. 3) PYINT(J) = ZYINT_XIN(KNIN) + FC2DC2DHM1(PXIN(KNIN),PYIN(KNIN),PYINPP(KNIN), &
              & PXIN(KNIN-1),PYIN(KNIN-1),PYINPP(KNIN-1),-H,PXOUT(J)-PXIN(KNIN))
          CASE (-3)
            PY(J) = FCCCC0(PYIN(KNIN),PYIN(KNIN-1),PYIN(KNIN-2),PYIN(KNIN-3), &
              & PXIN(KNIN),PXIN(KNIN-1),PXIN(KNIN-2),PXIN(KNIN-3),PXOUT(J))
            IF (KOPTDER .GE. 1) PYP(J) = FCCCC1(PYIN(KNIN),PYIN(KNIN-1),PYIN(KNIN-2),PYIN(KNIN-3), &
              & PXIN(KNIN),PXIN(KNIN-1),PXIN(KNIN-2),PXIN(KNIN-3),PXOUT(J))
            IF (KOPTDER .GE. 2) PYPP(J) = FCCCC2(PYIN(KNIN),PYIN(KNIN-1),PYIN(KNIN-2),PYIN(KNIN-3), &
              & PXIN(KNIN),PXIN(KNIN-1),PXIN(KNIN-2),PXIN(KNIN-3),PXOUT(J))
            IF (KOPTDER .GE. 3) PYINT(J) = ZYINT_XIN(KNIN) + FCCCCM1(PYIN(KNIN),PYIN(KNIN-1),PYIN(KNIN-2),PYIN(KNIN-3), &
              & PXIN(KNIN),PXIN(KNIN-1),PXIN(KNIN-2),PXIN(KNIN-3),PXOUT(J))
          END SELECT
        ELSE
          PRINT *,'OPTION  KOPTXPOL= ',KOPTXPOL,' NOT YET DEFINED'
          STOP 'KOPTXPOL 4'
        END IF
      END IF
    END IF
    !
100 END DO
  RETURN
END SUBROUTINE SPLIBNDA
