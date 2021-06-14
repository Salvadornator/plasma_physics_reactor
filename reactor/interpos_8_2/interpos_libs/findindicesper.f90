      SUBROUTINE FINDINDICESPER(PXIN,KNIN,PXOUT,KNOUT,IIN_XOUT,iper,ZXSHIFT,PERIOD)
      USE prec_rkind
      implicit none
      REAL(RKIND) :: PXIN(KNIN), PXOUT(KNOUT)
      REAL(RKIND) :: XA(KNIN+1), PERIOD
      INTEGER :: KNIN, KNOUT, iin_xout(KNOUT)
!
      REAL(RKIND) :: ZXSHIFT(KNOUT)
      INTEGER :: IPER(KNOUT), ISSHIFTED, i, j, iprevious
!
!   FIND MATCHING INTERVALS
!   assumes PXIN AND PXOUT MONOTONICAALY INCREASING
!-----------------------------------------------------------------------
!
!
      XA(1:KNIN) = PXIN
      XA(KNIN+1) = PXIN(1) + PERIOD
      IPREVIOUS = 1
      DO J=1,KNOUT
         !     SHIFT PXOUT POINTS WITHIN [PXIN(1),PXIN(KNIN)] WITH IPER * PERIOD
         IF (PXOUT(J) .LT. PXIN(1)) THEN
            ISSHIFTED = -1 
            IPER(J) = INT((PXIN(1)-PXOUT(J))/PERIOD) + 1
            ZXSHIFT(J) = PXOUT(J) + IPER(J)*PERIOD
         ELSE IF (PXOUT(J) .GT. PXIN(1)+PERIOD) THEN
            ISSHIFTED = 1
            IPER(J) = - (INT((PXOUT(J)-PXIN(1)-PERIOD)/PERIOD) + 1)
            ZXSHIFT(J) = PXOUT(J) + IPER(J)*PERIOD
         ELSE
            ISSHIFTED = 0
            IPER(J) = 0
            ZXSHIFT(J) = PXOUT(J)
         ENDIF
         ! FIND INTERVAL IN PXIN WHERE PXOUT+-IPER(J)*PERIOD BELONGS, START FROM PREVIOUS INTERVAL, EXCEPT IF THERE WAS A JUMP BACK
         IF (ZXSHIFT(J) .LT. PXIN(IPREVIOUS)) IPREVIOUS = 1
         DO I=IPREVIOUS,KNIN
           IF ((ZXSHIFT(J) .GE. XA(I)) .AND. (ZXSHIFT(J) .LE. XA(I+1))) THEN
             IIN_XOUT(J) = I
             IPREVIOUS = I
             EXIT
           END IF
         END DO
      END DO
    END SUBROUTINE findindicesper
