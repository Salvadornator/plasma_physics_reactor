SUBROUTINE CBSPLGNP(PXIN,PYIN,PYINNEW,PYINPP,KNIN,PXOUT,PYOUT, &
  &  PYOUTP,PYOUTPP,PYOUTINT,KNOUT,KOPT,PSIG,PERIOD,PXEXP0,PXEXPDL,KFLAG)
  !     =================================================================
  !
  !   ASSUME PERIODIC BOUNDARY CONDITION WITH Y(PXIN(1)+PERIOD) = PYIN(1)
  !   BUT PYIN(KNIN).NE. PYIN(1), SO DON-T INCLUDE EXTRA PERIODIC POINT
  !
  !     NOTE: THIS ROUTINE INCLUDES THE STANDARD CUBIC SPLINE IF PSIG=0:
  !           THEN PYINNEW IS NOT USED AND MDAMAT=3 IS SUFFICIENT
  !           (=> PYINNEW(1) OR PYINNEW=PYIN IS OK)
  !
  !     Interpolate (pxin,pyin) on (pxout,pyout) using
  !     Hirshman fitted cubic spline with ptaus value or
  !     standard cubic spline if PSIG=0
  !
  !     KOPT = 0: ONLY INTERPOLATE FUNCTION INTO PYOUT
  !     KOPT = 1: INTERPOLATE FUNCTION INTO PYOUT AND 1ST DER. INTO PYOUTP
  !     KOPT = 2: AS KOPT=1 PLUS 2ND DER. INTO PYOUTPP
  !     KOPT = 3: AS KOPT=2 PLUS PRIMITIVE(Y SPLINE) INTO PYOUTINT
  !
  !     SEE COMMENTS FOR ROUTINE CBFITBND FOR MORE INFORMATION
  !
  !     IF LAPACK ROUTINES NOT AVAILABLE, USE spgbtrf_s.f
  !     (LAPACK SOURCE COPIED FROM NETLIB.ORG)
  !
  !     PXIN(KNIN)    : INPUT ABSCISSA (GIVEN DATA)
  !     PYIN(KNIN)    : INPUT VALUES OF FUNCTION AT PXIN(I),I=1,KNIN
  !     PYINNEW(KNIN)   : IF PSIG.NE.0, THEN PYINNEW CONTAINS ON OUTPUT THE NEW VALUES
  !     .         OF THE FUNCTION AT PXIN(I) FOR THE CUBIC SPLINE FIT
  !     PYINPP(KNIN+1)  : SECOND DER. OF THE CUBIC SPLINE FIT FOR (PXIN,PYIN) IF PSIG=0 OR
  !     .         ON (PXIN,PYINNEW) OTHERWISE
  !     KNIN    : NUMBER OF INPUT POINTS
  !     PXOUT(KNOUT)   : X VALUES AT WHICH THE FUNCTION HAS TO BE INTERPOLATED (INPUT)
  !     PYOUT(KNOUT)   : INTERPOLATED VALUES AT PXOUT(I),I=1,KNOUT (OUTPUT)
  !     PYOUTP(KNOUT)  : INTERPOLATED VALUES OF 1ST DER. OF FUNCTIONS AT PXOUT(I) (OUTPUT, IF KOPT.GE.1)
  !     PYOUTPP(KNOUT) : INTERPOLATED VALUES OF 2ND DER. OF FUNCTIONS AT PXOUT(I) (OUTPUT, IF KOPT.GE.2)
  !     PYOUTINT(KNOUT): INTEGRAL OF YOUT AT PXOUT(I) (OUTPUT, IF KOPT.EQ.3)
  !     KNOUT   : NUMBER OF POINTS FOR OUTPUT
  !     KOPT    : SEE ABOVE
  !     PSIG   : WEIGHT OF SECOND DERIVATIVE IN THE CHI**2 TO BE MINIMIZED. PSIG=0 GIVES THE
  !     .         STANDARD CUBIC SPLINE. LARGER VALUES OF PSIG WILL SMOOTH MORE THE 2ND DER.
  !     PXEXP0  : PSIG IS WEIGHTED BY AN EXP(-((X-PXEXP0)/PXEXPDL)**2)
  !     PXEXPDL : IF PXEXP0 NOT IN [PXIN(1),PXIN(1)+PERIOD], EXP() IGNORED AND PSIG=CST
  !     .         (SEE ROUTINE CBFITBND BELOW)
  !     KFLAG   : ERROR FLAG: IF NOT 0, THERE IS A PROBLEM
  !
  !-----------------------------------------------------------------------
  USE prec_rkind
  implicit none
  REAL(RKIND) :: EPSILON
  PARAMETER(EPSILON = 1.0E-10)
  !
  REAL(RKIND) :: PXIN(KNIN), PYIN(KNIN) &
       &  ,PXOUT(KNOUT), PSIG(KNIN)
  REAL(RKIND) :: PYINNEW(KNIN), PYINPP(KNIN) &
       &  ,PYOUT(KNOUT), PYOUTP(KNOUT), PYOUTPP(KNOUT), PYOUTINT(KNOUT)
  !
  REAL(RKIND) :: PERIOD, PXEXP0, PXEXPDL
  INTEGER :: KOPT, KNIN, KNOUT
  INTEGER :: KFLAG
  !
  REAL(RKIND) :: ZY, ZYP, ZYPP, ZYINT, ZDXMIN, ZTAUS
  INTEGER :: I
  !-----------------------------------------------------------------------
  !
  !     0. CHECK INPUT CONSISTENCY
  !
  KFLAG=1
  ZTAUS = PSIG(1)
!!$  IF (ZTAUS .LT. 0) THEN
!!$    ZDXMIN = MINVAL(PXIN(2:KNIN)-PXIN(1:KNIN-1))
!!$    ZTAUS = ABS(ZTAUS) * ZDXMIN**3
!!$    PRINT *,'% TAUS CHANGED TO DEFAULT VALUE = ',ZTAUS
!!$  END IF
  !
  !   PXIN in ASCENDING ORDER
  !
  DO i=1,KNIN-1
    if (PXIN(i) .GE. PXIN(i+1)) then
      print *,' xin not in ascending order:'
      print *,' xin(',i,')= ',PXIN(i),'   >=   xin(',i+1,')= ', &
           &      PXIN(i+1)
      KFLAG = 2
      RETURN
    endif
  END DO
  !
  CALL CBFITPER(PXIN,PYIN,PYINNEW,KNIN,PYINPP,PSIG,PERIOD,PXEXP0,PXEXPDL)
  !
  !L    2. COMPUTE INTERPOLATED VALUE AT EACH PXOUT
  !
  KFLAG=2
  IF (ZTAUS .EQ. 0.0) THEN
    CALL SPLIPERA(PXIN,PYIN   ,PYINPP,KNIN,PXOUT,KNOUT,PYOUT,PYOUTP,PYOUTPP,PYOUTINT, &
      &      PERIOD,KOPT)
  ELSE
    CALL SPLIPERA(PXIN,PYINNEW,PYINPP,KNIN,PXOUT,KNOUT,PYOUT,PYOUTP,PYOUTPP,PYOUTINT, &
      &      PERIOD,KOPT)
  ENDIF
  !
  KFLAG=0
  RETURN
END SUBROUTINE CBSPLGNP
