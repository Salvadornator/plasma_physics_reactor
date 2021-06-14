MODULE prec_rkind
  !
  !   Precision for real and complex
  !
  INTEGER, PARAMETER :: RKIND = SELECTED_REAL_KIND(15,300)
  INTEGER, PARAMETER :: CKIND = RKIND
  INTEGER, PARAMETER :: ITM_I4 = SELECTED_INT_KIND (9)        ! Integer*4
  INTEGER, PARAMETER :: ITM_I8 = SELECTED_INT_KIND (18)       ! Integer*8
  !
END MODULE prec_rkind
