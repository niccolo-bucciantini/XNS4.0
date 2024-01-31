MODULE ROTATION
  !   ============================================================
  !   Purpose : This module contains subroutines for a differential rotator
  !   ============================================================
  
  USE SYSTEMXNS, ONLY: DIFFERENTIAL,OMGSPACE
  USE SYSTEMXNS, ONLY: JCONSTLAW,JCMODLAW,URYULAW3,URYULAW4
  USE SYSTEMXNS, ONLY: PROTDIFF,A2VALUE,RMVALUE
  USE SYSTEMXNS, ONLY: OMG,OMGMAX,NR,NTH
  USE SYSTEMXNS, ONLY: PSI,PSL,PSS,R,TH
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: CHECKROTDIFF,OMEGAVALUE,OMEGA3LVALUE
  
  REAL :: JMAX,XVAL,YVAL
  
CONTAINS
  
! ********************************************************
! ********************************************************
  
SUBROUTINE CHECKROTDIFF()
  !  ============================================================
  !  Purpose : check for consistency
  !  ============================================================
  LOGICAL,DIMENSION(:),ALLOCATABLE :: MASK_OS,MASK_JS
  CHARACTER(len=2000), ALLOCATABLE :: STR_OS(:),STR_JS(:),VAR(:)
  INTEGER,PARAMETER :: N_OS=2
  INTEGER,PARAMETER :: N_JS=2
  INTEGER :: TOT_OS,TOT_JS,TOT
  REAL :: CHECK_LIMIT,MVAL
  
  !=== Check for logical parameters ===
  ALLOCATE(MASK_OS(N_OS),MASK_JS(N_JS))
  ALLOCATE(STR_OS(N_OS),STR_JS(N_JS),VAR(1))
  
  MASK_OS = (/JCONSTLAW,JCMODLAW/)
  STR_OS = [CHARACTER(len=2000) :: 'JCONSTLAW','JCMODLAW']
  
  MASK_JS = (/URYULAW3,URYULAW4/)
  STR_JS = [CHARACTER(len=2000) :: 'URYULAW3','URYULAW4']

  TOT_OS = COUNT(MASK_OS)
  TOT_JS = COUNT(MASK_JS)
  TOT = TOT_OS+TOT_JS

  IF(TOT.NE.1) THEN
    WRITE(6,*)''
    WRITE(6,*)'Logical parameters for the rotational profile are conflicting.'
    WRITE(6,*)'One and only one must be set True for differential rotators.'
    STOP
  ENDIF
  IF(OMGSPACE .AND. TOT_JS.NE.0) THEN
    WRITE(6,*)''
    WRITE(6,*)'Logical parameters for the rotational profile are conflicting.'
    !VAR=STR_JS(FINDLOC(MASK_JS,.true.))
    WRITE(6,*)'If -OMEGASPACE- is .true., -',TRIM(VAR(1)),'- cannot be .true.'
    STOP
  ENDIF
  IF(.NOT.OMGSPACE .AND. TOT_OS.NE.0) THEN
    WRITE(6,*)''
    WRITE(6,*)'Logical parameters for the rotational profile are conflicting.'
    !VAR=STR_OS(FINDLOC(MASK_OS,.true.))
    WRITE(6,*)'If -OMEGASPACE- is .true., -',TRIM(VAR(1)),'- cannot be .true.'
    STOP
  ENDIF

  DEALLOCATE(MASK_OS,MASK_JS,STR_OS,STR_JS,VAR)

  !=== Check for numerical parameters ===
  !Omega-space
  IF(OMGSPACE) THEN
    IF(A2VALUE.LT.0) THEN
      WRITE(6,*)''
      WRITE(6,*)'WARNING: if -OMGSPACE- is True, -A2VALUE- cannot be less than zero.'
      STOP
    ENDIF
    IF(JCMODLAW .AND. PROTDIFF.LT.1) THEN
      WRITE(6,*)''
      WRITE(6,*)'WARNING: if -JCMODLAW- is True, -PROTDIFF- cannot be less than one.'
      STOP
    ENDIF
  ENDIF

  !j-space
  IF(.NOT.OMGSPACE) THEN
    IF(OMGMAX.LE.OMG .OR. RMVALUE.LE.0) THEN
      WRITE(6,*)''
      WRITE(6,*)'WARNING: if -OMGSPACE- is False, -OMGMAX- cannot be less than or equal to OMG'
      WRITE(6,*)'and -RMVALUE- cannot be less than or equal to zero.'
      STOP
    ENDIF
    MVAL=OMGMAX/(OMG + 1.D-10)
    CHECK_LIMIT=4.*PROTDIFF/(PROTDIFF+1.)**2
    IF(URYULAW3 .AND. (PROTDIFF.LT.1 .OR. MVAL.LE.CHECK_LIMIT)) THEN
      WRITE(6,*)''
      WRITE(6,*)'WARNING: if -URYULAW3- is True, -PROTDIFF- cannot be less than one'
      WRITE(6,*)'and OMGMAX/OMG must be greater than 4*PROTDIFF/(PROTDIFF+1)**2.'
      STOP
    ENDIF
  ENDIF

END SUBROUTINE CHECKROTDIFF

! ********************************************************
! ********************************************************

SUBROUTINE OMEGAVALUE(BETALOC,RSTARLOC,OMEGALOC)
  !   ============================================================
  !   Purpose : evaluate OMEGALOC
  !   ============================================================
  REAL :: RSTARLOC,BETALOC,OMEGALOC
  REAL :: FO,DFO,DOMEGALOC

  REAL,PARAMETER :: EPS1=1.E-16,TOL=1.E-12
  INTEGER :: ITER,ITER_MAX=1000

  OMEGALOC=OMG

  IF(DIFFERENTIAL)THEN
    !IF(URYULAW3 .OR. URYULAW4) THEN
    IF(.NOT.OMGSPACE) THEN
      CALL PARS_VALUE_JS
      OMEGALOC=(OMG+OMGMAX)/2.
    ELSE
      OMEGALOC=OMG/2.
    ENDIF

    DO ITER=1,ITER_MAX
      IF (OMGSPACE) THEN
        CALL FODFO_OS(BETALOC,RSTARLOC,OMEGALOC,FO,DFO)
      ELSE
        CALL FODFO_JS(BETALOC,RSTARLOC,OMEGALOC,FO,DFO)
      ENDIF

      DOMEGALOC = -FO/DFO

      IF (ABS(DOMEGALOC)<TOL*OMEGALOC) EXIT

      OMEGALOC=OMEGALOC+DOMEGALOC
      OMEGALOC=MAX(OMEGALOC,0.)
    END DO
  END IF

END SUBROUTINE OMEGAVALUE

! ********************************************************
! ********************************************************

SUBROUTINE OMEGA3LVALUE(BETALOC,RSTARLOC,OMEGALOC,A3L)
  !   ============================================================
  !   Purpose : evaluate OMEGALOC and A3L
  !   ============================================================
  REAL :: RSTARLOC,BETALOC
  REAL :: OMEGALOC,A3L

  !Evaluate OMEGALOC
  CALL OMEGAVALUE(BETALOC,RSTARLOC,OMEGALOC)

  !Evaluate A3L
  IF(DIFFERENTIAL) THEN
    IF (OMGSPACE) THEN
      CALL A3L_OS(OMEGALOC,A3L)
    ELSE
      CALL A3L_JS(BETALOC,RSTARLOC,OMEGALOC,A3L)
    ENDIF
  ELSE
    A3L=0.
  ENDIF

END SUBROUTINE OMEGA3LVALUE

! ********************************************************
! ********************************************************
! OMEGA-SPACE (_OS)
! ********************************************************
! ********************************************************

SUBROUTINE FODFO_OS(BETALOC,RSTARLOC,OMEGALOC,FO,DFO)
  !   ============================================================
  !   Purpose : evaluate FO and DFO in the Omega-space for a differential rotator
  !   ============================================================
  REAL :: RSTARLOC,BETALOC,OMEGALOC
  REAL :: FO,DFO
  REAL :: C1,C2,C3,C4,C5,C6,C7,C8

  !j-constant law
  IF(JCONSTLAW) THEN
    C1=A2VALUE*RSTARLOC
    C2=A2VALUE*RSTARLOC*(2*BETALOC-OMG)
    C3=A2VALUE*(BETALOC**2*RSTARLOC-1.-2*BETALOC*RSTARLOC*OMG)-RSTARLOC
    C4=A2VALUE*(OMG-BETALOC**2*OMG*RSTARLOC)-BETALOC*RSTARLOC

    FO=C1*OMEGALOC**3+C2*OMEGALOC**2+C3*OMEGALOC+C4
    DFO=3.*C1*OMEGALOC**2+2.*C2*OMEGALOC+C3
  ENDIF

  !Modified j-constant law
  IF(JCMODLAW) THEN
    C1=A2VALUE*RSTARLOC
    C2=2.*A2VALUE*RSTARLOC*BETALOC
    C3=-A2VALUE*RSTARLOC*OMG**PROTDIFF
    C4=A2VALUE*(BETALOC**2*RSTARLOC-1.)-RSTARLOC
    C5=-2*A2VALUE*BETALOC*RSTARLOC*OMG**PROTDIFF
    C6=A2VALUE*(OMG**PROTDIFF)*(1.-RSTARLOC*BETALOC**2)
    C7=-BETALOC*RSTARLOC

    FO=C1*OMEGALOC**3+C2*OMEGALOC**2+C3*OMEGALOC**(3.-PROTDIFF)+C4*OMEGALOC+&
      C5*OMEGALOC**(2.-PROTDIFF)+C6*OMEGALOC**(1.-PROTDIFF)+C7
    DFO=3.*C1*OMEGALOC**2+2.*C2*OMEGALOC+(3.-PROTDIFF)*C3*OMEGALOC**(2.-PROTDIFF)+&
      C4+(2.-PROTDIFF)*C5*OMEGALOC**(1.-PROTDIFF)+(1.-PROTDIFF)*C6/(OMEGALOC**PROTDIFF)
  ENDIF

END SUBROUTINE FODFO_OS

! ********************************************************
! ********************************************************

SUBROUTINE A3L_OS(OMEGALOC,A3L)
  !   ============================================================
  !   Purpose : evaluate A3L in the Omega-space for a differential rotator
  !   ============================================================
  REAL :: OMEGALOC,A3L
  REAL :: C1,C2,C3,C4

  !J-constant law
  IF(JCONSTLAW) THEN
    A3L=0.5*A2VALUE*(OMEGALOC-OMG)**2
  ENDIF

  !Modified j-constant law
  IF(JCMODLAW) THEN
    C1=(2.-PROTDIFF)*OMEGALOC**2
    C2=PROTDIFF*OMG**2
    C3=-2.*OMEGALOC**(2.-PROTDIFF)*OMG**PROTDIFF
    C4=2.*(2.-PROTDIFF)

    A3L=A2VALUE*(C1+C2+C3)/C4
  ENDIF

END SUBROUTINE A3L_OS


! ********************************************************
! ********************************************************
! J-SPACE (_JS)
! ********************************************************
! ********************************************************

SUBROUTINE PARS_VALUE_JS()
  !   ============================================================
  !   Purpose : evaluate JMAX, XVAL and YVAL see  Franceschetti et al 2022
  !   ============================================================
  INTEGER :: IM,IT=NTH/2+1
  REAL    :: BETAMAX,RSMAX,MVAL

  !Value of JMAX
  IM=MINLOC(ABS(R(1:NR)-RMVALUE),DIM=1)
  RSMAX=(R(IM)*SIN(TH(IT))*PSI(IT,IM)**3/PSL(IT,IM))**2
  BETAMAX=PSS(IT,IM)
  JMAX=RSMAX*(OMGMAX+BETAMAX)/(1.-RSMAX*(OMGMAX+BETAMAX)**2)

  !Value of MVAL
  MVAL=OMGMAX/(OMG+1.D-10)

  !Value of XVAL and YVAL (A8, A9, A18 in Franceschetti et al 2022)
  IF(URYULAW3) THEN
    !Uryu law with 3 parameters
    XVAL=(MVAL*(PROTDIFF+1.))**2-4.*PROTDIFF*MVAL !delta
    XVAL=1.+(SQRT(XVAL)-MVAL*(PROTDIFF+1.))/(2.*PROTDIFF)
    YVAL=-1.+MVAL/(1.-XVAL)
  ENDIF

  IF(URYULAW4) THEN
    !Uryu law with 4 parameters (p=1,q=3)
    XVAL=((MVAL-1.)/(3.*MVAL))**(1./4.)
    YVAL=4.*(MVAL-1.)/3.
  ENDIF

END SUBROUTINE PARS_VALUE_JS

! ********************************************************
! ********************************************************

SUBROUTINE FODFO_JS(BETALOC,RSTARLOC,OMEGALOC,FO,DFO)
  !   ============================================================
  !   Purpose : evaluate FO and DFO in the J-space for a differential rotator see  Franceschetti et al 2022
  !   ============================================================
  REAL :: RSTARLOC,BETALOC,OMEGALOC
  REAL :: FO,DFO
  REAL :: JLOC,DJLOC,CJ1,CJ2,CJ3
  REAL :: C1,C2,C3,C4,C5

  !Values of JLOC and DJLOC
  CJ1=OMEGALOC+BETALOC
  CJ2=RSTARLOC*CJ1
  CJ3=1.-RSTARLOC*CJ1**2
  JLOC=CJ2/CJ3
  DJLOC=(RSTARLOC+CJ2**2)/CJ3**2

  !Uryu law with 3 parameters (A2 in Franceschetti et al 2022)
  IF(URYULAW3) THEN
    C1=JLOC/JMAX
    C2=1.+YVAL*C1**PROTDIFF
    C3=1.-XVAL*C1
    C4=PROTDIFF*YVAL*C1**(PROTDIFF-1.)

    FO=OMEGALOC-OMG*C2*C3
    DFO=1.-(C4*C3-XVAL*C2)*OMG*DJLOC/JMAX
  ENDIF

  !Uryu law with 4 parameters (p=1,q=3) (A3 in Franceschetti et al 2022)
  IF(URYULAW4) THEN
    C1=JLOC/JMAX
    C2=1.+YVAL*C1
    C3=1.+(XVAL*C1)**4
    C4=4.*OMEGALOC*(XVAL**4)*(C1**3)
    C5=OMG*YVAL

    FO=OMEGALOC*C3-OMG*C2
    DFO=C3+(C4-C5)*DJLOC/JMAX
  ENDIF

END SUBROUTINE FODFO_JS

! ********************************************************
! ********************************************************

SUBROUTINE A3L_JS(BETALOC,RSTARLOC,OMEGALOC,A3L)
  !   ============================================================
  !   Purpose : evaluate A3L in the J-space for a differential rotator
  !   ============================================================
  REAL :: RSTARLOC,BETALOC,OMEGALOC,A3L
  REAL :: JLOC,CJ1,CJ2,CJ3
  REAL :: C1,C2,C3,C4,C5,C6

  !Value of JLOC
  CJ1=OMEGALOC+BETALOC
  CJ2=RSTARLOC*CJ1
  CJ3=1.-RSTARLOC*CJ1**2
  JLOC=CJ2/CJ3

  !Uryu law with 3 parameters
  IF(URYULAW3) THEN
    C1=JLOC/JMAX
    C2=YVAL*C1**PROTDIFF
    C3=XVAL*C1
    C4=PROTDIFF*C2/(PROTDIFF+1.)
    C5=(PROTDIFF+1.)*C2*C3/(PROTDIFF+2.)

    A3L=(C3/2.-C4+C5)*JLOC*OMG
  ENDIF

  !Uryu law with 4 parameters (p=1,q=3)
  IF(URYULAW4) THEN
    C1=OMG*JMAX/(4.*XVAL)
    C2=XVAL*JLOC/JMAX
    C3=2.*YVAL*ATAN(C2**2)/XVAL
    C4=ATAN(1.+SQRT(2.)*C2)
    C5=ATAN(1.-SQRT(2.)*C2)
    C6=ATANH(SQRT(2.)*C2/(1.+C2**2))

    A3L=-JLOC*OMEGALOC+C1*(C3+SQRT(2.)*(C4-C5+C6))
  ENDIF

END SUBROUTINE A3L_JS

! ********************************************************
! ********************************************************
END MODULE ROTATION
