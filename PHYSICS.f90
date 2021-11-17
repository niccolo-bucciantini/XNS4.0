MODULE PHYSICS
  !  ============================================================
  !  Purpose : storage module for grid-building related functions
  !  and for EoS related functions, and other auxiliary functions.
  !  ============================================================  
  
  USE SYSTEMXNS
  IMPLICIT NONE
  
  REAL :: RMAXINV,RMIDINV
  
CONTAINS

! **************************************************************
! **************************************************************

SUBROUTINE GRIDBUILD(R,DR)
  ! === Builds the grid and the quantities for the DGTSV solve ===
  ! === Options for uniform radial grid and regular+stretched ===
  IMPLICIT NONE
  REAL :: DOM  ! Radial domain of regular grid
  REAL :: ASTRREQ1,ASTRREQ2,LHS,RHS,ASTRR,ASTRR1,EX  ! For finding the stretching parameters
  
  REAL,DIMENSION(0:NR+1) :: R,DR
  INTEGER :: I,N
  INTEGER,PARAMETER :: MAXITER=100
  
  IF(.NOT. STRETCH) THEN
     DOM=RMAX-RMIN
     DO I=1,NR
        R(I)=DOM*(I-0.5)/NR
     END DO
     R(0)=-R(1)
     R(NR+1)=DOM*(NR+0.5)/NR
     DO I=1,NR
        !DR(I)=R(I+1)-R(I)   !!! CHANGE THIS
        DR(I)=R(I)-R(I-1)
     END DO
     !DR(NR+1)=DR(NR)
     DR(NR+1)=R(NR+1)-R(NR)
     DR(0)=DR(1)
     DO I=1,NR
        DRM(I)=R(I)-R(I-1)
        DRP(I)=R(I+1)-R(I)
     END DO
     DRM(0)=DRM(1)
     DRP(0)=R(1)-R(0)
     DRM(NR+1)=R(NR+1)-R(NR)
     DRP(NR+1)=DRP(NR)
     
     IF((DR(1) .GT. MINRESREG) .AND. VERBOSE)THEN
        WRITE(6,*)'WARNING: LOW RADIAL RESOLUTION: DR(1)=',DR(1),' => DECREASE INTERNAL RADIAL STEP'
     END IF
  ELSE ! Finds the STRR needed to reach RMAXSTR with a Newton-Raphson
     DOM=RREG-RMIN
     RHS=(RMAXSTR-RREG)*NRREG/RREG
     EX=NR-NRREG
     ! Second order analytic guess on ASTRR
     ASTRR=(3.*(EX - EX**2 + SQRT(-13.*EX**2 + 18.*EX**3 - 5.*EX**4 + 16.*EX*RHS - 24.*EX**2*RHS + &
          &8.*EX**3*RHS)/SQRT(3.)))/(2.*(2.*EX - 3.*EX**2 + EX**3))
     DO I=0,MAXITER
        LHS=-(1.-(1.+ASTRR)**(NR-NRREG))/ASTRR
        ASTRREQ1=LHS-RHS
        ASTRR1=ASTRR
        
        ASTRR=ASTRR*1.00001
        LHS=-(1.-(1.+ASTRR)**(NR-NRREG))/ASTRR
        ASTRREQ2=LHS-RHS
        ASTRR=ASTRR-ASTRREQ2*(ASTRR-ASTRR1)/(ASTRREQ2-ASTRREQ1)
        IF(ABS(ASTRR1-ASTRR) .LE. 1.0E-5)EXIT
     END DO
     STRR=1.+ASTRR
     DO I=1,NR
        IF(I .LE. NRREG) THEN
           R(I)=DOM*(I-0.5)/NRREG
        ELSE
           R(I)=R(I-1)+RREG/NRREG*STRR**(I-NRREG-1)
        END IF
     END DO
     IF(ABS(RMAXSTR-R(NR))/RMAXSTR .GE. 1.0E-2)THEN
        WRITE(6,*)'WARNING: MAXIMUM RADIUS IS LARGER THAN THE CHOSEN VALUE => INCREASE NR'
     END IF
     R(0)=-R(1)
     R(NR+1)=R(NR)+RREG/NRREG*STRR**(NR-NRREG)
     !(NR+1)=R(NR+1)*1.1
     DO I=1,NR
        !DR(I)=R(I+1)-R(I)    !!! CHANGE THIS
        DR(I)=R(I)-R(I-1)
     END DO
     !DR(NR+1)=DR(NR)
     DR(NR+1)=R(NR+1)-R(NR)
     DR(0)=DR(1)
     DO I=1,NR
        DRM(I)=R(I)-R(I-1)
        DRP(I)=R(I+1)-R(I)
     END DO
     DRM(0)=DRM(1)
     DRP(0)=R(1)-R(0)
     DRM(NR+1)=R(NR+1)-R(NR)
     DRP(NR+1)=DRP(NR)
     
     IF((DR(1) .GE. MINRESSTR) .AND. VERBOSE)THEN
        WRITE(6,*)'WARNING: LOW RADIAL RESOLUTION: DR(1)=',DR(1),' => DECREASE INTERNAL RADIAL STEP'
     END IF
  END IF
  
  ! Sefine the position of the medium radius of the domain (n4eded for mass convergence in TOVINI).
  IF(.NOT. STRETCH) THEN
     MIDGRID = NR/2.
  ELSE
     IF (RREG .lt. RMAX*0.5) THEN
        MIDGRID = INT( NRREG + LOG(1 + NRREG*(RMAX*0.5-RREG)*(STRR-1.)/(STRR*RREG))/LOG(STRR) )
     ELSE
        MIDGRID = INT( NRREG*RMAX*0.5/RREG)
     ENDIF
  END IF
  
  RMIDINV=1./R(MIDGRID)
  RMAXINV=1./R(NR)
  
  ! === Quantities for the DGTSV solve ===
  DO I=1,NR
     !DRM(I)=R(I)-R(I-1)
     !DRP(I)=R(I+1)-R(I)
     ACOEFF(I)=2.*(R(I)-DRP(I))/(R(I)*DRM(I)*(DRP(I)+DRM(I)))
     BCOEFF(I)=-2.*(R(I)+DRM(I)-DRP(I))/(R(I)*DRM(I)*DRP(I))
     CCOEFF(I)=2.*(R(I)+DRM(I))/(R(I)*DRP(I)*(DRP(I)+DRM(I)))
     IF(I .NE. NR)THEN
        DUMAT(I)=CCOEFF(I)
     END IF
     IF(I .NE. 1)THEN
        DLMAT(I-1)=ACOEFF(I)
     END IF
     IF(I==1)THEN
        DMAT(I)=ACOEFF(1)+BCOEFF(1)
     ELSE IF(I==NR)THEN
        DMAT(I)=BCOEFF(NR)+CCOEFF(NR)*R(NR)/R(NR+1)
     ELSE
        DMAT(I)=BCOEFF(I)
     END IF
  END DO
END SUBROUTINE GRIDBUILD

! **************************************************************
! **************************************************************

SUBROUTINE FUNCD_STRETCH(X,Y,DY,RHOVAR)
!   ============================================================
!   Purpose : Used by rtsafe to derive stretching factor
!             for the grid
!   ============================================================

  REAL :: X,Y,YN,Z,ZN,W,WN,DY,RHOVAR
  INTEGER :: NN

  NN=(NR-NRREG)

  Y=X*(1-X**NN)/(1-X)-NRREG*(RMAX-RREG)/(RREG-RMIN)
  DY= (NN*X**(NN+1.) - (1+NN)*X**NN +1.)/(X-1.)**2.

  RETURN

END SUBROUTINE FUNCD_STRETCH

! **************************************************************
! **************************************************************

  ! ********************************************************
  ! ********************************************************
  
  ! =============== READ EOS TABLE =============================
  SUBROUTINE EOSTABLEREAD

    INTEGER :: I,NPT
    ! Check if the array size are the same of the original resampled-table by python
    OPEN(12,FILE=FILEEOS)
    READ(12,*)NPT
    IF(NPT .NE. NPTRHO)THEN
       WRITE(6,*)'NPT is ',NPT,' not equal to NPTRHO'
       STOP
    END IF

    ! Read the bounds of the thermodynamical variables and the powerlaw for extrapolation
    READ(12,*)RHOTABMIN,RHOTABMAX
    READ(12,*)PRSTABMIN,PRSTABMAX,PRSTABINDL,PRSTABINDU
    READ(12,*)EINTABMIN,EINTABMAX,EINTABINDL,EINTABINDU
    READ(12,*)ENTTABMIN,ENTTABMAX,ENTTABINDL,ENTTABINDU

    ! Read the tables
    DO I=1,NPTRHO
       READ(12,*)RHOVEC1(I),PRSVEC1(I),EINVEC1(I),ENTVEC1(I),PRSIND1(I),EININD1(I),ENTIND1(I)
    END DO
    DO I=1,NPTRHO
       READ(12,*)PRSVEC2(I),RHOVEC2(I),RHOIND2(I)
    END DO
    DO I=1,NPTRHO
       READ(12,*)ENTVEC3(I),RHOVEC3(I),RHOIND3(I)
    END DO

    CLOSE(12)

  END SUBROUTINE EOSTABLEREAD
  ! ============================================================

  ! =============== RHO TO EOS =================================
  SUBROUTINE RHO2EOS(RHO,PRS,EIN,ENT)

    INTEGER :: IL
    REAL :: RHO,PRS,EIN,ENT
    REAL :: RHOLOG,PRSLOG,EINLOG,ENTLOG,RHOPT

    RHOLOG = LOG10(RHO)
    RHOPT = (RHOLOG-RHOTABMIN)/(RHOTABMAX-RHOTABMIN)
    IL = INT(RHOPT*(NPTRHO-1))

    IF(IL .GE. 0 .AND. IL .LE. 999)THEN
       PRSLOG = PRSVEC1(IL+1) +(RHOLOG -RHOVEC1(IL+1))*PRSIND1(IL+1)
       EINLOG = EINVEC1(IL+1) +(RHOLOG -RHOVEC1(IL+1))*EININD1(IL+1)
       ENTLOG = ENTVEC1(IL+1) +(RHOLOG -RHOVEC1(IL+1))*ENTIND1(IL+1)
    END IF
    IF(IL .LT. 0)THEN
       PRSLOG = PRSVEC1(1) +(RHOLOG -RHOVEC1(1))*PRSIND1(1)
       EINLOG = EINVEC1(1) +(RHOLOG -RHOVEC1(1))*EININD1(1)
       ENTLOG = ENTVEC1(1) +(RHOLOG -RHOVEC1(1))*ENTIND1(1)
    END IF
    IF(IL .GE. 1000)THEN
       PRSLOG = PRSVEC1(1000) +(RHOLOG -RHOVEC1(1000))*PRSIND1(1000)
       EINLOG = EINVEC1(1000) +(RHOLOG -RHOVEC1(1000))*EININD1(1000)
       ENTLOG = ENTVEC1(1000) +(RHOLOG -RHOVEC1(1000))*ENTIND1(1000)
    END IF
    PRS=10**PRSLOG
    EIN=10**EINLOG
    ENT=10**ENTLOG

  END SUBROUTINE RHO2EOS
  ! ============================================================

  ! =============== PRS TO EOS =================================
  SUBROUTINE PRS2EOS(PRS,RHO)

    INTEGER :: IL
    REAL :: RHO,PRS
    REAL :: RHOLOG,PRSLOG

    PRS = MAX(PRS,1.E-18)
    PRSLOG = LOG10(PRS)
    IL = INT((PRSLOG-PRSTABMIN)/(PRSTABMAX-PRSTABMIN)*(NPTRHO-1))
    IF(IL .GE. 0 .AND. IL .LE. 999)THEN
       RHOLOG = RHOVEC2(IL+1) +(PRSLOG -PRSVEC2(IL+1))*RHOIND2(IL+1)
    END IF
    IF(IL .LT. 0)THEN
     RHOLOG = RHOVEC2(1) +(PRSLOG -PRSVEC2(1))*RHOIND2(1)
    END IF
    IF(IL .GE. 1000)THEN
      RHOLOG = RHOVEC2(1000) +(PRSLOG -PRSVEC2(1000))/PRSIND1(1000)
    END IF
    RHO=10**RHOLOG

  END SUBROUTINE PRS2EOS
  ! ============================================================

   ! =============== ENT TO EOS =================================
  SUBROUTINE ENT2EOS(ENT,RHO)

    INTEGER :: IL
    REAL :: RHO,PRS,ENT
    REAL :: RHOLOG,PRSLOG,ENTLOG

    ENTLOG = LOG10(ENT)
    IL = INT((ENTLOG-ENTTABMIN)/(ENTTABMAX-ENTTABMIN)*(NPTRHO-1))
    IF(IL .GE. 0 .AND. IL .LE. 999)THEN
       RHOLOG = RHOVEC3(IL+1) +(ENTLOG -ENTVEC3(IL+1))*RHOIND3(IL+1)
    END IF
    IF(IL .LT. 0)THEN
     RHOLOG = RHOVEC3(1) +(ENTLOG -ENTVEC3(1))*RHOIND3(1)
    END IF
    IF(IL .GE. 1000)THEN
      RHOLOG = RHOVEC3(1000) +(ENTLOG -ENTVEC3(1000))/ENTIND1(1000)
    END IF
    RHO=10**RHOLOG

  END SUBROUTINE ENT2EOS
  ! ============================================================


  SUBROUTINE EOS(P,RHO,ENE,CHIS)

    ! === Gives the rest-mass density (RHO) and thermal energy (ENE and E) as functions of the pressure (P) (polytropic EoS) ===
    ! === EN is the ln(enthalpy)

    IMPLICIT NONE
    REAL :: RHOJ,PJ,EJ,ENJ
    REAL :: P,RHO,ENE,A,EN,E
    REAL,OPTIONAL :: CHIS
    LOGICAL :: EOSJOR

    IF(.NOT. PRESENT(CHIS)) EOSJOR=.TRUE.    ! If the optional variable CHIS is not present in the CALL EOS
    IF(PRESENT(CHIS)) EOSJOR=.FALSE.  ! If the optional variable CHIS is not present in the CALL EOS

    IF(GR)THEN
       A=1.
    ELSE
       IF(EOSJOR)THEN
          A=1.
       else
          A=EXP(ALPHA0*(CHIS-CHIINF)+0.5*BETA0*(CHIS-CHIINF)**2)
       endif
    END IF

    IF(P .GT. 1.D-18)THEN
       IF(EOSJOR)THEN ! Jordan frame EoS
          IF(EOSINT)THEN
             CALL PRS2EOS(P,RHO)
             CALL RHO2EOS(RHO,P,E,EN)
             ENE=E
          ELSE
             RHO=(P/K1)**(1./GAMMA)
             ENE=P/(GAMMA-1.)
          END IF
       ELSE  ! Einstein frame EoS
          IF(EOSINT)THEN
             PJ = P/A**4
             CALL PRS2EOS(PJ,RHOJ)
             CALL RHO2EOS(RHOJ,PJ,EJ,ENJ)
             RHO = RHOJ*A**4
             ENE = EJ*A**4
          ELSE
             RHO=(P/K1)**(1./GAMMA)*A**((4.*GAMMA-4.)/GAMMA)
             ENE=P/(GAMMA-1.)
          END IF
        END IF
       ELSE

          P   = 1.D-18
!           RHO = 1.D-15
! 	   CALL PRS2EOS(P,RHO)
		  RHO = (P/K1)**(1./GAMMA)*A**((4.*GAMMA-4.)/GAMMA)
          ENE  = 1.D-18
       END IF


       RETURN

     END SUBROUTINE EOS

SUBROUTINE FUNCD_EOS(X,Y,DY,RHOVAR)
!   ============================================================
!   Purpose : Used by rtsafe to derive central pressure
!             given the central density
!   ============================================================

  REAL :: X,Y,YN,Z,ZN,DY,RHOVAR,XN
  INTEGER :: NN


  CALL EOS(X,Y,Z)
  Y=Y-RHOVAR

  XN=1.0001*X
  CALL EOS(XN,YN,ZN)
  YN=YN-RHOVAR

  DY=(YN-Y)/(XN-X)

  RETURN

END SUBROUTINE FUNCD_EOS

! ***************************************************************************
! **************************************************************************

END MODULE PHYSICS
