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

SUBROUTINE GRIDBUILD(R,DR) !!! QUI
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

END MODULE PHYSICS
