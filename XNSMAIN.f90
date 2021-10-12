SUBROUTINE XNSMAIN(RHOVAR,QUCNV)
!   ============================================================
!   Purpose : This program solve the axisymmetric eq. for XCFC
!             in spherical coordinates, by Legendre Poly in Theta
!             direct solution in R. Uses XCFC.
!             This programs compute the stable equlibrium for
!             rotating and non-rotating NS, in XCFC.
!   ============================================================
!
! Input parameters
! RHOVAR = the desired central density
!
! Output parameter
! QUCNV = computed quantity of interest
!

! Parameters for the grid - set in module system
! R = radial grid points (+ boundaries) ; DR = increments
! TH = angular grid points (+ boundaries); DTH = increments
! XX = angular cos(th) points (+ boundaries)
!
! Initial TOV solution - set in module system - computed by tovini
! rhotv,prtv,etv   = density, pressure and energy density of the 1D-TOV solution
!
! RHOSRC,PSRC,ESRC = density, pressure and energy density of TOV on the 2D grid
! 				NB: during the initialization, ESRC does not contain the pressure, while it does contain the pressure
! 				    during the iterative process (ie, it beacomes the total energy rho*h)
! VPHI,VR,VTH      = azimuthal, radial, and theta component of the velocity (contrav)
! USRC,DSRC        = matter/em fields conserved variables (energy and density) corresponding to the initial condition
! USRCX            = scalar field conserved variable (energy) corresponding to the initial condition
! S1SRC,S2SRC,S3SRC= azimuthal, radial, and theta component of momentum of matter/em fields
! S1SRCX,S2SRCX,S3SRCX= azimuthal, radial, and theta component of momentum of the scalar field
!
! Computed quantities
! PSI = conformal factor - set in module system
! PSL = lapse - set in module system
! PSS = shift-phi - set in module system
!
!  Program by N. Bucciantini 2009 -  All rights reserved.
!   ============================================================


  USE SYSTEMXNS
  USE FUNCTIONS
  USE ROTATION, ONLY: OMEGAVALUE
  IMPLICIT NONE

  REAL :: QUCNV, RHOVAR

  REAL, DIMENSION(:), ALLOCATABLE :: RHOCVEC
  REAL, DIMENSION(:), ALLOCATABLE :: RHOCVECTEMP

  REAL :: B2,E2,BETALOC,RSTARLOC
  REAL :: OMEGALOC,OMGN,ALPH2GPP,MM,M0
  REAL :: B3S3,BRSR,BTSR,BYSR,EL2,EPHISR,EPOLROLD,EPOLRSR,EPOLTOLD,EPOLTSR,ERSR,ETSR,EYSR,GLF2
  REAL :: RHOCMAX,RHOCMIN,RNEWT,SSSC,V2
  INTEGER :: I,NCONV=0

  INTEGER, DIMENSION(1:NTH,1:NR) :: RHONINT, PNINT, PSIINT, V3NINT, ALPHAINT,B3NINT, PSSINT,LOGRINT

  INTEGER :: ILOOP
  LOGICAL :: EXT

  CHARACTER(len=200) :: delcommand

  ! Set initial density for the sequential model
   RHOINISEQ = RHOVAR
   RHOVAR = RHOVAR*MBARYONFC
  ! ............
  ! Check consistency
  ! ............

  IF(NTH .EQ.1)THEN
     IF(MLS .GT. 0)THEN
        WRITE(6,*)'ERROR: MLS > 0 IN 1D'
        STOP
     END IF
     IF(NGQ .GT. 1)THEN
        WRITE(6,*)'ERROR: NGQ > 1 IN 1D'
        STOP
     END IF
  END IF
  IF(MLS .GT. NTH)THEN
     WRITE(6,*)'ERROR: MLS > NTH'
     STOP
  END IF
  IF(NGQ .GT. NTH)THEN
     WRITE(6,*)'ERROR: NGQ > NTH'
     STOP
  END IF
  IF(MLS .GT. NGQ)THEN
     WRITE(6,*)'ERROR: MLS > NGQ'
     STOP
  END IF
  IF(EOSINT .AND. CTP)THEN
  	 WRITE(6,*)'ERROR: CTP can only be used with analytic polytropic EoS'
  	 STOP
  END IF

  ! ............
  ! Define grid
  ! ............

  ! Sets the radial grid
  CALL GRIDBUILD(R,DR)
  ! Set the angular theta grid (uniform for simplicity) and the cos grid
  ! and the boundary values for the cos grid out of [-1,1]
  ! Beware the Cos boundary are not trigonometrically coreect, but defined to ensure correct interpolation
  DO I=1,NTH
     TH(I)=PI*(I-0.5)/NTH
     XX(I)=COS(TH(I))
  END DO
  DO I=1,2
     TH(1-I)=-TH(I)
     XX(1-I)=+2.-XX(I)
     TH(NTH+I)=2*PI-TH(NTH+1-I)
     XX(NTH+I)=-2.-XX(NTH+1-I)
  END DO
  DO I=1,NTH+1
     DTH(I)=TH(I)-TH(I-1)
  END DO

  ! Write the grid
  IF((WRT.OR.(WRTF.AND.IDAT)).AND.WGRID)THEN
     OPEN(12,FILE=trim(adjustl(subdirpath))//'Grid.dat')
     WRITE(12,*)NTH,NR,NRREG,STRR
     DO IX=1,NTH
        WRITE(12,*)TH(IX)
     END DO
     DO IZ=1,NR
        WRITE(12,*)R(IZ)
     END DO
     CLOSE(12)
  END IF


  IF(EOSINT)THEN
	  CALL EOSTABLEREAD
  END IF

  IF(EOSINT .AND. (RHOVAR .GT. 10**RHOTABMAX))THEN
    WRITE(6,*)''
	WRITE(6,*)'WARNING: CHOSEN CENTRAL DENSITY'
	WRITE(6,*)'IS BEYOND THE MAXIMUM TABULATED VALUE OF ',10**RHOTABMAX
	WRITE(6,*)'POWER-LAW EXTRAPOLATION IS BEING USED'
	WRITE(6,*)''
  ENDIF
! ..........................
  ! Define initial conditions
  ! ..........................

  RHOCMAX=1. ;  RHOCMIN=0. !Used to check convergence on the central density

  ! Use TOV solution as inital data
  !   IF(VERBOSE) WRITE(6,*)'Computing TOV'

  CALL TOVINIMOD(RHOVAR)
!  DO I=1,NR
!     DR(I)=R(I)-R(I-1)
!  END DO
!  DR(NR+1)=DR(NR) !!! QUI
  ! Initialize the primitive variables in the domain + Atmosph.
  DO IX=1,nth
  	 WSURF(IX)=ISUR
     DO IZ=1,nr
        ! (Jordan) energy, density and pressure of the fluid
        RHOSRC(IX,IZ)= MAX(RHOTV(IZ),1.e-15)
        PSRC(IX,IZ)  = MAX(PRTV(IZ),1.e-18)
        ! ESRC(IX,IZ)  = RHOSRC(IX,IZ)+MAX(ETV(IZ),PSRC(IX,IZ))+PSRC(IX,IZ)
        ESRC(IX,IZ)  = RHOSRC(IX,IZ)+ETV(IZ)+PSRC(IX,IZ)
        ! (Jordan) contravariant velocities v^i (just a trial setup)
        VPHI(IX,IZ)  = 0.
        VTH(IX,IZ)   = 0.
        VR(IX,IZ)    = 0.
        ! (Jordan) contravariant magnetic field B^i (just a trial setup)
        BPHI(IX,IZ)  = 0.
        BPOLR(IX,IZ) = 0.
        BPOLT(IX,IZ) = 0.
        APHI(IX,IZ)  = 0.
        ! (Jordan) covariant electric field E_i (just a trial setup)
        EPHI(IX,IZ)  = 0.
        EPOLR(IX,IZ) = 0.
        EPOLT(IX,IZ) = 0.
        ATIM(IX,IZ)  = 0.
        ! (Einstein) metric terms
        PSI(IX,IZ)  = EXP(MU(IZ)/4.)
        PSL(IX,IZ)  = EXP(MU(IZ)/4.)*EXP(NU(IZ)/2.)
        PSS(IX,IZ)  = 0.
        PSSR(IX,IZ) = 0.
        PSST(IX,IZ) = 0.
        ! (Einstein) scalar field chi and its derivatives	P,Q_mu
        CHI(IX,IZ) = CHITV(IZ)
        PSCAL(IX,IZ) = 0.
        QSCALR(IX,IZ) = DCHITV(IZ)
        QSCALT(IX,IZ) = 0.
        QSCALP(IX,IZ) = 0.
        ! (Jordan) scalar field
        ASCAL(IX,IZ) = EXP(ALPHA0*(CHI(IX,IZ)-CHIINF)+0.5*BETA0*(CHI(IX,IZ)-CHIINF)**2)
        !Stellar interior
        IF(RHOSRC(ix,iz) .GT. 1.e-6)THEN
           BETALOC=0.
           RSTARLOC=(R(IZ)*SIN(TH(IX))*EXP(MU(IZ)/2.-NU(IZ)/2.))**2 ! Same in J and E-frames
           CALL OMEGAVALUE(BETALOC,RSTARLOC,OMEGALOC)
           OMGN=OMEGALOC
           VPHI(IX,IZ) = OMGN/(ASCAL(IX,IZ)*EXP(NU(IZ)/2.))
           VTH(IX,IZ)  = 0.0
           VR(IX,IZ)   = 0.0
           ALPH2GPP    = (R(IZ)*SIN(TH(IX))*EXP(MU(IZ)/2.+NU(IZ)/2.))**2
           IF(IMAG.AND.ITOR)THEN
              ! BPHI(IX,IZ)  = BCOEF*((RHOSRC(IX,IZ)+GAMMA/(GAMMA-1)*PSRC(IX,IZ))*ASCAL(IX,IZ)**4*ALPH2GPP)**MAGIND &
!                    / ALPH2GPP * EXP(NU(IZ)/2.)/ASCAL(IX,IZ)**3
				BPHI(IX,IZ)  = BCOEF*((ESRC(IX,IZ))*ASCAL(IX,IZ)**4*ALPH2GPP)**MAGIND &
                   / ALPH2GPP * EXP(NU(IZ)/2.)/ASCAL(IX,IZ)**3
           END IF
        END IF
     END DO
  END DO



  ! Write the source mass and energy term (used for comparing with TOV)
  IF(CHUP.AND.WRT) THEN
     OPEN(12,FILE=trim(adjustl(subdirpath))//'Source.dat') !FOLDER
     WRITE(12,*)NTH,NR
     DO IX=1,NTH
        DO IZ=1,NR
           WRITE(12,*)RHOSRC(IX,IZ),PSRC(IX,IZ),ESRC(IX,IZ)
        END DO
     END DO
     CLOSE(12)
  END IF

  ! Initialize the (Jordan) conserved variables in the domain: sqrt(gamma)*(E,D,S_j)=sqrt(f)*(hat(E),hat(D),hat(S_j))
  ! Assume that the (Einstein) metric is the one given by the STT TOV (correct only for 0 velocity)									!!!STT!!!
  DO IX=1,NTH
    DO IZ=1,NR
      B2=ASCAL(IX,IZ)**2*PSI(IX,IZ)**4*(BPHI(IX,IZ)*BPHI(IX,IZ)*R(IZ)**2*SIN(TH(IX))**2 + &
           BPOLR(IX,IZ)*BPOLR(IX,IZ) + BPOLT(IX,IZ)*BPOLT(IX,IZ)*R(IZ)**2)

      E2=ASCAL(IX,IZ)**(-2)/(PSI(IX,IZ)**4)*(EPHI(IX,IZ)*EPHI(IX,IZ)/R(IZ)**2./SIN(TH(IX))**2.+ &
           EPOLR(IX,IZ)*EPOLR(IX,IZ) + EPOLT(IX,IZ)*EPOLT(IX,IZ)/R(IZ)**2.)

      GLF2=1./(1.-ASCAL(IX,IZ)**2*PSI(IX,IZ)**4*(VR(IX,IZ)*VR(IX,IZ)+VTH(IX,IZ)*VTH(IX,IZ)*R(IZ)**2 + &
           VPHI(IX,IZ)*VPHI(IX,IZ)*R(IZ)**2*SIN(TH(IX))**2))

      ! USRC=sqrt(gamma)*E=sqrt(f)*hat(E)
      USRC(IX,IZ) =((ESRC(IX,IZ))*GLF2-PSRC(IX,IZ)+B2/2.+E2/2.)* &
      	    ASCAL(IX,IZ)**3*PSI(IX,IZ)**6*R(IZ)**2*SIN(TH(IX))
      ! DSRC=sqrt(gamma)*D=sqrt(f)*hat(D) with D=RHOSRC*GFL2
      DSRC(IX,IZ) =RHOSRC(IX,IZ)*ASCAL(IX,IZ)**3*PSI(IX,IZ)**6*R(IZ)**2*SIN(TH(IX))*SQRT(GLF2)
      ! COVARIANT CONSERVED MOMENTUM (SQRT(GAMMA)*S_J = SQRT(GAMMA) GAMMA_JI S^I)
      ! - USE TOV METRIC AS FIRST APPROX
      S1SRC(IX,IZ)=ASCAL(IX,IZ)**5*((ESRC(IX,IZ))*VR(IX,IZ)*GLF2 &
           *PSI(IX,IZ)**6*R(IZ)**2*SIN(TH(IX))*PSI(IX,IZ)**4)
      S2SRC(IX,IZ)=ASCAL(IX,IZ)**5*((ESRC(IX,IZ))*VTH(IX,IZ)*GLF2 &
           *PSI(IX,IZ)**6*R(IZ)**2*SIN(TH(IX))*R(IZ)**2*PSI(IX,IZ)**4)
      S3SRC(IX,IZ)=ASCAL(IX,IZ)**5*((ESRC(IX,IZ))*VPHI(IX,IZ)*GLF2 + &
           (EPOLR(IX,IZ)*BPOLT(IX,IZ)*PSI(IX,IZ)**4*R(IZ)**2.- &
           EPOLT(IX,IZ)*BPOLR(IX,IZ)*PSI(IX,IZ)**4)/(ASCAL(IX,IZ)*PSI(IX,IZ)**6*R(IZ)**2*SIN(TH(IX)))) &
           *PSI(IX,IZ)**6*R(IZ)**2*SIN(TH(IX))*R(IZ)**2*SIN(TH(IX))**2*PSI(IX,IZ)**4
    END DO
  END DO
!   WRITE(6,*)'PRIMA',RHOSRC(30,40),DSRC(30,40),ASCAL(30,40),PSI(30,40)**6.
  ! Initialize the (Einstein) conserved variables for the scalar field																!
  DO IX=1,NTH
     DO IZ=1,NR
        ! Computes the (Einstein) derivatives of the scalar field chi
        QSCAL2(IX,IZ) = 1./PSI(IX,IZ)**4*(QSCALR(IX,IZ)**2+QSCALT(IX,IZ)**2/R(IZ)**2+QSCALP(IX,IZ)**2/(R(IZ)**2*SIN(TH(IX))**2))
        ! USRCX=sqrt(gamma)*(P^2+Q^2)/(8*Pi)
        USRCX(IX,IZ) = (PSCAL(IX,IZ)**2+QSCAL2(IX,IZ))/(8.*PI)*PSI(IX,IZ)**6*R(IZ)**2*SIN(TH(IX))
        ! SiSRCX=(SQRT(GAMMA)*S(s)_i/(4*Pi)
        S1SRCX(IX,IZ)=1./(4.*PI)*PSCAL(IX,IZ)*QSCALR(IX,IZ)*PSI(IX,IZ)**6*R(IZ)**2*SIN(TH(IX))
        S2SRCX(IX,IZ)=1./(4.*PI)*PSCAL(IX,IZ)*QSCALT(IX,IZ)*PSI(IX,IZ)**6*R(IZ)**2*SIN(TH(IX))
        S3SRCX(IX,IZ)=1./(4.*PI)*PSCAL(IX,IZ)*QSCALP(IX,IZ)*PSI(IX,IZ)**6*R(IZ)**2*SIN(TH(IX))
     END DO
  END DO


  ! ------------------------------------------------------------------------------
  ! ------- Begin iteration loop to converge to model ----------------------------
  ! ------------------------------------------------------------------------------

  WRT=.TRUE.
  !If only final files are desired (WRTF=.TRUE.), set WRT to false
  !WRT will become true only when convergence is reached
  IF(WRTF) WRT=.FALSE.
  EXT=.FALSE.
  ILOOP=1

  ALLOCATE(RHOCVEC(MAXLOOP))

  DO !ILOOP=1,MAXLOOP
     IF(.NOT. MPICODE)THEN
        WRITE(6,*)''
        WRITE(6,*)'Loop number',ILOOP,'out of',MAXLOOP
     END IF
     IF((ILOOP .EQ. MAXLOOP) .AND. (WRTF.AND.IDAT))THEN
        WRT=.TRUE.
     END IF
     ! write(6,*)rhonew(1,1)
     ! ------------------------------------------------------------------------------
     ! ------- X-vector Equation Solver ---------------------------------------------
     ! ------------------------------------------------------------------------------

     ! Solve for the x-vector in phi

     ! Initialize the (Einstein) phi-xshift coefficient in the domain - 8*Pi*f^phiphi*(hat(S)_phi+hat(Ss)_phi) 						!!!STT!!!
     DO IX=1,NTH
        DO IZ=1,NR
           ES3SRC(IX,IZ)=8.*PI*(S3SRC(IX,IZ)+S3SRCX(IX,IZ))/ &
                (R(IZ)**2*SIN(TH(IX)))/(R(IZ)**2*SIN(TH(IX))**2)
        END DO
     END DO

     ! Solve phi-component of vector poisson
     CALL SHIFTPHI(ES3SRC)
     IF(VERBOSE) WRITE(6,*)'X-Shift phi .... done'

     ! Write the XShift-phi component and its source (for test)
     IF(CHUP.AND.WRT)THEN
        OPEN(12,FILE=trim(adjustl(subdirpath))//'XShiftphi.dat') !FOLDER
        WRITE(12,*)NTH,NR
        DO IX=1,NTH
           DO IZ=1,NR
              WRITE(12,*) PSS(IX,IZ),PSS(IX,IZ)/R(IZ)/SIN(TH(IX)),ES3SRC(IX,IZ)
           END DO
        END DO
        CLOSE(12)
     END IF

     ! Write new source term for the phi-Shift Equation (X-vector contravariant)
     DO IX=1,NTH
        DO IZ=1,NR
           ES3SRC(IX,IZ)=PSS(IX,IZ)/R(IZ)/SIN(TH(IX))
        END DO
     END DO

     ! Initialize the (Einstein) pol-xshift coefficient in the domain
     DO IX=1,NTH
        DO IZ=1,NR
           ES1SRC(IX,IZ)=8.*PI*(S1SRC(IX,IZ)+S1SRCX(IX,IZ))/(R(IZ)**2*SIN(TH(IX)))
           ES2SRC(IX,IZ)=8.*PI*(S2SRC(IX,IZ)+S2SRCX(IX,IZ))/(R(IZ)**2*SIN(TH(IX)))/(R(IZ)**2)
        END DO
     END DO

     ! Solve poloidal-component of vector poisson (not needed for Equilibrium solution)
     ! call shiftpol(es1src,es2src)
     PSSR(:,:)=0.
     PSST(:,:)=0.

     ! Write new source term for the poloidal-Shift Equation (X-vector contravariant)
     DO IX=1,NTH
        DO IZ=1,NR
           ES1SRC(IX,IZ)=PSSR(IX,IZ)
           ES2SRC(IX,IZ)=PSST(IX,IZ)/R(IZ)
        END DO
     END DO

     ! Initialize the (Einstein) conserved variables for the scalar field
     DO IX=1,NTH
        DO IZ=1,NR
           ! Computes the (Einstein) derivatives of the scalar field chi
           QSCAL2(IX,IZ) = 1./PSI(IX,IZ)**4*(QSCALR(IX,IZ)**2+QSCALT(IX,IZ)**2/R(IZ)**2+QSCALP(IX,IZ)**2/(R(IZ)**2*SIN(TH(IX))**2))
           ! USRCX=sqrt(gamma)*(P^2+Q^2)/(8*Pi)
           USRCX(IX,IZ) = (PSCAL(IX,IZ)**2+QSCAL2(IX,IZ))/(8.*PI)*PSI(IX,IZ)**6*R(IZ)**2*SIN(TH(IX))
           ! SiSRCX=(SQRT(GAMMA)*S(s)_i/(4*Pi)
           S1SRCX(IX,IZ)=1./(4.*PI)*PSCAL(IX,IZ)*QSCALR(IX,IZ)*PSI(IX,IZ)**6*R(IZ)**2*SIN(TH(IX))
           S2SRCX(IX,IZ)=1./(4.*PI)*PSCAL(IX,IZ)*QSCALT(IX,IZ)*PSI(IX,IZ)**6*R(IZ)**2*SIN(TH(IX))
           S3SRCX(IX,IZ)=1./(4.*PI)*PSCAL(IX,IZ)*QSCALP(IX,IZ)*PSI(IX,IZ)**6*R(IZ)**2*SIN(TH(IX))
        END DO
     END DO

     ! ------------------------------------------------------------------------------
     ! ------- Conformal factor Psi Equation Solver ---------------------------------
     ! ------------------------------------------------------------------------------

     ! Initialize the (Einstein) conformal energy source the domain - -2*Pi*(hat(E)+hat(Es))											!!!STT!!!
     DO IX=1,NTH
        DO IZ=1,NR
           ECSRC(IX,IZ)=-2.*PI*(ASCAL(IX,IZ)*USRC(IX,IZ)+USRCX(IX,IZ))/(R(IZ)**2*SIN(TH(IX)))
        END DO
     END DO

     ! Initialize the (Einstein) conformal curvature source the domain A^ij*A_ij / 8
     CALL CURV1

     IF(DEBUG .AND. (.NOT. MPICODE))WRITE(6,*)'before PSI: ',PSI(1,1)**4.,EXP(MU(1))
!      WRITE(6,*)'before PSI: ',PSI(1,1),EXP(MU(1)/4.)
!      WRITE(6,*)'source PSI: ',ASCAL(1,1),USRC(1,1),USRCX(1,1)
     ! Solve poisson equation by recursion
     PSI0=PSI
     WRITE(6,*)'CONF CHECK 1: ',(PSI(NTH/2,NR-2)-1)*R(NR-2)*2
     CALL CONFORMAL(ECSRC,CURVC)

     ! Conformal factor Psi
     WRITE(6,*)'CONF CHECK 2: ',PSI(NTH/2,NR-2)*R(NR-2)*2
     DO IX=1,NTH
        DO IZ=1,NR
           PSI(IX,IZ)=1.+PSI(IX,IZ)
        END DO
     END DO
     PSI = QFACTORCONF*PSI + (1.-QFACTORCONF)*PSI0
     IF(DEBUG .AND. (.NOT. MPICODE))WRITE(6,*)'after PSI:  ',PSI(1,1)**4.,EXP(MU(1))
!      WRITE(6,*)'after PSI:  ',PSI(1,1),EXP(MU(1)/4.)

     ! Write the conformal factor and its sources
     IF(CHUP.AND.WRT)THEN
        OPEN(12,FILE=trim(adjustl(subdirpath))//'Conformal.dat') !FOLDER
        WRITE(12,*)NTH,NR
        DO IX=1,NTH
           DO IZ=1,NR
              WRITE(12,*)PSI(IX,IZ),ECSRC(IX,IZ),CURVC(IX,IZ)
           END DO
        END DO
        CLOSE(12)
     END IF

! 	WRITE(6,*)'For eos check (conformal)',PSI(30,40)**6.,EXP(3.*MU(40)/2.),RHOSRC(30,40)

     ! ------------------------------------------------------------------------------
     ! ------- Derivation of primitive variables ------------------------------------
     ! ------------------------------------------------------------------------------

     ! Once the conformal factor is known the cons_to_prim can be executed

     ! Simplified version just for TOV static case (assume S_i = 0)
     ! To be changed with the true cons_to_prim

     ! If purely toroidal
     IF(.NOT.(IPOL.OR.ITWT))THEN
        DO IX=1,NTH
           DO IZ=1,NR
              ! (Jordan) fluid/em variables
              USR=USRC(IX,IZ)/(R(IZ)**2*SIN(TH(IX)))/(ASCAL(IX,IZ)**3*PSI(IX,IZ)**6)      ! E
              DSR=DSRC(IX,IZ)/(R(IZ)**2*SIN(TH(IX)))/(ASCAL(IX,IZ)**3*PSI(IX,IZ)**6)      ! D
              SYSR=S3SRC(IX,IZ)/(R(IZ)**2*SIN(TH(IX)))/(ASCAL(IX,IZ)**3*PSI(IX,IZ)**6)    ! S_phi
              SY2SR=SYSR**2/(R(IZ)**2*SIN(TH(IX))**2)/(ASCAL(IX,IZ)**2*PSI(IX,IZ)**4)     ! S_phi*S^phi
              B2YSR=BPHI(IX,IZ)**2*R(IZ)**2*SIN(TH(IX))**2*ASCAL(IX,IZ)**2*PSI(IX,IZ)**4  ! B^2=B_phi*B^phi
              B3S3=BPHI(IX,IZ)*S3SRC(IX,IZ)                            		    ! B^phi*sqrt(f)*hat(S_phi)

              ! (Einstein) scalar field variables
              USRX=USRCX(IX,IZ)/(R(IZ)**2*SIN(TH(IX)))/PSI(IX,IZ)**6                      ! Es
              SYSRX=S3SRCX(IX,IZ)/(R(IZ)**2*SIN(TH(IX)))/PSI(IX,IZ)**6                    ! Ss_phi
              SY2SRX=SYSRX**2/(R(IZ)**2*SIN(TH(IX))**2)/PSI(IX,IZ)**4			    ! Ss_phi*Ss^phi

              ! IF(IX==30 .AND. IZ==40)THEN
! 				WRITE(6,*)'before CTP',BPHI(IX,IZ),RHOSRC(IX,IZ),PSRC(IX,IZ),SSS(IX,IZ)
! 			  ENDIF

			  IF(CTP)THEN
! 			  	  IF(IX==10 .AND. IZ==40)WRITE(6,*)'BEFORE CTP',USR,DSR
				  CALL CONS_TO_PRIM(USR,SY2SR,DSR,B2YSR,DENSC,PRESSC,VELOCC,BMAG,SSSC)
! 				  IF(IX==10 .AND. IZ==76)WRITE(6,*)'AFTER CTP',USR,DSR,SSSC
! 				  IF(IX==10 .AND. IZ==76)STOP
! 				  ESRC(IX,IZ)   = USRC(IX,IZ)/(R(IZ)**2*SIN(TH(IX)))/(ASCAL(IX,IZ)**3*PSI(IX,IZ)**6)	! E
				  RHOSRC(IX,IZ) = DENSC                							! rho=D/GammaL
				  PSRC(IX,IZ)   = PRESSC
				  VPHI(IX,IZ)   = VELOCC/(R(IZ)*SIN(TH(IX))*ASCAL(IX,IZ)*PSI(IX,IZ)**2)
				  BPHI(IX,IZ)   = BMAG/(R(IZ)*SIN(TH(IX))*ASCAL(IX,IZ)*PSI(IX,IZ)**2)
				  SSS(IX,IZ)    = SSSC	! S=tr(S^munu) if purely toroidal
			  ELSE
				  V2=(VPHI(IX,IZ)*(R(IZ)*SIN(TH(IX))*ASCAL(IX,IZ)*PSI(IX,IZ)**2))**2
				  B2=(BPHI(IX,IZ)*(R(IZ)*SIN(TH(IX))*ASCAL(IX,IZ)*PSI(IX,IZ)**2))**2
				  SSS(IX,IZ)=(ESRC(IX,IZ))*V2/(1.-V2)-B2+3*(PSRC(IX,IZ)+0.5*B2)
			  ENDIF
! 			IF(IX==30 .AND. IZ==40)THEN
! 				WRITE(6,*)'SSSC',IX,IZ,RHOSRC(IX,IZ),DENSC
! 			ENDIF
			! IF(IX==30 .AND. IZ==40)THEN
! 				WRITE(6,*)'after CTP',BPHI(IX,IZ),RHOSRC(IX,IZ),PSRC(IX,IZ),SSS(IX,IZ)
! 			ENDIF
           END DO
        END DO
     END IF

     ! If poloidal or twisted-torus
     IF(IPOL.OR.ITWT)THEN
        DO IX=1,NTH
           DO IZ=1,NR
              ! (Jordan) fluid/em variables
              USR=USRC(IX,IZ)/(R(IZ)**2*SIN(TH(IX)))/(ASCAL(IX,IZ)**3*PSI(IX,IZ)**6)		 ! E
              DSR=DSRC(IX,IZ)/(R(IZ)**2*SIN(TH(IX)))/(ASCAL(IX,IZ)**3*PSI(IX,IZ)**6)		 ! D
              SYSR=S3SRC(IX,IZ)/(R(IZ)**2*SIN(TH(IX)))/(ASCAL(IX,IZ)**3*PSI(IX,IZ)**6)    ! S_phi
              SY2SR=SYSR**2/(R(IZ)**2*SIN(TH(IX))**2)/(ASCAL(IX,IZ)**2*PSI(IX,IZ)**4)     ! S_phi*S^phi
              BYSR=BPHI(IX,IZ)  ! B^phi
              BRSR=BPOLR(IX,IZ) ! B^r
              BTSR=BPOLT(IX,IZ)	! B^th
              EPHISR=EPHI(IX,IZ)		 ! E_phi
              EPOLRSR=EPOLR(IX,IZ)			 ! E_r
              EPOLTSR=EPOLT(IX,IZ)					 ! E_th
              EPOLROLD=EPOLRSR/(ASCAL(IX,IZ)**2*PSI(IX,IZ)**4)            ! E_phi
              EPOLTOLD=EPOLTSR/R(IZ)**2/(ASCAL(IX,IZ)**2*PSI(IX,IZ)**4)	! E_th

              ! (Einstein) scalar field variables
              USRX=USRCX(IX,IZ)/(R(IZ)**2*SIN(TH(IX)))/PSI(IX,IZ)**6						 ! Es
              SYSRX=S3SRCX(IX,IZ)/(R(IZ)**2*SIN(TH(IX)))/PSI(IX,IZ)**6					 ! Ss_phi
              SY2SRX=SYSRX**2/(R(IZ)**2*SIN(TH(IX))**2)/PSI(IX,IZ)**4			     	 ! Ss_phi*Ss^phi

              CALL COVTERM(IX,IZ)     	! Computes the covariant terms of the metric tensor in the E-frame
              IF(CTP) THEN			! if FALSE the code doesn't use CONS_TO_PRIM
!               	 IF((IZ .EQ. 10) .AND. (IX .EQ. 100))WRITE(6,*)'CTP1',USR,DSR
                 CALL CONS_TO_PRIM_POL(USR,SYSR,DSR,BYSR,BRSR,BTSR,EPHISR,EPOLRSR,EPOLTSR,DENSC,PRESSC,VELOCC,SSSC,IX,IZ)  	!!!STT!!!
! 	 			 IF((IZ .EQ. 10) .AND. (IX .EQ. 100))WRITE(6,*)'CTP2',DENSC,PRESSC,USR,DSR,RHOSRC(IX,IZ),DSRC(IX,IZ),VELOCC
	 			 RHOSRC(IX,IZ)=DENSC
                 PSRC(IX,IZ)=PRESSC
                 VPHI(IX,IZ)=VELOCC/(R(IZ)*SIN(TH(IX))*ASCAL(IX,IZ)*PSI(IX,IZ)**2)
                 SSS(IX,IZ)=SSSC
              ELSE
                 V2=(VPHI(IX,IZ)*(R(IZ)*SIN(TH(IX))*ASCAL(IX,IZ)*PSI(IX,IZ)**2))**2.	     ! v^2=v^phi*v_phi								!!!STT!!!
                 B2=ASCAL(IX,IZ)**2*(BYSR*BYSR*GCOVP+BRSR*BRSR*GCOVR+BTSR*BTSR*GCOVT)      ! B^2=B^i*B_i
                 EL2=ASCAL(IX,IZ)**(-2)*(EYSR*EYSR/GCOVP+ERSR*ERSR/GCOVR+ETSR*ETSR/GCOVT)  ! E^2=E^i*E_i
                 SSS(IX,IZ)=(ESRC(IX,IZ))*V2/(1.-V2)-B2-EL2+3*(PSRC(IX,IZ)+0.5*B2+0.5*EL2) ! S=tr(S^munu)	if poloidal or twt
              ENDIF
           END DO
        END DO
     END IF

     ! Write the primitive variables - to compared with TOV
     IF(CHUP.AND.WRT)THEN
        OPEN(12,FILE=trim(adjustl(subdirpath))//'Primitive.dat') !FOLDER
        WRITE(12,*)NTH,NR
        DO IX=1,NTH
           DO IZ=1,NR
              WRITE(12,*)RHOSRC(IX,IZ),PSRC(IX,IZ),ESRC(IX,IZ),VPHI(IX,IZ),BPHI(IX,IZ)
           END DO
        END DO
        CLOSE(12)
        IF(IPOL.OR.ITWT)THEN
           OPEN(12,FILE=trim(adjustl(subdirpath))//'Primitive_mag.dat') !FOLDER
           WRITE(12,*)NTH,NR
           DO IX=1,NTH
              DO IZ=1,NR
                 WRITE(12,*)BPHI(IX,IZ),BPOLR(IX,IZ),BPOLT(IX,IZ)
              END DO
           END DO
           CLOSE(12)
        END IF
     END IF

     ! ------------------------------------------------------------------------------
     ! ------- Lapse equation Equation Solver ---------------------------------
     ! ------------------------------------------------------------------------------

     ! Solve for the lapse*Psi -1.

     ! Initialize the (Einstein) lapse energy source the domain - 2*Pi*(hat(E)+hat(Es)+2*hat(S)+2*hat(Ss))/PSI^2						!!!STT!!!
     DO IX=1,NTH
        DO IZ=1,NR
           ELSRC(IX,IZ)=2.*PI*((ASCAL(IX,IZ)*USRC(IX,IZ)+USRCX(IX,IZ))/(R(IZ)**2*SIN(TH(IX)))+ &
                2*(ASCAL(IX,IZ)**4*SSS(IX,IZ)*PSI(IX,IZ)**6-USRCX(IX,IZ)/(R(IZ)**2*SIN(TH(IX)))))/PSI(IX,IZ)**2						! scrivere generalizzato
        END DO
     END DO

     IF(DEBUG .AND. (.NOT. MPICODE))WRITE(6,*)'before PSL: ',PSL(1,1),EXP(MU(1)/4.)*EXP(NU(1)/2.),CURVC(1,1),ELSRC(1,1)
     WRITE(6,*)'before PSL: ',PSL(1,1),EXP(MU(1)/4.)*EXP(NU(1)/2.),ELSRC(1,1),CURVC(1,1)
!     WRITE(6,*)'source PSL: ',ELSRC(1,1),ELSRC(1,50),ELSRC(1,100),CURVC(1,1),CURVC(1,50),CURVC(1,100)
! 	WRITE(6,*)'source PSL: ',ASCAL(1,1),USRC(1,1),USRCX(1,1),SSS(1,1),PSI(1,1)
     ! Solve poisson equation by recursion

     CALL LAPSE(ELSRC,CURVC)

     DO IX=1,NTH
        DO IZ=1,NR
!if(isnan(PSL(ix,iz))) then
!print*, 'PSL NAN in', ix,iz
!print*, 'XNSMAIN'
!stop
!endif
           PSL(IX,IZ)=1.+PSL(IX,IZ)
        END DO
     END DO

     IF(DEBUG .AND. (.NOT. MPICODE))WRITE(6,*)'after PSL:   ',PSL(1,1),EXP(MU(1)/4.)*EXP(NU(1)/2.)
     WRITE(6,*)'after PSL:   ',PSL(1,1),EXP(MU(1)/4.)*EXP(NU(1)/2.),PSI(1,1),EXP(MU(1)/4.)

     ! Write the lapse*Psi and its energy-source
     IF(CHUP.AND.WRT)THEN
        OPEN(12,FILE=trim(adjustl(subdirpath))//'Lapse.dat')
        WRITE(12,*)NTH,NR
        DO IX=1,NTH
           DO IZ=1,NR
              WRITE(12,*)PSL(IX,IZ),ELSRC(IX,IZ)
           END DO
        END DO
        CLOSE(12)
     END IF

     ! ------------------------------------------------------------------------------
     ! ------- Shift-vector Equation Solver -----------------------------------------
     ! ------------------------------------------------------------------------------

     ! Solve for the shift in phi

     ! Initialize the shift source in the domain
     CALL CURV2

     IF(DEBUG .AND. (.NOT. MPICODE))WRITE(6,*)'before PSS: ',PSS(1,1)

     ! Solve phi-component of vector poisson
     CALL SHIFTPHI(CURVP)
     IF(VERBOSE) WRITE(6,*)'Shift phi ... done'
     IF(DEBUG .AND. (.NOT. MPICODE))WRITE(6,*)'after PSS:  ',PSS(1,1)

     ! Phi component of the shift vector - contravariant
     DO IX=1,NTH
        DO IZ=1,NR
           PSS(IX,IZ)=PSS(IX,IZ)/R(IZ)/SIN(TH(IX))
        END DO
     END DO

     ! Write the Shift-phi component and its source (for test)
     IF(CHUP.AND.WRT)THEN
        OPEN(12,FILE=trim(adjustl(subdirpath))//'Shiftphi.dat')
        WRITE(12,*)NTH,NR
        DO IX=1,NTH
           DO IZ=1,NR
              WRITE(12,*)PSS(IX,IZ),CURVP(IX,IZ)
           END DO
        END DO
        CLOSE(12)
     END IF

     PSSR(:,:)=0.
     PSST(:,:)=0.

     ! ------------------------------------------------------------------------------
     ! -------------------- Scalar field chi Equation Solver ------------------------
     ! ------------------------------------------------------------------------------

     DO IX=1,NTH
        DO IZ=1,NR
           IF(IZ .LE. WSURF(IX))THEN
           	  TRACEM(IX,IZ)=4.*PSRC(IX,IZ)-(ESRC(IX,IZ))	! ESRC = rho*h
           ELSE
           	  TRACEM(IX,IZ)=0
           ENDIF
        END DO
     END DO

     !WRITE(6,*)'TRACEM',PSRC(100,50),TRACEM(100,50),ESRC(100,50),RHOSRC(100,50)

!      WRITE(6,*)'TRACEM',TRACEM(100,1),PSRC(100,1),ESRC(100,1),RHOSRC(100,1)
!      TRACEM(100,50),TRACEM(100,100)

     IF(DEBUG .AND. (.NOT. MPICODE))WRITE(6,*)'before CHI:  ',chi(1,1),chitv(1)
    ! WRITE(6,*)'before CHI:  ',chi(1,1),chitv(1)
     CALL CHISOL(TRACEM)
    ! WRITE(6,*)'after CHI:   ',chi(1,1),chitv(1)
     IF(DEBUG .AND. (.NOT. MPICODE))WRITE(6,*)'after CHI:   ',chi(1,1),chitv(1)

     ! Write the scalar field
     IF(CHUP.AND.WRT)THEN
        OPEN(12,FILE=trim(adjustl(subdirpath))//'Chi.dat') !FOLDER
        WRITE(12,*)NTH,NR
        DO IX=1,NTH
           DO IZ=1,NR
              WRITE(12,*)CHI(IX,IZ),TRACEM(IX,IZ)
           END DO
        END DO
        CLOSE(12)
     END IF

     ! ------------------------------------------------------------------------------
     ! ------- Hydrostatic Equilibrium  -----------------------------------------
     ! ------------------------------------------------------------------------------

     ! Solve the hydrostatic equlibrium in the newly computed CFC metric
     ENDID=.FALSE.
     CALL HYDROEQ(RHOVAR,ILOOP)
     IF(ENDID)THEN
  	GOTO 200
     ENDIF

     IF(DEBUG .AND. (.NOT. MPICODE))WRITE(6,*)'RHONEW:     ',RHONEW(1,1),rhotv(1)

     IF(XNSERR.NE.0) RETURN
     IF(VERBOSE) WRITE(6,*)'Hydro equilib. ... done'

     ! Re-Initialize the conserved variables in the domain
     ! Assume that metric is the one coputed previously and new density/pressure

     IF(DEBUG .AND. (.NOT. MPICODE))WRITE(6,*)'before B2:  ',BPHI(1,1)*BPHI(1,1)*R(1)**2*SIN(TH(1))**2*ASCAL(1,1)**2*PSI(1,1)**4

!      WRITE(6,*)'BEFORE B',BPHI(30,40)
     IF(.NOT.(IPOL.OR.ITWT))THEN
        DO IX=1,NTH
           DO IZ=1,NR
              ESRC(IX,IZ)=RHONEW(IX,IZ)+ENEW(IX,IZ)+PNEW(IX,IZ)
              RHOSRC(IX,IZ)=RHONEW(IX,IZ)
              PSRC(IX,IZ)=PNEW(IX,IZ)
              VPHI(IX,IZ)=V3NEW(IX,IZ)
              BPHI(IX,IZ)=B3NEW(IX,IZ)
              B2  = BPHI(IX,IZ)*BPHI(IX,IZ)*R(IZ)**2*SIN(TH(IX))**2*ASCAL(IX,IZ)**2*PSI(IX,IZ)**4
              GLF2=1./(1.-VPHI(IX,IZ)*VPHI(IX,IZ)*R(IZ)**2*SIN(TH(IX))**2*ASCAL(IX,IZ)**2*PSI(IX,IZ)**4)

              USRC(IX,IZ) =(ESRC(IX,IZ)*GLF2-PNEW(IX,IZ)+B2/2.)*ASCAL(IX,IZ)**3*PSI(IX,IZ)**6*R(IZ)**2*SIN(TH(IX))

              DSRC(IX,IZ) =RHOSRC(IX,IZ)*ASCAL(IX,IZ)**3*PSI(IX,IZ)**6*R(IZ)**2*SIN(TH(IX))*SQRT(GLF2)

              ! Covariant conserved momentum (sqrt(gamma)*S_j = sqrt(gamma) gamma_ji S^i)
              ! - use TOV metric as first approx
              S1SRC(IX,IZ)=ASCAL(IX,IZ)**5*(ESRC(IX,IZ))*VR(IX,IZ)*GLF2 &
                   *PSI(IX,IZ)**6*R(IZ)**2*SIN(TH(IX))*PSI(IX,IZ)**4
              S2SRC(IX,IZ)=ASCAL(IX,IZ)**5*(ESRC(IX,IZ))*VTH(IX,IZ)*GLF2 &
                   *PSI(IX,IZ)**6*R(IZ)**2*SIN(TH(IX))*R(IZ)**2*PSI(IX,IZ)**4
              S3SRC(IX,IZ)=ASCAL(IX,IZ)**5*(ESRC(IX,IZ))*VPHI(IX,IZ)*GLF2 &
                   *PSI(IX,IZ)**6*R(IZ)**2*SIN(TH(IX))*R(IZ)**2*SIN(TH(IX))**2*PSI(IX,IZ)**4
           END DO
        END DO
     END IF
!      WRITE(6,*)'AFTER B',BPHI(30,40)
     IF(IPOL.OR.ITWT)THEN
        DO IX=1,NTH
           DO IZ=1,NR
				ESRC(IX,IZ)=RHONEW(IX,IZ)+ENEW(IX,IZ)+PNEW(IX,IZ)																			!!!STT!!!
				RHOSRC(IX,IZ)=RHONEW(IX,IZ)
				PSRC(IX,IZ)=PNEW(IX,IZ)
				VPHI(IX,IZ)=V3NEW(IX,IZ)
				BPHI(IX,IZ)=B3NEW(IX,IZ)
				EPHI(IX,IZ)=E3NEW(IX,IZ)
				B2=ASCAL(IX,IZ)**2*(BPHI(IX,IZ)*BPHI(IX,IZ)*R(IZ)**2.*SIN(TH(IX))**2.*PSI(IX,IZ)**4. + &
				   BPOLR(IX,IZ)*BPOLR(IX,IZ)*PSI(IX,IZ)**4. + &
				   BPOLT(IX,IZ)*BPOLT(IX,IZ)*R(IZ)**2*PSI(IX,IZ)**4.)
				E2=ASCAL(IX,IZ)**(-2)*(EPHI(IX,IZ)*EPHI(IX,IZ)/R(IZ)**2./SIN(TH(IX))**2./PSI(IX,IZ)**4. + &
				   EPOLR(IX,IZ)*EPOLR(IX,IZ)/PSI(IX,IZ)**4. + &
				   EPOLT(IX,IZ)*EPOLT(IX,IZ)/R(IZ)**2./PSI(IX,IZ)**4.)
				GLF2=1./(1.-VPHI(IX,IZ)*VPHI(IX,IZ)*R(IZ)**2.*SIN(TH(IX))**2*ASCAL(IX,IZ)**2*PSI(IX,IZ)**4)
				USRC(IX,IZ) =(ESRC(IX,IZ)*GLF2-PNEW(IX,IZ)+B2/2.+E2/2.)*ASCAL(IX,IZ)**3*PSI(IX,IZ)**6.*R(IZ)**2.*SIN(TH(IX))
				DSRC(IX,IZ) =RHOSRC(IX,IZ)*SQRT(GLF2)*ASCAL(IX,IZ)**3*PSI(IX,IZ)**6.*R(IZ)**2.*SIN(TH(IX))
				! Covariant conserved momentum (sqrt(gamma)*S_j = sqrt(gamma) gamma_ji S^i)
				! - use TOV metric as first approx
				S1SRC(IX,IZ)=ASCAL(IX,IZ)**5*(ESRC(IX,IZ))*VR(IX,IZ)*GLF2 &
					 *PSI(IX,IZ)**6*R(IZ)**2.*SIN(TH(IX))*PSI(IX,IZ)**4
				S2SRC(IX,IZ)=ASCAL(IX,IZ)**5*(ESRC(IX,IZ))*VTH(IX,IZ)*GLF2 &
					 *PSI(IX,IZ)**6*R(IZ)**2.*SIN(TH(IX))*R(IZ)**2*PSI(IX,IZ)**4
				S3SRC(IX,IZ)=ASCAL(IX,IZ)**5*(ESRC(IX,IZ)*VPHI(IX,IZ)*GLF2 + &
					 (EPOLR(IX,IZ)*BPOLT(IX,IZ)*PSI(IX,IZ)**4.*R(IZ)**2.- &
					 EPOLT(IX,IZ)*BPOLR(IX,IZ)*PSI(IX,IZ)**4.)/(ASCAL(IX,IZ)*PSI(IX,IZ)**6.*R(IZ)**2.*SIN(TH(IX)))) &
					 *PSI(IX,IZ)**6*R(IZ)**2*SIN(TH(IX)) &
					 *R(IZ)**2*SIN(TH(IX))**2*PSI(IX,IZ)**4
      END DO
    END DO
  END IF

  IF(DEBUG .AND. (.NOT. MPICODE))WRITE(6,*)'after B2:   ',BPHI(1,1)*BPHI(1,1)*R(1)**2*SIN(TH(1))**2*ASCAL(1,1)**2*PSI(1,1)**4

  ! Write the  mass and pressure term (used for comparing Hydrost equil with TOV)

  IF(WRT)THEN
    OPEN(12,FILE=trim(adjustl(subdirpath))//'Hydroeq.dat')
    WRITE(12,*)NTH,NR,OMG
    DO IX=1,NTH
      DO IZ=1,NR
        WRITE(12,*)RHONEW(IX,IZ)/MBARYONFC,PNEW(IX,IZ),PSI(IX,IZ),V3NEW(IX,IZ),&
             PSL(IX,IZ)/PSI(IX,IZ),PSS(IX,IZ),CHI(IX,IZ),QSCALR(IX,IZ),QSCALT(IX,IZ)
      END DO
    END DO
    CLOSE(12)

    IF(IMAG)THEN
      OPEN(12,FILE=trim(adjustl(subdirpath))//'Hydroeq_mag.dat')
      WRITE(12,*)NTH,NR
      DO IX=1,NTH
        DO IZ=1,NR
        WRITE(12,*) B3NEW(IX,IZ),BPOLR(IX,IZ),BPOLT(IX,IZ),APHI(IX,IZ), &
             E3NEW(IX,IZ),EPOLR(IX,IZ),EPOLT(IX,IZ),ATIM(IX,IZ),&
             JPHI(IX,IZ),JRR(IX,IZ),JTH(IX,IZ)
        END DO
      END DO
      CLOSE(12)
    END IF

    IF(IMAG.AND.WRTF.AND.(OMG.NE.0.))THEN
      OPEN(12,FILE=trim(adjustl(subdirpath))//'Mxwll_test.dat')
      WRITE(12,*) NTH,NR,KBPOL
      DO IX=1,NTH
        DO IZ=1,NR
        WRITE(12,*) RHOEMXL(IX,IZ),JPHIMXL(IX,IZ),ATIMIN(IX,IZ),ATIMOUT(IX,IZ),ATIMARM(IX,IZ),&
              OMGMET(IX,IZ),GAMLOC(IX,IZ),DRGAML(IX,IZ),DTGAML(IX,IZ)
        END DO
      END DO
      CLOSE(12)
    ENDIF
  END IF

  ! Write the central density (this is used to check the convergence)
  RNEWT=(RHONEW(1,1)*R(2)**2-RHONEW(1,2)*R(1)**2)/(R(2)**2-R(1)**2)
  IF(VERBOSE .AND. (.NOT. MPICODE)) WRITE(6,'(A10,I3,A20,E12.7)')'Step = ',ILOOP,' Central density = ',RNEWT !real(rhonew(1,1))

  !-------------------------------------------
  !----------Check convergence ---------------
  !-------------------------------------------
  IF(.NOT. MPICODE)THEN
	  WRITE(6,*)'RHO CENTRAL',RNEWT
  ENDIF
  RHOCVEC(ILOOP)=RNEWT
  IF(EXT)THEN
    IF(WRTF.AND.IDAT) WRT=.FALSE.
    EXIT
  END IF

  OPEN(12,FILE=trim(adjustl(subdirpath))//'Rhovec.dat',access='append')
        WRITE(12,*) ILOOP,RNEWT,PSI(1,1),PSL(1,1)/PSI(1,1),CHI(1,1),&
           &ASCAL(1,1)**2*(BPHI(1,1)*BPHI(1,1)*R(IZ)**2.*SIN(TH(IX))**2.*PSI(1,1)**4. + &
           &BPOLR(1,1)*BPOLR(1,1)*PSI(1,1)**4.+BPOLT(1,1)*BPOLT(1,1)*R(IZ)**2*PSI(1,1)**4.)
  CLOSE(12)

  IF(.NOT. MPICODE)THEN
	  WRITE(6,*)'METRIC CONVERGENCE',PSL(50,1),PSI(50,1)
  ENDIF
  ! The convergence is oscillatory so one needs to find its bounds
  ! Check convergence on the central density
  !WARNING: in some cases of damped relaxation (i.e. QAPHI.NE.1 and QFACTOR.NE.1)
  !        the solution may require additional interative steps with QAPHI and QFACTOR
  !        manually set to 1. In this case, increase the number in the following IF statement.
  IF(ILOOP .GT. 5)THEN
!     IF((RHOCVEC(ILOOP-1)-RHOCVEC(ILOOP-2))*(RHOCVEC(ILOOP)-RHOCVEC(ILOOP-1)) .LT. 0)THEN 			! if the solution is oscillating
!       IF(RHOCVEC(ILOOP-1) .GT. RHOCVEC(ILOOP)) RHOCMAX= RHOCVEC(ILOOP-1)
!       IF(RHOCVEC(ILOOP-1) .LT. RHOCVEC(ILOOP)) RHOCMIN= RHOCVEC(ILOOP-1)
!       IF((.NOT. MPICODE))WRITE(6,*)'Oscillatory convergence check: reached',ABS(RHOCVEC(ILOOP)-RHOCVEC(ILOOP-1)),',needed',OSCCONV
!    	IF((.NOT. MPICODE))WRITE(6,*)'METRIC CONVERGENCE',PSL(1,1),PSI(1,1)
!     ELSE  ! if the solution is not oscillating
!       IF((.NOT. MPICODE))WRITE(6,*)'Monotonous convergence check: reached',ABS(RHOCVEC(ILOOP)-RHOCVEC(ILOOP-1)),&
!       	&',needed',MONCONV
		!IF((.NOT. MPICODE))WRITE(6,*)'METRIC CONVERGENCE',PSL(50,1),PSI(50,1)
!     END IF

!     IF(ABS(RHOCMAX-RHOCMIN).LT. OSCCONV) THEN 					 ! convergence if the solution is oscillating
!       IF(WRTF.AND.IDAT) WRT=.TRUE.
!       EXT=.TRUE.
!     END IF
! 	WRITE(6,*)'CONV: ',ABS(PSI(50,1)-PSI0(50,1))/PSI0(50,1)
    IF(ABS(PSI(50,1)-PSI0(50,1))/PSI0(50,1) .LT. MONCONV) THEN     ! convergence if the solution is not oscillating
    	NCONV=NCONV+1
    ELSE
    	NCONV=0
    ENDIF
    IF(NCONV==5)THEN ! Checks how many consecutive times the solution hasn't changed
      IF(WRTF.AND.IDAT) WRT=.TRUE.
      EXT=.TRUE.
    END IF
  ENDIF
  IF(ILOOP .EQ. MAXLOOP)THEN
  	WRITE(6,*)'WARNING: convergence to desired level has not been obtained.'
  ENDIF

  ! IF(ILOOP .EQ. MAXLOOP)THEN
!   	ALLOCATE(RHOCVECTEMP(3*MAXLOOP))
!   	RHOCVECTEMP(1:MAXLOOP)=RHOCVEC
!     DEALLOCATE(RHOCVEC)
!     CALL MOVE_ALLOC(RHOCVECTEMP,RHOCVEC)
!   ENDIF


  IF(ILOOP .EQ. MAXLOOP)EXIT
  ILOOP=ILOOP+1
  !End convergence loop
END DO

  ! ------------------------------------------------------------------------------
  ! ------- CALCULATE STELLAR QUANTITIES  ----------------------------------------
  ! ------------------------------------------------------------------------------

  OPEN(13,FILE=trim(adjustl(subdirpath))//'Surf.dat')
  DO IX=1,NTH
  	WRITE(13,*) R(WSURF(IX)+1)
  END DO
  CLOSE(13)

  CALL QUANTITIES (RNEWT,RHOVAR, MM, M0, ILOOP)

#ifdef NWTRPS
  IF(QUOC .EQ. 0) QUCNV=RNEWT
  IF(QUOC .EQ. 1) QUCNV=MM
  IF(QUOC .EQ. 2) QUCNV=M0
#endif

200 END SUBROUTINE XNSMAIN


! ********************************************************
! ********************************************************

SUBROUTINE CONFORMAL(RHOS,RHOC) !!! QUI
!   ============================================================
!   Purpose : This subroutine solve the axisymmetric eq. for Conformal
!             factor in spherical coorinates, by Legendre Poly in Theta
!             direct solution in R (solve for conf-fact -1.)
!   ============================================================
!
! Parameters for the grid - set in module system
!
! R = radial grid points (+ boundaries); DR = increments
! TH = angular grid points (+ boundaries)
! XX = agular cos(th) points (+ boundaries)
!
! Source - in input
! RHOS = source (total energy in conservative form)
! RHOC = source (curvature term)
!
! Conformal factor in output - set in module system
! PSI = potential on the original grid -1.
!
! Parameters of Legendre series
!
! NGQ = number of point in the Gauss quadrature
! MLS = number of Legendre polinomial for expansion in theta
! XGQ = gauss quadrature points
! WGQ = gauss quadrature weights
! Y   = function interpolated in gauss quadrature points
! PN  = Legendre polyn da 0-MLS evaluated in gauss quadrature points
! PD  = Deriv. Legendre polyn da 0-MLS evaluated in gauss quadrature points
! COEF  = Coefficient of the legendre expansion series of the function
! TAB = radial array of coeff of legendre expansion
!
! Parameters of radial matrix inversion
!
! IL = order of the polyn.
! DC = main diagonal
! DU = upper diagonal
! DL = lower diagonal
! PHI = potential in series of Legendre polyn.
!
!   ============================================================

  USE SYSTEMXNS
  USE FUNCTIONS

  INTEGER,PARAMETER :: MAXIT = 200 !25
  INTEGER :: I
  REAL,DIMENSION(-1:NTH+2,0:NR+1):: RHOS,RHO,RHOC

  ! Array for the Legendre expansion
  REAL,DIMENSION(NGQ) :: XGQ,WGQ,Y
  REAL,DIMENSION(0:MLS,NGQ) :: PN,PD
  REAL,DIMENSION(0:MLS) :: COEF,PDX,PNX
  REAL,DIMENSION(0:NR+1,0:MLS) :: TAB

  ! Array for the radial solution
  INTEGER :: IL,ITER,INFO
  REAL :: A1,A2,A3,A4,A5,A6,DCI,DCF,ERROR
  REAL,DIMENSION(NR) :: DC,DC1,DCP,SS1
  REAL,DIMENSION(NR-1) :: DL1,DU1,DLP,DUP
  REAL,DIMENSION(1:NR,0:MLS) :: PHI
  REAL,DIMENSION(1:NTH,1:NR) :: PSIOLD

  !....................
  ! Loop for convergence on the non linear term
  ! ...................

  ! Initialize the potential -1.
  DO IX=1,NTH
    DO IZ=1,NR
!       PSI(IX,IZ)=0.
        PSI(IX,IZ)=PSI(IX,IZ)-1.
    END DO
  END DO

  ! Iterate on the potential
  DO ITER=1,MAXIT

  ! Initialize the nonlinear source in the domain
  DO IX=1,NTH
    DO IZ=1,NR
      RHO(IX,IZ) = RHOS(IX,IZ)/(1.+PSI(IX,IZ)) + RHOC(IX,IZ)/(1.+PSI(IX,IZ))**7.
    END DO
  END DO

  ! ...................
  ! Solve the problem in theta with series expansion in legendre polynomials
  ! This has to be repeated at all R
  ! ...................

  ! Compute the Gauss-quadrature points XGQ and weights WGQ
  CALL LEGZO(NGQ,XGQ,WGQ)

  ! Compute the matrix of Weighted Legendre Polyn in the quadrature points
  DO I=1,NGQ
    CALL LPN(MLS,XGQ(I),PN(0:MLS,I),PD(0:MLS,I))
    PN(0:MLS,I)=PN(0:MLS,I)*WGQ(I)
  END DO

  ! Fill boundary values for rho in theta at fixed r
  DO IZ=1,NR
    DO I=1,2
      RHO(1-I,IZ)=RHO(I,IZ)
      RHO(NTH+I,IZ)=RHO(NTH+1-I,IZ)
    END DO

    ! Interpolate the disctretized function on the quadrature points
    CALL POLINT(XX(-1:NTH+2),RHO(-1:NTH+2,IZ),NTH+4,XGQ,Y,NGQ)

    ! Evaluate the coefficient of the Legendre poly. expansion
    DO I=0,MLS
      COEF(I)=(2.*I+1.)/2.*DOT_PRODUCT(PN(I,1:NGQ),Y(1:NGQ))
    END DO

    ! Compute the matrix of legendre coefficients
    TAB(IZ,0:MLS)=COEF(0:MLS)
  END DO

  ! ...................
  ! Solve the problem in r for each Legendre polynomial
  ! This must be loopen on all the Legendre coeff rows
  ! Use direct inversion of the tridiag matrix
  ! ...................

  ! Set the matrix elements (case L=0)
  ! Other cases need to add a term to the main diagonal and modify inner and outer
  ! BC for the first and last diag element

  I=1
  A1= 2./DR(I)/(DR(I)+DR(I+1))
  A2= 2./DR(I+1)/(DR(I)+DR(I+1))
  A3=-2./DR(I)/DR(I+1)
  A4=-DR(I+1)/DR(I)/(DR(I)+DR(I+1))
  A5= DR(I)/DR(I+1)/(DR(I)+DR(I+1))
  A6= (DR(I+1)-DR(I))/DR(I)/DR(I+1)

  DCI= A3+2.*A6/R(I)
  DUP(I)= A2+2.*A5/R(I)
  DO I=2,NR-1
    A1= 2./DR(I)/(DR(I)+DR(I+1))
    A2= 2./DR(I+1)/(DR(I)+DR(I+1))
    A3=-2./DR(I)/DR(I+1)
    A4=-DR(I+1)/DR(I)/(DR(I)+DR(I+1))
    A5= DR(I)/DR(I+1)/(DR(I)+DR(I+1))
    A6= (DR(I+1)-DR(I))/DR(I)/DR(I+1)
    DCP(I)  = A3+2.*A6/R(I)
    DLP(I-1)  = A1+2.*A4/R(I)
    DUP(I)    = A2+2.*A5/R(I)
  END DO
  I=NR
  A1= 2./DR(I)/(DR(I)+DR(I+1))
  A2= 2./DR(I+1)/(DR(I)+DR(I+1))
  A3=-2./DR(I)/DR(I+1)
  A4=-DR(I+1)/DR(I)/(DR(I)+DR(I+1))
  A5= DR(I)/DR(I+1)/(DR(I)+DR(I+1))
  A6= (DR(I+1)-DR(I))/DR(I)/DR(I+1)
  DLP(I-1)= A1+2.*A4/R(I)
  DCF  = A3+2.*A6/R(I)

  DO IL=0,MLS
    ! BC at inner radius ( Origin )
    ! For L .ne. 0 the function goes to 0, parity depends on L
    ! We assume odd parity because for NS  L=odd term are 0.
    A1= 2./DR(1)/(DR(1)+DR(1+1))
    A4=-DR(1+1)/DR(1)/(DR(1)+DR(1+1))
    IF(IL .EQ.0)THEN
      DCP(1)=DCI+(A1+2.*A4/R(1))
    ELSE
      DCP(1)=DCI-(A1+2.*A4/R(1))
    END IF

    ! BC at outer radius
    ! Various multipole must dacay as R^(L+1)
    ! The monopole must decay to 1
    A2= 2./DR(NR+1)/(DR(NR)+DR(NR+1))
    A5= DR(NR)/DR(NR+1)/(DR(NR)+DR(NR+1))

    DCP(NR)=DCF+(A2+2.*A5/R(NR))*(R(NR)/(R(NR)+DR(NR+1)))**(IL+1.)

    ! Set main diagonal

    DO I=1,NR
      DC(I)  = DCP(I) - IL*(IL+1.)/R(I)**2.
    END DO

    SS1(1:NR)=TAB(1:NR,IL)

    DL1=DLP
    DC1=DC
    DU1=DUP

    CALL DGTSV(NR,DL1,DC1,DU1,SS1,NR,INFO)

    PHI(1:NR,IL)=SS1(1:NR)
    !IF(IL==0)THEN
    !  WRITE(6,*)'INTERNAL :',PHI(NR,IL)*2*R(NR)
    !  WRITE(6,*)'INTERNAL :',R(NR-2:NR),DR(NR-1:NR+1)
    !ENDIF
  END DO
!STOP
  ! Compute the matrix of Legendre Polyn in the grid points
  DO IX=1,NTH
    CALL LPN(MLS,XX(IX),PNX(0:MLS),PDX(0:MLS))
    ! Compute the potential on the grid
    DO IZ=1,NR
      PSIOLD(IX,IZ)=PSI(IX,IZ)
      PSI(IX,IZ)=DOT_PRODUCT(PNX(0:MLS),PHI(IZ,0:MLS))
      PSI(IX,IZ)=QFACTORMETRIC*PSI(IX,IZ) + (1.-QFACTORMETRIC)*PSIOLD(IX,IZ)
    END DO
  END DO

  ERROR=MAXVAL(ABS(PSIOLD(1:NTH,1:NR)-PSI(1:NTH,1:NR)))
!   write(6,*)psi(1,1)
  IF(VERBOSE .AND. (.NOT. MPICODE)) WRITE(6,*)'Max err Conf = ',error
  IF(ERROR .LT. TOLCONV) EXIT
  END DO

END SUBROUTINE  CONFORMAL


! ********************************************************
! ********************************************************

SUBROUTINE LAPSE(RHOS,RHOC) !!! QUI

!   ============================================================
!   Purpose : This subroutine solve the axisymmetric eq. for Lapse
!             factor in spherical coorinates, by Legendre Poly in Theta
!             direct solution in R (solve for lapse -1.)
!   ============================================================
!
! Parameters for the grid - set in module system
!
! R = radial grid points (+ boundaries); DR = increments
! TH = angular grid points (+ boundaries)
! XX = agular cos(th) points (+ boundaries)
!
! Source - in input
! RHOS = source (total energy in conservative form)
! RHOC = source (curvature term)
!
! Conformal factor in output - set in module system
! PSL = lapse on the original grid -1.
!
! Parameters of Legendre series
!
! NGQ = number of point in the Gauss quadrature
! MLS = number of Legendre polinomial for expansion in theta
! XGQ = gauss quadrature points
! WGQ = gauss quadrature weights
! Y   = function interpolated in gauss quadrature points
! PN  = Legendre polyn da 0-MLS evaluated in gauss quadrature points
! PD  = Deriv. Legendre polyn da 0-MLS evaluated in gauss quadrature points
! COEF  = Coefficient of the legendre expansion series of the function
! TAB = radial array of coeff of legendre expansion
!
! Parameters of radial matrix inversion
!
! IL = order of the polyn.
! DC = main diagonal
! DU = upper diagonal
! DL = lower diagonal
! PHI = potential in series of Legendre polyn.
!
!   ============================================================

  USE SYSTEMXNS
  USE FUNCTIONS

  INTEGER,PARAMETER :: MAXIT = 200!25
  INTEGER :: I,INFO
  REAL,DIMENSION(-1:NTH+2,0:NR+1):: RHOS,RHO,RHOC

  ! Array for the Legendre expansion
  REAL,DIMENSION(NGQ) :: XGQ,WGQ,Y
  REAL,DIMENSION(0:MLS,NGQ) :: PN,PD
  REAL,DIMENSION(0:MLS) :: COEF,PDX,PNX
  REAL,DIMENSION(0:NR+1,0:MLS) :: TAB

  ! Array for the radial solution
  INTEGER :: IL,ITER
  REAL :: A1,A2,A3,A4,A5,A6,DCI,DCF,ERROR
  REAL,DIMENSION(NR) :: DC,DC1,DCP,SS1
  REAL,DIMENSION(NR-1) :: DL1,DU1,DLP,DUP
  REAL,DIMENSION(1:NR,0:MLS) :: PHI
  REAL,DIMENSION(1:NTH,1:NR) :: PSIOLD

  !....................
  ! Loop for convergence on the non linear term
  ! ...................
  ! Initialize the lapse -1.
  DO IX=1,NTH
    DO IZ=1,NR
       PSL(IX,IZ)=0.
!      PSL(IX,IZ)=PSL(IX,IZ)-1.
    END DO
  END DO

  ! Iterate on the potential
  DO ITER=1,MAXIT
  ! Initialize the nonlinear source in the domain
  DO IX=1,NTH
    DO IZ=1,NR
      RHO(IX,IZ) = RHOS(IX,IZ)*(1.+PSL(IX,IZ)) - 7.*RHOC(IX,IZ)*(1.+PSL(IX,IZ))/PSI(IX,IZ)**8
    END DO
  END DO

!print*, 'OOOOOOKKKKK -- 3'

!   write(6,*)'ciao'

  ! ...................
  ! Solve the problem in theta with series expansion in legendre polynomials
  ! This has to be repeated at all R
  ! ...................

  ! Compute the Gauss-quadrature points XGQ and weights WGQ
  CALL LEGZO(NGQ,XGQ,WGQ)

!PRINT*, 'OOOOOOOOOOOOOOKKKKKKKKKK -- 4'

  ! Compute the matrix of Weighted Legendre Polyn in the quadrature points
  DO I=1,NGQ
    CALL LPN(MLS,XGQ(I),PN(0:MLS,I),PD(0:MLS,I))
    PN(0:MLS,I)=PN(0:MLS,I)*WGQ(I)
  END DO

!PRINT*, 'OOOOOOOOOKKKKKK -- 5'

  DO IZ=1,NR
    ! Fill boundary values for rho in theta at fixed r
    DO I=1,2
      RHO(1-I,IZ)=RHO(I,IZ)
      RHO(NTH+I,IZ)=RHO(NTH+1-I,IZ)
    END DO

!PRINT*, 'OOOOOOOOOOOOKKKKKKKK -- 6'

    ! Interpolate the disctretized function on the quadrature points
    CALL POLINT(XX(-1:NTH+2),RHO(-1:NTH+2,IZ),NTH+4,XGQ,Y,NGQ)

!PRINT*, 'OOOOOOOOOOOOOOOOKKKKKKKKKK -- 7'

    ! Evaluate the coefficient of the Legendre poly. expansion
    DO I=0,MLS
      COEF(I)=(2.*I+1.)/2.*DOT_PRODUCT(PN(I,1:NGQ),Y(1:NGQ))
      !if(abs(dot_product(PN(I,1:NGQ),Y(1:NGQ)))>HUGE(dbl_prec_var)) then
       !      print*, 'DOT PRODUCT'
        !     do ix=1,NGQ
         !    if(abs(y(ix))>HUGE(dbl_prec_var)) then
          !           print*,'Y',Y(ix),ix,i
          !   endif
          !   enddo
          !   stop
     !endif
      !if(isnan(coef(i))) then
       !       print*,'COEFFFFF',I
        !      stop
      !endif
    END DO

    ! Compute the matrix of legendre coefficients
    TAB(IZ,0:MLS)=COEF(0:MLS)
  END DO

  ! ...................
  ! Solve the problem in r for each Legendre polynomial
  ! This must be loopen on all the Legendre coeff rows
  ! Use direct inversion of the tridiag matrix
  ! ...................

  ! Set the matrix elements (case L=0)
  ! Other cases need to add a term to the main diagonal and modify inner and outer
  ! BC for the first and last diag element

  I=1
  A1= 2./DR(I)/(DR(I)+DR(I+1))
  A2= 2./DR(I+1)/(DR(I)+DR(I+1))
  A3=-2./DR(I)/DR(I+1)
  A4=-DR(I+1)/DR(I)/(DR(I)+DR(I+1))
  A5= DR(I)/DR(I+1)/(DR(I)+DR(I+1))
  A6= (DR(I+1)-DR(I))/DR(I)/DR(I+1)

  DCI   = A3+2.*A6/R(I)
  DUP(I)= A2+2.*A5/R(I)

  DO I=2,NR-1
    A1= 2./DR(I)/(DR(I)+DR(I+1))
    A2= 2./DR(I+1)/(DR(I)+DR(I+1))
    A3=-2./DR(I)/DR(I+1)
    A4=-DR(I+1)/DR(I)/(DR(I)+DR(I+1))
    A5= DR(I)/DR(I+1)/(DR(I)+DR(I+1))
    A6= (DR(I+1)-DR(I))/DR(I)/DR(I+1)
    DCP(I)   = A3+2.*A6/R(I)
    DLP(I-1) = A1+2.*A4/R(I)
    DUP(I)   = A2+2.*A5/R(I)
  END DO

  I=NR
  A1= 2./DR(I)/(DR(I)+DR(I+1))
  A2= 2./DR(I+1)/(DR(I)+DR(I+1))
  A3=-2./DR(I)/DR(I+1)
  A4=-DR(I+1)/DR(I)/(DR(I)+DR(I+1))
  A5= DR(I)/DR(I+1)/(DR(I)+DR(I+1))
  A6= (DR(I+1)-DR(I))/DR(I)/DR(I+1)
  DLP(I-1) = A1+2.*A4/R(I)
  DCF      = A3+2.*A6/R(I)

  DO IL=0,MLS
    ! BC at inner radius ( Origin )
    ! For L .ne. 0 the function goes to 0, parity depends on L
    ! We assume odd parity because for NS  L=odd term are 0.
    A1= 2./DR(1)/(DR(1)+DR(1+1))
    A4=-DR(1+1)/DR(1)/(DR(1)+DR(1+1))
    IF(IL .EQ.0)THEN
      DCP(1)=DCI+(A1+2.*A4/R(1))
    ELSE
      DCP(1)=DCI-(A1+2.*A4/R(1))
    END IF

    ! BC at outer radius
    ! Various multipole must dacay as R^(L+1)
    ! The monopole must decay to 1
    A2= 2./DR(NR+1)/(DR(NR)+DR(NR+1))
    A5= DR(NR)/DR(NR+1)/(DR(NR)+DR(NR+1))

    DCP(NR)=DCF+(A2+2.*A5/R(NR))*(R(NR)/(R(NR)+DR(NR+1)))**(IL+1.)

    ! Set main diagonal
    DO I=1,NR
      DC(I)  = DCP(I) - IL*(IL+1.)/R(I)**2.
    END DO

    SS1(1:NR)=TAB(1:NR,IL)

    DL1=DLP
    DC1=DC
    DU1=DUP

    CALL DGTSV(NR,DL1,DC1,DU1,SS1,NR,INFO)

    PHI(1:NR,IL)=SS1(1:NR)

  END DO

  ! Compute the matrix of Legendre Polyn in the grid points
  DO IX=1,NTH
    CALL LPN(MLS,XX(IX),PNX(0:MLS),PDX(0:MLS))
    ! Compute the potential on the grid
    DO IZ=1,NR
      PSIOLD(IX,IZ)=PSL(IX,IZ)
      PSL(IX,IZ)=DOT_PRODUCT(PNX(0:MLS),PHI(IZ,0:MLS))
      PSL(IX,IZ)=QFACTORMETRIC*PSL(IX,IZ) + (1.-QFACTORMETRIC)*PSIOLD(IX,IZ)
      if(isnan(phi(ix,iz))) stop 'PHIIII'
    END DO
  END DO
!   WRITE(6,*)'LAPSE ',1.+PSIOLD(1,1),1.+PSL(1,1),RHOS(1,1),RHOC(1,1)
  ERROR=MAXVAL(ABS(PSIOLD(1:NTH,1:NR)-PSL(1:NTH,1:NR)))
  IF(VERBOSE .AND. (.NOT. MPICODE)) WRITE(6,*)'Max err Lapse= ',error
  IF(ERROR .LT. TOLCONV)EXIT
  END DO

END SUBROUTINE LAPSE

! ********************************************************
! ********************************************************

SUBROUTINE SHIFTPHI(RHOS)!!! QUI

!   ============================================================
!   Purpose : This subroutine solve the axisymmetric eq. for the
!    phi- component of the vector poisson equation, in spherical
!    coordinates using coordinate quantities
!   ============================================================
!
! Parameters for the grid - set in module system
!
! R = radial grid points (+ boundaries); DR = increments
! TH = angular grid points (+ boundaries)
! XX = agular cos(th) points (+ boundaries)
!
! Source - in input
! RHOS = source (total energy in conservative form)
!
! Shift vector in output - set in module system
! PSS = shift on the original grid (in coord quant)
!
! Parameters of Legendre series
!
! NGQ = number of point in the Gauss quadrature
! MLS = number of Legendre polinomial for expansion in theta
! XGQ = gauss quadrature points
! WGQ = gauss quadrature weights
! Y   = function interpolated in gauss quadrature points
! PN  = Legendre polyn da 0-MLS evaluated in gauss quadrature points
! PD  = Deriv. Legendre polyn da 0-MLS evaluated in gauss quadrature points
! COEF  = Coefficient of the legendre expansion series of the function
! TAB = radial array of coeff of legendre expansion
!
! Parameters of radial matrix inversion
!
! IL = order of the polyn.
! DC = main diagonal
! DU = upper diagonal
! DL = lower diagonal
! PHI = potential in series of Legendre polyn.
!
!   ============================================================

  USE SYSTEMXNS
  USE FUNCTIONS

  INTEGER,PARAMETER :: MAXIT = 20
  INTEGER :: I,INFO
  REAL,DIMENSION(-1:NTH+2,0:NR+1):: RHOS,RHO

  ! Array for the Legendre expansion
  REAL,DIMENSION(NGQ) :: XGQ,WGQ,Y
  REAL,DIMENSION(0:MLS,NGQ) :: PN,PD
  REAL,DIMENSION(0:MLS) :: COEF,PDX,PNX
  REAL,DIMENSION(0:NR+1,0:MLS) :: TAB

  ! Array for the radial solution
  INTEGER :: IL
  REAL :: A1,A2,A3,A4,A5,A6,DCI,DCF,TMP
  REAL,DIMENSION(NR) :: DC,DC1,DCP,SS1
  REAL,DIMENSION(NR-1) :: DL1,DU1,DLP,DUP
  REAL,DIMENSION(1:NR,0:MLS) :: PHI

  !....................
  ! Solve first the separate equation for the phi component
  ! ...................

  ! Initialize the source in the domain - initialize to coord quantities
  DO IX=1,NTH
    DO IZ=1,NR
      RHO(IX,IZ) = RHOS(IX,IZ)*R(IZ)*SIN(TH(IX))
    END DO
  END DO

  ! ...................
  ! Solve the problem in theta with series expansion in legendre polynomials
  ! This has to be repeated at all R
  ! ...................

  ! Compute the Gauss-quadrature points XGQ and weights WGQ
  CALL LEGZO(NGQ,XGQ,WGQ)

  ! Compute the matrix of Weighted Spherical Harmonics Deriv in the quadrature points
  DO I=1,NGQ
    CALL LPN(MLS,XGQ(I),PN(0:MLS,I),PD(0:MLS,I))
    PN(0:MLS,I)=-PD(0:MLS,I)*WGQ(I)*SQRT(1.-XGQ(I)**2)
  END DO

  COEF(0:MLS)=0.
  DO IZ=1,NR
    ! Fill boundary values for rho in theta at fixed r
    DO I=1,2
      RHO(1-I,IZ)=-RHO(I,IZ)
      RHO(NTH+I,IZ)=-RHO(NTH+1-I,IZ)
    END DO

    ! Interpolate the disctretized function on the quadrature points
    CALL POLINT(XX(-1:NTH+2),RHO(-1:NTH+2,IZ),NTH+4,XGQ,Y,NGQ)

    ! Evaluate the coefficient of the Spheric Harm. expansion (L=0 => 0)
    COEF(0)=0.
    TAB(IZ,0)=0.
    DO I=1,MLS
      TMP=SQRT(Pi*(2.*I+1.))*DOT_PRODUCT(PN(I,1:NGQ),Y(1:NGQ))/(I*(I+1.))
      COEF(I)=TMP
    END DO

    ! Compute the matrix of legendre coefficients
    TAB(IZ,0:MLS)=COEF(0:MLS)
  END DO

  ! ...................
  ! Solve the problem in r for each Legendre polynomial
  ! This must be looped on all the Legendre coeff rows
  ! Use direct inversion of the tridiag matrix
  ! ...................

  ! Set the matrix elements (case L=0)
  ! Other cases need to add a term to the main diagonal and modify inner and outer
  ! BC for the first and last diag element

  I=1
  A1= 2./DR(I)/(DR(I)+DR(I+1))
  A2= 2./DR(I+1)/(DR(I)+DR(I+1))
  A3=-2./DR(I)/DR(I+1)
  A4=-DR(I+1)/DR(I)/(DR(I)+DR(I+1))
  A5= DR(I)/DR(I+1)/(DR(I)+DR(I+1))
  A6= (DR(I+1)-DR(I))/DR(I)/DR(I+1)

  DCI= A3+2.*A6/R(I)
  DUP(I)= A2+2.*A5/R(I)
  DO I=2,NR-1
    A1= 2./DR(I)/(DR(I)+DR(I+1))
    A2= 2./DR(I+1)/(DR(I)+DR(I+1))
    A3=-2./DR(I)/DR(I+1)
    A4=-DR(I+1)/DR(I)/(DR(I)+DR(I+1))
    A5= DR(I)/DR(I+1)/(DR(I)+DR(I+1))
    A6= (DR(I+1)-DR(I))/DR(I)/DR(I+1)
    DCP(I)  = A3+2.*A6/R(I)
    DLP(I-1)  = A1+2.*A4/R(I)
    DUP(I)    = A2+2.*A5/R(I)
  END DO
  I=NR
  A1= 2./DR(I)/(DR(I)+DR(I+1))
  A2= 2./DR(I+1)/(DR(I)+DR(I+1))
  A3=-2./DR(I)/DR(I+1)
  A4=-DR(I+1)/DR(I)/(DR(I)+DR(I+1))
  A5= DR(I)/DR(I+1)/(DR(I)+DR(I+1))
  A6= (DR(I+1)-DR(I))/DR(I)/DR(I+1)
  DLP(I-1)= A1+2.*A4/R(I)
  DCF  = A3+2.*A6/R(I)

  ! Parity imposes that there is not L=0 term (we work on coord quantity)
  DO IL=1,MLS
    ! BC at inner radius ( Origin )
    ! The function goes to 0, parity depends on L
    ! We assume odd parity because for NS  L=odd term are 0.
    A1= 2./DR(1)/(DR(1)+DR(1+1))
    A4=-DR(1+1)/DR(1)/(DR(1)+DR(1+1))

    DCP(1)=DCI-(A1+2.*A4/R(1))

    ! BC at outer radius
    ! Various multipole must dacay as R^(L+1)
    A2= 2./DR(NR+1)/(DR(NR)+DR(NR+1))
    A5= DR(NR)/DR(NR+1)/(DR(NR)+DR(NR+1))

    DCP(NR)=DCF+(A2+2.*A5/R(NR))*(R(NR)/(R(NR)+DR(NR+1)))**(IL+1.)

    ! Set main diagonal
    DO I=1,NR
      DC(I)  = DCP(I) - IL*(IL+1.)/R(I)**2.
    END DO

    SS1(1:NR)=TAB(1:NR,IL)

    DL1=DLP
    DC1=DC
    DU1=DUP

    CALL DGTSV(NR,DL1,DC1,DU1,SS1,NR,INFO)

    PHI(1:NR,IL)=SS1(1:NR)

  END DO

  ! Compute the matrix of Legendre Polyn in the grid points
  DO IX=1,NTH
    CALL LPN(MLS,XX(IX),PNX(0:MLS),PDX(0:MLS))
    DO I=1,MLS
      PDX(I)=SQRT((2*I+1.)/4./Pi)*PDX(I)
    END DO
    ! Compute the potential on the grid
    DO IZ=1,NR
      PSS(IX,IZ)=-1.*DOT_PRODUCT(PDX(1:MLS)*SIN(TH(IX)),PHI(IZ,1:MLS))
    END DO
  END DO

END SUBROUTINE  SHIFTPHI

! ********************************************************
! ********************************************************

SUBROUTINE CURV1 !!! QUI

!   ============================================================
!   Purpose : This subroutine computes the term A^ij * A_ij that
!             is used as curvature source term in the routines
!             for the solution of conformal factor and lapse
!   ============================================================
!
! Parameters for the grid - set in module system
!
! R = radial grid points (+ boundaries); DR = increments
! TH = angular grid points (+ boundaries)
! XX = agular cos(th) points (+ boundaries)
!
! Output
! CURVC = curvature term corresponding to  A^ij * A_ij
!
!   ============================================================

  USE SYSTEMXNS
  USE FUNCTIONS

  REAL,DIMENSION(0:NTH+1,0:NR+1) :: XR
  REAL,DIMENSION(0:NTH+1,0:NR+1) :: XT
  REAL,DIMENSION(0:NTH+1,0:NR+1) :: XP

  REAL :: A1,A2,A3,B1,B2,B3
  REAL :: DXRDR,DXRDT,DXTDR,DXTDT,DXPDR,DXPDT
  REAL :: SIN2,CTG


  ! Component of the X-vector (contravariant)
  DO IX=1,NTH
    DO IZ=1,NR
      XR(IX,IZ)=ES1SRC(IX,IZ)
      XT(IX,IZ)=ES2SRC(IX,IZ)
      XP(IX,IZ)=ES3SRC(IX,IZ)
    END DO
  END DO

  ! Assume parity on axis
  DO IZ=1,NR
    XR(0,IZ) =  XR(1,IZ)
    XT(0,IZ) = -XT(1,IZ)
    XP(0,IZ) =  XP(1,IZ)
    XR(NTH+1,IZ) =  XR(NTH,IZ)
    XT(NTH+1,IZ) = -XT(NTH,IZ)
    XP(NTH+1,IZ) =  XP(NTH,IZ)
  END DO

  ! Assume parity at center
  ! Assume that the shift is smooth at the outer boundary.
  DO IX=1,NTH
   !XR(IX,0) = -XR(IX,1)
   !XT(IX,0) =  XT(IX,1)
   !XP(IX,0) =  XP(IX,1)
    XR(IX,0) = -XR(NTH-IX+1,1)
    XT(IX,0) =  XT(NTH-IX+1,1)
    XP(IX,0) =  XP(NTH-IX+1,1)
    XR(IX,NR+1) =  XR(IX,NR) + DR(NR)/DR(NR-1)*(XR(IX,NR)-XR(IX,NR-1))
    XT(IX,NR+1) =  XT(IX,NR) + DR(NR)/DR(NR-1)*(XT(IX,NR)-XT(IX,NR-1))
    XP(IX,NR+1) =  XP(IX,NR) + DR(NR)/DR(NR-1)*(XP(IX,NR)-XP(IX,NR-1))
  END DO

 ! Compute the various derivative terms
  DO IX=1,NTH
    DO IZ=1,NR
      A1=-DR(IZ+1)/DR(IZ)/(DR(IZ)+DR(IZ+1))
      A2= (DR(IZ+1)-DR(IZ))/DR(IZ)/DR(IZ+1)
      A3= DR(IZ)/DR(IZ+1)/(DR(IZ)+DR(IZ+1))
      B1=-DTH(IX+1)/DTH(IX)/(DTH(IX)+DTH(IX+1))
      B2= (DTH(IX+1)-DTH(IX))/DTH(IX)/DTH(IX+1)
      B3= DTH(IX)/DTH(IX+1)/(DTH(IX)+DTH(IX+1))

      DXRDR = A1*XR(IX,IZ-1)+A3*XR(IX,IZ+1)+A2*XR(IX,IZ)
      DXTDR = A1*XT(IX,IZ-1)+A3*XT(IX,IZ+1)+A2*XT(IX,IZ)
      DXPDR = A1*XP(IX,IZ-1)+A3*XP(IX,IZ+1)+A2*XP(IX,IZ)
      DXRDT = B1*XR(IX-1,IZ)+B3*XR(IX+1,IZ)+B2*XR(IX,IZ)
      DXTDT = B1*XT(IX-1,IZ)+B3*XT(IX+1,IZ)+B2*XT(IX,IZ)
      DXPDT = B1*XP(IX-1,IZ)+B3*XP(IX+1,IZ)+B2*XP(IX,IZ)

      CTG=COS(TH(IX))/SIN(TH(IX))
      SIN2=SIN(TH(IX))**2

      CURVC(IX,IZ) = -1./8.*2./3.*(4.*(DXRDR**2+DXTDT**2+(XR(IX,IZ)/R(IZ))**2+   &
           (CTG*XT(IX,IZ))**2-2*XR(IX,IZ)*DXRDR/R(IZ) &
           +XR(IX,IZ)/R(IZ)*DXTDT-DXRDR*DXTDT-                            &
           CTG*XT(IX,IZ)*(DXTDT+DXRDR-XR(IX,IZ)/R(IZ))) + &
           3.*(R(IZ)**2*SIN2*DXPDR**2+SIN2*DXPDT**2   &
           +R(IZ)**2*DXTDR**2+DXRDT**2/R(IZ)**2+2.*DXTDR*DXRDT))
     END DO
  END DO

END SUBROUTINE CURV1

! ********************************************************
! ********************************************************

SUBROUTINE CURV2 !!! QUI

!   ============================================================
!   Purpose : This subroutine computes the term A^ij * A_ij that
!             is used as curvature source term in the routines
!             for the solution of conformal factor and lapse
!   ============================================================
!
! Parameters for the grid - set in module system
!
! R = radial grid points (+ boundaries); DR = increments
! TH = angular grid points (+ boundaries)
! XX = agular cos(th) points (+ boundaries)
!
! Output
! CURV1 = curvature term corresponding to  A^ij * A_ij
!
!   ============================================================

  USE SYSTEMXNS
  USE FUNCTIONS

  REAL,DIMENSION(0:NTH+1,0:NR+1) :: XR
  REAL,DIMENSION(0:NTH+1,0:NR+1) :: XT
  REAL,DIMENSION(0:NTH+1,0:NR+1) :: XP
  REAL,DIMENSION(0:NTH+1,0:NR+1) :: XC

  REAL :: A1,A2,A3,B1,B2,B3
  REAL :: DXRDR,DXRDT,DXTDR,DXTDT,DXPDR,DXPDT,DXCDR,DXCDT
  REAL :: SIN2
  REAL :: ARR,ATT,APP,ART,ARP,ATP,CTG

  DO IX=1,NTH
    DO IZ=1,NR
      XR(IX,IZ)=ES1SRC(IX,IZ)
      XT(IX,IZ)=ES2SRC(IX,IZ)
      XP(IX,IZ)=ES3SRC(IX,IZ)
      XC(IX,IZ)=2.*PSL(IX,IZ)/PSI(IX,IZ)**7
    END DO
  END DO

  ! Assume parity on axis
  DO IZ=1,NR
    XR(0,IZ) =  XR(1,IZ)
    XT(0,IZ) = -XT(1,IZ)
    XP(0,IZ) =  XP(1,IZ)
    XC(0,IZ) =  XC(1,IZ)
    XR(NTH+1,IZ) =  XR(NTH,IZ)
    XT(NTH+1,IZ) = -XT(NTH,IZ)
    XP(NTH+1,IZ) =  XP(NTH,IZ)
    XC(NTH+1,IZ) =  XC(NTH,IZ)
  END DO

  ! Assume parity at center
  ! Assume that the shift is smooth at the outer boundary.
  DO IX=1,NTH
   !  XR(IX,0) = -XR(IX,1)
   !  XT(IX,0) =  XT(IX,1)
   !  XP(IX,0) =  XP(IX,1)
   !  XC(IX,0) =  XC(IX,1)
    XR(IX,0) = -XR(NTH-IX+1,1)
    XT(IX,0) =  XT(NTH-IX+1,1)
    XP(IX,0) =  XP(NTH-IX+1,1)
    XC(IX,0) =  XC(NTH-IX+1,1)
    XR(IX,NR+1) =  XR(IX,NR) + DR(NR)/DR(NR-1)*(XR(IX,NR)-XR(IX,NR-1))
    XT(IX,NR+1) =  XT(IX,NR) + DR(NR)/DR(NR-1)*(XT(IX,NR)-XT(IX,NR-1))
    XP(IX,NR+1) =  XP(IX,NR) + DR(NR)/DR(NR-1)*(XP(IX,NR)-XP(IX,NR-1))
    XC(IX,NR+1) =  XC(IX,NR) + DR(NR)/DR(NR-1)*(XC(IX,NR)-XC(IX,NR-1))
  END DO

 ! Compute the various derivative terms
  DO IX=1,NTH
    DO IZ=1,NR
      A1=-DR(IZ+1)/DR(IZ)/(DR(IZ)+DR(IZ+1))
      A2= (DR(IZ+1)-DR(IZ))/DR(IZ)/DR(IZ+1)
      A3= DR(IZ)/DR(IZ+1)/(DR(IZ)+DR(IZ+1))
      B1=-DTH(IX+1)/DTH(IX)/(DTH(IX)+DTH(IX+1))
      B2= (DTH(IX+1)-DTH(IX))/DTH(IX)/DTH(IX+1)
      B3= DTH(IX)/DTH(IX+1)/(DTH(IX)+DTH(IX+1))

      DXRDR = A1*XR(IX,IZ-1)+A3*XR(IX,IZ+1)+A2*XR(IX,IZ)
      DXTDR = A1*XT(IX,IZ-1)+A3*XT(IX,IZ+1)+A2*XT(IX,IZ)
      DXPDR = A1*XP(IX,IZ-1)+A3*XP(IX,IZ+1)+A2*XP(IX,IZ)
      DXCDR = A1*XC(IX,IZ-1)+A3*XC(IX,IZ+1)+A2*XC(IX,IZ)
      DXRDT = B1*XR(IX-1,IZ)+B3*XR(IX+1,IZ)+B2*XR(IX,IZ)
      DXTDT = B1*XT(IX-1,IZ)+B3*XT(IX+1,IZ)+B2*XT(IX,IZ)
      DXPDT = B1*XP(IX-1,IZ)+B3*XP(IX+1,IZ)+B2*XP(IX,IZ)
      DXCDT = B1*XC(IX-1,IZ)+B3*XC(IX+1,IZ)+B2*XC(IX,IZ)

      CTG=COS(TH(IX))/SIN(TH(IX))
      SIN2=SIN(TH(IX))**2

      ARR = 2.*DXRDR-2./3.*(DXRDR+DXTDT+2.*XR(IX,IZ)/R(IZ)+CTG*XT(IX,IZ))
      ATT = 2.*DXTDT/R(IZ)**2+2.*XR(IX,IZ)/R(IZ)**3-   &
            2./3./R(IZ)**2*(DXRDR+DXTDT+2.*XR(IX,IZ)/R(IZ)+CTG*XT(IX,IZ))
      APP = 2.*XR(IX,IZ)/R(IZ)**3/SIN2+2*CTG*XT(IX,IZ)/R(IZ)**2/SIN2+ &
            2./3./R(IZ)**2/SIN2*(DXRDR+DXTDT+2.*XR(IX,IZ)/R(IZ)+CTG*XT(IX,IZ))
      ART = DXTDR+DXRDT/R(IZ)**2
      ARP = DXPDR
      ATP = DXPDT/R(IZ)**2

      CURVR(ix,iz) = XC(IX,IZ)*8.*Pi*S1SRC(ix,iz)/(R(iz)**2*SIN(TH(IX)))
      CURVT(ix,iz) = XC(IX,IZ)*8.*Pi*S2SRC(ix,iz)/(R(iz)**2*SIN(TH(IX)))/(R(iz)**2)
      CURVP(ix,iz) = XC(IX,IZ)*8.*Pi*S3SRC(ix,iz)/(R(iz)**2*SIN(TH(IX)))/(R(iz)**2*SIN(TH(IX))**2)

      CURVR(IX,IZ) = CURVR(IX,IZ) + ART*DXCDT + ARR*DXCDR
      CURVT(IX,IZ) = CURVT(IX,IZ) + ATT*DXCDT + ART*DXCDR
      CURVP(IX,IZ) = CURVP(IX,IZ) + ATP*DXCDT + ARP*DXCDR
    END DO
  END DO

END SUBROUTINE CURV2

! ***************************************************************************
! ***************************************************************************

SUBROUTINE SOURCECHI(DRMETTERM1,DTMETTERM1) !!! QUI
!---------------------------------------------------------
! Compute the metric terms and metric derivatives used in
! the source terms of the scalar field equation
! METERM1(:,:)  = Log(Alpha*Psi^2)
! DRMETTERM1, DTMETTERM1 -> derivatives of METERM1
!--------------------------------------------------------

  USE SYSTEMXNS
  IMPLICIT NONE

  REAL :: A1,A2,A3,B1,B2,B3
  REAL,DIMENSION(0:NTH+1,0:NR+1) :: METERM1,DRMETTERM1,DTMETTERM1

  DO IX=1,NTH
    DO IZ=1,NR
      METERM1(IX,IZ) = LOG(PSL(IX,IZ)*PSI(IX,IZ))
    END DO
  END DO

  !..................
  ! Compute the metric source terms for the vector potential
  !..................
  ! Assume parity on axis
  DO IZ=1,NR
    METERM1(0,IZ)     =  METERM1(1,IZ)
    METERM1(NTH+1,IZ) =  METERM1(NTH,IZ)
  END DO

  ! Assume parity at center
  ! Assume that the metric source term is smooth at the outer boundary.
  DO IX=1,NTH
    METERM1(IX,0)    = METERM1(NTH-IX+1,1)
    METERM1(IX,NR+1) = METERM1(IX,NR) + DR(NR)/DR(NR-1)*(METERM1(IX,NR)-METERM1(IX,NR-1)) !
  END DO

  DO IX=1,NTH
    DO IZ=1,NR
      A1=-DR(IZ+1)/DR(IZ)/(DR(IZ)+DR(IZ+1))
      A2= (DR(IZ+1)-DR(IZ))/DR(IZ)/DR(IZ+1)
      A3= DR(IZ)/DR(IZ+1)/(DR(IZ)+DR(IZ+1))
      B1=-DTH(IX+1)/DTH(IX)/(DTH(IX)+DTH(IX+1))
      B2= (DTH(IX+1)-DTH(IX))/DTH(IX)/DTH(IX+1)
      B3= DTH(IX)/DTH(IX+1)/(DTH(IX)+DTH(IX+1))

      DRMETTERM1(IX,IZ) = A1*METERM1(IX,IZ-1)+A3*METERM1(IX,IZ+1)+A2*METERM1(IX,IZ)
      DTMETTERM1(IX,IZ) = B1*METERM1(IX-1,IZ)+B3*METERM1(IX+1,IZ)+B2*METERM1(IX,IZ)

    END DO
  END DO
! write(6,*)'metterm ',drmetterm1(2,50),drmetterm1(199,50)

END SUBROUTINE SOURCECHI

! ***************************************************************************
! ***************************************************************************

SUBROUTINE SOLVECHI(RHO) !!! QUI
  !  ============================================================
  !  Purpose : this subroutine solves for the scalar field
  !     equation. The procedure and notation are similar to the
  !     ones in SOLVEATIME subroutine, used to solve for the
  !     electric potential ATIME in HYDROEQ.f90
  !
  !  RHO = Source term
  !  ============================================================

  USE SYSTEMXNS
  USE FUNCTIONS
  IMPLICIT NONE

  REAL, DIMENSION(-1:NTH+2,0:NR+1) :: RHO
!   REAL, DIMENSION(0:NTH,0:NR) :: METERM1,DRMETTERM1,DTMETTERM1

  INTEGER :: I, INFO
  ! Array for the Legendre expansion
  REAL,DIMENSION(NGQ) :: XGQ,WGQ,Y
  REAL,DIMENSION(0:MLS,NGQ) :: PN,PD
  REAL,DIMENSION(0:MLS) :: COEF,PDX,PNX
  REAL,DIMENSION(0:NR+1,0:MLS) :: TAB
  ! Array for the radial solution
  INTEGER :: IL
  REAL :: A1,A2,A3,A4,A5,A6,DCI,DCF,TMP
  REAL,DIMENSION(NR) :: DC,DC1,DCP,SS1
  REAL,DIMENSION(NR-1) :: DL1,DU1,DLP,DUP
  REAL,DIMENSION(1:NR,0:MLS) :: PHI

  ! Compute the Gauss-quadrature points XGQ and weights WGQ
  CALL LEGZO(NGQ,XGQ,WGQ)

  ! Compute the matrix of Weighted Legendre Polyn in the quadrature points
  DO I=1,NGQ
    CALL LPN(MLS,XGQ(I),PN(0:MLS,I),PD(0:MLS,I))
    PN(0:MLS,I)=PN(0:MLS,I)*WGQ(I)
  END DO

  DO IZ=1,NR
  ! Fill boundary values for rho in theta at fixed r
    DO I=1,2
      RHO(1-I,IZ)=RHO(I,IZ)
      RHO(NTH+I,IZ)=RHO(NTH+1-I,IZ)
    END DO

    ! Interpolate the disctretized function on the quadrature points
    CALL POLINT(XX(-1:NTH+2),RHO(-1:NTH+2,IZ),NTH+4,XGQ,Y,NGQ)

    ! Evaluate the coefficient of the Legendre poly. expansion
    DO I=0,MLS
      COEF(I)=(2.*I+1.)/2.*DOT_PRODUCT(PN(I,1:NGQ),Y(1:NGQ))
    END DO
    ! Compute the matrix of legendre coefficients
    TAB(IZ,0:MLS)=COEF(0:MLS)
  END DO

  ! ...................
  ! Solve the problem in r for each Legendre polynomial
  ! This must be loopen on all the Legendre coeff rows
  ! Use direct inversion of the tridiag matrix
  ! ...................

  ! Set the matrix elements (case L=0)
  ! Other cases need to add a term to the main diagonal and modify inner and outer
  ! BC for the first and last diag element

  I=1
  A1= 2./DR(I)/(DR(I)+DR(I+1))
  A2= 2./DR(I+1)/(DR(I)+DR(I+1))
  A3=-2./DR(I)/DR(I+1)
  A4=-DR(I+1)/DR(I)/(DR(I)+DR(I+1))
  A5= DR(I)/DR(I+1)/(DR(I)+DR(I+1))
  A6= (DR(I+1)-DR(I))/DR(I)/DR(I+1)

  DCI= A3+2.*A6/R(I)
  DUP(I)= A2+2.*A5/R(I)
  DO I=2,NR-1
    A1= 2./DR(I)/(DR(I)+DR(I+1))
    A2= 2./DR(I+1)/(DR(I)+DR(I+1))
    A3=-2./DR(I)/DR(I+1)
    A4=-DR(I+1)/DR(I)/(DR(I)+DR(I+1))
    A5= DR(I)/DR(I+1)/(DR(I)+DR(I+1))
    A6= (DR(I+1)-DR(I))/DR(I)/DR(I+1)
    DCP(I)  = A3+2.*A6/R(I)
    DLP(I-1)  = A1+2.*A4/R(I)
    DUP(I)    = A2+2.*A5/R(I)
  END DO
  I=NR
  A1= 2./DR(I)/(DR(I)+DR(I+1))
  A2= 2./DR(I+1)/(DR(I)+DR(I+1))
  A3=-2./DR(I)/DR(I+1)
  A4=-DR(I+1)/DR(I)/(DR(I)+DR(I+1))
  A5= DR(I)/DR(I+1)/(DR(I)+DR(I+1))
  A6= (DR(I+1)-DR(I))/DR(I)/DR(I+1)
  DLP(I-1)= A1+2.*A4/R(I)
  DCF  = A3+2.*A6/R(I)

  DO IL=0,MLS
    ! BC at inner radius ( Origin )
    ! For L .ne. 0 the function goes to 0, parity depends on L
    ! We assume odd parity because for NS  L=odd term are 0.
     A1= 2./DR(1)/(DR(1)+DR(1+1))
     A4=-DR(1+1)/DR(1)/(DR(1)+DR(1+1))
     IF(IL .EQ.0)THEN
       DCP(1)=DCI+(A1+2.*A4/R(1))
     ELSE
       DCP(1)=DCI-(A1+2.*A4/R(1))
     END IF

     ! BC at outer radius
     ! Various multipole must dacay as R^(L+1)
     ! The monopole must decay to 1
     A2= 2./DR(NR+1)/(DR(NR)+DR(NR+1))
     A5= DR(NR)/DR(NR+1)/(DR(NR)+DR(NR+1))

     DCP(NR)=DCF+(A2+2.*A5/R(NR))*(R(NR)/(R(NR)+DR(NR+1)))**(IL+1.)

     ! Set main diagonal
     DO I=1,NR
       DC(I)  = DCP(I) - IL*(IL+1.)/R(I)**2.
     END DO

     SS1(1:NR)=TAB(1:NR,IL)

     DL1=DLP
     DC1=DC
     DU1=DUP

     CALL DGTSV(NR,DL1,DC1,DU1,SS1,NR,INFO)
     PHI(1:NR,IL)=SS1(1:NR)		! A_l(r) coefficients in eq. 72, BD 2011
  END DO

  ! Compute the matrix of Legendre Polyn in the grid points
  DO IX=1,NTH

    CALL LPN(MLS,XX(IX),PNX(0:MLS),PDX(0:MLS))		! P_l(cos theta) functions in eq. 73 BD 2011

    ! Compute the scalar field on the grid
    DO IZ=1,NR
      CHI(IX,IZ)=DOT_PRODUCT(PNX(0:MLS),PHI(IZ,0:MLS))
    END DO
  END DO


END SUBROUTINE SOLVECHI

! ********************************************************
! ********************************************************

SUBROUTINE CHISOL(TRACEMM) !!! QUI

!  =============================================================
!  Purpose: this subroutine solves the scalar field equation.
!
!  CHI     = Einstein frame scalar field
!  CHIIIN   = CHI inside the star
!  CHIIOUT  = CHI outside the star
!  RHO     = source term of the scalar field equation
!  CHIIM    = maximum value of CHI in the domain
!  CC1     = monopolar content of CHI
!  =============================================================

  USE SYSTEMXNS
  USE FUNCTIONS
  IMPLICIT NONE

  INTEGER, PARAMETER :: MAXIT = 100 !1500
  INTEGER :: I, INFO,ILOOP
	REAL, DIMENSION(-1:NTH+2,0:NR+1) :: DCHIDR,DCHIDT!,DRMETTERM2,DTMETTERM2
  REAL, DIMENSION(0:NTH+1,0:NR+1) :: DRMETTERM2,DTMETTERM2
  REAL, DIMENSION(-1:NTH+2,0:NR+1) :: RHO,CHIOLD,RHOSCAL,CHIOLD2
  REAL, DIMENSION(1:NTH,1:NR) :: TRACEMM

  ! Array for the Legendre expansion
  REAL,DIMENSION(NGQ) :: XGQ,WGQ,Y
  REAL,DIMENSION(0:MLS,NGQ) :: PN,PD
  REAL,DIMENSION(0:MLS) :: COEF,PDX,PNX
  REAL,DIMENSION(0:NR+1,0:MLS) :: TAB

  ! Array for the radial solution
  REAL,DIMENSION(NR) :: DC,DC1,DCP,SS1
  REAL,DIMENSION(NR-1) :: DL1,DU1,DLP,DUP
  REAL :: A1,A2,A3,A4,A5,A6,B1,B2,B3
  REAL :: DCI,DCF,ERROR,TMP
  REAL :: CHIM,CC1,CONSTANT,CONSTANTTMP
  REAL :: TERM1,TERMa,TERMb,TERM5
  REAL :: TERM2,TERM2a,TERM2b
  REAL :: TERM3,TERM3a,TERM3b
  REAL :: TERM4,TERM4a,TERM4b
  INTEGER :: IL,ITER,ITP


  CHIOLD2=CHI

  DO ITER=1, MAXIT
     !Impose boundary conditions on the scalar field
     !Assume parity on axis
     DO IZ=1,NR
        CHI(0,IZ) =  CHI(1,IZ)
        CHI(NTH+1,IZ) =  CHI(NTH,IZ)
     END DO

     !Assume parity at center and smoothness at the outer boundary
     DO IX=1,NTH
        CHI(IX,0) =  CHI(IX,1)
        CHI(IX,NR+1) =  CHI(IX,NR)*R(NR)/R(NR+1)
     END DO

     !Evaluate derivatives to be used in the source terms
     DO IX=1,NTH
        DO IZ=1,NR
           A1=-DR(IZ+1)/DR(IZ)/(DR(IZ)+DR(IZ+1))
           A2= (DR(IZ+1)-DR(IZ))/DR(IZ)/DR(IZ+1)
           A3= DR(IZ)/DR(IZ+1)/(DR(IZ)+DR(IZ+1))
           B1=-DTH(IX+1)/DTH(IX)/(DTH(IX)+DTH(IX+1))
           B2= (DTH(IX+1)-DTH(IX))/DTH(IX)/DTH(IX+1)
           B3= DTH(IX)/DTH(IX+1)/(DTH(IX)+DTH(IX+1))

           DCHIDR(IX,IZ) = A1*CHI(IX,IZ-1)+A3*CHI(IX,IZ+1)+A2*CHI(IX,IZ)
           DCHIDT(IX,IZ) = B1*CHI(IX-1,IZ)+B3*CHI(IX+1,IZ)+B2*CHI(IX,IZ)

        END DO
     END DO


     !------------------------------------------------------------------------------
     !Solves the scalar field equation
     !------------------------------------------------------------------------------

     !Initialize the charge density
     DO IX=1,NTH
        DO IZ=1,NR

           ASCAL(IX,IZ) = EXP(ALPHA0*(CHI(IX,IZ)-CHIINF)+0.5*BETA0*(CHI(IX,IZ)-CHIINF)**2)
           TERM1=-4.*PI*ASCAL(IX,IZ)**4*TRACEMM(IX,IZ)
           TERM2=ALPHA0+BETA0*(CHI(IX,IZ)-CHIINF)																! GENERALIZZARE
           RHOSCAL(IX,IZ)=TERM1*TERM2

           ! Impose vacuum outside the star
           ! IF(IZ .GT. WSURF(IX)) RHOSCAL(IX,IZ)=0.

        END DO
     END DO

     ! Computes the derivatives of the metric terms in the source term
     CALL SOURCECHI(DRMETTERM2,DTMETTERM2)

     !Initialize the source of Maxwell-Gauss equation
     DO IX=1,NTH
        DO IZ=1,NR
           TERM1= PSI(IX,IZ)**4*RHOSCAL(IX,IZ)

           TERM2a=DRMETTERM2(IX,IZ)*DCHIDR(IX,IZ)
           TERM2b=1./R(IZ)**2*DTMETTERM2(IX,IZ)*DCHIDT(IX,IZ)
           TERM2 = -(TERM2a + TERM2b)

           ! IF(IZ .GT. WSURF(IX)) TERM1=0. !Redundant
           RHO(IX,IZ)= TERM1 + TERM2

        END DO
     END DO

     !Solve the equation by using a semi-spectral method
     CHIOLD=CHI

     CALL SOLVECHI(RHO)

     CHIM=0
     DO IX=1,NTH
        DO IZ=1,NR
           CHIM=MAX(CHIM,ABS(CHI(IX,IZ)))
        END DO
     END DO

     !Check for convergence
     ERROR=MAXVAL(ABS(CHIOLD(1:NTH,1:NR)-CHI(1:NTH,1:NR)))
     IF(VERBOSE .AND. (.NOT. MPICODE) .AND. (.NOT. GR)) WRITE(6,*) ITER,'MaxChi & MaxErrChi =', CHIM,ERROR
     IF(ERROR .LT. TOLCHI) EXIT

  END DO

  CHI(1:NTH,1:NR)=(1.-QFACTORCHI)*CHIOLD2(1:NTH,1:NR)+QFACTORCHI*CHI(1:NTH,1:NR)

  CALL CHIDERIVS

  IF(WCONVC) THEN
     OPEN(12,FILE=TRIM(ADJUSTL(subdirpath))//'Chiconv.dat',access='append')
     WRITE(12,*) ILOOP,CHIM,ITER,ERROR
     CLOSE(12)
  END IF

END SUBROUTINE CHISOL

! ********************************************************
! ********************************************************

SUBROUTINE LEGZO(N,X,W)

! =========================================================
! Purpose : Compute the zeros of Legendre polynomial Pn(x)
!           in the interval [1,-1], and the corresponding
!           weighting coefficients for Gauss-Legendre
!           integration
! Input :   n    --- Order of the Legendre polynomial
! Output:   X(n) --- Zeros of the Legendre polynomial
!           W(n) --- Corresponding weighting coefficients
! =========================================================

  INTEGER :: N,NR,N0,I,K,J
  REAL :: PI,Z,Z0,P,F0,F1,PF,PD,FD,Q,WP,GD
  REAL,DIMENSION(N):: X,W


  X(1)=0.
  W(1)=2.
  IF (N .GT. 1)THEN
  X=0.
  W=0.
  N0=(N+1)/2
  DO  NR=1,N0
    PI=ACOS(-1.D0)
    Z=COS(PI*(NR-0.25D0)/N)

    DO WHILE(ABS(Z-Z0) .GT. ABS(Z)*1.D-14)
      Z0=Z
      P=1.0D0
      DO  I=1,NR-1
        P=P*(Z-X(I))
      END DO
      F0=1.0D0
      IF ((NR.EQ.N0).AND.(N.NE.2*INT(N/2))) Z=0.0D0
      F1=Z
      DO  K=2,N
        PF=(2.0D0-1.0D0/K)*Z*F1-(1.0D0-1.0D0/K)*F0
        PD=K*(F1-Z*PF)/(1.0D0-Z*Z)
        F0=F1
        F1=PF
      END DO

      IF (Z.NE.0.0) THEN
        FD=PF/P
        Q=0.0D0
        DO  I=1,NR
          WP=1.0D0
          DO  J=1,NR
            IF (J.NE.I) WP=WP*(Z-X(J))
          END DO
          Q=Q+WP
        END DO
        GD=(PD-Q*FD)/P
        Z=Z-FD/GD
      END IF

      IF (Z.EQ.0.0)EXIT
    END DO

    X(NR)=Z
    X(N+1-NR)=-Z
    W(NR)=2.0D0/((1.0D0-Z*Z)*PD*PD)
    W(N+1-NR)=W(NR)

  END DO
  END IF

  RETURN

END SUBROUTINE LEGZO

! ***************************************************************************
! ***************************************************************************

SUBROUTINE LPN(N,X,PN,PD)

!  ===============================================
!  Purpose: Compute Legendre polynomials Pn(x)
!           and their derivatives Pn'(x)
!  Input :  x --- Argument of Pn(x)
!           n --- Degree of Pn(x) ( n = 0,1,...)
!  Output:  PN(n) --- Pn(x)
!           PD(n) --- Pn'(x)
!  ===============================================

  INTEGER :: K,N
  REAL :: P0,P1,PF,X
  REAL,DIMENSION(0:N):: PN,PD

  PN(0)=1.0D0
  PN(1)=X
  PD(0)=0.0D0
  PD(1)=1.0D0
  P0=1.0D0
  P1=X
  DO  K=2,N
   PF=(2.0D0*K-1.0D0)/K*X*P1-(K-1.0D0)/K*P0					! Bonnet's recursion formula
   PN(K)=PF
   IF (ABS(X).EQ.1.0D0) THEN
     PD(K)=0.5D0*X**(K+1)*K*(K+1.0D0)
   ELSE
     PD(K)=K*(P1-X*PF)/(1.0D0-X*X)
   ENDIF
   P0=P1
   P1=PF
  END DO

  RETURN
END SUBROUTINE LPN

! ***************************************************************************
! ***************************************************************************

! SUBROUTINE DGTSV(N,DL,D,DU,B,LDB,INFO)
!
! !  ============================================================
! !  -- LAPACK routine (version 3.2) --
! !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
! !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
! !     November 2006
! !  ============================================================
! !  Purpose : solves the equation  A*X = B,
! !            where A is an n by n tridiagonal matrix, by
! !            Gaussian elimination with partial pivoting.
! !
! !            Note that the equation  A'*X = B  may be solved by
! !            interchanging the order of the arguments DU and DL.
! !
! !  Arguments
! !  =========
! !
! !  N       (input) INTEGER
! !          The order of the matrix A.  N >= 0.
! !
! !
! !  DL      (input/output) DOUBLE PRECISION array, dimension (N-1)
! !          On entry, DL must contain the (n-1) sub-diagonal elements of
! !          A.
! !
! !          On exit, DL is overwritten by the (n-2) elements of the
! !          second super-diagonal of the upper triangular matrix U from
! !          the LU factorization of A, in DL(1), ..., DL(n-2).
! !
! !  D       (input/output) DOUBLE PRECISION array, dimension (N)
! !          On entry, D must contain the diagonal elements of A.
! !
! !          On exit, D is overwritten by the n diagonal elements of U.
! !
! !  DU      (input/output) DOUBLE PRECISION array, dimension (N-1)
! !          On entry, DU must contain the (n-1) super-diagonal elements
! !          of A.
! !
! !          On exit, DU is overwritten by the (n-1) elements of the first
! !          super-diagonal of U.
! !
! !  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
! !          On entry, the N  matrix of right hand side matrix B.
! !          On exit, if INFO = 0, the N  solution matrix X.
! !
! !  LDB     (input) INTEGER
! !          The leading dimension of the array B.  LDB >= max(1,N).
! !
! !  INFO    (output) INTEGER
! !          = 0: successful exit
! !          < 0: if INFO = -i, the i-th argument had an illegal value
! !          > 0: if INFO = i, U(i,i) is exactly zero, and the solution
! !               has not been computed.  The factorization has not been
! !               completed unless i = N.
! !
! !  =====================================================================
!
!   INTEGER :: INFO,LDB,N,I
!   REAL ::   B(LDB), D(N), DL(N-1), DU(N-1)
!   REAL,PARAMETER :: ZERO = 0.0D0
!   REAL :: FACT,TEMP
!
!   INFO = 0
!
!   IF(N.LE.0)THEN
!     INFO = -1
!     WRITE(6,*)' N <= 0 IN DGTSV'
!     STOP
!   ELSE IF(NRHS.LT.0)THEN
!     INFO = -2
!     WRITE(6,*)' RHSN < 0 IN DGTSV'
!     STOP
!   ELSE IF(LDB .LT. MAX(1,N))THEN
!     INFO = -7
!     WRITE(6,*)' LBD < MAX(1,N) IN DGTSV'
!     STOP
!   END IF
!
!   DO I = 1, N - 2
!     IF( ABS( D( I ) ).GE.ABS( DL( I ) ) ) THEN
!       ! No row interchange required
!       IF( D( I ).NE.ZERO ) THEN
!         FACT = DL( I ) / D( I )
!         D( I+1 ) = D( I+1 ) - FACT*DU( I )
!         B( I+1 ) = B( I+1 ) - FACT*B( I )
!       ELSE
!         INFO = I
!         RETURN
!       END IF
!       DL( I ) = ZERO
!     ELSE  ! Interchange rows I and I+1
!       FACT = D( I ) / DL( I )
!       D( I ) = DL( I )
!       TEMP = D( I+1 )
!       D( I+1 ) = DU( I ) - FACT*TEMP
!       DL( I ) = DU( I+1 )
!       DU( I+1 ) = -FACT*DL( I )
!       DU( I ) = TEMP
!       TEMP = B( I )
!       B( I ) = B( I+1 )
!       B( I+1 ) = TEMP - FACT*B( I+1 )
!     END IF
!   END DO
!
!   IF( N.GT.1 ) THEN
!     I = N - 1
!     IF( ABS( D( I ) ).GE.ABS( DL( I ) ) ) THEN
!       IF( D( I ).NE.ZERO ) THEN
!         FACT = DL( I ) / D( I )
!         D( I+1 ) = D( I+1 ) - FACT*DU( I )
!         B( I+1 ) = B( I+1 ) - FACT*B( I )
!       ELSE
!         INFO = I
!         RETURN
!       END IF
!     ELSE
!       FACT = D( I ) / DL( I )
!       D( I ) = DL( I )
!       TEMP = D( I+1 )
!       D( I+1 ) = DU( I ) - FACT*TEMP
!       DU( I ) = TEMP
!       TEMP = B( I )
!       B( I ) = B( I+1 )
!       B( I+1 ) = TEMP - FACT*B( I+1 )
!     END IF
!   END IF
!   IF( D( N ).EQ.ZERO ) THEN
!     INFO = N
!     RETURN
!   END IF
!
!   ! Back solve with the matrix U from the factorization.
!   B( N ) = B( N ) / D( N )
!   IF( N.GT.1 )THEN
!     B( N-1 ) = ( B( N-1 )-DU( N-1 )*B( N ) ) / D( N-1 )
!   END IF
!   DO I = N - 2, 1, -1
!     B( I ) = ( B( I )-DU( I )*B( I+1 )-DL( I )* &
!              B( I+2 ) ) / D( I )
!   END DO
!
!   RETURN
!
! END SUBROUTINE DGTSV

! ***************************************************************************
! ***************************************************************************

SUBROUTINE POLINT(XA,YA,N,X,Y,M)

!  ===============================================
!  Purpose:    Polynomial 2nd order Interpolation
!  Input  :    XA --- Table of abcissas  (N)
!              YA --- Table of ordinates (N)
!              N  --- Number of points
!              X  --- Interpolating abscissa value    (M)
!  Output :    Y  --- Returned estimation of function for X  (M)
!
!  ===============================================

  INTEGER,PARAMETER :: NMAX=10
  REAL,DIMENSION(N) :: XA,YA
  REAL,DIMENSION(M) :: X,Y
  INTEGER,DIMENSION(M) :: JJ
  REAL :: D1C,D1R,D1L,D2R,D2L
  INTEGER :: N,I,M,JIN,II ,J

  JIN=1
  DO II=1,M
    I=II
    J=JIN
    DO WHILE((X(I)-XA(J))*(X(I)-XA(J+1)) .GT. 0)
      J=J+1
    END DO
    JJ(I)=J
    JIN=J
  END DO

  DO II=1,M
    I=II
    D1C=(YA(JJ(I)+1)-YA(JJ(I)))/(XA(JJ(I)+1)-XA(JJ(I)))
    D1R=(YA(JJ(I)+2)-YA(JJ(I)))/(XA(JJ(I)+2)-XA(JJ(I)))
    D1L=(YA(JJ(I)+1)-YA(JJ(I)-1))/(XA(JJ(I)+1)-XA(JJ(I)-1))

    D2L=4*(YA(JJ(I)+1)+YA(JJ(I)-1)-2.*YA(JJ(I)))/(XA(JJ(I)+1)-XA(JJ(I)-1))**2.
    D2R=4*(YA(JJ(I)+2)+YA(JJ(I))-2.*YA(JJ(I)+1))/(XA(JJ(I)+2)-XA(JJ(I)))**2.

!PRINT*, 'POLINT --3',II
!PRINT*,D2L,D2R

    IF(D2L*D2R .GT. 0)THEN
!            print*, 'POLINT -- 3A',II
      Y(I)=YA(JJ(I))+D1L*(X(I)-XA(JJ(I)))+0.5*D2L*(X(I)-XA(JJ(I)))**2.
      !print*,'POLINT -- 4',II
      Y(I)=Y(I)+YA(JJ(I)+1)+D1R*(X(I)-XA(JJ(I)+1))+0.5*D2R*(X(I)-XA(JJ(I)+1))**2.
      !PRINT*, 'POLINT -- 5',II
      Y(I)=0.5*Y(I)
    ELSE
 !           PRINT*, 'POLINT -- 3B',II
      Y(I)=YA(JJ(I))+D1C*(X(I)-XA(JJ(I)))
      !PRINT*,'POLINT -- 6',II
    END IF
  END DO

  RETURN

END SUBROUTINE POLINT
