MODULE SYSTEMXNS


!   ============================================================
!   Purpose : This module contains common block array definition
!             and the physisc of the problem  - EoS
!
!   All parameters to be changed by user should be set here.
!
!   ============================================================
!   ============================================================

  REAL :: PI = 3.1415926535897932D+0 ! Definition of pi

! ====================================================================================
! Parameters for convergence
! ====================================================================================

  INTEGER,PARAMETER :: NVALUE  = 100   ! Iteration for the Newton scheme in XNS
  INTEGER :: MAXLOOP = 100      ! Max Number of loops for convergence in XNSMAIN(1 = TOV solution)
  REAL,PARAMETER :: OSCCONV = 1.0E-9   ! Desired precision in case of oscillatory convergence
  REAL,PARAMETER :: MONCONV = 1.0E-8   ! Desired precision in case of monotonous convergence
! ====================================================================================
! Grid parameters
! ====================================================================================

  INTEGER,PARAMETER :: NR    = 250 ! Radial mesh points
  INTEGER,PARAMETER :: NTH   = 100 ! Angular mesh points (1 = 1D TOV solution)
  INTEGER :: IX,IZ
  REAL,DIMENSION(0:NR+1)   :: R,DR,DRM,DRP
  REAL :: DOM
  REAL,DIMENSION(-1:NTH+2) :: XX,TH,DTH
  REAL,PARAMETER :: RMIN = 0.         ! Center
  REAL,PARAMETER :: RMAX = 20.        ! Outer radius for the regular grid if stretch=.FALSE.
  REAL :: REQMAX = 12.!15.!12.!6.9000000000000000     ! Truncation radius for equilibrium solution (to avoid explosion)
  REAL,PARAMETER :: MINRESREG = 0.125 ! Minimum resolution of the grid (if uniform)
  REAL,PARAMETER :: RINI = 1.0000000000000001E-005		  	! Radius used for the expansion of the TOV equations

  INTEGER,PARAMETER :: NRREG = 600       ! Points for the regular grid
  REAL,PARAMETER :: MINRESSTR = 8.D-2	! Minimum resolution of the grid (if stretched)
  REAL,PARAMETER :: RREG = 10.          ! Radius for the regular Grid if Stretch=.TRUE.
  REAL :: STRR                          !Stretching Ratio
  INTEGER :: MIDGRID
  REAL,PARAMETER :: RMAXSTR = 100.00000000000000         ! Radius for the total grid if stretch=.TRUE.

  REAL,DIMENSION(0:NR) :: ACOEFF,BCOEFF,CCOEFF   ! Coefficients of the tridiagonal matrix A used in DGTSV
  INTEGER,PARAMETER :: LDB=NR,NMAT=NR				 ! Dimensions of the elements of the matrix equation
  REAL :: BMAT(NR), DMAT(NR), DLMAT(NR-1), DUMAT(NR-1)		 ! Diagonals of the A matrix and source vector B
  REAL :: DS(NR), DSL(NR-1), DSU(NR-1)						 ! Dummy diagonals
  REAL :: RHOSURF = 1.0E-08    ! Density at which the surface is set
! ====================================================================================
! Scalar field parameters
! ====================================================================================

  REAL :: ALPHA0 = 0.0!-2.0E-004 			! Alpha0 parameter for the scalar field
  REAL :: BETA0 = 0.0!-5.5   		   	! Beta0 parameter for the scalar field
  REAL,PARAMETER :: CHIINF = 0.0000000000000000			! Cosmological value of the scalar field

! ====================================================================================
! EoS/fluid parameters
! ====================================================================================

  REAL,PARAMETER :: RHOINI = 1.28e-3!3.81e-4!2.6158e-3!1.28E-03!1.4090909090909089EE-03! Central density in the Jordan frame (beware scheme converges to QUCONV)
  REAL,PARAMETER :: MBARYONFC = 1.0!0.86 ! Ratio between tabulated reduced baryon mass and true baryon mass
  REAL :: RHOINISEQ   ! Initial central density for sequances
  REAL,PARAMETER :: RHOCENTSTART = 5.0D-4 	! Initial central density used as starting point for the mass-rho diagram
  REAL,PARAMETER :: RHOCENTEND = 8.0D-3  	! Final central density of the mass-rho diagram
  REAL,PARAMETER :: K1 = 110.00000000000000 !100.         ! Politropic coefficient
  REAL,PARAMETER :: GAMMA = 2.0000000000000000        ! Politropic exponent

  LOGICAL :: EOSINT = .FALSE. ! If true use an interpolated/tabulated EoS
  CHARACTER(LEN=30) :: FILEEOS = 'XXX_resampled.dat'!'PL2_resampled.dat'
  INTEGER,PARAMETER :: NPTRHO = 1000 ! Points of the tabulated EoS from resample.py
  REAL :: RHOTABMIN,RHOTABMAX,PRSTABMIN,PRSTABMAX,EINTABMIN,EINTABMAX,ENTTABMIN,ENTTABMAX  ! Tabulated maxima and minima of rho,prs,ein,ent
  REAL :: PRSTABINDL,PRSTABINDU,EINTABINDL,EINTABINDU,ENTTABINDL,ENTTABINDU ! Tabulated maxima and minima power-law index of of prs,ein,ent
  REAL,DIMENSION(NPTRHO) :: RHOVEC1,PRSVEC1,EINVEC1,ENTVEC1,PRSIND1,EININD1,ENTIND1  ! From primary table
  REAL,DIMENSION(NPTRHO) :: RHOVEC2,PRSVEC2,RHOIND2 ! From secondary table
  REAL,DIMENSION(NPTRHO) :: ENTVEC3,RHOVEC3,RHOIND3 ! From thirtiary table

! ====================================================================================
! Logical variables
! ====================================================================================

  LOGICAL :: GR = .TRUE.	! False = scalar-tensor theory; true = general relativity. Must be false here
  LOGICAL,PARAMETER :: STRETCH = .False. ! If false the grid is uniform; if true it is uniform up to RREG and logarithmically stretched beyond
  LOGICAL :: DEBUG = .FALSE. 		! Set to true to print things for debugging purposes
  LOGICAL :: SINGLESURF = .TRUE.	! If true the MPSAMPLING.f90 program will build the sampling only for the RHOINI central density
  LOGICAL :: COUNTDOWN = .TRUE.	        ! Set to true to print the remaining number of configurations during the sequence computation
  LOGICAL :: LOGFILE = .TRUE.		! Set to true to save a log file
  LOGICAL :: ANALYTIC = .FALSE.		! Set to true if you want to solve for phi using the analytic equation D(D(phi))=...
  LOGICAL :: EOSJOR = .FALSE.		! Used to change the frame (Jor. or Ein.) in which the EoS is computed (LEAVE FALSE HERE)
  LOGICAL :: CONVHELP = .FALSE.         ! If true, it activates an option for RHOCENT in HYDROEQ.f90 to help achieve convergence
  					! WARNING: if true, the final central density will be slightly different than the chosen RHOINI

! ====================================================================================
! Flags for writing output files. Put false for speed test.
! ====================================================================================

! Will always wrtite output for the final step of convergence loop
  LOGICAL :: WRT = .TRUE.      ! Write all files
  LOGICAL :: WRTF = .FALSE.     ! Write only final files
  LOGICAL :: VERBOSE = .FALSE.  ! Print all infos
  LOGICAL :: CHUP = .TRUE.     ! Write all files (including the metric test)
  LOGICAL :: WGRID = .TRUE.     ! Write the grid file
  LOGICAL :: WCONVA = .FALSE.    ! Write file to check convergence of Aphi and Atim
  LOGICAL :: WCONVC = .FALSE.    ! Write file to check convergence of chi
  LOGICAL :: IDAT               ! Write infos on the convergence for QUOC

! ====================================================================================
! Parameters for  the overall Newton-Raphson scheme (if precompiled)
! ====================================================================================

#ifdef NWTRPS
  INTEGER,PARAMETER :: QUOC = 0          !Quantity of convergence for sequences
                                       !QUOC=0 => central density
                                       !QuOC=1 => gravitational mass
                                       !QUOC=2 => baryonic mass
  REAL,PARAMETER :: QUCONV = 0.!RHOINI  !Value of quantity to convergence
#endif

! ====================================================================================
! Physics parameters - Rotation
! ====================================================================================

  REAL,PARAMETER :: OMG = 2.22e-2!2.633e-2!1.97549e-2!2.633e-2!2.22E-2!5.0E-3       ! Central Rotation rate
  LOGICAL        :: DIFFERENTIAL = .TRUE.       ! Differential rotation option
  LOGICAL        :: OMGSPACE = .FALSE.           ! If .true. j(Omega), else Omega(j)
  LOGICAL        :: JCONSTLAW = .FALSE.           ! Rotation law: j constant A^2*(omega_c-omega)
  LOGICAL        :: JCMODLAW = .FALSE.           ! Rotation law: modified j constant A^2*omega*[(omega_c/omega)^p-1]
  LOGICAL        :: URYULAW3 = .FALSE.           ! Rotation law: Uryu with 3 Parameters
  LOGICAL        :: URYULAW4 = .TRUE.           ! Rotation law: Uryu with 4 Parameters
  REAL,PARAMETER :: PROTDIFF = 3./2.           ! Rotation index if URYULAW3=.true.
  REAL,PARAMETER :: A2VALUE = 70.0000000000000000           ! Differential rotation coeff
  REAL,PARAMETER :: OMGMAX = 2.*OMG           !Omega_max/Omega_c
  REAL,PARAMETER :: RMVALUE = 2.5!6.54           !R_max (Omega=Omega_max)


! ====================================================================================
! Physics - Magnetic Fields
! ====================================================================================

  LOGICAL :: IMAG = .FALSE.    ! Magnetized cases
  LOGICAL :: ITOR = .FALSE.   ! Purely toroidal B-field
  LOGICAL :: IPOL = .FALSE.   ! Purely poloidal B-field
  LOGICAL :: ITWT = .FALSE.    ! Mixed B-field

! ====================================================================================
! Physics - purely TOROIDAL B-field only!
! ====================================================================================

  REAL :: BCOEF  = 0.0     !Toroidal magnetization constant
  REAL,PARAMETER :: MAGIND = 1.0000000000000000      !Toroidal magnetization index (>=1)

! ====================================================================================
! Physics - purely POLOIDAL B-field only!
! ====================================================================================

  REAL :: KBPOL = 0.0                ! Magnetic coefficient
  REAL,PARAMETER :: NPOL  = 0.0000000000000000          ! Magnetic powerlaw index
  REAL,PARAMETER :: CSI = 0.0000000000000000            ! Current coefficient
  LOGICAL,PARAMETER :: QNULL = .TRUE.    ! Logical parameter for the star charge (IF OMG.NE.0)
  LOGICAL,PARAMETER :: CTP = .FALSE.     ! Use cons to prim routines

! ====================================================================================
! Physics - TWISTED-TORUS magnetic field
! ====================================================================================

  REAL :: KBTT = 0.0       ! Magnetic coefficient
  REAL,PARAMETER :: ATWT = 1.4E-3        ! TT magnetic coefficient
  REAL,PARAMETER :: ZETA = 0.0000000000000000            ! TT magnetic index
  REAL,PARAMETER :: CUT = 4.0000000000000000              ! Maximum distance for twisted magnetosphere

! ====================================================================================
! Parameters to help or stabilize convergence
! ====================================================================================

  REAL,PARAMETER    :: CONVF =  1.d-4     ! Convergence of the newton scheme in XNS
  REAL    :: QFACTOR = 0.7!0.5!0.85       ! Damping factor for press the convergence loop  qnew = qfactro*qnew + (1-qfactor)*qold
  REAL    :: QFACTORMETRIC = 0.5!0.35 ! Damping factor for psi/psl for the internal convergence loop  qnew = qfactro*qnew + (1-qfactor)*qold
  REAL    :: QFACTORCONF = 0.5    ! Damping factor for psi for the external convergence loop  qnew = qfactro*qnew + (1-qfactor)*qold
  REAL    :: QFACTORCHI = 0.45    ! Damping factor for chi for the convergence loop  qnew = qfactro*qnew + (1-qfactor)*qold
  REAL,PARAMETER    :: QAPHI = 1.!0.5         ! Damping factor for the convergence of the vector potential
  REAL,PARAMETER    :: EPS = 9.9999999999999995E-008

! ====================================================================================
! Array for the Legendre expansion
! ====================================================================================

 INTEGER,PARAMETER :: MLS = 20 !20!8       ! Number of Legendre polinomia for expansion in theta (0 = 1D)
  INTEGER,PARAMETER :: NGQ = 50 !50      ! Number of point in the Gauss quadrature (1 = 1D)
  INTEGER,PARAMETER :: MLSL = 20     ! Number of Legendre polinomia used to solve Laplace equation
  REAL,PARAMETER    :: TOLCONV = 1.D-10 ! Convergence for the lapse and conformal factor iterative solvers
  REAL,PARAMETER    :: TOLCHI = 1.D-10 	! Convergence for the scalar field iterative solvers

! ====================================================================================
! Parameters of the fluid field and metric
! ====================================================================================

  REAL,DIMENSION(1:NTH,1:NR) :: PSI     ! Conformal factor
  REAL,DIMENSION(1:NTH,1:NR) :: PSI0     ! Conformal factor (old)
  REAL,DIMENSION(1:NTH,1:NR) :: PSL     ! Lapse * conformal factor
  REAL,DIMENSION(1:NTH,1:NR) :: PSS     ! Phi - shift
  REAL,DIMENSION(1:NTH,1:NR) :: PSSR    ! R - shift
  REAL,DIMENSION(1:NTH,1:NR) :: PSST    ! Theta - shift

  REAL,DIMENSION(-1:NTH+2,0:NR+1):: RHOSRC,ESRC,PSRC,VPHI,VR,VTH,BPHI,EPHI,SSS
  REAL,DIMENSION(-1:NTH+2,0:NR+1):: USRC,DSRC,S3SRC,S1SRC,S2SRC
  REAL,DIMENSION(-1:NTH+2,0:NR+1):: USRCX,S3SRCX,S1SRCX,S2SRCX
  REAL,DIMENSION(-1:NTH+2,0:NR+1):: ECSRC,ELSRC,ES1SRC,ES2SRC,ES3SRC
  REAL,DIMENSION(-1:NTH+2,0:NR+1):: CURVC,CURVR,CURVT,CURVP

! ====================================================================================
! For the initial TOV solution
! ====================================================================================

  REAL,DIMENSION(1:NR) :: RHOTV,PRTV,ETV,RHOTVJOR,ETVJOR,PRTVJOR
  REAL,DIMENSION(0:NR+1) :: CHITV,NU,MU,DCHITV,DDCHITV
  INTEGER,PARAMETER :: RELIT = 100								 ! Number of iterations used in the RELAXLOOP and in the PHILOOP
  REAL,PARAMETER :: CONV2 = 1.0D-6								 ! Convergence parameter for the relaxation loop
  INTEGER,PARAMETER :: MAXITM0 = 30       	! Maximum number of iterations for the GR shooting
  INTEGER,PARAMETER :: MAXSTEP = 90     	! Maximum number of iterations for the STT shooting
  REAL,PARAMETER :: CONV = 1.0D-5			! Convergence parameter for the GR shooting
  REAL,PARAMETER :: DELTMU0 = 0.0001		! Delta mu0 used in the main shooting
  REAL :: MMID = 1.9			! Initial guess on the masses at MIDGRID and NR
  REAL :: M = 1.9
  REAL :: MUIN = 0.802!1.0e-2!0.802!1.0000000000000000E-002	        ! Initial guess on the central mu
  REAL,PARAMETER :: QRELAX = 0.29999999999999999								 ! Used in the filter to force convergence of CHI

! ====================================================================================
! Paramaters of the Equilibrium solution via bernoulli integral
! ====================================================================================

  REAL,DIMENSION(1:NTH,1:NR) :: PNEW,RHONEW,ENEW,V3NEW,B3NEW,E3NEW,TRACEM

  INTEGER,DIMENSION(1:NTH) :: WSURF,CUTSURF    !Integer vectors for surfaces
  REAL,DIMENSION(-1:NTH+2) :: ELSURF           !Fitted surface
  REAL :: NSPE

! ====================================================================================
! Scalar field variables
! ====================================================================================

  REAL,DIMENSION(-1:NTH+2,0:NR+1) :: CHI,PSCAL,QSCALTIM,QSCALR,QSCALT,QSCALP,QSCAL2,ASCAL
  REAL :: CHIS
  REAL :: USR,SYSR,SY2SR,DSR,DENSC,PRESSC,VELOCC,BMAG,B2YSR
  REAL :: USRX,SYSRX,SY2SRX
  REAL,DIMENSION(-1:NTH+2,0:NR+1) :: CHIIN,CHIOUT

! ====================================================================================
! Parameters for the magnetic vector potential
! ====================================================================================

  REAL,DIMENSION(1:NTH,1:NR) :: RHOTERM,OMTERM,VMETERM,VLOC
  REAL,DIMENSION(1:NTH,1:NR) :: JPHIMXL,RHOEMXL
  REAL,DIMENSION(1:NTH,1:NR) :: DRMETTERM,DTMETTERM
  REAL,DIMENSION(1:NTH,1:NR) :: DROMGM,DTOMGM
  REAL,DIMENSION(1:NTH,1:NR) :: DRGAML,DTGAML
  REAL,DIMENSION(1:NTH,1:NR) :: BPOLT,BPOLR,EPOLT,EPOLR,ATHETA
  REAL,DIMENSION(1:NTH,1:NR) :: JPHI, JTH, JRR

  REAL,DIMENSION(-1:NTH+2,0:NR+1) :: METERM,GAMLOC,LOGGAM,OMGMET
  REAL,DIMENSION(-1:NTH+2,0:NR+1) :: APHI,ATIM,ATIMARM,ATIMIN,ATIMOUT

  REAL :: APHIMAX,APM,RFCUT

! ====================================================================================
! Parameters for covterm
! ====================================================================================

  REAL    :: ALPHA,GP,GM,GCOVR,GCOVT,GCOVP
  INTEGER :: ISUR,ELOOP
  REAL    :: REQCHECK

! ====================================================================================
! Stuff for parallel-MPI runs (for parameters study)
! ====================================================================================

  CHARACTER(len=100) :: subdirpath
  INTEGER :: XNSERR
  LOGICAL :: ENDID
#ifdef MPIV
	LOGICAL :: MPICODE = .TRUE.
#else
	LOGICAL :: MPICODE = .FALSE.
#endif

! ====================================================================================
! Indices/stuff
! ====================================================================================

  INTEGER :: ILOC
  INTEGER :: K,L,STEP,RHOSTEP = 0
  REAL :: START,FINISH
  CHARACTER(8)  :: DATE
  CHARACTER(10) :: TIME
  CHARACTER(5)  :: ZONE
  INTEGER,DIMENSION(8) :: VALUES

! ====================================================================================
! Used in the sampling VAR_SAMPL.f90
! ====================================================================================

  INTEGER,PARAMETER :: JRHOMAX = 9     		! Number of configurations of the MPSAMPLING sampling
  INTEGER,PARAMETER :: NUMCYCLES = 100  	! Number of samples of CHI0 and mu0 in MPSAMPLING
  INTEGER,PARAMETER :: ABCYCLES = 30  		! Number of samples of alpha0 and beta0 in ABSAMPLING
  REAL :: CHI0START = 0.2               	! Ranges of the sampling of the CHI0-mu0 space in MPSAMPLING
  REAL :: CHI0END = -0.2
  REAL,PARAMETER :: MU0START = 2
  REAL,PARAMETER :: MU0END = 0.
  REAL,PARAMETER :: BETA0START = -2.			! Ranges of the sampling of the alpha0-beta0 space in ABSAMPLING
  REAL,PARAMETER :: BETA0END = -6.
  REAL,PARAMETER :: ALPHA0START = -0.001
  REAL,PARAMETER :: ALPHA0END = -0.01


CONTAINS

! ********************************************************
! ********************************************************

  ! =============== READ EOS TABLE =============================
  SUBROUTINE EOSTABLEREAD

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
!     write(6,*)IL
    IF(IL .GE. 0 .AND. IL .LE. 999)THEN
       PSRLOG = PRSVEC1(IL+1) +(RHOLOG -RHOVEC1(IL+1))*PRSIND1(IL+1)
       EINLOG = EINVEC1(IL+1) +(RHOLOG -RHOVEC1(IL+1))*EININD1(IL+1)
       ENTLOG = ENTVEC1(IL+1) +(RHOLOG -RHOVEC1(IL+1))*ENTIND1(IL+1)
    END IF
    IF(IL .LT. 0)THEN
       PSRLOG = PRSVEC1(1) +(RHOLOG -RHOVEC1(1))*PRSIND1(1)
       EINLOG = EINVEC1(1) +(RHOLOG -RHOVEC1(1))*EININD1(1)
       ENTLOG = ENTVEC1(1) +(RHOLOG -RHOVEC1(1))*ENTIND1(1)
    END IF
    IF(IL .GE. 1000)THEN
       PSRLOG = PRSVEC1(1000) +(RHOLOG -RHOVEC1(1000))*PRSIND1(1000)
       EINLOG = EINVEC1(1000) +(RHOLOG -RHOVEC1(1000))*EININD1(1000)
       ENTLOG = ENTVEC1(1000) +(RHOLOG -RHOVEC1(1000))*ENTIND1(1000)
    END IF
    PRS=10**PSRLOG
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
    REAL :: RHO,PRS
    REAL :: RHOLOG,PRSLOG

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
       !write(6,*)'eos',p,rho,en

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
!  write(6,*)'eos1',rhovar,y,x
  XN=1.0001*X
  CALL EOS(XN,YN,ZN)
  YN=YN-RHOVAR

  DY=(YN-Y)/(XN-X)
 ! write(6,*)'eos',rhovar

  RETURN

END SUBROUTINE FUNCD_EOS

! ********************************************************

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

! ***********************************************************************
! **************************** NR routines ******************************

REAL FUNCTION RTSAFEG(X0,X1,X2,XACC,RHOVAR,FUNCD)

  INTERFACE
    SUBROUTINE FUNCD(a,b,c,d)
      real :: a,b,c,d
    END SUBROUTINE FUNCD
  END INTERFACE
  real :: X0,X1,X2,XACC,RHOVAR,FL,DF,DX,DXOLD,F,FH,XL,XH,TEMP,SWAP,XN
  integer :: j


  INTEGER,PARAMETER :: MAXIT=2000

 ! write(6,*)'rtsafeg',x1,fl,df,rhovar

  CALL FUNCD(X1,FL,DF,RHOVAR)
!   write(6,*)'rtsafeg'

  IF(FL.EQ.0.) THEN
     RTSAFEG=X1
     RETURN
  ENDIF
  CALL FUNCD(X2,FH,DF,RHOVAR)

  IF(FH.EQ.0.) THEN
     RTSAFEG=X2
     RETURN
  ENDIF
  IF(FL*FH.GT.0.)THEN
     WRITE(6,*)'FL ',FL ,'FH ',FH
     !STOP 'Root must be bracketed'
     XACC=-1.
     RETURN
  ENDIF
  IF(FL.LT.0.)THEN
     XL=X1
     XH=X2
  ELSE
     XH=X1
     XL=X2
     SWAP=FL
     FL=FH
     FH=SWAP
  ENDIF
  RTSAFEG=X0
  DXOLD=ABS(X2-X1)
  DX=DXOLD
  CALL FUNCD(RTSAFEG,F,DF,RHOVAR)
  DO 11 J=1,MAXIT
     IF(((RTSAFEG-XH)*DF-F)*((RTSAFEG-XL)*DF-F).GE.0. &
          .OR. ABS(2.*F).GT.ABS(DXOLD*DF) ) THEN
        DXOLD=DX
        DX=0.5*(XH-XL)
        RTSAFEG=XL+DX
        IF(XL.EQ.RTSAFEG)RETURN
     ELSE
        DXOLD=DX
        DX=F/DF
        TEMP=RTSAFEG
        RTSAFEG=RTSAFEG-DX
        IF(TEMP.EQ.RTSAFEG)RETURN
     ENDIF
     IF(ABS(DX).LT.XACC) RETURN
     CALL FUNCD(RTSAFEG,F,DF,RHOVAR)
     IF(F.LT.0.) THEN
        XL=RTSAFEG
        FL=F
     ELSE
        XH=RTSAFEG
        FH=F
     ENDIF
11    CONTINUE
     WRITE(6,*) 'RTSAFE exceeding maximum iterations'

     XACC=+1
     RETURN
END FUNCTION RTSAFEG

! **************************************************************************

SUBROUTINE XNS2ECHO_OUT()

  open(10,file='XNS2ECHO_init.dat',form='unformatted')
  print*, NTH,NR
  write(10) NTH
  write(10) NR
  write(10) RMIN
  write(10) RMAX
  write(10) RHONEW/MBARYONFC
  write(10) V3NEW
  write(10) PNEW
  write(10) EPOLR
  write(10) EPOLT
  write(10) BPOLR
  write(10) BPOLT
  write(10) B3NEW
  write(10) PSI
  write(10) PSL/PSI
  write(10) PSS
  close(10)

END SUBROUTINE XNS2ECHO_OUT


END MODULE SYSTEMXNS
