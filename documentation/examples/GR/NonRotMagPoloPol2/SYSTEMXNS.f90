MODULE SYSTEMXNS
  !   ============================================================
  !   Purpose : This module contains common block array definition
  !             and parameters for the simulations
  !
  !   All parameters to be changed by user should be set here.
  !
  !   ============================================================
  !   ============================================================
  
  REAL :: PI = 3.1415926535897932D+0 ! Definition of pi

  ! ====================================================================================
  ! Parameters for convergence of the Newton-Raphson scheme of XNS
  ! ====================================================================================
  
  INTEGER,PARAMETER :: NVALUE  = 100   ! Iteration for the Newton scheme in XNS
  REAL,PARAMETER    :: CONVF =  1.d-4  ! Convergence of the newton scheme in XNS

  ! ====================================================================================
  ! Parameters for convergence of XNSMAIN over a model
  ! ====================================================================================
  
  INTEGER :: MAXLOOP = 100                 ! Max Number of loops for convergence in XNSMAIN (1 = TOV solution)
  LOGICAL,PARAMETER :: CONVRHO = .FALSE.   ! Set .TRUE. if convergence should be checked on density - otherwise checked on lapse 
  REAL,PARAMETER :: MONCONV = 1.0E-8       ! Desired precision in case of monotonous convergence (absolute on RHOC, relative on PSI)
  REAL,PARAMETER :: OSCCONV = 1.0E-9       ! Desired precision in case of oscillatory convergence (absolute on RHOC)
  
  ! ====================================================================================
  ! Grid parameters
  ! ====================================================================================

  LOGICAL,PARAMETER :: STRETCH = .TRUE. ! If false the grid is uniform; if true it is uniform up to RREG and logarithmically stretched beyond
  
  INTEGER,PARAMETER :: NR    = 900 ! Radial mesh points
  INTEGER,PARAMETER :: NTH   = 100 ! Angular mesh points (1 = 1D TOV solution)
  INTEGER,PARAMETER :: NRREG = 600       ! Points for the regular grid if stretch=.TRUE.

  REAL,PARAMETER :: RMIN = 0.         ! Center
  REAL,PARAMETER :: RMAX = 100.       ! Outer radius for the regular grid if stretch=.FALSE.
  REAL,PARAMETER :: RMAXSTR = 100.0   ! Outer radius for the total grid if stretch=.TRUE.
  REAL,PARAMETER :: RREG = 13.        ! Radius for the regular Grid if Stretch=.TRUE.
  
  REAL :: REQMAX = 13.50               ! Truncation radius for equilibrium solution (to avoid explosion)
  REAL,PARAMETER :: MINRESREG = 0.125 ! Minimum resolution of the grid (if uniform)
  REAL,PARAMETER :: MINRESSTR = 8.D-2 ! Minimum resolution of the grid (if stretched)
  REAL,PARAMETER :: RINI = 1.E-005    ! Radius used for the expansion of the TOV equations
  
  REAL :: RHOSURF = 1.0E-08             ! Density at which the surface is set

  REAL :: STRR                          !Stretching Ratio
  REAL :: DOM
  INTEGER :: MIDGRID
  INTEGER :: IX,IZ
  REAL,DIMENSION(0:NR+1)   :: R,DR,DRM,DRP  ! Radial grid arrays
  REAL,DIMENSION(-1:NTH+2) :: XX,TH,DTH     ! Angular grid arrays
  
  REAL,DIMENSION(0:NR) :: ACOEFF,BCOEFF,CCOEFF          ! Coefficients of the tridiagonal matrix A used in DGTSV
  INTEGER,PARAMETER :: LDB=NR,NMAT=NR	                ! Dimensions of the elements of the matrix equation
  REAL :: BMAT(NR), DMAT(NR), DLMAT(NR-1), DUMAT(NR-1)	! Diagonals of the A matrix and source vector B
  REAL :: DS(NR), DSL(NR-1), DSU(NR-1)			! Dummy diagonals
  
  ! ====================================================================================
  ! Gravitational Theory  parameters
  ! ====================================================================================

  LOGICAL :: GR = .TRUE.	        ! False = scalar-tensor theory; true = general relativity. Must be false here
  
  REAL :: ALPHA0 = -0.0E-004 	! Alpha0 parameter for the scalar field (def = -2.0E-004)
  REAL :: BETA0  = -0.0   		! Beta0 parameter for the scalar field (def = -6.0)
  REAL,PARAMETER :: CHIINF = 0.0        ! Cosmological value of the scalar field (beware not fully tested)

  ! ====================================================================================
  ! EoS/fluid parameters
  ! ====================================================================================

  REAL,PARAMETER :: RHOINI = 8.3e-4 ! Central density in the Jordan frame (beware scheme converges to QUCONV)
  REAL,PARAMETER :: MBARYONFC = 1.0  ! Ratio between tabulated reduced baryon mass and true baryon mass
  LOGICAL,PARAMETER :: VACUUM = .FALSE. ! Set to zero the physical source terms in the conf-lapse solver outside the star 
  
  REAL,PARAMETER :: K1 = 110.0           ! Politropic coefficient
  REAL,PARAMETER :: GAMMA = 2.0          ! Politropic exponent

  REAL :: RHOINISEQ   ! Initial central density for sequences
  REAL,PARAMETER :: RHOCENTSTART = 5.0D-4 	! Initial central density used as starting point for the mass-rho diagram
  REAL,PARAMETER :: RHOCENTEND = 8.0D-3  	! Final central density of the mass-rho diagram

  LOGICAL,PARAMETER :: CTP = .FALSE.     ! Use cons to prim routines
  LOGICAL :: EOSINT = .FALSE. ! If true use an interpolated/tabulated EoS
  CHARACTER(LEN=30) :: FILEEOS = 'XXX_resampled.dat'!'XXX_resampled.dat'
  INTEGER,PARAMETER :: NPTRHO = 1000 ! Points of the tabulated EoS from resample.py
  LOGICAL :: EOSJOR = .FALSE.		! Used to change the frame (Jor. or Ein.) in which the EoS is computed (LEAVE FALSE HERE)
  
  REAL :: RHOTABMIN,RHOTABMAX,PRSTABMIN,PRSTABMAX,EINTABMIN,EINTABMAX,ENTTABMIN,ENTTABMAX  ! Tabulated maxima and minima of rho,prs,ein,ent
  REAL :: PRSTABINDL,PRSTABINDU,EINTABINDL,EINTABINDU,ENTTABINDL,ENTTABINDU ! Tabulated maxima and minima power-law index of of prs,ein,ent
  REAL,DIMENSION(NPTRHO) :: RHOVEC1,PRSVEC1,EINVEC1,ENTVEC1,PRSIND1,EININD1,ENTIND1  ! From primary table
  REAL,DIMENSION(NPTRHO) :: RHOVEC2,PRSVEC2,RHOIND2 ! From secondary table
  REAL,DIMENSION(NPTRHO) :: ENTVEC3,RHOVEC3,RHOIND3 ! From thirtiary table

  ! ====================================================================================
  ! Physics parameters - Rotation
  ! ====================================================================================

  LOGICAL        :: DIFFERENTIAL = .FALSE.   ! Differential rotation option
  LOGICAL        :: OMGSPACE     = .TRUE.    ! If .true. j(Omega), else Omega(j)
  LOGICAL        :: JCONSTLAW    = .FALSE.   ! Rotation law: j constant A^2*(omega_c-omega)
  LOGICAL        :: JCMODLAW     = .FALSE.   ! Rotation law: modified j constant A^2*omega*[(omega_c/omega)^p-1]
  LOGICAL        :: URYULAW3     = .FALSE.   ! Rotation law: Uryu with 3 Parameters
  LOGICAL        :: URYULAW4     = .FALSE.    ! Rotation law: Uryu with 4 Parameters
  REAL,PARAMETER :: OMG = 0.00               ! Central Rotation rate
  REAL,PARAMETER :: PROTDIFF = 3./2.         ! Rotation index if URYULAW3=.true.
  REAL,PARAMETER :: A2VALUE = 0.0            ! Differential rotation coeff
  REAL,PARAMETER :: OMGMAX = 2.*OMG          ! Omega_max/Omega_c
  REAL,PARAMETER :: RMVALUE = 2.5!6.54       ! R_max (Omega=Omega_max)

  ! ====================================================================================
  ! Physics - Magnetic Fields
  ! ====================================================================================
  
  LOGICAL :: IMAG = .TRUE.    ! Magnetized cases
  LOGICAL :: ITOR = .FALSE.   ! Purely toroidal B-field
  LOGICAL :: IPOL = .TRUE.   ! Purely poloidal B-field
  LOGICAL :: ITWT = .FALSE.    ! Mixed B-field
  
  ! ====================================================================================
  ! Physics - purely TOROIDAL B-field only!
  ! ====================================================================================
  
  REAL :: BCOEF  = 0.0           !Toroidal magnetization constant
  REAL,PARAMETER :: MAGIND = 1.0   !Toroidal magnetization index (>=1)
  
  ! ====================================================================================
  ! Physics - purely POLOIDAL B-field only!
  ! ====================================================================================
  
  REAL :: KBPOL = 0.34135                    ! Magnetic coefficient
  REAL,PARAMETER :: NPOL  = 0.0          ! Magnetic powerlaw index
  REAL,PARAMETER :: CSI = 0.0            ! Current coefficient
  LOGICAL,PARAMETER :: QNULL = .TRUE.    ! Logical parameter for the star charge (IF OMG.NE.0)
  
  ! ====================================================================================
  ! Physics - TWISTED-TORUS magnetic field
  ! ====================================================================================
  
  REAL :: KBTT = 0.0                  ! Magnetic coefficient
  REAL,PARAMETER :: ATWT = 0.0     ! TT magnetic coefficient
  REAL,PARAMETER :: ZETA = 0.0        ! TT magnetic index
  REAL,PARAMETER :: CUT = 0.0         ! Maximum distance for twisted magnetosphere
  
  ! ====================================================================================
  ! For the initial TOV solution
  ! ====================================================================================

  INTEGER,PARAMETER :: RELIT = 100	! Number of iterations used in the RELAXLOOP 
  REAL,PARAMETER :: CONVT = 1.0D-6	! Convergence parameter for the relaxation loop (also on  CHITVLOOP)
  REAL,PARAMETER :: CONV = 1.0D-5	! Convergence parameter for the GR shooting
  INTEGER,PARAMETER :: MAXSTEPTV = 100    ! Maximum number of iterations for the GR shooting at fixed scalalr field in STEPLOOP
  INTEGER,PARAMETER :: MAXSTEPCH = 100    ! Maximum number of iterations for the STT shooting at fixed metric in CHITVLOOP
  REAL,PARAMETER :: DELTMU0 = 0.0001	! Delta mu0 used in the main shooting
  
  REAL :: MMID = 1.3			! Initial guess on the masses at MIDGRID and NR (will be overwritten)
  REAL :: M = 1.3
  REAL :: MUIN = 1.0e-2	                ! Initial guess on the central mu
  
  REAL,PARAMETER :: QRELAX = 0.3    ! Used in the filter to force convergence of CHI
  LOGICAL :: ANALYTIC = .TRUE.      ! Set to true if you want to solve for phi using the analytic equation D(D(phi))=...
  
  REAL,DIMENSION(1:NR) :: RHOTV,PRTV,ETV,RHOTVJOR,ETVJOR,PRTVJOR
  REAL,DIMENSION(0:NR+1) :: CHITV,NU,MU,DCHITV,DDCHITV
  
  ! ====================================================================================
  ! Logical variables
  ! ====================================================================================
  
  LOGICAL :: DEBUG = .FALSE. 		! Set to true to print things for debugging purposes
  LOGICAL :: SINGLESURF = .TRUE.	! If true the MPSAMPLING.f90 program will build the sampling only for the RHOINI central density
  LOGICAL :: COUNTDOWN = .TRUE.	        ! Set to true to print the remaining number of configurations during the sequence computation
  LOGICAL :: CONVHELP = .FALSE.         ! If true, it activates an option for RHOCENT in HYDROEQ.f90 to help achieve convergence
  					! WARNING: if true, the final central density will be slightly different than the chosen RHOINI

  ! ====================================================================================
  ! Flags for writing output files. Put false for speed test.
  ! ====================================================================================
  
  ! Will always wrtite output for the final step of convergence loop
  LOGICAL :: LOGFILE = .TRUE.	! Set to true to save a log file
  LOGICAL :: WRT = .TRUE.       ! Write all files
  LOGICAL :: WRTF = .FALSE.     ! Write only final files
  LOGICAL :: VERBOSE = .FALSE.  ! Print all infos
  LOGICAL :: CHUP = .FALSE.     ! Write all files (including the metric test)
  LOGICAL :: WGRID = .TRUE.     ! Write the grid file
  LOGICAL :: WCONVA = .FALSE.   ! Write file to check convergence of Aphi and Atim
  LOGICAL :: WCONVC = .FALSE.   ! Write file to check convergence of chi
  LOGICAL :: IDAT               ! Write infos on the convergence for QUOC

  ! ====================================================================================
  ! Parameters for  the overall Newton-Raphson scheme (if precompiled)
  ! ====================================================================================
  
#ifdef NWTRPS
  INTEGER,PARAMETER :: QUOC = 0        !Quantity of convergence for sequences
          !QUOC=0 => central density
          !QUOC=1 => gravitational mass
          !QUOC=2 => baryonic mass
  REAL,PARAMETER :: QUCONV = 0.!RHOINI  !Value of quantity to convergence
#endif
    
  ! ====================================================================================
  ! Parameters to help or stabilize convergence
  ! ====================================================================================
  
  REAL :: QFACTOR = 0.85       ! Damping factor for press the convergence loop  qnew = qfactro*qnew + (1-qfactor)*qold
  REAL :: QFACTORMETRIC = 0.35          ! Damping factor for psi/psl for the internal convergence loop  qnew = qfactor*qnew + (1-qfactor)*qold
  REAL :: QFACTORCONF = 0.35             ! Damping factor for psi for the external convergence loop  qnew = qfactor*qnew + (1-qfactor)*qold
  REAL :: QFACTORCHI = 0.45             ! Damping factor for chi for the convergence loop  qnew = qfactro*qnew + (1-qfactor)*qold
  REAL,PARAMETER    :: QAPHI = 0.5      ! Damping factor for the convergence of the vector potential
  REAL,PARAMETER    :: EPS = 1.E-007
  
  ! ====================================================================================
  ! Parameters for the Legendre expansion & elliptic solvers
  ! ====================================================================================
  
  INTEGER,PARAMETER :: MLS = 20 !20!8       ! Number of Legendre polinomia for expansion in theta (0 = 1D)
  INTEGER,PARAMETER :: NGQ = 50 !50      ! Number of point in the Gauss quadrature (1 = 1D)
  INTEGER,PARAMETER :: MLSL = 20     ! Number of Legendre polinomia used to solve Laplace equation
  INTEGER,PARAMETER :: MLST = 20     ! Number of Legendre polinomia used to solve Maxwell-Gauss equation
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
  
END MODULE SYSTEMXNS
