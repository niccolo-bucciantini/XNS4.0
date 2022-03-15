SUBROUTINE TOVINIMOD(RHOVAR)
  ! ==============================================================================
  !   Purpose : This program solves the TOV equations in isotropic coordinates
  !             used to provide initial conditions for the XNS program
  !             Uses a relaxation method to find the convergent solution
  !             Metric: ds^2 = -exp(nu)*dt^2 + exp(mu)*( dr^2 + r^2*dOmega ) in the Einstein frame
  !             Coupling function: A(r)=exp(alpha0*(chi(r)-chi_inf)+0.5*beta0*(chi(r)-chi_inf)^2)
  !             The parameters for the shooting are the central values of mu(0) anche chi(0)
  !             First: shoot for mu(0) at given chi(0) to get fluid structure
  !             Second: on the given fluid structure compute chi(r) from second order eq with tridiag-inversion
  !             Third: reiterate first step on new chi(r) and proceed up to convergence
  ! ==============================================================================
  !
  ! In input:
  ! RHOVAR = central density in adimesional units in the Jordan frame
  !
  ! In output:
  ! RHOTV,PRTV,ETV = density,pressure and energy density of the TOV solution in the Jordan frame
  ! MU,NU = metric exponents of the TOV metric solution in the Einstein frame
  ! CHITV = scalar field in the Einstein frame
  !
  ! Output files:
  ! PROFILES.dat -> contains the radial profiles of some quantities
  ! DEBUG.dat (if DEBUG flag is true) -> contains data for debugging purposes (set by user if needed)
  !
  ! Grid: Radial, with NR = number of radial gridpoints, and either uniform or logarithmically stretched
  ! R = radial grid points (+ boundaries)
  ! DRP = increments R(I+1) - R(I)
  ! DRM = increments R(I) - R(I-1)
  ! DMAT,DUMAT,DLMAT = coefficients of the tridiagonal matrix for the Poisson Solver (from GRIDBUILD)
  !
  !
  ! From SYSTEM:
  ! GR = logical paramater (.True. = GR; .False. = STT)
  ! ALPHA0, BETA0 = parameters of the coupling function in STT
  ! CHIINF = cosmological value of the scalar field in STT
  ! EOSINT = logical parameter for EoS type (.True. = analytical plytropic; .False. = tabulated) 
  ! K1 = polytropic coefficient (for analytical polytropic EoS)
  ! GAMMA =  polytropic exponent (for analytical polytropic EoS)
  ! ANALYTIC = if True direct integration of scalalr potential from TOV equation, otherwise use the one at previous step.
  !
  ! Internal:
  ! DCHITV,DDCHITV = first and second derivative of the scalar field
  ! RHOTV,PRTV,ETV = density,pressure and energy density of the TOV solution in the Einstein frame
  ! RHOTVJOR,PRTVJOT,ETVJOR = density,pressure and energy density of the TOV solution in the Jordan frame
  ! CHITV0, MU0 = central values for the shooting
  ! COEFF,COEFFN = function temr for the shooting on MU0 and CHITV0 (zero at convergnce)
  ! RHOSURFTOV, PSURF = threshold values to set the NS surface 
  !
  ! Output files:
  ! TOVINIMOD_PROFILES.dat -> contains the radial profiles of some quantities
  ! TOVINIMOD_DEBUG.dat (if DEBUG flag is true) -> contains data for debugging purposes
  ! ==============================================================

  USE SYSTEMXNS
  USE PHYSICS
  USE FUNCTIONS
  IMPLICIT NONE

  REAL,PARAMETER :: RHOSURFTOV=1.0E-8,PSURF=1.0D-15
  REAL :: RHOVAR
  REAL :: RHOCENT,PCENT,ECENT,ENCENT,ESURF
  REAL :: MU0,CHITV0,DDCHITV0,DMU0,DCHITV0,MU0OLD,CHITV0OLD,CHITVTEMP,CHITVTEMP1,NUP,MUP
  REAL :: MOLD,MJOR,DIFFM2,DIFFM1,DM,MK,MKJOR,MSCAL
  REAL :: DCOEFF,COEFF,COEFFN,FPREV,DMU,MUPREV,COEFFPRV
  REAL :: CC,CC2,SCALCH1,SCALCH2
  
  INTEGER :: I,J,INFO
  REAL :: RHOX, EX, XSUR
  REAL,DIMENSION(6) :: Y,YIN
  
  ! Check for computational time
  CALL CPU_TIME(START)

  ! In GR set the STT coupling parameters to zero  
  IF(GR)THEN
     ALPHA0 = 0.
     BETA0 = 0.
  END IF
  
  ! ====== Set Initial Conditions =================
  MU0 = MUIN             ! Central value of mu
  RHOCENT = RHOVAR       ! Central value of rho
  ! Sets pressure  and total energy in Jordan frame
  IF(EOSINT)THEN
     CALL RHO2EOS(RHOCENT,PCENT,ECENT,ENCENT)  
  ELSE
     PCENT = K1*RHOCENT**GAMMA         
     ECENT = RHOCENT+PCENT/(GAMMA-1.)   
  END IF
  ! Set initial scalar field to zero and its first and second derivatives (+BC) in the Einstein frame
  CHITV0=0.     ! Central value of the scalar field
  DDCHITV0 = 0. ! Central value of the second drivative of the scalar field
  DO I=0,NR+1
     CHITV(I) = 0.
  END DO
  DO I=1,NR
     DCHITV(I) = (DRM(I)**2*CHITV(I+1)-DRP(I)**2*CHITV(I-1)-(DRM(I)**2-DRP(I)**2)*CHITV(I))/&
          &(DRP(I)*DRM(I)*(DRP(I)+DRM(I)))
     DDCHITV(I) = 2.*(DRM(I)*CHITV(I+1)+DRP(I)*CHITV(I-1)-(DRM(I)+DRP(I))*CHITV(I))/&
          &(DRM(I)*DRP(I)*(DRP(I)+DRM(I)))
  END DO
  DCHITV(0)     = -DCHITV(1)
  DCHITV(NR+1)  =  DCHITV(NR)*(R(NR)/R(NR+1))**2  
  DDCHITV(0)    =  DDCHITV(1)
  DDCHITV(NR+1) =  DDCHITV(NR)*(R(NR)/R(NR+1))**3
  
  ! ====== Begin serching for a solution using double shooting method ==============
  ! ====== Outer loop is on the Central Value of the Scalar Field ==================
  ! ====== Inner loop is on the Central Value of Mu in the Einstein Frame ==========
  
  ! Iterate over the Central Value of the Scalar Field
  CHITVTEMP = CHITV0 ! Set the value for convergence on CHITV0
  RELAXLOOP : DO K=1,RELIT
     IF(VERBOSE .AND. (.NOT. MPICODE))WRITE(6,*)K

     ! STEPLOOP: solves the STT TOV with a fixed CHITV (given by CHITVLOOP) until convergence of M ====
     FPREV = 0. ! Initialize to zero the function for the shooting
     STEPLOOP: DO STEP=1,MAXSTEPTV
        ! Recompute central pressure and total energy for a given central density and scalar field
        ! In the Einstein Frame 
        IF(EOSINT)THEN
           CALL RHO2EOS(RHOCENT,PCENT,ECENT,ENCENT)
           PCENT = PCENT*EXP(4*ALPHA0*(CHITV0-CHIINF)+2*BETA0*(CHITV0-CHIINF)**2)
           ECENT = (ECENT+RHOCENT)*EXP(4*ALPHA0*(CHITV0-CHIINF)+2*BETA0*(CHITV0-CHIINF)**2)
        ELSE
           PCENT=K1*RHOCENT**GAMMA*EXP(4*ALPHA0*(CHITV0-CHIINF)+2*BETA0*(CHITV0-CHIINF)**2)
           ECENT=RHOCENT*EXP(4*ALPHA0*(CHITV0-CHIINF)+2*BETA0*(CHITV0-CHIINF)**2)+PCENT/(GAMMA-1.)
        ENDIF
        ! From the center compute the vaues in RINI using Taylor Expansion of STT-TOV Eqs.
        ! We start with an initial guess for MU = MU0
        CALL EXPANSION(MU0,PCENT,ECENT,RINI,Y(1:6),CHITV0,DDCHITV0) 
        ! Using Runge-Kutta 4th compute values of all fields in the first radial point R(1)
        YIN(1:6)=Y(1:6)
        ILOC=0
        CALL RK4(YIN,6,RINI,R(1)-RINI,Y)
        ! If analytic then TOV provides also the solution for the scalar field
        ! Then reset the scalar field
        IF(ANALYTIC)THEN
           CHITV(1) = Y(5)
           DCHITV(1) = Y(6)
        END IF
        ! From the new pressure compute the density and total energy in the Einstein Frame
        CALL EOS(Y(3),RHOX,EX,CHITV(1))
        RHOTV(1)=RHOX
        ETV(1)=EX
        PRTV(1)=Y(3)
        MU(1)=Y(1)
        NU(1)=Y(4)        
        
        ! Integrate the STT-TOV equantions
        DO I=1,NR-1
           YIN(1:4)=Y(1:4)
           ! Set the scalalr field to the one previously computed 
           YIN(5)=CHITV(I)
           YIN(6)=DCHITV(I)
           ILOC=I
           ! Integrate from R(i) to R(I+1) in the Einstein Frame
           CALL RK4(YIN,6,R(I),DRP(I),Y)
           ! If analytic then TOV provides also the solution for the scalar field
           ! Then reset the scalar field
           IF(ANALYTIC)THEN
              CHITV(I+1) = Y(5)
              DCHITV(I+1) = Y(6)
           END IF
           ! From the new pressure compute the density and total energy in the Einstein Frame
           CALL EOS(Y(3),RHOX,EX,CHITV(I+1))
           RHOTV(I+1)=RHOX
           ETV(I+1)=EX
           ! Check if density below limit for stellar surface condition, if so set stellar surface radius
           IF(RHOX .GT. RHOSURFTOV)THEN
              XSUR=R(I+1)
              ISUR=I+1
           END IF
           PRTV(I+1)=Y(3)
           MU(I+1)=Y(1)
           NU(I+1)=Y(4)
           
        END DO
        ! Compute the asympthotic NU trend and its ratio with respect to the Scalar Field derivative
        NUP=(DRM(NR-1)**2*NU(NR)-DRP(NR-1)**2*NU(NR-2)-(DRM(NR-1)**2-DRP(NR-1)**2)*NU(NR-1))/&
             &(DRP(NR-1)*DRM(NR-1)*(DRP(NR-1)+DRM(NR-1)))
        ! Compute the Scalar Charge as ratio of the chi and nu derivatives
        CC=DCHITV(NR-1)/NUP
        
        ! Solves equation for M at NR and MIDGRID
        CC2=CC**2
        CALL MASSFIND(CC2,MU(MIDGRID),R(MIDGRID),MMID)
        CALL MASSFIND(CC2,MU(NR),R(NR),M)
        ! If the solution is correct the Mass shopuld be the same
        COEFF=MMID-M
        ! If confervence exit loop
        IF(ABS(COEFF) .LT. CONV )THEN
           IF(VERBOSE .AND. (.NOT. MPICODE))THEN
              WRITE(6,*)'EXITED STEPLOOP => MASSES HAVE CONVERGED',M
           END IF
           EXIT STEPLOOP
        END IF
        ! If too many step mass has not converged 
        IF(STEP .EQ. MAXSTEPTV)THEN
           IF(VERBOSE .AND. (.NOT. MPICODE))THEN
              WRITE(6,*)'EXITED STEPLOOP => MASSES HAVE NOT CONVERGED',M,MMID,DIFFM2,DIFFM1,DM
           END IF
           WRITE(6,*)'EXITED STEPLOOP => MASSES HAVE NOT CONVERGED',M,MMID
           EXIT STEPLOOP
        END IF
        ! If surface too close to outer boundary
        IF(ISUR+10 .GT. MIDGRID-1)THEN
           WRITE(6,*)"ERROR: from TOVINI, NS surf too close to outer boundary"
           WRITE(6,*)"ISUR+10="," ",ISUR+10,">","MIDGRID-1=",MIDGRID-1
           WRITE(6,*)"R_SUR=","   ",R(ISUR),"R_MID=",R(MIDGRID)
           STOP
        END IF
        
        ! If the convergence has not being achieved we need to find a new value of MU0
        ! If the function has changed sign then the best option is to use a secant
        ! If the function has changed sign then the best option is a Newton Raphson
        COEFFPRV = COEFF
        IF (COEFFPRV*FPREV .LT. 0.) THEN 
           DCOEFF=(COEFF-FPREV)/(MU0-MUPREV)
           MUPREV = MU0         ! Set previous MU0
           DMU  = -COEFF/DCOEFF ! Set increment by secant
           FPREV = COEFF        ! Set previous function for next loop
        ELSE IF (COEFFPRV*FPREV .GE. 0.) THEN
           MUPREV = MU0         ! Set previous MU0
           FPREV = COEFF        ! Set previous function for next loop
           ! Compute DMU with NewtRaph
           MU0=MU0*(1.+DELTMU0)
           ! From the center compute the vaues in RINI using Taylor Expansion of STT-TOV Eqs.
           ! We start with an initial guess for MU = MU0
           CALL EXPANSION(MU0,PCENT,ECENT,RINI,Y(1:6),CHITV0,DDCHITV0)
           YIN(1:6)=Y(1:6)
           ILOC=0
           ! Using Runge-Kutta 4th compute values of all fields in the first radial point R(1)
           CALL RK4(YIN,6,RINI,R(1)-RINI,Y)
           ! If analytic then TOV provides also the solution for the scalar field
           ! Then reset the scalar field
           IF(ANALYTIC)THEN
              CHITV(1) = Y(5)
              DCHITV(1) = Y(6)
           END IF
           ! From the new pressure compute the density and total energy in the Einstein Frame
           CALL EOS(Y(3),RHOX,EX,CHITV(1))
           RHOTV(1)=RHOX
           ETV(1)=EX
           PRTV(1)=Y(3)
           MU(1)=Y(1)
           NU(1)=Y(4)
          
           !write(6,*)
           ! Integrate the STT-TOV equantions
           DO I=1,NR-1
              YIN(1:4)=Y(1:4)
              ! Set the scalalr field to the one previously computed 
              YIN(5)=CHITV(I)
              YIN(6)=DCHITV(I)
              ILOC=I
              ! Integrate from R(i) to R(I+1) in the Einstein Frame
              CALL RK4(YIN,6,R(I),DRP(I),Y)
              ! If analytic then TOV provides also the solution for the scalar field
              ! Then reset the scalar field
              IF(ANALYTIC)THEN
                 CHITV(I+1) = Y(5)
                 DCHITV(I+1) = Y(6)
              END IF
              ! From the new pressure compute the density and total energy in the Einstein Frame
              CALL EOS(Y(3),RHOX,EX,CHITV(I+1))
              RHOTV(I+1)=RHOX
              ETV(I+1)=EX
              ! Check if density below limit for stellar surface condition, if so set stellar surface radius
              IF(RHOX .GT. RHOSURFTOV)THEN
                 XSUR=R(I+1)
                 ISUR=I+1
              END IF
              PRTV(I+1)=Y(3)
              MU(I+1)=Y(1)
              NU(I+1)=Y(4)
             
           END DO
           ! Compute the asympthotic NU trend and its ratio with respect to the Scalar Field derivative
           NUP=(DRM(NR-1)**2*NU(NR)-DRP(NR-1)**2*NU(NR-2)-(DRM(NR-1)**2-DRP(NR-1)**2)*NU(NR-1))/&
                &(DRP(NR-1)*DRM(NR-1)*(DRP(NR-1)+DRM(NR-1)))
           ! Compute the Scalar Charge as ratio of the chi and nu derivatives
           CC=DCHITV(NR-1)/NUP
           
           ! Solves equation for M at NR and MIDGRID
           CC2=CC**2
           CALL MASSFIND(CC2,MU(MIDGRID),R(MIDGRID),MMID)
           CALL MASSFIND(CC2,MU(NR),R(NR),M)
           ! Recmpute the new function in the new MU0
           COEFFN=MMID-M
           ! Compute the DMU increment with NewtRaph
           DCOEFF=(COEFFN-COEFF)/MU0/DELTMU0
           DMU  = -COEFFN/DCOEFF
        END IF
        ! Update MU0 for reiterating the MU0 convergence
        MU0=MAX(MU0+DMU,0.01)
         
     END DO STEPLOOP
    
     ! Set the BD for MU and NU
     MU(0)=(MU(2)*R(1)**2-MU0*(R(1)**2-R(2)**2))/R(2)**2 ! Quadratic interpolation to get MU in the Origin
     NU(0)=NU(1)
     MU(NR+1)=MU(NR)*R(NR)/R(NR+1)
     NU(NR+1)=NU(NR)*R(NR)/R(NR+1)
     
     ! Transform fluid quantities (density, pressure, energy) from Einstein to Jordan Frame for use in CHI solver
     PRTVJOR(1:NR)=PRTV(1:NR)/EXP(4.*ALPHA0*(CHITV(1:NR)-CHIINF)+2.*BETA0*(CHITV(1:NR)-CHIINF)**2) 
     ETVJOR(1:NR)=ETV(1:NR)/EXP(4.*ALPHA0*(CHITV(1:NR)-CHIINF)+2.*BETA0*(CHITV(1:NR)-CHIINF)**2)
     RHOTVJOR(1:NR)=RHOTV(1:NR)/EXP(4.*ALPHA0*(CHITV(1:NR)-CHIINF)+2.*BETA0*(CHITV(1:NR)-CHIINF)**2)

     ! At this point we have computed the correct TOV solution in STT for a given value of the Scalar Field in the
     ! origin CHITV0 (or for a given derivative of the scalar DRCITV(I) if not ANALYTIC).
     ! We need to recompute the Scalalr Field on top of this new matter-energy distribution to get a new
     ! value for CHITV0 (or a new radial derivative DCHITV(I)).
     ! This is done by solving the second order equation for CHI using tridiagonal-inversion
     ! Given that the Eq for CHI is non-linear we need to iterate
     
     CHITVLOOP : DO J=1,MAXSTEPCH
        ! Set the BC
        CHITV(0)=CHITV(1)
        CHITV(NR+1) =CHITV(NR)*R(NR)/R(NR+1)
        ! Compute the source term for the Scalar Field Equation to be used in tridiagonal solver
        DO I=1,NR
           MUP=(DRM(I)**2*MU(I+1)-DRP(I)**2*MU(I-1)-(DRM(I)**2-DRP(I)**2)*MU(I))/(DRP(I)*DRM(I)*(DRP(I)+DRM(I)))
           NUP=(DRM(I)**2*NU(I+1)-DRP(I)**2*NU(I-1)-(DRM(I)**2-DRP(I)**2)*NU(I))/(DRP(I)*DRM(I)*(DRP(I)+DRM(I)))
           DCHITV(I)=(DRM(I)**2*CHITV(I+1)-DRP(I)**2*CHITV(I-1)-(DRM(I)**2-DRP(I)**2)*CHITV(I))/&
                 (DRP(I)*DRM(I)*(DRP(I)+DRM(I)))
           BMAT(I)=-0.5*(NUP+MUP)*DCHITV(I)-4.*PI*(ALPHA0+BETA0*(CHITV(I)-CHIINF))*EXP(MU(I))* &
                (3.*PRTVJOR(I)-ETVJOR(I)-RHOTVJOR(I))* &
                EXP(4.*ALPHA0*(CHITV(I)-CHIINF)+2.*BETA0*(CHITV(I)-CHIINF)**2)
        END DO
        ! Set the lower-upper-central diagonal for thidiuagona matrix - DMAT, DLMAT, DUMAT from GRIDBUILD
        DS=DMAT  ! DS, DSL, DSU are dummy arrays used because DGTSV modifies the diagonals in output
        DSL=DLMAT
        DSU=DUMAT
        CALL DGTSV(NMAT,DSL,DS,DSU,BMAT,LDB,INFO)

        IF(INFO .NE. 0 .AND. (.NOT. MPICODE))WRITE(6,*)"INFO : ",INFO
        ! Upadate CHITV using also the old CHITV. QRELAX used to force convergence.
        DO I=1,NR
           CHITV(I)=(1.-QRELAX)*BMAT(I)+QRELAX*CHITV(I)    
        END DO
        ! Check for converence, starting only from second iteration.
        ! Convergence is checked on the central value of the Scalar Field
        IF((J .NE. 1) .AND. (ABS((CHITV(2)*R(1)**2-CHITV(1)*R(2)**2)/(R(1)**2-R(2)**2)-CHITVTEMP1) .LT. CONVT))THEN
           IF(VERBOSE .AND. (.NOT. MPICODE))THEN
              WRITE(6,*)'EXITED CHITVLOOP => CHITV HAS CONVERGED'
           END IF
           EXIT CHITVLOOP
        END IF
        ! If not converged
        IF(J==MAXSTEPCH)THEN 
           IF(VERBOSE .AND. J==MAXSTEPCH .AND. (.NOT. MPICODE))THEN
              WRITE(6,*)'EXITED CHITVLOOP WITHOUT CONVERGENCE: CHITV0 CONVERGED TO A PRECISION OF ',&
                   &ABS((CHITV(2)*R(1)**2-CHITV(1)*R(2)**2)/(R(1)**2-R(2)**2)-CHITVTEMP1)
           END IF
           WRITE(6,*)'EXITED CHITVLOOP WITHOUT CONVERGENCE: CHITV0 CONVERGED TO A PRECISION OF ',&
                &ABS((CHITV(2)*R(1)**2-CHITV(1)*R(2)**2)/(R(1)**2-R(2)**2)-CHITVTEMP1)
        END IF
        ! Reset the Central Value for the Scalar Field
        CHITVTEMP1=(CHITV(2)*R(1)**2-CHITV(1)*R(2)**2)/(R(1)**2-R(2)**2)
        ! Compute the BC on the derivative
        DCHITV(0)    = -DCHITV(1)
        DCHITV(NR+1) =  DCHITV(NR)*(R(NR)/R(NR+1))**2
        
     END DO CHITVLOOP
     
     ! Compute the derivative of the new scalar field
     CHITV(0)=CHITV(1)
     CHITV(NR+1) =CHITV(NR)*R(NR)/R(NR+1)
     DO I=1,NR
        DCHITV(I)  = (DRM(I)**2*CHITV(I+1)-DRP(I)**2*CHITV(I-1)-(DRM(I)**2-DRP(I)**2)*CHITV(I))/&
             & (DRP(I)*DRM(I)*(DRP(I)+DRM(I)))
        DDCHITV(I) = 2.*(DRM(I)*CHITV(I+1)+DRP(I)*CHITV(I-1)-(DRM(I)+DRP(I))*CHITV(I))/&
             &(DRM(I)*DRP(I)*(DRP(I)+DRM(I)))
     END DO
     DCHITV(0)    = -DCHITV(1)
     DCHITV(NR+1) =  DCHITV(NR)*(R(NR)/R(NR+1))**2
     DDCHITV(0)    = DDCHITV(1)
     DDCHITV(NR+1) =  DDCHITV(NR)*(R(NR)/R(NR+1))**3
     ! Compute the new value of the Scalar field in the center with Quadratic Interpolation
     CHITV0=(CHITV(2)*R(1)**2-CHITV(1)*R(2)**2)/(R(1)**2-R(2)**2)
     ! And its second derivative in the center
     DDCHITV0=2.*(CHITV(1)-CHITV0)/R(1)**2

     ! Check the convergence on CHITV0 with respect to the previous one
     COEFFN=ABS(CHITV0-CHITVTEMP)
     ! If converged
     IF((K .NE. 1) .AND. (ABS(CHITV0-CHITVTEMP) .LT. CONVT))THEN
        IF(VERBOSE .AND. (.NOT. MPICODE))THEN
           WRITE(6,*)'EXITED RELAXLOOP => RELAXATION HAS CONVERGED'
        END IF
        EXIT RELAXLOOP
     END IF
     ! If not converged
     IF(K==RELIT)THEN
        IF(VERBOSE .AND. (.NOT. MPICODE))THEN
           WRITE(6,*)'EXITED RELAXLOOP WITHOUT CONVERGENCE: CHITV0 CONVERGED TO A PRECISION OF ',&
             &ABS(CHITV0-CHITVTEMP)
        END IF
        WRITE(6,*)'EXITED RELAXLOOP WITHOUT CONVERGENCE: CHITV0 CONVERGED TO A PRECISION OF ',&
             &ABS(CHITV0-CHITVTEMP)
     END IF
     ! Reset for next loop
     CHITVTEMP=CHITV0
     
  END DO RELAXLOOP

  ! Add offset to rescale the NU metric coefficient
  DO I=0,NR+1
     NU(I)=NU(I)-NU(NR)+2.*LOG((2*R(NR)-(M+CC2*M**2*1./R(NR+1)-CC2*M**3/6.*1./R(NR+1)**2+&
          &CC2*(1.+3.*CC2)*M**4/12.*1./R(NR+1)**3-M**5/120.*CC2*(3.+11.*CC2)*1./R(NR+1)**4))/&
          &(2*R(NR)+M-CC2*M**2*1./R(NR+1)-CC2*M**3/6.*1./R(NR+1)**2-CC2*(1.+3.*CC2)*M**4/12./&
          &R(NR+1)**3-M**5/120.*CC2*(3.+11.*CC2)*1./R(NR+1)**4))
  END DO

  ! Compute the Global Quantities of the Model
  ! Calculates the Komar mass in the Einsteih frames by integrating over the star
  MK=0.
  DO I=1,ISUR
     MK=MK+4.*PI*(ETV(I)+RHOTV(I)+3.*PRTV(I))*EXP(NU(I)/2.)*EXP(3.*MU(I)/2.)*R(I)**2*DR(I)
  END DO
  ! Compute the scalar mass and the Komar mass in Jprdn Frame
  MSCAL=0.5*R(NR)*(CHITV(NR)-CHIINF)
  MKJOR=MK-2.*ALPHA0*MSCAL
  ! Calculates the ADM mass in the Jordan frame
  CALL MASSFIND(CC2,MU(NR)+2.*ALPHA0*(CHITV(NR)-CHIINF)+BETA0*(CHITV(NR)-CHIINF)**2,R(NR),MJOR)
  ! Compute the Scalar charge integrating over the star
  SCALCH1=0.
  DO I=1,ISUR+10
     SCALCH1=SCALCH1+2.*PI*((ALPHA0+BETA0*(CHITV(I)-CHIINF))*(3.*PRTV(I)-ETV(I)-RHOTV(I))*&
          &EXP(3.*MU(I)/2.))*R(I)**2*DR(I)
  END DO

  ! =========================== Write the outcome of STT-TOV  ====================
  ! Write to screen 
  IF(.NOT. MPICODE)THEN
     WRITE(6,*)''
     WRITE(6,*)'==================== TOV SOLUTION IN STT ===================='
     WRITE(6,*)''
     WRITE(6,*)'Central density (J)      ','    ',RHOCENT/MBARYONFC	! Final central density	
     WRITE(6,*)'Central density (E)      ','    ',RHOCENT*EXP(4*ALPHA0*(CHITV0-CHIINF)+2*BETA0*(CHITV0-CHIINF)**2)/MBARYONFC
     WRITE(6,*)'Central scalar field (E) ','    ',CHITV0! Central scalar field and mu
     WRITE(6,*)'Central mu (E)           ','    ',MU0
     WRITE(6,*)'Central exp(mu) (E)      ','    ',EXP(MU0)
     WRITE(6,*)'ADM mass (E)             ','    ',M	
     WRITE(6,*)'ADM mass (J)             ','    ',MJOR
     WRITE(6,*)'Komar mass (E)           ','    ',MK
     WRITE(6,*)'Komar mass (J)           ','    ',MKJOR
     WRITE(6,*)'Scalar mass              ','    ',MSCAL
     WRITE(6,*)'Scalar charge (E)        ','    ',CC
     WRITE(6,*)'Iso. Radius (E)          ','    ',R(ISUR)                  ! Stellar radius in Isotropic coordinates
     WRITE(6,*)'Sch. Radius (E)          ','    ',R(ISUR)*EXP(MU(ISUR)/2.) ! Stellar radius in Schwarzschild coordinates
     WRITE(6,*)'Isurf                    ','    ',ISUR
     WRITE(6,*)'Midgrid                  ','    ',MIDGRID
     WRITE(6,*)'Grid radius (E)          ','    ',R(NR)	! To check that isurf<midgrid
     WRITE(6,*)'COEFF (mass)             ','    ',COEFF	! Convergence quantities
     WRITE(6,*)'COEFFN (scalar field)    ','    ',COEFFN
     WRITE(6,*)''
     WRITE(6,*)'============================================================='
     WRITE(6,*)''
  END IF

  ! =========================== Prepare output for XNSMAIN  ====================  
  ! Converts PRTV, ETV and RHOTV to the Jordan frame, ready to be used by XNS
  PRTV(1:NR)=PRTV(1:NR)/EXP(4.*ALPHA0*(CHITV(1:NR)-CHIINF)+2.*BETA0*(CHITV(1:NR)-CHIINF)**2)
  ETV(1:NR)=ETV(1:NR)/EXP(4.*ALPHA0*(CHITV(1:NR)-CHIINF)+2.*BETA0*(CHITV(1:NR)-CHIINF)**2)
  RHOTV(1:NR)=RHOTV(1:NR)/EXP(4.*ALPHA0*(CHITV(1:NR)-CHIINF)+2.*BETA0*(CHITV(1:NR)-CHIINF)**2)
  ! Overwrite with original mu0 value for the log file
  MU0=MUIN

  IF(DEBUG)THEN
     OPEN(14,FILE='TOVINIMOD_DEBUG.dat')
     WRITE(6,*)'Here go the infos for debug'
     CLOSE(14)
  END IF
  
  ! Write STT-TOV solution to file
  OPEN(13,FILE='TOVINIMOD_PROFILES.dat')
  DO IZ=1,NR
     WRITE(13,*)IZ,RHOTV(IZ)/MBARYONFC,PRTV(IZ),NU(IZ),MU(IZ),CHITV(IZ)
  END DO
  CLOSE(13)

  ! Find computational time for TOVINI
  CALL CPU_TIME(FINISH)
  IF(.NOT. MPICODE)PRINT '("Elapsed CPU time = ",F6.3," seconds")',FINISH-START
  
END SUBROUTINE TOVINIMOD

! **************************************************************
! **************************************************************
  
SUBROUTINE EXPANSION(MU0,PCENT,ECENT,RINIT,YINI,CHITV0,DDCHITV0)
  ! ==============================================================
  ! Purpose: Expansion of the TOV eqs around the origin, to start the integration
  !
  ! MU0 = central value of the metric term MU
  ! PCENT, ECENT = pressure and total energy at the center in Einstein frame
  ! CHTV0, DDCHITV0 = scalalr field and its second derivative in the center
  ! RINIT = radius to which the expansion is required
  ! YINI = values of the quantities expanded in RINIT
  ! ==============================================================

  USE SYSTEMXNS, ONLY : PI,ALPHA0,BETA0,CHIINF,ANALYTIC,GR
  IMPLICIT NONE
  REAL,DIMENSION(6) :: YINI
  REAL :: MU0,CHITV0,DDCHITV0,PCENT,ECENT
  REAL :: RINIT
  REAL :: P0,E0,DDMU0,DDP0,DDNU0
  
  E0 = ECENT
  P0 = PCENT
  ! Check is GR or STT
  IF(GR)THEN
     DDP0 = -2.*PI/3.*EXP(MU0)*(E0**2+3*P0**2+4*E0*P0)
     DDCHITV0 = 0.
  ELSE
     DDP0 = -2.*PI/3.*EXP(MU0)*(E0**2*(1+(ALPHA0+BETA0*(CHITV0-CHIINF))**2)+3*P0**2*(1+3*(ALPHA0+BETA0*(CHITV0-CHIINF))**2) &
          +2*E0*P0*(2-3*(ALPHA0+BETA0*(CHITV0-CHIINF))**2))
     ! If ANALYTIC then recompute the second derivative of the scalar field
     IF(ANALYTIC)THEN
        DDCHITV0 = EXP(MU0)/6.*4*PI*(ALPHA0+BETA0*(CHITV0-CHIINF))*(E0-3*P0)
     END IF
  END IF
  DDMU0 = -4./3.*EXP(MU0)*PI*E0
  DDNU0 = 4./3.*EXP(MU0)*PI*(E0+3*P0)

  ! Compute values in RINI using the Taylor Expansion 
  YINI(1) = MU0 + DDMU0*RINIT**2
  YINI(2) = 2*RINIT*DDMU0
  YINI(3) = P0 + DDP0*RINIT**2
  YINI(4) = DDNU0*RINIT**2
  YINI(5) = CHITV0 + 0.5*DDCHITV0*RINIT**2
  YINI(6) = RINIT*DDCHITV0
  
  RETURN
  
END SUBROUTINE EXPANSION

! **************************************************************
! **************************************************************

SUBROUTINE TOVEQS(R,Y,DY)
  ! ==============================================================
  ! Purpose: solve the STT-TOV Equations in Isotropic Coordinates
  ! Works in the Einstein Frame.
  !
  ! ILOC = Index of the local radial position: R = R(ILOC)
  ! R = Starting Radius
  ! Y(1) = mu
  ! Y(2) = mu'
  ! Y(3) = p
  ! Y(4) = nu
  ! Y(5) = CHI
  ! Y(6) = CHI'
  ! ==============================================================

  USE SYSTEMXNS, ONLY : GR,ANALYTIC,PI,ALPHA0,BETA0,CHIINF,DCHITV,DDCHITV,ILOC
  USE PHYSICS, ONLY : EOS
  IMPLICIT NONE
  !LOGICAL :: GR
  REAL,DIMENSION(6) :: Y,DY
  REAL :: R,PRSX,RHOX,ENX,CHIX,ETOT
  
  ! Given a value of the Pressure and Scalar Field in the Einstein frame
  ! derive the internal energy and density in the Einstein frame
  PRSX = Y(3)
  CHIX = Y(5)
  CALL EOS(PRSX,RHOX,ENX,CHIX)
  ! Total Enengy in the Einstein Frame
  ETOT=RHOX+ENX
  
  IF(GR)THEN ! GR TOV
     DY(1) = Y(2)
     DY(2) = -8.*PI*ETOT*EXP(Y(1)) -2./R*Y(2) -0.25*Y(2)**2
     DY(3) = -(8.*PI*Y(3)*EXP(Y(1)) -0.25*Y(2)**2 -Y(2)/R)/(Y(2)/2.+1./R)*(ETOT+Y(3))/2.
     DY(4) = +(8.*PI*Y(3)*EXP(Y(1)) -0.25*Y(2)**2 -Y(2)/R)/(Y(2)/2.+1./R)
     DY(5) = 0.
     DY(6) = 0.
  ELSE       ! STT TOV
     DY(1) = Y(2)
     DY(2) = -8.*PI*ETOT*EXP(Y(1)) -2./R*Y(2) -0.25*Y(2)**2 -Y(6)**2
     DY(3) = -(8.*PI*Y(3)*EXP(Y(1)) -0.25*Y(2)**2 -Y(2)/R +Y(6)**2)/(Y(2)/2.+1./R)*(ETOT+Y(3))/2.&
          +(ALPHA0+BETA0*(Y(5)-CHIINF))*(3*Y(3)-ETOT)*Y(6)
     DY(4) = +(8.*PI*Y(3)*EXP(Y(1)) -0.25*Y(2)**2 -Y(2)/R +Y(6)**2)/(Y(2)/2.+1./R)
     IF (ANALYTIC)THEN !Integrate also for the scalalr filed
        DY(5) = Y(6)
        DY(6) = -Y(6)*(DY(4)/2.+Y(2)/2.+2./R) -EXP(Y(1))*4*PI*(ALPHA0+BETA0*(Y(5)-CHIINF))*(3*Y(3)-ETOT)
     ELSE ! Use as local approximantion the average between R(ILOC) and  R(ILOC+1)
        DY(5) = 0.5*(DCHITV(ILOC) + DCHITV(ILOC+1))
        DY(6) = 0.5*(DDCHITV(ILOC) + DDCHITV(ILOC+1))
     ENDIF
  END IF
  
  RETURN
  
END SUBROUTINE TOVEQS

! **************************************************************
! **************************************************************

SUBROUTINE RK4(Y,N,RI,H,YOUT)
  ! ==============================================================
  ! Purpose: solve the TOV-ODE from RI to RI + H using
  ! 4th-order Runge-Kutta Integrator 
  !
  ! N = number of equations
  ! RI = starting position
  ! Y(N) = starting values in RI
  ! H = increment
  ! YOUT(N) = final values in RI + H
  ! ==============================================================

  IMPLICIT NONE
  INTEGER :: N
  REAL::  H,RI,Y(N),YOUT(N)
  INTEGER :: I
  REAL :: H2,DYDR1(N),Y1(N),DYDR2(N),Y2(N),DYDR3(N),Y3(N),DYDR4(N)
  
  H2=H/2.    
  CALL TOVEQS(RI,Y,DYDR1)
  Y1(1:N)=Y(1:N)+H2*DYDR1(1:N)
  CALL TOVEQS(RI+H2,Y1,DYDR2)
  Y2(1:N)=Y(1:N)+H2*DYDR2(1:N)
  CALL TOVEQS(RI+H2,Y2,DYDR3)
  Y3(1:N)=Y(1:N)+H*DYDR3(1:N)
  CALL TOVEQS(RI+H,Y3,DYDR4)
  YOUT(1:N)=Y(1:N)+(DYDR1(1:N)+2*DYDR2(1:N)+2*DYDR3(1:N)+DYDR4(1:N))*H/6.
  
  RETURN
  
END SUBROUTINE RK4

! **************************************************************
! **************************************************************

SUBROUTINE MASSFIND(CC2,MUREF,RREF,MOUT)
  ! ==============================================================
  ! Purpose: Find the Mass at a given radius from the metric terms
  ! and the scalar charge, using vacuum expansion for exp(mu)
  !
  ! CC2 = square of the scalar charge
  ! RREF = radius
  ! MUREF = value of MU in RRREF
  ! MUOT = Mass
  ! ==============================================================
  
  USE SYSTEMXNS, ONLY : CONV
  IMPLICIT NONE
  REAL :: CC2,MUREF,RREF,MOUT
  REAL :: A1,A2,A3,A4,A5,DIFFM,DDIFFM,RREF1
  INTEGER,PARAMETER :: MAXITER = 100
  INTEGER :: J
  
  RREF1 = 1./RREF
  ! Coefficients of the Vacuum Expansion up to 5th-order
  A1 = 1.
  A2 = -CC2*RREF1
  A3 = -CC2/6.*RREF1**2
  A4 = -CC2*(1.+3.*CC2)/12.*RREF1**3
  A5 = -CC2*(3.+11.*CC2)*RREF1**4/120.
  !Solve by Newt-Raph. the Equation for the Vacuum Expansion of exp(mu)
  MOUT=1.4 ! Initial Guess
  DO J=1,MAXITER
     DIFFM = MOUT*(A1 + A2*MOUT + A3*MOUT**2 + A4*MOUT**3. + A5*MOUT**4) - 2.*RREF*(EXP(MUREF/4.)-1)
     IF(ABS(DIFFM) .LT. CONV)EXIT
     DDIFFM=(A1 + 2*A2*MOUT + 3*A3*MOUT**2 + 4*A4*MOUT**3. + 5*A5*MOUT**4)
     MOUT=MOUT-DIFFM/DDIFFM
  END DO

  RETURN
  
END SUBROUTINE MASSFIND

! **************************************************************
! **************************************************************
