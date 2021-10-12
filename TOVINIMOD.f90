! PROGRAM TEST
!
! 	USE VAR_GENERAL
! 	USE VAR_RELAX
! 	USE VAR_SHOOT
! 	USE FUNCTIONS
! 	IMPLICIT NONE
!
! 	CALL TOVINISTT(RHOCENTCHOICE)
!
! END PROGRAM TEST

SUBROUTINE TOVINIMOD(RHOVAR) !!! QUI

! ==============================================================
!   Purpose : This program solves the TOV equations in isotropic coordinates
!             used to provide initial conditions for the XNS program
!             Uses a relaxation method to find the convergent solution
!             EoS: polytropic in the Jordan frame
!             Metric: ds^2 = -exp(nu)*dt^2 + exp(mu)*( dr^2 + r^2*dOmega ) in the Einstein frame
!			  Coupling function: A(CHITV)=exp(alpha0*(CHITV-CHIINF)+0.5*beta0*(CHITV-CHIINF)^2)
! ==============================================================
!
! Grid:
! Radial, with NR = number of radial gridpoints, and either uniform or logarithmically stretched
! R = radial grid points (+ boundaries)
! DR = increments
!
! In input:
! RHOVAR = central density in adimasional units in the Jordan frame
! ALPHA0, BETA0 = parameters of the coupling function
! CHIINF = cosmological value of the scalar field
! K1 = polytropic coefficient
! GAMMA =  polytropic exponent
!
! In output:
! RHOTV,PRTV,ETV = density,pressure and energy density of the TOV solution in the Einstein frame
! MU,NU = metric exponents of the TOV metric solution in the Einstein frame
! CHITV = scalar field in the Einstein frame
!
! Output files:
! PROFILES.dat -> contains the radial profiles of some quantities
! DEBUG.dat (if DEBUG flag is true) -> contains data for debugging purposes
! ==============================================================

USE SYSTEMXNS
USE FUNCTIONS
IMPLICIT NONE

REAL :: RHOCENT,PCENT,ECENT,ENCENT,RHOSURFTOV=1.0E-8,ESURF,PSURF=1.0D-15
REAL :: RHOVAR
REAL :: MU0,CHITV0,DDCHITV0,DMU0,DCHITV0,MU0OLD,CHITV0OLD,CHITVTEMP,CHITVTEMP1,NUP,MUP
REAL :: MOLD,MJOR,DIFFM2,DIFFM1,DM,MK,MKJOR,MSCAL
REAL :: DCOEFF,COEFF,COEFFN,FPREV,DMU,MUPREV
REAL :: CC,CC2,SCALCH1,SCALCH2

INTEGER :: I,J,INFO
REAL :: RHOX, EX, XSUR
REAL,DIMENSION(6) :: Y,YIN

! REAL,DIMENSION(0:NR) :: ACOEFF,BCOEFF,CCOEFF   ! Coefficients of the tridiagonal matrix A used in DGTSV
! INTEGER :: LDB=NR,N=NR									 ! Dimensions of the elements of the matrix equation
! REAL :: B(NR), D(NR), DL(NR-1), DU(NR-1)					 ! Diagonals of the A matrix and source vector B
! REAL :: DS(NR), DSL(NR-1), DSU(NR-1)						 ! Dummy diagonals
REAL :: FACT,TEMP											 ! Used in DGTSV


CALL CPU_TIME(START)

IF(GR)THEN
	ALPHA0=0.
	BETA0=0.
END IF

! MMID=1.9			! Initial guess on the masses at MIDGRID and NR
! M=1.9
! MUIN=0.01			! Initial guess on the central mu
! CHITVINI=0.001		! Initial guess on the central CHITV

! ============= Grid setting =============


! ========================================

MU0=MUIN
RHOCENT=RHOVAR
WRITE(6,*)DRM,DRP
IF(EOSINT)THEN
   CALL RHO2EOS(RHOCENT,PCENT,ECENT,ENCENT)
ELSE
   PCENT=K1*RHOCENT**GAMMA          ! Sets pressure in jordan frame
   ECENT=RHOCENT+PCENT/(GAMMA-1.) 	 ! Sets total energy in jordan frame
END IF
DO I=0,NR+1
	CHITV(I)=0.
END DO
CHITV0=0.

DO I=1,NR
	DCHITV(I)=(DRM(I)**2*CHITV(I+1)-DRP(I)**2*CHITV(I-1)-(DRM(I)**2-DRP(I)**2)*CHITV(I))/&
		&(DRP(I)*DRM(I)*(DRP(I)+DRM(I)))
	DDCHITV(I)=2.*(DRM(I)*CHITV(I+1)+DRP(I)*CHITV(I-1)-(DRM(I)+DRP(I))*CHITV(I))/&
		&(DRM(I)*DRP(I)*(DRP(I)+DRM(I)))
END DO

DCHITV(0)    = -DCHITV(1)
DCHITV(NR+1) =  DCHITV(NR)*(R(NR)/R(NR+1))**2
DDCHITV(0)    = DDCHITV(1)
DDCHITV(NR+1) =  DDCHITV(NR)*(R(NR)/R(NR+1))**3

CHITV0   = 0.
DCHITV0  = 0.
DDCHITV0 = 0.

! ==== RELAXLOOP: relaxation loop until convergence of CHITV0 ====

RELAXLOOP : DO K=1,RELIT

! ==== STEPLOOP: solves the STT TOV with a fixed CHITV (given by CHITVLOOP) until convergence of M ====

   IF(VERBOSE .AND. (.NOT. MPICODE))WRITE(6,*)K

   FPREV = 0.
   STEPLOOP: DO STEP=1,MAXSTEP

      IF(EOSINT)THEN
         CALL RHO2EOS(RHOCENT,PCENT,ECENT,ENCENT)
         PCENT = PCENT*EXP(4*ALPHA0*(CHITV0-CHIINF)+2*BETA0*(CHITV0-CHIINF)**2)
         ECENT = (ECENT+RHOCENT)*EXP(4*ALPHA0*(CHITV0-CHIINF)+2*BETA0*(CHITV0-CHIINF)**2)
      ELSE
        PCENT=K1*RHOCENT**GAMMA*EXP(4*ALPHA0*(CHITV0-CHIINF)+2*BETA0*(CHITV0-CHIINF)**2)
        ECENT=RHOCENT*EXP(4*ALPHA0*(CHITV0-CHIINF)+2*BETA0*(CHITV0-CHIINF)**2)+PCENT/(GAMMA-1.)
!         CALL RHO2EOS(RHOCENT,PCENT,ECENT,ENCENT)
!          PCENT = PCENT*EXP(4*ALPHA0*(CHITV0-CHIINF)+2*BETA0*(CHITV0-CHIINF)**2)
!          ECENT = (ECENT+RHOCENT)*EXP(4*ALPHA0*(CHITV0-CHIINF)+2*BETA0*(CHITV0-CHIINF)**2)
      ENDIF
      CALL EXPANSION(MU0,PCENT,ECENT,RINI,Y(1:6),CHITV0)
!       print(Y(1))
      YIN(1:6)=Y(1:6)
      ILOC=0
!       write(6,*)yin
      CALL RK4(YIN,6,RINI,R(1)-RINI,Y)
!       write(6,*)y
      CALL EOS(Y(3),RHOX,EX,CHITV(1))
! 	  write(6,*)y(3),rhox,ex
! 	stop
      RHOTV(1)=RHOX
      ETV(1)=EX
      PRTV(1)=Y(3)
      MU(1)=Y(1)
      NU(1)=Y(4)

      IF(ANALYTIC)THEN
         CHITV(1) = Y(5)
         DCHITV(1) = Y(6)
      END IF


      ! TOV integration

      DO I=1,NR-1
         YIN(1:4)=Y(1:4)
         YIN(5)=CHITV(I)
         YIN(6)=DCHITV(I)
         ILOC=I

         CALL RK4(YIN,6,R(I),DRP(I),Y)
         CALL EOS(Y(3),RHOX,EX,CHITV(I+1))
         RHOTV(I+1)=RHOX
         ETV(I+1)=EX
!          write(6,*)RHOX,I
         IF(RHOX .GT. RHOSURFTOV)THEN
            XSUR=R(I+1)
            ISUR=I+1
         END IF
         PRTV(I+1)=Y(3)
         MU(I+1)=Y(1)
         NU(I+1)=Y(4)
         IF(ANALYTIC)THEN
            CHITV(I) = Y(5)
            DCHITV(I) = Y(6)
         END IF
      END DO

      NUP=(DRM(NR-1)**2*NU(NR)-DRP(NR-1)**2*NU(NR-2)-(DRM(NR-1)**2-DRP(NR-1)**2)*NU(NR-1))/&
           &(DRP(NR-1)*DRM(NR-1)*(DRP(NR-1)+DRM(NR-1)))
      CC=DCHITV(NR-1)/NUP
      CC2=CC**2

      ! Solves equation for M at NR and MIDGRID
      CALL MASSFIND(CC2,MU(MIDGRID),R(MIDGRID),MMID)
      CALL MASSFIND(CC2,MU(NR),R(NR),M)

      COEFF=MMID-M

      IF(ABS(COEFF) .LT. CONV )THEN
         IF(VERBOSE .AND. (.NOT. MPICODE))THEN
            WRITE(6,*)'EXITED STEPLOOP => MASSES HAVE CONVERGED',M
         END IF
         EXIT STEPLOOP
      END IF
      IF(STEP .EQ. MAXSTEP)THEN
         IF(VERBOSE .AND. (.NOT. MPICODE))THEN
            WRITE(6,*)'EXITED STEPLOOP => MASSES HAVE NOT CONVERGED',M,MMID,DIFFM2,DIFFM1,DM
         END IF
         EXIT STEPLOOP
      END IF
      IF(ISUR+10 .GT. MIDGRID-1)THEN
         WRITE(6,*)"ERROR: from TOVINI, NS surf too close to outer boundary"
         WRITE(6,*)"ISUR+10="," ",ISUR+10,">","MIDGRID-1=",MIDGRID-1
         WRITE(6,*)"R_SUR=","   ",R(ISUR),"R_MID=",R(MIDGRID)
         STOP
      END IF

      ! Second guess of mu0
      IF (COEFF*FPREV .LT. 0.) THEN ! IF THE FUNCTION CHANGE SIGN USE A SECANT SCHEME FOR STABILITY
         DCOEFF=(COEFF-FPREV)/(MU0-MUPREV)
         DMU  = -COEFF/DCOEFF
!           WRITE(6,*)mu0,COEFF,muprev,FPREV
         FPREV = COEFF
         MUPREV = MU0

      ELSE IF (COEFF*FPREV .GE. 0.) THEN ! IF THE FUNCTION NOT CHANGE SIGN USE A NEWTRAP FOR FASTNESS
         FPREV = COEFF
         MUPREV = MU0
         MU0=MU0*(1.+DELTMU0)

         CALL EXPANSION(MU0,PCENT,ECENT,RINI,Y(1:6),CHITV0)

         YIN(1:6)=Y(1:6)
         ILOC=0

         CALL RK4(YIN,6,RINI,R(1)-RINI,Y)

         CALL EOS(Y(3),RHOX,EX,CHITV(1))

         RHOTV(1)=RHOX
         ETV(1)=EX
         PRTV(1)=Y(3)
         MU(1)=Y(1)
         NU(1)=Y(4)

         IF(ANALYTIC)THEN
            CHITV(1) = Y(5)
            DCHITV(1) = Y(6)
         END IF


         ! TOV integration

         DO I=1,NR-1
            YIN(1:4)=Y(1:4)
            YIN(5)=CHITV(I)
            YIN(6)=DCHITV(I)
            ILOC=I

            CALL RK4(YIN,6,R(I),DRP(I),Y)
            CALL EOS(Y(3),RHOX,EX,CHITV(I+1))
            RHOTV(I+1)=RHOX
            ETV(I+1)=EX
            IF(RHOX .GT. RHOSURFTOV)THEN
               XSUR=R(I+1)
               ISUR=I+1
            END IF

            ! IF(I==222 .AND. STEP .ge. 0)THEN
            !              WRITE(6,*)YIN
            ! 			 WRITE(6,*)Y
            ! 			 WRITE(6,*)RHOTV(I),RHOTV(I+1)
            ! 		! 	 STOP
            ! 		 ENDIF
            PRTV(I+1)=Y(3)
            MU(I+1)=Y(1)
            NU(I+1)=Y(4)
            IF(ANALYTIC)THEN
               CHITV(I) = Y(5)
               DCHITV(I) = Y(6)
            END IF
         END DO

         NUP=(DRM(NR-1)**2*NU(NR)-DRP(NR-1)**2*NU(NR-2)-(DRM(NR-1)**2-DRP(NR-1)**2)*NU(NR-1))/&
              &(DRP(NR-1)*DRM(NR-1)*(DRP(NR-1)+DRM(NR-1)))
         CC=DCHITV(NR-1)/NUP
         CC2=CC**2

         ! Solves equation for M at NR and MIDGRID
         CALL MASSFIND(CC2,MU(MIDGRID),R(MIDGRID),MMID)
         CALL MASSFIND(CC2,MU(NR),R(NR),M)

         COEFFN=MMID-M
         DCOEFF=(COEFFN-COEFF)/MU0/DELTMU0
         DMU  = -COEFFN/DCOEFF
      END IF
      MU0=MAX(MU0+DMU,0.01)

   END DO STEPLOOP

   MU(0)=(MU(2)*R(1)**2-MU0*(R(1)**2-R(2)**2))/R(2)**2		! Quadratic interpolation
   NU(0)=NU(1)
   MU(NR+1)=MU(NR)*R(NR)/R(NR+1)
   NU(NR+1)=NU(NR)*R(NR)/R(NR+1)

   PRTVJOR(1:NR)=PRTV(1:NR)/EXP(4.*ALPHA0*(CHITV(1:NR)-CHIINF)+2.*BETA0*(CHITV(1:NR)-CHIINF)**2) ! Transforms PRTV, ETV and RHOTV to the Jordan frame
   ETVJOR(1:NR)=ETV(1:NR)/EXP(4.*ALPHA0*(CHITV(1:NR)-CHIINF)+2.*BETA0*(CHITV(1:NR)-CHIINF)**2)
   RHOTVJOR(1:NR)=RHOTV(1:NR)/EXP(4.*ALPHA0*(CHITV(1:NR)-CHIINF)+2.*BETA0*(CHITV(1:NR)-CHIINF)**2)

   ! ==== CHITVLOOP: solves the equation for CHITV until convergence on CHITV0 ====

   CHITVLOOP : DO J=1,RELIT

      CHITV(0)=CHITV(1)
      CHITV(NR+1) =CHITV(NR)*R(NR)/R(NR+1)

      DO I=1,NR

         MUP=(DRM(I)**2*MU(I+1)-DRP(I)**2*MU(I-1)-(DRM(I)**2-DRP(I)**2)*MU(I))/(DRP(I)*DRM(I)*(DRP(I)+DRM(I)))
         NUP=(DRM(I)**2*NU(I+1)-DRP(I)**2*NU(I-1)-(DRM(I)**2-DRP(I)**2)*NU(I))/(DRP(I)*DRM(I)*(DRP(I)+DRM(I)))
         DCHITV(I)=(DRM(I)**2*CHITV(I+1)-DRP(I)**2*CHITV(I-1)-(DRM(I)**2-DRP(I)**2)*CHITV(I))/&
              &(DRP(I)*DRM(I)*(DRP(I)+DRM(I)))
         BMAT(I)=-0.5*(NUP+MUP)*DCHITV(I)-4.*PI*(ALPHA0+BETA0*CHITV(I))*EXP(MU(I))*(3.*PRTVJOR(I)-ETVJOR(I)-RHOTVJOR(I))*&
              &EXP(4.*ALPHA0*(CHITV(I)-CHIINF)+2.*BETA0*(CHITV(I)-CHIINF)**2)
      END DO
      DS=DMAT  ! DS, DSL, DSU are dummy arrays used because DGTSV modifies the diagonals in output
      DSL=DLMAT
      DSU=DUMAT
      CALL DGTSV(NMAT,DSL,DS,DSU,BMAT,LDB,INFO)

      IF(INFO .NE. 0 .AND. (.NOT. MPICODE))WRITE(6,*)"INFO : ",INFO

      DO I=1,NR
         CHITV(I)=(1.-QRELAX)*BMAT(I)+QRELAX*CHITV(I)    ! Filter that uses the updates CHITV using also the old CHITV; used to force convergence
      END DO

      IF((J .NE. 1) .AND. (ABS((CHITV(2)*R(1)**2-CHITV(1)*R(2)**2)/(R(1)**2-R(2)**2)-CHITVTEMP1) .LT. CONV2))THEN
         IF(VERBOSE .AND. (.NOT. MPICODE))THEN
            WRITE(6,*)'EXITED CHITVLOOP => CHITV HAS CONVERGED'
         END IF
         EXIT CHITVLOOP
      END IF
      IF(VERBOSE .AND. J==RELIT .AND. (.NOT. MPICODE))THEN
         WRITE(6,*)'EXITED CHITVLOOP WITHOUT CONVERGENCE: CHITV0 CONVERGED TO A PRECISION OF ',&
              &ABS((CHITV(2)*R(1)**2-CHITV(1)*R(2)**2)/(R(1)**2-R(2)**2)-CHITVTEMP1)
      END IF

      CHITVTEMP1=(CHITV(2)*R(1)**2-CHITV(1)*R(2)**2)/(R(1)**2-R(2)**2)

      DCHITV(0)    = -DCHITV(1)
      DCHITV(NR+1) =  DCHITV(NR)*(R(NR)/R(NR+1))**2

   END DO CHITVLOOP

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

   CHITV0=(CHITV(2)*R(1)**2-CHITV(1)*R(2)**2)/(R(1)**2-R(2)**2)	! Quadratic interpolation

!!! PERCHÉ IL 2 IN DCHITV0? !!!

   DCHITV0=2.*(CHITV(1)-CHITV0)/R(1)
   DDCHITV0=2.*(CHITV(1)-CHITV0)/R(1)**2

   ! if(debug)then
   ! 		write(6,*) "CHITV0 = ", CHITV0, "mu0 = ", mu0
   ! 	end if

   CHITV(0)=CHITV(1)
   CHITV(NR+1) =CHITV(NR)*R(NR)/R(NR+1)
   COEFFN=ABS(CHITV0-CHITVTEMP)

   IF((K .NE. 1) .AND. (ABS(CHITV0-CHITVTEMP) .LT. CONV2))THEN
      IF(VERBOSE .AND. (.NOT. MPICODE))THEN
         WRITE(6,*)'EXITED RELAXLOOP => RELAXATION HAS CONVERGED'
      END IF
      EXIT RELAXLOOP
   END IF
   IF(VERBOSE .AND. K==RELIT .AND. (.NOT. MPICODE))THEN
      WRITE(6,*)'EXITED RELAXLOOP WITHOUT CONVERGENCE: CHITV0 CONVERGED TO A PRECISION OF ',&
           &ABS(CHITV0-CHITVTEMP)
   END IF
   CHITVTEMP=CHITV0

END DO RELAXLOOP

DO I=0,NR+1
   NU(I)=NU(I)-NU(NR)+2.*LOG((2*R(NR)-(M+CC2*M**2*1./R(NR+1)-CC2*M**3/6.*1./R(NR+1)**2+&
        &CC2*(1.+3.*CC2)*M**4/12.*1./R(NR+1)**3-M**5/120.*CC2*(3.+11.*CC2)*1./R(NR+1)**4))/&
        &(2*R(NR)+M-CC2*M**2*1./R(NR+1)-CC2*M**3/6.*1./R(NR+1)**2-CC2*(1.+3.*CC2)*M**4/12./&
        &R(NR+1)**3-M**5/120.*CC2*(3.+11.*CC2)*1./R(NR+1)**4))
END DO

! Calculates the Komar mass in both frames

MK=0.

DO I=1,ISUR
   MK=MK+4.*PI*(ETV(I)+RHOTV(I)+3.*PRTV(I))*EXP(NU(I)/2.)*EXP(3.*MU(I)/2.)*R(I)**2*DR(I)
END DO

! Calculates the ADM (M) mass in the Jordan frame

CALL MASSFIND(CC2,MU(NR)+2.*ALPHA0*(CHITV(NR)-CHIINF)+BETA0*(CHITV(NR)-CHIINF)**2,R(NR),MJOR)

MSCAL=0.5*R(NR)*(CHITV(NR)-CHIINF)
MKJOR=MK-2.*ALPHA0*MSCAL

SCALCH1=0.
! SCALCH2=0.

DO I=1,ISUR+10
   SCALCH1=SCALCH1+2.*PI*((ALPHA0+BETA0*(CHITV(I)-CHIINF))*(3.*PRTV(I)-ETV(I)-RHOTV(I))*&
        &EXP(3.*MU(I)/2.))*R(I)**2*DR(I)
END DO

DO I=1,NR
   NUP=(DRM(I)**2*NU(I+1)-DRP(I)**2*NU(I-1)-(DRM(I)**2-DRP(I)**2)*NU(I))/(DRP(I)*DRM(I)*(DRP(I)+DRM(I)))
   ! 	SCALCH2=SCALCH2+EXP(MU(I)/2.)/2.*NUP*DCHITV(I)*R(I)**2*DR(I)
END DO


IF(.NOT. MPICODE)THEN
   WRITE(6,*)''
   WRITE(6,*)'==================== TOV SOLUTION IN STT ===================='
   WRITE(6,*)''
   WRITE(6,*)'Central density (J)      ','    ',RHOCENT/MBARYONFC															! Final central density
   WRITE(6,*)'Central density (E)      ','    ',RHOCENT*EXP(4*ALPHA0*(CHITV0-CHIINF)+2*BETA0*(CHITV0-CHIINF)**2)/MBARYONFC
   WRITE(6,*)'Central scalar field (E) ','    ',CHITV0															! Central scalar field and mu
   WRITE(6,*)'Central mu (E)           ','    ',MU0
   WRITE(6,*)'Central exp(mu) (E)      ','    ',EXP(MU0)
   WRITE(6,*)'ADM mass (E)             ','    ',M																! Masses and scalar charge
   WRITE(6,*)'ADM mass (J)             ','    ',MJOR
   WRITE(6,*)'Komar mass (E)           ','    ',MK
   WRITE(6,*)'Komar mass (J)           ','    ',MKJOR
   WRITE(6,*)'Scalar mass              ','    ',MSCAL
   WRITE(6,*)'Scalar charge (E)        ','    ',CC
   WRITE(6,*)'Iso. Radius (E)          ','    ',R(ISUR)															! Stellar radius in Isotropic coordinates
   WRITE(6,*)'Sch. Radius (E)          ','    ',R(ISUR)*EXP(MU(ISUR)/2.)        								! Stellar radius in Schwarzschild coordinates
   WRITE(6,*)'Isurf                    ','    ',ISUR
   WRITE(6,*)'Midgrid                  ','    ',MIDGRID
   WRITE(6,*)'Grid radius (E)          ','    ',R(NR)			 												! To check that isurf<midgrid
   WRITE(6,*)'COEFF (mass)             ','    ',COEFF															! Convergence quantities
   WRITE(6,*)'COEFFN (scalar field)    ','    ',COEFFN
   WRITE(6,*)''
   WRITE(6,*)'============================================================='
   WRITE(6,*)''
END IF

! Converts PRTV, ETV and RHOTV to the Jordan frame, ready to be used by XNS
PRTV(1:NR)=PRTV(1:NR)/EXP(4.*ALPHA0*(CHITV(1:NR)-CHIINF)+2.*BETA0*(CHITV(1:NR)-CHIINF)**2)
ETV(1:NR)=ETV(1:NR)/EXP(4.*ALPHA0*(CHITV(1:NR)-CHIINF)+2.*BETA0*(CHITV(1:NR)-CHIINF)**2)
RHOTV(1:NR)=RHOTV(1:NR)/EXP(4.*ALPHA0*(CHITV(1:NR)-CHIINF)+2.*BETA0*(CHITV(1:NR)-CHIINF)**2)
! Saves a log file
MU0=MUIN

! if(debug)then
! 	open(14,file='TOVINIMOD_DEBUG.dat')
! end if
OPEN(13,FILE='TOVINIMOD_PROFILES.dat')
DO IZ=1,NR
	WRITE(13,*)IZ,RHOTV(IZ)/MBARYONFC,PRTV(IZ),NU(IZ),MU(IZ),CHITV(IZ)
END DO
CLOSE(13)

! if(debug)then
! 	close(14)
! end if


CALL CPU_TIME(FINISH)
IF(.NOT. MPICODE)PRINT '("Elapsed CPU time = ",F6.3," seconds")',FINISH-START

END SUBROUTINE TOVINIMOD
