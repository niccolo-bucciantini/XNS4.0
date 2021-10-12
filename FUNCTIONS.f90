! Defines the functions used in MPSAMPLING.f90, SEQUENCE.f90 and SINGLERHO.f90
MODULE FUNCTIONS
	!pippo pluto
	USE SYSTEMXNS
	IMPLICIT NONE

	REAL :: RMAXINV,RMIDINV

	CONTAINS

		! SUBROUTINE EOSEIN(P,RHO,EN,CHI)
!
! 		! === Gives the rest-mass density (RHO) and thermal energy (EN) as functions of the pressure (P) (polytropic EoS) in the Einstein frame ===
!
! 			IMPLICIT NONE
! 			REAL :: P,RHO,EN,CHI,ASCAL
!
! 			! Check if TES or GR
! 			IF(GR)THEN
! 				ASCAL=1.
! 			ELSE
! 				ASCAL=EXP(ALPHA0*(CHI-CHIINF)+0.5*BETA0*(CHI-CHIINF)**2)
! 			END IF
!
! 			IF(P .GT. 1.D-15)THEN
! 				RHO=(P/K1)**(1./GAMMA)*ASCAL**((4.*GAMMA-4.)/GAMMA)
! 			ELSE
! 				P=1.D-18
! 				RHO=1.D-18     ! usare questo per domini più estesi?
! 			!	RHO=(P/K1)**(1./GAMMA)*ASCAL**((4.*GAMMA-4.)/GAMMA)
! 			END IF
! 			EN=P/(GAMMA-1.)
!
! 			RETURN
!
! 		END SUBROUTINE EOSEIN

		SUBROUTINE EXPANSION(MU0,PCENT,ECENT,RINI,YINI,CHITV0)

		! === Expansion of the TOV eqs around the origin, to start the integration ===

			IMPLICIT NONE
			REAL,DIMENSION(6) :: YINI
			REAL :: MU0,CHITV0,PCENT,ECENT
    		REAL :: RINI,P0,E0,DDMU0,DDP0,DDNU0,DDCHITV0

			E0=ECENT
			P0=PCENT

			IF(GR)THEN
				DDP0=-2.*PI/3.*EXP(MU0)*(E0**2+3*P0**2+4*E0*P0)
				DDCHITV0=0.
			ELSE
				DDP0=-2.*PI/3.*EXP(MU0)*(E0**2*(1+(ALPHA0+BETA0*(CHITV0-CHIINF))**2)+3*P0**2*(1+3*(ALPHA0+BETA0*(CHITV0-CHIINF))**2)&
					&+2*E0*P0*(2-3*(ALPHA0+BETA0*(CHITV0-CHIINF))**2))
				IF(ANALYTIC)THEN
					DDCHITV0=EXP(MU0)/6.*4*PI*(ALPHA0+BETA0*(CHITV0-CHIINF))*(E0-3*P0)
				END IF
			END IF
			DDMU0=-4./3.*EXP(MU0)*PI*E0
			DDNU0=4./3.*EXP(MU0)*PI*(E0+3*P0)

			YINI(1) = MU0+DDMU0*RINI**2
			YINI(2) = 2*RINI*DDMU0
			YINI(3) = P0+DDP0*RINI**2
			YINI(4) = DDNU0*RINI**2
			YINI(5)	= CHITV0+0.5*DDCHITV0*RINI**2
			YINI(6)	= RINI*DDCHITV0

			RETURN

		END SUBROUTINE EXPANSION


		SUBROUTINE TOVEQS(R,Y,DY)

			! === TOV eqs in isotropic coordinates ===
			! 	  Y(1): mu
			! 	  Y(2): mu'
			! 	  Y(3): p
			!	  Y(4): nu
			!	  Y(5): CHI
			!	  Y(6): CHI'
			! ========================================

			IMPLICIT NONE
			LOGICAL :: GR
			REAL,DIMENSION(6) :: Y,DY
			REAL :: R,RHOX,ENX,ETOT

			CALL EOS(Y(3),RHOX,ENX,Y(5))
! 			write(6,*)'tov',Y(3),RHOX,ENX
			ETOT=RHOX+ENX

			IF(GR)THEN ! GR TOV
			   DY(1) = Y(2)
			   DY(2) = -8.*PI*ETOT*EXP(Y(1))-2./R*Y(2)-0.25*Y(2)**2
			   DY(3) = -(8.*PI*Y(3)*EXP(Y(1))-0.25*Y(2)**2-1./R*Y(2))/(Y(2)/2.+1./R)*(ETOT+Y(3))/2.
			   DY(4) = (8.*PI*Y(3)*EXP(Y(1))-0.25*Y(2)**2-1./R*Y(2))/(Y(2)/2.+1./R)
			   DY(5) = 0.
			   DY(6) = 0.
			ELSE       ! STT TOV
			   DY(1) = Y(2)
			   DY(2) = -8.*PI*ETOT*EXP(Y(1))-2./R*Y(2)-0.25*Y(2)**2-Y(6)**2
			   DY(3) = -(8.*PI*Y(3)*EXP(Y(1))-0.25*Y(2)**2-1./R*Y(2)+Y(6)**2)/(Y(2)/2.+1./R)*&
					&(ETOT+Y(3))/2.+(ALPHA0+BETA0*(Y(5)-CHIINF))*(3*Y(3)-ETOT)*Y(6)
			   DY(4) = (8.*PI*Y(3)*EXP(Y(1))-0.25*Y(2)**2-1./R*Y(2)+Y(6)**2)/(Y(2)/2.+1./R)
			   IF (ANALYTIC)THEN
				  DY(5) = Y(6)
				  DY(6) = -Y(6)*(DY(4)/2.+Y(2)/2.+2./R)-EXP(Y(1))*4*PI*(ALPHA0+BETA0*(Y(5)-CHIINF))*(3*Y(3)-ETOT)
			   ELSE
				  DY(5) = 0.5*(DCHITV(ILOC) + DCHITV(ILOC+1))
				  DY(6) = 0.5*(DDCHITV(ILOC) + DDCHITV(ILOC+1))
			   ENDIF
			END IF

			RETURN

		  END SUBROUTINE TOVEQS


		SUBROUTINE RK4(Y,N,RI,H,YOUT)

		! === 4th order Runge-Kutta integrator ===

			IMPLICIT NONE
			INTEGER N,NMAX
			REAL H,RI,DYDR(N),Y(N),YOUT(N)
			PARAMETER (NMAX=6)
			INTEGER I
			REAL H6,HH,RH,DYM(NMAX),DYT(NMAX),YT(NMAX)

			HH=H*0.5
			H6=H/6.
			RH=RI+HH

			CALL TOVEQS(RI,Y,DYDR)
			DO I=1,N
			   YT(I)=Y(I)+HH*DYDR(I)
			END DO
! 			WRITE(6,*)'RK4',DYDR(3)
			CALL TOVEQS(RH,YT,DYT)
			DO I=1,N
			   YT(I)=Y(I)+HH*DYT(I)
			END DO
			CALL TOVEQS(RH,YT,DYM)
			DO I=1,N
			   YT(I)=Y(I)+H*DYM(I)
			   DYM(I)=DYT(I)+DYM(I)
			END DO
			CALL TOVEQS(RI+H,YT,DYT)
			DO I=1,N
			   YOUT(I)=Y(I)+H6*(DYDR(I)+DYT(I)+2.*DYM(I))
			END DO

			RETURN

		END SUBROUTINE RK4

		SUBROUTINE MASSFIND(CC2,MUREF,RREF,MOUT)

		! === Solves for the M mass ===

			IMPLICIT NONE
			REAL :: CC2,MUREF,RREF,MOUT,RREF1
			REAL :: A1,A2,A3,A4,A5,DIFFM,DDIFFM
			INTEGER,PARAMETER :: MAXITER = 100
			INTEGER :: J

			RREF1 = 1./RREF

			A1 = 1.
			A2 = -CC2*RREF1
			A3 = -CC2/6.*RREF1**2
			A4 = -CC2*(1.+3.*CC2)/12.*RREF1**3
			A5 = -CC2*(3.+11.*CC2)*RREF1**4/120.

			MOUT=1.4
			DO J=1,MAXITER
			   DIFFM = MOUT*(A1 + A2*MOUT + A3*MOUT**2 + A4*MOUT**3. + A5*MOUT**4) - 2.*RREF*(EXP(MUREF/4.)-1)
			   IF(ABS(DIFFM) .LT. CONV)EXIT
			   DDIFFM=(A1 + 2*A2*MOUT + 3*A3*MOUT**2 + 4*A4*MOUT**3. + 5*A5*MOUT**4)
			   MOUT=MOUT-DIFFM/DDIFFM
			END DO

		END SUBROUTINE MASSFIND

		SUBROUTINE CHIDERIVS

			IMPLICIT NONE
			INTEGER :: IX,IZ
			REAL :: A1,A2,A3,B1,B2,B3

			DO IZ=1,NR
			  CHI(0,IZ) =  CHI(1,IZ)
			  CHI(NTH+1,IZ) =  CHI(NTH,IZ)
			END DO

		    !Assume parity at center and smoothness at the outer boundary
		   	DO IX=1,NTH
			  CHI(IX,0) =  CHI(IX,1)
			  CHI(IX,NR+1) =  CHI(IX,NR)*R(NR)/R(NR+1)
			END DO

			!Evaluate Q_r and Q_th
			DO IX=1,NTH
			  DO IZ=1,NR
				A1=-DR(IZ+1)/DR(IZ)/(DR(IZ)+DR(IZ+1))
				A2= (DR(IZ+1)-DR(IZ))/DR(IZ)/DR(IZ+1)
				A3= DR(IZ)/DR(IZ+1)/(DR(IZ)+DR(IZ+1))
				B1=-DTH(IX+1)/DTH(IX)/(DTH(IX)+DTH(IX+1))
				B2= (DTH(IX+1)-DTH(IX))/DTH(IX)/DTH(IX+1)
				B3= DTH(IX)/DTH(IX+1)/(DTH(IX)+DTH(IX+1))

				QSCALR(IX,IZ) = A1*CHI(IX,IZ-1)+A3*CHI(IX,IZ+1)+A2*CHI(IX,IZ)
				QSCALT(IX,IZ) = B1*CHI(IX-1,IZ)+B3*CHI(IX+1,IZ)+B2*CHI(IX,IZ)
				ASCAL(IX,IZ) = EXP(ALPHA0*(CHI(IX,IZ)-CHIINF)+0.5*BETA0*(CHI(IX,IZ)-CHIINF)**2)

			  END DO
			END DO

 END SUBROUTINE CHIDERIVS

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

  SUBROUTINE DGTSV(N,DL,D,DU,B,LDB,INFO)

		!  ============================================================
		!  -- LAPACK routine (version 3.2) --
		!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
		!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
		!     November 2006
		!  ============================================================
		!  Purpose : solves the equation  A*X = B,
		!            where A is an n by n tridiagonal matrix, by
		!            Gaussian elimination with partial pivoting.
		!
		!            Note that the equation  A'*X = B  may be solved by
		!            interchanging the order of the arguments DU and DL.
		!
		!  Arguments
		!  =========
		!
		!  N       (input) INTEGER
		!          The order of the matrix A.  N >= 0.
		!
		!
		!  DL      (input/output) DOUBLE PRECISION array, dimension (N-1)
		!          On entry, DL must contain the (n-1) sub-diagonal elements of
		!          A.
		!
		!          On exit, DL is overwritten by the (n-2) elements of the
		!          second super-diagonal of the upper triangular matrix U from
		!          the LU factorization of A, in DL(1), ..., DL(n-2).
		!
		!  D       (input/output) DOUBLE PRECISION array, dimension (N)
		!          On entry, D must contain the diagonal elements of A.
		!
		!          On exit, D is overwritten by the n diagonal elements of U.
		!
		!  DU      (input/output) DOUBLE PRECISION array, dimension (N-1)
		!          On entry, DU must contain the (n-1) super-diagonal elements
		!          of A.
		!
		!          On exit, DU is overwritten by the (n-1) elements of the first
		!          super-diagonal of U.
		!
		!  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
		!          On entry, the N  matrix of right hand side matrix B.
		!          On exit, if INFO = 0, the N  solution matrix X.
		!
		!  LDB     (input) INTEGER
		!          The leading dimension of the array B.  LDB >= max(1,N).
		!
		!  INFO    (output) INTEGER
		!          = 0: successful exit
		!          < 0: if INFO = -i, the i-th argument had an illegal value
		!          > 0: if INFO = i, U(i,i) is exactly zero, and the solution
		!               has not been computed.  The factorization has not been
		!               completed unless i = N.
		!
		!  =====================================================================

			  INTEGER :: INFO,LDB,N,I,NRHS
			  REAL ::   B(LDB), D(N), DL(N-1), DU(N-1)
			  REAL :: FACT,TEMP

			  INFO = 0
			  IF(N.LE.0)THEN
				INFO = -1
				WRITE(6,*)' N <= 0 IN DGTSV'
				STOP
			  ELSE IF(NRHS.LT.0)THEN
				INFO = -2
				WRITE(6,*)' RHSN < 0 IN DGTSV'
				STOP
			  ELSE IF(LDB .LT. MAX(1,N))THEN
				INFO = -7
				WRITE(6,*)' LBD < MAX(1,N) IN DGTSV'
				STOP
			  END IF

			  DO I = 1, N - 2
				IF( ABS( D( I ) ).GE.ABS( DL( I ) ) ) THEN
				  ! No row interchange required
				  IF( D( I ).NE. 0. ) THEN
					FACT = DL( I ) / D( I )
					D( I+1 ) = D( I+1 ) - FACT*DU( I )
					B( I+1 ) = B( I+1 ) - FACT*B( I )
				  ELSE
					INFO = I
					RETURN
				  END IF
				  DL( I ) = 0.
				ELSE  ! Interchange rows I and I+1
				  FACT = D( I ) / DL( I )
				  D( I ) = DL( I )
				  TEMP = D( I+1 )
				  D( I+1 ) = DU( I ) - FACT*TEMP
				  DL( I ) = DU( I+1 )
				  DU( I+1 ) = -FACT*DL( I )
				  DU( I ) = TEMP
				  TEMP = B( I )
				  B( I ) = B( I+1 )
				  B( I+1 ) = TEMP - FACT*B( I+1 )
				END IF
			  END DO

			  IF( N.GT.1 ) THEN
				I = N - 1
				IF( ABS( D( I ) ).GE.ABS( DL( I ) ) ) THEN
				  IF( D( I ).NE.0. ) THEN
					FACT = DL( I ) / D( I )
					D( I+1 ) = D( I+1 ) - FACT*DU( I )
					B( I+1 ) = B( I+1 ) - FACT*B( I )
				  ELSE
					INFO = I
					RETURN
				  END IF
				ELSE
				  FACT = D( I ) / DL( I )
				  D( I ) = DL( I )
				  TEMP = D( I+1 )
				  D( I+1 ) = DU( I ) - FACT*TEMP
				  DU( I ) = TEMP
				  TEMP = B( I )
				  B( I ) = B( I+1 )
				  B( I+1 ) = TEMP - FACT*B( I+1 )
				END IF
			  END IF
			  IF( D( N ).EQ.0. ) THEN
				INFO = N
				RETURN
			  END IF

			  ! Back solve with the matrix U from the factorization.
			  B( N ) = B( N ) / D( N )
			  IF( N.GT.1 )THEN
				B( N-1 ) = ( B( N-1 )-DU( N-1 )*B( N ) ) / D( N-1 )
			  END IF
			  DO I = N - 2, 1, -1
				B( I ) = ( B( I )-DU( I )*B( I+1 )-DL( I )* &
						 B( I+2 ) ) / D( I )
			  END DO

			  RETURN

			END SUBROUTINE DGTSV



END MODULE FUNCTIONS
