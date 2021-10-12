SUBROUTINE HYDROEQ(RHOVAR,ILOOP)
  ! ciao
  !  ============================================================
  !  Purpose : This program solve for hydrostatic equlibrium in
  !            a given CFC metric -
  !            Work  for a single politropic equation
  !  ============================================================
  !
  !  Parameters for the source
  !
  !  R = radial grid points (+ boundaries)
  !  TH = angular grid points (+ boundaries)
  !  XX = angular cos(th) points (+ boundaries)
  !
  !  Parameters
  !  rhocent,pcent,hcent = central density, pressure and hentalphy
  !  NB: EN (or ENEW) is the thermal energy; the total energy is rho*h=rhonew+enew+pnew
  !  ============================================================

  USE SYSTEMXNS
  IMPLICIT NONE
 
  REAL    :: RHOVAR
  INTEGER :: ILOOP

  REAL :: BCENT,PCENT,ENTCENT
  REAL :: TOL,X0,X1,X2
  REAL :: RHOCENT,ECENT,HCENT,RHOTMP
  REAL :: ALF1,ALF2,LPS1,LPS2,ASCAL1,ASCAL2,ALF0,LPS0,ASCAL0,NU0,NUC																!!!STT!!!
  REAL :: BPOT


  INTEGER :: ICENT

  !Derive the central value of the lapse function via interpolation
  ICENT=1
  IF(NTH .GT. 1)THEN
     ICENT=NTH/2
  END IF

  ! Compute metric and scalar factor at the center
  ALF1=PSI(ICENT,1)
  ALF2=PSI(ICENT,2)
  LPS1=PSL(ICENT,1)
  LPS2=PSL(ICENT,2)
  ASCAL1=ASCAL(ICENT,1)
  ASCAL2=ASCAL(ICENT,2)
  ALF0=(ALF1*R(2)**2-ALF2*R(1)**2)/(R(2)**2-R(1)**2)
  LPS0=(LPS1*R(2)**2-LPS2*R(1)**2)/(R(2)**2-R(1)**2)
  NU0=LOG(LPS0/ALF0)		! nu0=log(alpha) -> diverso da nu=log(alpha^2)
  ASCAL0=(ASCAL1*R(2)**2-ASCAL2*R(1)**2)/(R(2)**2-R(1)**2)
  NUC=(NU(1)*R(2)**2-NU(2)*R(1)**2)/(R(2)**2-R(1)**2)	! Computed on TOVINI

  ! Choose central density
  IF(CONVHELP)THEN
     RHOCENT=RHOVAR*EXP(NUC/2.)/EXP(NU0)
  ELSE
     RHOCENT=RHOVAR
  END IF

  ! Derive central pressure
  TOL=1.E-8
  X0=1.
  X1=1.E-8
  X2=10.*RHOCENT

  ! This is the definition of the central enthalpy
  ! Beware: not associated with RHOCENT or RHOVAR for ease of convergence
  IF(EOSINT)THEN
     ! Different choices depending on stability. First one preferred.
     ! WRITE(6,*)'RHO CENTRAL 1',RHOCENT,RHOVAR
     RHOTMP=2*RHOVAR-RHOCENT
     CALL RHO2EOS(RHOTMP,PCENT,ECENT,ENTCENT)
     ! CALL RHO2EOS(RHOCENT,PCENT,ECENT,ENTCENT)

     BCENT=ENTCENT+NU0+LOG(ASCAL0)
     ! WRITE(6,*)'BCENT T',ENTCENT
  ELSE
     ! Different choices depending on stability. First one preferred.
     PCENT=K1*RHOVAR**GAMMA
     ! 	PCENT=K1*RHOCENT**GAMMA

     ! Derive central energy density
     ECENT=RHOCENT+PCENT/(GAMMA-1.)

     !Computing the bernoulli integral for equlibrium
     HCENT=1.+GAMMA/(GAMMA-1.)*PCENT/RHOCENT
     BCENT=LOG(HCENT)+NU0+LOG(ASCAL0)
     ! WRITE(6,*)'BCENT F',LOG(HCENT)
  ENDIF

  ! Different choices depending on stability. First one preferred.
  ! 	PCENT = RTSAFEG(X0,X1,X2,TOL,RHOVAR,FUNCD_EOS)
  !   PCENT = RTSAFEG(X0,X1,X2,TOL,RHOCENT,FUNCD_EOS)
  !   WRITE(6,*)'PCENT',PCENT,K1*RHOVAR**GAMMA

  IF(tol .EQ. -1.) THEN
     WRITE(6,*)"ERROR: from HYDROEQ, RTSAFEG failure, root must be bracketed",RHOCENT,NUC,NU0
#ifdef MPIV
     ENDID=.TRUE.
     if(VERBOSE)print*,'setting endid to true',endid
     GOTO 100
#else
     STOP
#endif

  END IF

  ! Compute the hydro and magentic variables at all radii
  IF(.NOT.IMAG)     CALL HYDROVAR(BCENT,PCENT)
  IF(IMAG.AND.ITOR) CALL HYDROVAR_TOR(BCENT,PCENT)
  IF(IMAG.AND.IPOL) CALL HYDROVAR_POL(BCENT,PCENT,ILOOP)
  IF(IMAG.AND.ITWT) CALL HYDROVAR_POL(BCENT,PCENT,ILOOP)

100 END SUBROUTINE HYDROEQ

  !*********************************************************
  !*********************************************************

SUBROUTINE HYDROVAR(BCENT,PCENT)
  !---------------------------------------------------------
  ! Computation of hydrodynamic quantities (RHONEW,PNEW,ENEW
  ! V3NEW) for unmagnetized models
  !---------------------------------------------------------

  USE SYSTEMXNS
  USE FUNCTIONS
  USE ROTATION, ONLY: OMEGA3LVALUE
  IMPLICIT NONE
  REAL, INTENT(IN) :: BCENT, PCENT
  REAL :: P, RHO, EN, POLD, PNEWLOC
  REAL :: NULOC, POR, HENT, ENT
  REAL :: BETALOC,OMGN,OMEGALOC
  REAL :: RSTARLOC,V3LOC,V3L,A3L

  LOGICAL :: SURF
!   write(6,*)'before',rhonew(1,1)
  DO IX=1,NTH
     SURF=.TRUE.
     DO IZ=1,NR
        !Compute the local equilibrium hentalpy
        BETALOC=PSS(IX,IZ)		! beta^phi
        RSTARLOC=(R(IZ)*SIN(TH(IX))*PSI(IX,IZ)**3/PSL(IX,IZ))**2 ! R^2/alpha^2
        CALL OMEGA3LVALUE(BETALOC,RSTARLOC,OMEGALOC,A3L)
        OMGN=OMEGALOC		! Omega
        NULOC=LOG(PSL(IX,IZ)/PSI(IX,IZ))        ! nu=log(alpha)
        V3LOC =MIN((OMGN+PSS(IX,IZ))*R(IZ)*SIN(TH(IX))/EXP(NULOC)*PSI(IX,IZ)**2,0.9999999)	! sqrt(v_phi*v^phi)/alpha
        V3L=0.5*LOG(1.-V3LOC**2)		! -log(GammaL)
        HENT=EXP(BCENT-LOG(ASCAL(IX,IZ))-NULOC-V3L+A3L)	! Enthalpy from the Bernoulli int.
        HENT=MAX(HENT,1+5.0E-16)
        IF(EOSINT)THEN
           CALL ENT2EOS(LOG(HENT),RHO)
           CALL RHO2EOS(RHO,P,EN,ENT)
           POR=P/RHO
        ELSE
        								! p
!         PNEW(IX,IZ)=MAX(QFACTOR*P+(1.-QFACTOR)*PSRC(IX,IZ),1.E-15)		! p modified with damping factor Q
!         P=MIN(PNEW(IX,IZ),PCENT)

           POR=(HENT-1.)*(GAMMA-1.)/GAMMA		! p/rho
           P=((POR**GAMMA)/K1)**(1./(GAMMA-1.))
           CALL EOS(P,RHO,EN)
        ENDIF
        !Compute the local equilibrium  v^phi and B^phi
        V3NEW(IX,IZ)=MAX(QFACTOR*(OMGN+PSS(IX,IZ))/(ASCAL(IX,IZ)*EXP(NULOC)) +(1.-QFACTOR)*VPHI(IX,IZ),0.)! v^phi modified with damping factor Q !!!STT!!!
        B3NEW(IX,IZ)=0.
        ! IF(IX==NTH/2 .AND. IZ .LE. 230)THEN
        !    WRITE(6,*)'NEW 1',IZ,HENT,RHO,P
        ! ENDIF
        !Maximum radius for the stability of the convergence scheme (to be changes id needed)
        IF((R(iz) .GT. REQMAX) .OR. (RHO .LE. 1.0E-12))THEN
           POR=0.
!            if(ix==1 .and. iz==1)write(6,*)'ciao',rho
        END IF

        !Compute local equilibrium pressure, density and energy density (reset velocity outside star)
        IF(POR .GT. 0 .AND. SURF)THEN
           IF(EOSINT)THEN
              ! PNEWLOC=MAX(QFACTOR*P+(1.-QFACTOR)*PSRC(IX,IZ),1.E-18)
              PNEWLOC=MAX(EXP(QFACTOR*LOG(P)+(1.-QFACTOR)*LOG(PSRC(IX,IZ))),1.E-18)
              CALL PRS2EOS(PNEWLOC,RHO)
              CALL RHO2EOS(RHO,P,EN,ENT)
              RHONEW(IX,IZ)=MAX(RHO,1.E-15)
              ENEW(IX,IZ)=MAX(EN,1.E-18)
              PNEW(IX,IZ)=MAX(P,1.E-18)
              ! IF(IX==NTH/2 .AND. IZ .LE. 230)THEN
              ! WRITE(6,*)'NEW 2',IZ,RHONEW(IX,IZ),PNEW(IX,IZ)
              ! ENDIF
           ELSE
              P=((POR**GAMMA)/K1)**(1./(GAMMA-1.))				! p
              PNEW(IX,IZ)=MAX(QFACTOR*P+(1.-QFACTOR)*PSRC(IX,IZ),1.E-15)	! p modified with damping factor Q
              P=PNEW(IX,IZ)!MIN(PNEW(IX,IZ),PCENT)
              CALL EOS(P,RHO,EN)
              RHONEW(IX,IZ)=MAX(RHO,1.E-12)
              ENEW(IX,IZ)=MAX(EN,1.E-12)
           ENDIF
           IF(P .LT. 1.E-12)THEN
              V3NEW(IX,IZ)=0.
           END IF
           B3NEW(IX,IZ)=0.
           IF(RHONEW(IX,IZ) .GT. RHOSURF)THEN
              WSURF(IX)=IZ
           END IF
        ELSE
           SURF=.FALSE.
           P=1.E-18
           P=MAX(EXP(QFACTOR*LOG(P)+(1.-QFACTOR)*LOG(PSRC(IX,IZ))),1.E-18)
           IF(EOSINT)THEN
              CALL PRS2EOS(P,RHO)
              CALL RHO2EOS(RHO,P,EN,ENT)
              RHONEW(IX,IZ)=MAX(RHO,1.E-15)
              ENEW(IX,IZ)=MAX(EN,1.E-18)
              PNEW(IX,IZ)=1.E-18
           ELSE
              PNEW(IX,IZ)=1.E-15
              CALL EOS(P,RHO,EN)
              RHONEW(IX,IZ)=MAX(RHO,1.E-15)
              ENEW(IX,IZ)=MAX(EN,1.E-15)
           ENDIF
           V3NEW(IX,IZ)=0.
           B3NEW(IX,IZ)=0.
        END IF
     END DO
  END DO
! write(6,*)'after',rhonew(1,1)
END SUBROUTINE HYDROVAR

!*********************************************************
!*********************************************************

SUBROUTINE HYDROVAR_TOR(BCENT,PCENT)
  !---------------------------------------------------------
  !  Computation of hydrodynamic quantities (RHONEW,PNEW,ENEW
  !  V3NEW,B3NEW) for models with purely toroidal magnetic
  !  field.
  !---------------------------------------------------------

  USE SYSTEMXNS
  USE FUNCTIONS
  USE ROTATION, ONLY: OMEGA3LVALUE
  IMPLICIT NONE
  REAL, INTENT(IN) :: BCENT, PCENT
  REAL :: P, RHO, EN, POLD, B3, FF
  REAL :: DER, DPN, DFF, DHGUESS
  REAL :: DEGUESS,DPGUESS,DRHOGUESS, DLOGHGUESS
  REAL :: RHOGUESS, HGUESS, EGUESS, PGUESS, CHIGUESS, LOGHGUESS
  REAL :: ALPH2GPP, NULOC, POR, HENTMAG, PNEWLOC, ENTLOC
  REAL :: BETALOC,OMGN,OMEGALOC
  REAL :: RSTARLOC,V3LOC,V3L,A3L

  INTEGER :: ITER,ITER_MAX=100
  REAL,PARAMETER :: TOLL=1.E-12

  REAL :: REQMAXOLD

  LOGICAL :: SURF

  REQMAXOLD=REQMAX  ! METTERE FUORI DAL LOOP
  ! WRITE(6,*)'BCENT',BCENT,PCENT
  DO IX=1,NTH
     SURF=.TRUE.
     DO IZ=1,NR
        !Compute the local equilibrium magnetized-hentalpy function
        BETALOC=PSS(IX,IZ)
        RSTARLOC=(R(IZ)*SIN(TH(IX))*PSI(IX,IZ)**3/PSL(IX,IZ))**2
        CALL OMEGA3LVALUE(BETALOC,RSTARLOC,OMEGALOC,A3L)
        OMGN=OMEGALOC

        NULOC=LOG(PSL(IX,IZ)/PSI(IX,IZ))
        V3LOC=MIN((OMGN+PSS(IX,IZ))*R(IZ)*SIN(TH(IX))/EXP(NULOC)*PSI(IX,IZ)**2,0.9999999)
        V3L=0.5*LOG(1.-V3LOC**2)
        HENTMAG=BCENT-LOG(ASCAL(IX,IZ))-NULOC-V3L+A3L

        !Compute local equilibrium pressure, density and energy density (reset velocity outside star)
        PGUESS=MAX(PSRC(IX,IZ),1.E-18) !!!!!!!! PROVARE A LEVARE IL MAX

        ALPH2GPP=(PSL(IX,IZ)*PSI(IX,IZ)*R(IZ)*SIN(TH(IX)))**2	! alpha^2*g_phiphi

        DO ITER=1,ITER_MAX
           IF(EOSINT)THEN
              CALL PRS2EOS(PGUESS,RHOGUESS)
              CALL RHO2EOS(RHOGUESS,PGUESS,EGUESS,LOGHGUESS)
              ! IF((IZ .EQ. 75) .AND. (IX .EQ. 100))WRITE(6,*)'1 EOSINT T',PGUESS,RHOGUESS,LOGHGUESS,BCENT,NULOC,PSI(IX,IZ),PSL(IX,IZ)
           ELSE
              CALL EOS(PGUESS,RHOGUESS,EGUESS)
              LOGHGUESS=LOG(1+GAMMA/(GAMMA-1)*PGUESS/RHOGUESS)
              ! IF((IZ .EQ. 75) .AND. (IX .EQ. 100))WRITE(6,*)'1 EOSINT F',PGUESS,RHOGUESS,LOGHGUESS,BCENT,NULOC,PSI(IX,IZ),PSL(IX,IZ)
           ENDIF
           ALPH2GPP=(PSL(IX,IZ)*PSI(IX,IZ)*R(IZ)*SIN(TH(IX)))**2	! alpha^2*g_phiphi
           FF=LOGHGUESS+BCOEF**2*(MAGIND/(2*MAGIND-1))*(RHOGUESS*(EXP(LOGHGUESS))*ASCAL(IX,IZ)**4*ALPH2GPP)**(2*MAGIND-1)-HENTMAG

           DPGUESS=MAX(PGUESS*(1.+1.E-6),1.0E-18)
           IF(EOSINT)THEN
              CALL PRS2EOS(DPGUESS,DRHOGUESS)
              CALL RHO2EOS(DRHOGUESS,DPGUESS,DEGUESS,DLOGHGUESS)
           ELSE
              CALL EOS(DPGUESS,DRHOGUESS,DEGUESS)
              DLOGHGUESS=LOG(1+GAMMA/(GAMMA-1)*DPGUESS/DRHOGUESS)
           ENDIF

           DFF=DLOGHGUESS+BCOEF**2*(MAGIND/(2*MAGIND-1))*(DRHOGUESS*EXP(DLOGHGUESS)*ASCAL(IX,IZ)**4*ALPH2GPP)**(2*MAGIND-1)-HENTMAG
           DER=(DFF-FF)/(DPGUESS-PGUESS)

           DPN=-FF/DER
           ! IF((IZ .EQ. 1) .AND. (IX .EQ. 50))WRITE(6,*)'PGUESS2',IX,IZ,PGUESS,RHOGUESS,LOGHGUESS
           ! IF((IZ .EQ. 1) .AND. (IX .EQ. 50))WRITE(6,*)'DPGUESS',IX,IZ,DPGUESS,DRHOGUESS,DLOGHGUESS,PGUESS+DPN
           IF (ABS(DPN)<TOLL*PGUESS) EXIT
           ! PGUESS=MAX(PGUESS+DPN,1.E-14)
           PGUESS=PGUESS+DPN
           B3=BCOEF*(EXP(LOGHGUESS)*RHOGUESS*ASCAL(IX,IZ)**4*ALPH2GPP)**MAGIND / ALPH2GPP *(PSL(IX,IZ)/PSI(IX,IZ))/ASCAL(IX,IZ)**3     ! B^phi !!!STT!!!

        END DO

        PGUESS=MAX(PGUESS,1.0E-18)
        ! Compute local equilibrium v^phi
        V3NEW(IX,IZ)=MAX(QFACTOR*(OMGN+PSS(IX,IZ))/(ASCAL(IX,IZ)*EXP(NULOC)) +(1.-QFACTOR)*VPHI(IX,IZ),0.)
        ! Maximum radius for the stability of the convergence scheme (to be changes id needed)
        POR=PGUESS/RHOGUESS
        ! IF((IZ .EQ. 75) .AND. (IX .EQ. 100))WRITE(6,*)'RHOGUESS',RHOGUESS,IX,IZ,PGUESS
        IF((R(iz) .GT. REQMAX) .AND. (REQMAX .LT. RMAXSTR/1.2))THEN
           IF(RHOGUESS .LT. 1.0E-8)THEN
              POR=0.
           ELSE
              WRITE(6,*)'ERROR: star extends beyond REQMAX'
              STOP
!               REQMAX = REQMAX*1.
           END IF
        END IF

        ! Compute local equilibrium pressure, density and energy density (reset velocity outside star)
        !IF(RHOGUESS .GT. 1.0E-8 .AND. SURF)THEN
        IF(POR .GT. 0 .AND. SURF)THEN
           IF(EOSINT)THEN
              P=PGUESS
              !PNEWLOC=MAX(QFACTOR*P+(1.-QFACTOR)*PSRC(IX,IZ),1.E-18)
              PNEWLOC=MAX(EXP(QFACTOR*LOG(P)+(1.-QFACTOR)*LOG(PSRC(IX,IZ))),1.E-18)
              CALL PRS2EOS(PNEWLOC,RHO)
              CALL RHO2EOS(RHO,P,EN,ENTLOC)
              RHONEW(IX,IZ)=MAX(RHO,1.E-15)
              ENEW(IX,IZ)=MAX(EN,1.E-18)
              PNEW(IX,IZ)=MAX(P,1.E-18)
           ELSE
              P=((POR**GAMMA)/K1)**(1./(GAMMA-1.))				! p
              PNEW(IX,IZ)=MAX(QFACTOR*P+(1.-QFACTOR)*PSRC(IX,IZ),1.E-15)	! p modified with damping factor Q
              P=PNEW(IX,IZ)!MIN(PNEW(IX,IZ),PCENT)
              CALL EOS(P,RHO,EN)
              RHONEW(IX,IZ)=MAX(RHO,1.E-12)
              ENEW(IX,IZ)=MAX(EN,1.E-12)
           ENDIF
           ! POLD=PSRC(IX,IZ)
           B3NEW(IX,IZ)=QFACTOR*B3+(1.-QFACTOR)*BPHI(IX,IZ)
           IF(P .LT. 1.E-12)THEN
              V3NEW(IX,IZ)=0.
              B3NEW(IX,IZ)=0.
           END IF
           IF(RHONEW(IX,IZ) .GT. RHOSURF)THEN
              WSURF(IX)=IZ
           END IF
        ELSE ! outside the star
           SURF=.FALSE.
           P=1.E-18
           P=MAX(EXP(QFACTOR*LOG(P)+(1.-QFACTOR)*LOG(PSRC(IX,IZ))),1.E-18)
           IF(EOSINT)THEN
              CALL PRS2EOS(P,RHO)
              CALL RHO2EOS(RHO,P,EN,ENTLOC)
              RHONEW(IX,IZ)=MAX(RHO,1.E-15)
              ENEW(IX,IZ)=MAX(EN,1.E-18)
              PNEW(IX,IZ)=1.E-18
           ELSE
              PNEW(IX,IZ)=1.E-15
              CALL EOS(P,RHO,EN)
              RHONEW(IX,IZ)=MAX(RHO,1.E-15)
              ENEW(IX,IZ)=MAX(EN,1.E-15)
           ENDIF
           V3NEW(IX,IZ)=0.
           B3NEW(IX,IZ)=0.
        END IF
        ! WRITE(6,*)'RHONEW ENEW',RHONEW(1,1),ENEW(1,1)
     END DO
  END DO

  REQMAX=REQMAXOLD

END SUBROUTINE HYDROVAR_TOR

!*********************************************************
!*********************************************************

SUBROUTINE HYDROVAR_POL(BCENT,PCENT,ILOOP)		 !!! QUI
!---------------------------------------------------------
!  Computation of hydrodynamic quantities (RHONEW,PNEW,ENEW
!  V3NEW,B3NEW,BPOLR,BPOLT,E3NEW,EPOLR,EPOLT) for models with
!  poloidal magnetic field.
!---------------------------------------------------------

  USE SYSTEMXNS
  USE FUNCTIONS
  USE ROTATION, ONLY: OMEGA3LVALUE
  IMPLICIT NONE
  REAL, INTENT(IN) :: BCENT,PCENT
  REAL :: P,RHO,EN,POLD,FF,ENT
  REAL :: DER,DPN,DFF,DHGUESS
  REAL :: DEGUESS,DPGUESS,DRHOGUESS
  REAL :: RHOGUESS,HGUESS,EGUESS,PGUESS
  REAL :: ALPH2GPP,NULOC,POR,HENT,BPOT,ENTLOC
  REAL :: BETALOC,OMGN,OMEGALOC,HENTMAG,PNEWLOC
  REAL :: RSTARLOC,V3LOC,V3L,A3L
  REAL :: A1,A2,A3,B1,B2,B3

  INTEGER :: ILOOP
  INTEGER :: ITER,ITER_MAX=1000
  REAL,PARAMETER :: TOLL=1.E-12

  LOGICAL :: SURF

  !In the case of rotating models compute the phi-component
  !and the t-component of the four-vector potential by solving
  !Maxwell Equations, otherwhise compute the phi-component
  !of the vector potential by solving the Grad-Shafranov Equation.
  CALL SOURCEPOT
  IF(OMG .EQ. 0.) THEN
    CALL VECPOTPHI
    ATIM=0.*APHI
  ELSE
    CALL MXWLSOL(ILOOP)
  ENDIF

  !Impose bc on the phi and time components of the vector potential (\tilde{A}_phi, A_t)
  !in order to evaluate numerical derivatives (assume axisymmetry)
  DO IZ=1,NR
    APHI(0,IZ)    =-APHI(1,IZ)
    APHI(NTH+1,IZ)=-APHI(NTH,IZ)
    ATIM(0,IZ)    = ATIM(1,IZ)
    ATIM(NTH+1,IZ)= ATIM(NTH,IZ)

    ATIMIN(0,IZ)     =  ATIMIN(1,IZ)
    ATIMIN(NTH+1,IZ) =  ATIMIN(NTH,IZ)

    ATIMOUT(0,IZ)     =  ATIMOUT(1,IZ)
    ATIMOUT(NTH+1,IZ) =  ATIMOUT(NTH,IZ)
  END DO
  DO IX=1,NTH
    APHI(IX,0)   =-APHI(NTH-IX+1,1)
    APHI(IX,NR+1)= APHI(IX,NR) + APHI(IX,NR)/DR(NR-1)*(APHI(IX,NR)-APHI(IX,NR-1))
    ATIM(IX,0)   = ATIM(NTH-IX+1,1)
    ATIM(IX,NR+1)= ATIM(IX,NR) + DR(NR)/DR(NR-1)*(ATIM(IX,NR)-ATIM(IX,NR-1))

    ATIMIN(IX,0) =  ATIMIN(NTH-IX+1,1)
    ATIMIN(IX,NR+1) =  ATIMIN(IX,NR) + DR(NR)/DR(NR-1)*(ATIMIN(IX,NR)-ATIMIN(IX,NR-1))

    ATIMOUT(IX,0) =  ATIMOUT(NTH-IX+1,1)
    ATIMOUT(IX,NR+1) =  ATIMOUT(IX,NR) + DR(NR)/DR(NR-1)*(ATIMOUT(IX,NR)-ATIMOUT(IX,NR-1))
  END DO

  !Compute the poloidal components of the magnetic field (B^r, B^th) and
  !of the electric field (E_r,E_th)
  DO IX=1,NTH
     DO IZ=1,NR
        A1=-DR(IZ+1)/DR(IZ)/(DR(IZ)+DR(IZ+1))
        A2= (DR(IZ+1)-DR(IZ))/DR(IZ)/DR(IZ+1)
        A3= DR(IZ)/DR(IZ+1)/(DR(IZ)+DR(IZ+1))
        B1=-DTH(IX+1)/DTH(IX)/(DTH(IX)+DTH(IX+1))
        B2= (DTH(IX+1)-DTH(IX))/DTH(IX)/DTH(IX+1)
        B3= DTH(IX)/DTH(IX+1)/(DTH(IX)+DTH(IX+1))

        BPOLT(IX,IZ) = -SIN(TH(IX))*((A1*APHI(IX,IZ-1)+A3*APHI(IX,IZ+1)+A2*APHI(IX,IZ))*R(IZ)+ &
             APHI(IX,IZ))/(ASCAL(IX,IZ)**3*PSI(IX,IZ)**6*R(IZ)**2*SIN(TH(IX)))
        BPOLR(IX,IZ) = ( R(IZ)*SIN(TH(IX))*(B1*APHI(IX-1,IZ)+B3*APHI(IX+1,IZ)+B2*APHI(IX,IZ)) + &
             R(IZ)*COS(TH(IX))*APHI(IX,IZ))/(ASCAL(IX,IZ)**3*PSI(IX,IZ)**6*R(IZ)**2*SIN(TH(IX)))

        IF(OMG .EQ. 0.) THEN
           EPOLR(IX,IZ)=0.
           EPOLT(IX,IZ)=0.
           E3NEW(IX,IZ)=0.
        ELSE
           IF( R(IZ) .LE. ELSURF(IX)) THEN   				! E_i=(omega-Omega)/alpha*grad(PSI)
              EPOLR(IX,IZ) = PSI(IX,IZ)/PSL(IX,IZ)/ASCAL(IX,IZ)*(A1*ATIMIN(IX,IZ-1)+A3*ATIMIN(IX,IZ+1)+A2*ATIMIN(IX,IZ) - &
                   ASCAL(IX,IZ)**3*PSI(IX,IZ)**6*R(IZ)**2*SIN(TH(IX))*OMGMET(IX,IZ)*BPOLT(IX,IZ))
              EPOLT(IX,IZ) = PSI(IX,IZ)/PSL(IX,IZ)/ASCAL(IX,IZ)*(B1*ATIMIN(IX-1,IZ)+B3*ATIMIN(IX+1,IZ)+B2*ATIMIN(IX,IZ) + &
                   ASCAL(IX,IZ)**3*PSI(IX,IZ)**6*R(IZ)**2*SIN(TH(IX))*OMGMET(IX,IZ)*BPOLR(IX,IZ))
              E3NEW(IX,IZ)=0.
           ELSE
              EPOLR(IX,IZ) = PSI(IX,IZ)/PSL(IX,IZ)/ASCAL(IX,IZ)*(A1*ATIMOUT(IX,IZ-1)+A3*ATIMOUT(IX,IZ+1)+A2*ATIMOUT(IX,IZ) - &
                   ASCAL(IX,IZ)**3*PSI(IX,IZ)**6*R(IZ)**2*SIN(TH(IX))*OMGMET(IX,IZ)*BPOLT(IX,IZ))
              EPOLT(IX,IZ) = PSI(IX,IZ)/PSL(IX,IZ)/ASCAL(IX,IZ)*(B1*ATIMOUT(IX-1,IZ)+B3*ATIMOUT(IX+1,IZ)+B2*ATIMOUT(IX,IZ) + &
                   ASCAL(IX,IZ)**3*PSI(IX,IZ)**6*R(IZ)**2*SIN(TH(IX))*OMGMET(IX,IZ)*BPOLR(IX,IZ))
              E3NEW(IX,IZ)=0.
           END IF
        END IF

        !The azimutal component of the magnetic field is provided by the prescription on the
        !current function (B3NEW to be changed if current distribution is changed)
        IF(IPOL) B3NEW=0.
        IF(ITWT) THEN
           IF(APM .EQ. 0.) APM=1.
           B3NEW(IX,IZ) = PSI(IX,IZ)/PSL(IX,IZ)*ATWT/(ZETA+1.)* &	! B^phi= I/(alpha*psi^4*r^2*sin(th)^2)
                MAX( APHI(IX,IZ)*R(IZ)*SIN(TH(IX))-APHIMAX,0.)*&
                (APHI(IX,IZ)*R(IZ)*SIN(TH(IX))-APHIMAX)**ZETA/&
                (R(IZ)**2*PSI(IX,IZ)**4*SIN(TH(IX))**2.)/APM**(ZETA+0.5)/ASCAL(IX,IZ)**3
        END IF

        !Evaluate the current density for non-rotating or rotating stars accordingly to
        !the prescribed current function
        IF(OMG .EQ. 0) THEN
           JPHI(IX,IZ)= ATWT**2./(ZETA+1.)/(R(IZ)**2*PSL(IX,IZ)**2.*PSI(IX,IZ)**2.*SIN(TH(IX))**2.)*&
                MAX(APHI(IX,IZ)*R(IZ)*SIN(TH(IX))-APHIMAX,0.)*&
                (APHI(IX,IZ)*R(IZ)*SIN(TH(IX))-APHIMAX)**(2.*ZETA)/APM**(2.*ZETA+1.)/ASCAL(IX,IZ)**4+&
                (ENEW(IX,IZ)+RHONEW(IX,IZ))*KBTT
           JTH(IX,IZ) = PSI(IX,IZ)/PSL(IX,IZ)*ATWT*BPOLT(IX,IZ)*MAX(APHI(IX,IZ)*R(IZ)*SIN(TH(IX))-APHIMAX,0.)*&
                (APHI(IX,IZ)*R(IZ)*SIN(TH(IX))-APHIMAX)**(ZETA-1.)/APM**(ZETA+0.5)/ASCAL(IX,IZ)
           JRR(IX,IZ) = PSI(IX,IZ)/PSL(IX,IZ)*ATWT*BPOLR(IX,IZ)*MAX(APHI(IX,IZ)*R(IZ)*SIN(TH(IX))-APHIMAX,0.)*&
                (APHI(IX,IZ)*R(IZ)*SIN(TH(IX))-APHIMAX)**(ZETA-1.)/APM**0.5/ASCAL(IX,IZ)
        ELSE
           JTH(IX,IZ)=0.
           JRR(IX,IZ)=0.
           JPHI(IX,IZ)=JPHIMXL(IX,IZ)
        ENDIF
     END DO
  END DO

  DO IZ=1,NR-1
     IF(B3NEW(NTH/2,IZ) .GT. 0. ) THEN	! RFCUT is the magnetosphere radius
        RFCUT = R(IZ)
        EXIT
     ENDIF
  END DO

  DO IX=1,NTH
     SURF=.TRUE.
     DO IZ=1,NR
        ! Compute the local equilibrium hentalpy function
        BETALOC=PSS(IX,IZ)
        RSTARLOC=(R(IZ)*SIN(TH(IX))*PSI(IX,IZ)**3/PSL(IX,IZ))**2
        CALL OMEGA3LVALUE(BETALOC,RSTARLOC,OMEGALOC,A3L)
        OMGN=OMEGALOC

        NULOC=LOG(PSL(IX,IZ)/PSI(IX,IZ))
        V3LOC =MIN((OMGN+PSS(IX,IZ))*R(IZ)*SIN(TH(IX))/EXP(NULOC)*PSI(IX,IZ)**2/ASCAL(IX,IZ),0.9999999)								!!!STT!!!
        V3L=0.5*LOG(1.-V3LOC**2)

        !Compute the magnetization function BPOT=-M(A_phi)
        IF(IPOL .AND. (OMG .EQ. 0) ) BPOT = -KBPOL*(R(IZ)*SIN(TH(IX))*APHI(IX,IZ)+0.5*CSI*(R(IZ)*SIN(TH(IX))*APHI(IX,IZ))**2.)
        IF(IPOL .AND. (OMG .GT. 0) ) BPOT = -KBPOL*(R(IZ)*SIN(TH(IX))*APHI(IX,IZ))
        IF(ITWT) BPOT = -KBTT*(R(IZ)*SIN(TH(IX))*APHI(IX,IZ))

        HENTMAG=BCENT-NULOC-LOG(ASCAL(IX,IZ))-V3L+A3L-BPOT
        HENTMAG=MAX(HENTMAG,1.0E-12)

        IF(EOSINT)THEN
           CALL ENT2EOS(HENTMAG,RHO)
           CALL RHO2EOS(RHO,P,EN,ENT)
           POR=P/RHO
           RHOGUESS=RHO
        ELSE
           HENT=MAX(EXP(HENTMAG),1.)
           POR=(HENT-1.)*(GAMMA-1.)/GAMMA
           P=((POR**GAMMA)/K1)**(1./(GAMMA-1.))
           RHOGUESS=(P/K1)**(1./GAMMA)
        ENDIF

        IF(RHOGUESS .LT. 1.0E-12)THEN
           POR=0.
        END IF
        ! Compute the local equilibrium  v^phi
        V3NEW(IX,IZ)=MAX(QFACTOR*(OMGN+PSS(IX,IZ))/EXP(NULOC)/ASCAL(IX,IZ)**2 +(1.-QFACTOR)*VPHI(IX,IZ),0.)
        ! Maximum radius for the stability of the convergence scheme (to be changes if needed)
        IF(R(IZ) .GT. REQMAX)THEN
           POR=0.
        END IF

        ! Compute local equilibrium pressure, density and energy density (reset velocity outside star)
        !IF(RHOGUESS .GT. 1.0E-8 .AND. SURF)THEN
        IF(POR .GT. 0 .AND. SURF)THEN
           !PNEWLOC=MAX(QFACTOR*P+(1.-QFACTOR)*PSRC(IX,IZ),1.E-18)
           PNEWLOC=MAX(EXP(QFACTOR*LOG(P)+(1.-QFACTOR)*LOG(PSRC(IX,IZ))),1.E-18)
           IF(EOSINT)THEN
              CALL PRS2EOS(PNEWLOC,RHO)
              CALL RHO2EOS(RHO,P,EN,ENT)
              RHONEW(IX,IZ)=MAX(RHO,1.E-15)
              ENEW(IX,IZ)=MAX(EN,1.E-18)
              PNEW(IX,IZ)=MAX(P,1.E-18)
           ELSE
              P=((POR**GAMMA)/K1)**(1./(GAMMA-1.))				! p
              PNEW(IX,IZ)=MAX(QFACTOR*P+(1.-QFACTOR)*PSRC(IX,IZ),1.E-15)! p modified with damping factor Q
              P = PNEW(IX,IZ)
              CALL EOS(P,RHO,EN)
              RHONEW(IX,IZ)=MAX(RHO,1.E-12)
              ENEW(IX,IZ)=MAX(EN,1.E-12)
           ENDIF
           IF(P .LT. 1.E-12)THEN
              V3NEW(IX,IZ)=0.
           END IF
           IF(RHONEW(IX,IZ) .GT. RHOSURF)THEN
              WSURF(IX)=IZ
           END IF
        ELSE ! outside the star
           SURF=.FALSE.
           P=1.E-18
           P=MAX(EXP(QFACTOR*LOG(P)+(1.-QFACTOR)*LOG(PSRC(IX,IZ))),1.E-18)
           IF(EOSINT)THEN
              CALL PRS2EOS(P,RHO)
              CALL RHO2EOS(RHO,P,EN,ENTLOC)
              RHONEW(IX,IZ)=MAX(RHO,1.E-15)
              ENEW(IX,IZ)=MAX(EN,1.E-18)
              PNEW(IX,IZ)=1.E-18
           ELSE
              PNEW(IX,IZ)=1.E-15
              CALL EOS(P,RHO,EN)
              RHONEW(IX,IZ)=MAX(RHO,1.E-15)
              ENEW(IX,IZ)=MAX(EN,1.E-15)
           ENDIF
           V3NEW(IX,IZ)=0.
        END IF
     END DO
  END DO

END SUBROUTINE HYDROVAR_POL

! ********************************************************
! ********************************************************

SUBROUTINE COVTERM(IXX,IZZ)

  USE SYSTEMXNS
  IMPLICIT NONE

  INTEGER :: IXX,IZZ

  ALPHA = PSL(IXX,IZZ)/PSI(IXX,IZZ)
  GP = PSI(IXX,IZZ)**6*R(IZZ)**2*SIN(TH(IXX))
  GM = 1./GP
  GCOVR = PSI(IXX,IZZ)**4
  GCOVT = PSI(IXX,IZZ)**4*R(IZZ)**2
  GCOVP = PSI(IXX,IZZ)**4*R(IZZ)**2*SIN(TH(IXX))**2

END SUBROUTINE COVTERM

! ********************************************************
! ********************************************************

SUBROUTINE CONS_TO_PRIM(U,SY2,D,BY2,DENS,PRESS,VELOC,BMG,SSC)

  USE SYSTEMXNS
  IMPLICIT NONE

  REAL :: U,SY2,D,DENS,PRESS,VELOC,BY2,B2,BMG,SSC
  REAL :: ET,S2,W,V2,RH,PG,FW,DV2,DRH,DPG,DFW,DW
  REAL :: CHII,DCHII,GLF,GAM,SB2

  REAL,PARAMETER :: EPS1=1.E-16,TOL=1.E-12
  INTEGER :: ITER,ITER_MAX=1000

  GAM=(GAMMA-1.)/GAMMA

  ET=U
  S2=SY2
  B2=BY2
  SB2=S2*B2

  ! Non magnetic case!
  IF(.NOT. IMAG)THEN
     ! Initial guess
     W=(SQRT(2.*(ET)+MAX(4.*(ET)**2-3.*(S2),0.)))/3.    ! COS'È W???
     DO ITER=1,ITER_MAX
        V2=(S2)/(W**2)
        RH=D*SQRT(MAX(EPS1,1.-V2))
        PG=GAM*(W*(1.-V2)-RH)
        FW=W-PG-ET

        DV2=-2.*(S2)/(W)**3
        DRH=-0.5*D**2/RH*DV2
        DPG=GAM*(1.-V2-W*DV2-DRH)
        DFW=1.-DPG

        DW=-FW/DFW

        IF (ABS(DW)<TOL*W) EXIT
        W=W+DW
    END DO

    PRESS=MAX(PG,EPS1)
    DENS=D*SQRT(MAX(EPS1,1.-V2))
    VELOC=SQRT(V2)
    BMG=SQRT(B2)
 END IF

 ! Magnetic case
 IF(IMAG.AND.ITOR)THEN
    ! Initial guess
    W=(SQRT(2.*MAX((ET-B2),1.e-3)+MAX(4.*(ET-B2)**2-3.*(S2-B2*(2.*ET-B2)),0.)))/3.
    DO ITER=1,ITER_MAX
       V2=(W*W*S2+(2.*W+B2)*SB2)/(W*W*(W+B2)**2)
       GLF=1./SQRT(1.-V2)
       RH=D/GLF
       CHII=(W-D*GLF)/GLF**2
       PG=GAM*CHII
       FW=W-ET-PG+.5*B2

       DV2=-2.*(S2+SB2*(3.*W*(W+B2)+B2**2)/W**3)/(W+B2)**3
       DCHII=1./GLF**2-W*DV2+0.5*GLF*D*DV2
       DPG=GAM*DCHII
       DFW=1.-DPG

       DW=-FW/DFW

       IF (ABS(DW)<TOL*W) EXIT
       W=W+DW
    END DO
    PRESS=MAX(PG,EPS1)
    DENS=D*SQRT(MAX(EPS1,1.-V2))
    VELOC=SQRT(V2)
    BMG=SQRT(B2)
 END IF

 IF(IMAG.AND.IPOL)THEN
    ! Initial guess
    W=(SQRT(2.*(ET)+MAX(4.*(ET)**2-3.*(S2),0.)))/3.

    DO ITER=1,ITER_MAX
       V2=(S2)/(W**2)
       RH=D*SQRT(MAX(EPS1,1.-V2))
       PG=GAM*(W*(1.-V2)-RH)
       FW=W-PG-ET

       DV2=-2.*(S2)/(W)**3
       DRH=-0.5*D**2/RH*DV2
       DPG=GAM*(1.-V2-W*DV2-DRH)
       DFW=1.-DPG

       DW=-FW/DFW

       IF (ABS(DW)<TOL*W) EXIT
       W=W+DW
    END DO

    PRESS=MAX(PG,EPS1)
    DENS=D*SQRT(MAX(EPS1,1.-V2))
    VELOC=SQRT(V2)
    BMG=SQRT(B2)
 END IF

 ! In the limit of purely toroidal magnetic field and velocity
 SSC=(DENS+GAMMA/(GAMMA-1)*PRESS)*V2/(1.-V2)-B2+3*(PRESS+0.5*B2)
 !   WRITE(6,*)'CTP UNMAG',V2,DENS,PRESS,W,ET,U,S2,D

END SUBROUTINE CONS_TO_PRIM

! ********************************************************
! ********************************************************

SUBROUTINE CONS_TO_PRIM_POL(U,SY,D,BY,BR,BT,EY,ER,ET,DENS,PRESS,VELOC,SSC,IXT,IZT)

  USE SYSTEMXNS
  USE FUNCTIONS
  IMPLICIT NONE

  REAL :: U,SY,D,BY,BR,BT,EY,ER,ET,DENS,PRESS,VELOC,SSC
  REAL :: S2,B2,SB2,GAM,ENT,EL2,PG,VYCOV
  REAL :: RHO,V2,W,C0,C2,C3,DW,DC3,DC2,DLOGW,WB,VB2,F,DF,DV2

  REAL,PARAMETER :: EPS1=1.E-16,TOL1=1.E-6,TOL2=1.E-10
  INTEGER :: ITERW,ITERV,ITERV_MAX=1000,IXt,IZt

  GAM=(GAMMA-1.)/GAMMA

  S2=SY*SY/GCOVP/ASCAL(IX,IZ)**2	         				 !Ok (poloidal altrimenti manca il termine -B(v scalar B))				!!!STT!!!
  B2=ASCAL(IX,IZ)**2*(BY*BY*GCOVP+BR*BR*GCOVR+BT*BT*GCOVT)	 !OK
  SB2=(SY*BY)**2.      						 				 ! ok, ma nel mio caso e' nullo


  !IF(ABS(SY).LE.TOL1)THEN
  IF(IZT .GT. WSURF(IXT)) THEN
    VELOC = 0.
    V2=0.
    VYCOV=0.
    EL2=(EY*EY/GCOVP+ER*ER/GCOVR+ET*ET/GCOVT)/ASCAL(IX,IZ)**2																		!!!STT!!!
    PRESS = EPS1 !(U-0.5*B2-D)*(GAMMA-1.)
    CALL EOS(PRESS,DENS,ENT)
!     WRITE(6,*)'INSIDE CTP2',D,RHO,V2
    !DENS  = D
  END IF

  !IF(ABS(SY).GT.TOL1)THEN
  IF(IZT .LE. WSURF(IXT)) THEN
    V2 = 0.5
    DO ITERV=1,ITERV_MAX
      RHO = D*SQRT(1.-V2)
      C3 = 1-GAM*(1.-V2)
      C2 = GAM*RHO+0.5*B2*(1+V2)-U
      C0 = -0.5*SB2
      W = MAX(MAX(-C2/C3,(-C0/C3)**(1./3.)),1.e-8)
      DO ITERW=1,100
         DW = -((C3*W+C2)*W**2+C0)/((3*C3*W+2*C2)*W)
         IF(ABS(DW/W).LT.TOL2) EXIT
         W = W+DW
      END DO
      DC3 = GAM
      DC2 = 0.5*(B2-GAM*RHO/(1.-V2))
      DLOGW = -(DC3*W+DC2)/(3*C3*W+2*C2)

      WB = W+B2
      VB2 = SB2/W**2

      F = WB**2*V2-(2*W+B2)*VB2-S2
      DF = WB*(WB+2*DLOGW*(W*V2+VB2))

      DV2 = -F/DF
      IF( ABS(DV2/V2) .LT. TOL2)EXIT
      V2 = MAX(MIN(V2+DV2,0.99),0.)
    END DO

    PG = GAM*(W*(1.-V2)-RHO)

    PRESS=MAX(PG,EPS1)
    DENS=D*SQRT(MAX(EPS1,1.-V2))
    VELOC=SQRT(V2)

    VYCOV = VELOC*SQRT(GCOVP)*ASCAL(IX,IZ)																							!!!STT!!!
    ET = -GM*(VYCOV*BR*GCOVR)/ASCAL(IX,IZ)
    ER =  GM*(VYCOV*BT*GCOVT)/ASCAL(IX,IZ)
    EL2 = ASCAL(IX,IZ)**2*(ET**2*GCOVT +ER**2*GCOVR)
!   	WRITE(6,*)'INSIDE CTP1',D,RHO,V2
  END IF

  SY=W*VYCOV+ASCAL(IX,IZ)**3*GP*(ER*BT-ET*BR) 																						!!!STT!!!
  D=RHO/SQRT(1-V2)

  ! Compute the source term for lapse equation
  SSC=(DENS+GAMMA/(GAMMA-1)*PRESS)*V2/(1.-V2)-B2-EL2+3*(PRESS+0.5*B2+0.5*EL2)

END SUBROUTINE CONS_TO_PRIM_POL

!*********************************************************
!*********************************************************

SUBROUTINE OMEGAVALUE(BETALOC,RSTARLOC,OMEGALOC)

  USE SYSTEMXNS
  IMPLICIT NONE

  REAL :: OMEGALOC,BETALOC,RSTARLOC
  REAL :: C1,C2,C3,C4,FO,DFO,DOMEGALOC

  REAL,PARAMETER :: EPS1=1.E-16,TOL=1.E-12
  INTEGER :: ITER,ITER_MAX=1000

  OMEGALOC=OMG

  IF(DIFFERENTIAL)THEN
    OMEGALOC=OMG/2.
    DO ITER=1,ITER_MAX
      C1=A2VALUE*RSTARLOC
      C2=A2VALUE*RSTARLOC*(2*BETALOC-OMG)
      C3=A2VALUE*(BETALOC**2*RSTARLOC-1.-2*BETALOC*RSTARLOC*OMG)-RSTARLOC
      C4=A2VALUE*(OMG-BETALOC**2*OMG*RSTARLOC)-BETALOC*RSTARLOC

      FO=C1*OMEGALOC**3 + C2*OMEGALOC**2 + C3*OMEGALOC + C4
      DFO=3*C1*OMEGALOC**2 + 2*C2*OMEGALOC + C3

      DOMEGALOC = -FO/DFO

      IF (ABS(DOMEGALOC)<TOL*OMEGALOC) EXIT

      OMEGALOC=OMEGALOC+DOMEGALOC

    END DO
  END IF

END SUBROUTINE OMEGAVALUE

! ********************************************************
! ********************************************************

SUBROUTINE QUANTITIES(RHOCENT,RHOVAR,MM,M0,ILOOP)

  USE SYSTEMXNS
  USE ROTATION, ONLY: OMEGAVALUE
  IMPLICIT NONE

  REAL :: RHOCENT,RHOVAR,MM,M0,MP,MPJ,SCHARGE,SQUAD
  REAL :: SINIX,DTT,DET
  REAL :: BETALOC,RSTARLOC,OMEGALOC,OMGN
  REAL :: RPOL,REQ,RATIO,RCIRC,REQKM,OMGEQ
  REAL :: GLF,V2,JJ,TT,VYCOV
  REAL :: IZZM,IXXM,IZZS,IXXS,IZZ,IXX,DEF,TQUAD,TMONO
  REAL :: AA1,AA2,AA3,BB1,BB2,BB3

  REAL :: DELTAPHI
  REAL :: BPOLO,BCENT,B2
  REAL :: B2MAX,B2POL,B2TOR,BFCG
  REAL :: BMAXG,BMAXGPOL,BMAXGTOR
  REAL :: B2MAXPOL,B2MAXTOR
  REAL :: EPOLO,ECHRG,ECHRGIN,E2,EEBB
  REAL :: JMPHI,JMPHIIN,JMPHIEXT,DIPM
  REAL :: HH,HHTOR,HHTOREXT,HHTORIN,HHPOL,FLUX,HHE

  REAL,DIMENSION(1:NTH,1:NR) :: B2VEC,DCHIDR,DCHIDT

  INTEGER :: ILOOP

  B2MAX=0.; B2MAXPOL=0.; B2MAXTOR=0.
  BCENT=0.;    BPOLO=0.; BPOLO=0.
  HHTOR=0.;  HHTORIN=0.; HHPOL=0.
  EPOLO=0.;    ECHRG=0.; ECHRGIN=0.
  JMPHI=0.;  JMPHIIN=0.; FLUX=0.
  MM=0.;  M0=0.;  MP=0.; MPJ=0.
  TT=0.;  JJ=0.;  HH=0.; HHE = 0.
  SCHARGE=0.; TQUAD=0.; TMONO=0.;
  SQUAD=0; IXXM=0.; IZZM=0.;
  IXXS=0.; IZZS=0.;IXX=0.; IZZ=0.;

  ! initialise the boundary conditions for chi
  CHI(1:NTH,0)=CHI(1:NTH,1)
  CHI(1:NTH,NR+1)=CHI(1:NTH,NR)
  CHI(0,0:NR+1)=CHI(1,0:NR+1)
  CHI(NTH+1,0:NR+1)=CHI(NTH,0:NR+1)

	WRITE(6,*)'WSURF',WSURF(1),WSURF(NTH/2)

  IF(.NOT.(IPOL.OR.ITWT))THEN
  DO IX=1,NTH
    SINIX=SIN(TH(IX))
    DTT=DTH(IX)
    IF(NTH .EQ. 1)THEN
      SINIX=2.
      DTT=1.
    END IF

    DO IZ=1,WSURF(IX)
       V2=V3NEW(IX,IZ)*V3NEW(IX,IZ)*R(IZ)**2*ASCAL(IX,IZ)**2*PSI(IX,IZ)**4*SINIX**2
       GLF=1./SQRT(1.-V2)
       B2=B3NEW(IX,IZ)*B3NEW(IX,IZ)*R(IZ)**2*ASCAL(IX,IZ)**2*PSI(IX,IZ)**4*SINIX**2
       B2MAX=MAX(B2,B2MAX)
       !IF(IX .eq. 1 .and. IZ .eq. 1) BCENT = sqrt(B2)
       !IF(IX .eq. 1 .AND. IZ .eq. WSURF(1) ) BPOLO = sqrt(B2)

       BETALOC=PSS(IX,IZ)
       RSTARLOC=(R(IZ)*SINIX*PSI(IX,IZ)**3/PSL(IX,IZ))**2
       CALL OMEGAVALUE(BETALOC,RSTARLOC,OMEGALOC)
       OMGN=OMEGALOC
       IF(IX .EQ. NTH/2+1 .AND. PNEW(IX,IZ) .GT. 1.E-15) THEN
          OMGEQ=OMGN
       END IF
       !!! QUI
       DET=PSI(IX,IZ)**6*R(IZ)**2*SINIX*DR(IZ)*DTT
       MM = MM+2*PI*ASCAL(IX,IZ)**4*(2.*PNEW(IX,IZ)+(ENEW(IX,IZ)+PNEW(IX,IZ)+RHONEW(IX,IZ))*GLF**2* &	! gravitational (Komar) mass M in the E frame
            &(1+V2-2*PSI(IX,IZ)**3/PSL(IX,IZ)*R(IZ)*SINIX*SQRT(V2)*PSS(IX,IZ))+B2)*DET*PSL(IX,IZ)/PSI(IX,IZ)
       ! Baryonic mass M0
       M0 = M0+2*PI*ASCAL(IX,IZ)**3*RHONEW(IX,IZ)/MBARYONFC*GLF*DET
       MP = MP+2*PI*ASCAL(IX,IZ)**4*(ENEW(IX,IZ)+RHONEW(IX,IZ))*GLF*DET	! proper mass Mp in the E frame
       MPJ = MPJ+2*PI*ASCAL(IX,IZ)**3*(ENEW(IX,IZ)+RHONEW(IX,IZ))*GLF*DET ! proper mass Mp in the J frame
!       IF(IX .LE. 2 .and. IZ .LE. 85)WRITE(6,*)'MM',IX,IZ,MM
       TT=TT+PI*ASCAL(IX,IZ)**4*(ENEW(IX,IZ)+PNEW(IX,IZ)+RHONEW(IX,IZ))*GLF**2*SQRT(V2)*OMGN*R(IZ)*SINIX*DET*PSI(IX,IZ)**2 ! rotational kinetic energy in the E frame
       IF(OMG .NE. 0)  JJ=JJ+2*PI*ASCAL(IX,IZ)**4*(ENEW(IX,IZ)+PNEW(IX,IZ)+&	! Komar angular momentum in the E frame
            &RHONEW(IX,IZ))*GLF**2*SQRT(V2)*R(IZ)*SINIX*DET*PSI(IX,IZ)**2

       HH = HH + PI*ASCAL(IX,IZ)**3*B2*DET  ! Magnetic energy H in the J frame
       HHE = HHE + PI*ASCAL(IX,IZ)**4*B2*DET  ! Magnetic energy H in the E frame
       FLUX = FLUX+ASCAL(IX,IZ)**2*SQRT(B2)*PSI(IX,IZ)**4.*R(IZ)*DR(IZ)*DTT      	! flux of B^phi through meridional plane

       IXXM=IXXM+PI*ASCAL(IX,IZ)**4*(ENEW(IX,IZ)+RHONEW(IX,IZ))*R(IZ)**4.*SINIX*(2.-SINIX**2.)*DR(IZ)*DTT 	! fluid's inertia moment I_xx in the E frame
       IZZM=IZZM+2.*PI*ASCAL(IX,IZ)**4*(ENEW(IX,IZ)+RHONEW(IX,IZ))*R(IZ)**4.*SINIX**3.*DR(IZ)*DTT          	! fluid's inertia moment I_zz in the E frame

       SCHARGE=SCHARGE+2.*PI*PSL(IX,IZ)/PSI(IX,IZ)*(ALPHA0+BETA0*(CHI(IX,IZ))-CHIINF)*ASCAL(IX,IZ)**4*&
            &(3.*PNEW(IX,IZ)-(ENEW(IX,IZ)+RHONEW(IX,IZ)))*DET 																					! Scalar charge Q_chi in the E frame

       TQUAD=TQUAD+2.*PI*(ALPHA0+BETA0*(CHI(IX,IZ))-CHIINF)*ASCAL(IX,IZ)**4*&			! quadrupole of the trace in the E frame (Newtonian)
            &(3.*PNEW(IX,IZ)-(ENEW(IX,IZ)+RHONEW(IX,IZ)))*(2.-3.*SINIX**2)*R(IZ)**4*SINIX*DR(IZ)*DTT
       TMONO=TMONO+2.*PI*(ALPHA0+BETA0*(CHI(IX,IZ))-CHIINF)*ASCAL(IX,IZ)**4*&		! normalization factor for TQUAD
            &(3.*PNEW(IX,IZ)-(ENEW(IX,IZ)+RHONEW(IX,IZ)))*R(IZ)**2*SINIX*DR(IZ)*DTT
       ! Quadrupole of the trace in the E frame (Newtonian)
       SQUAD=SQUAD+2.*PI*PSL(IX,IZ)/PSI(IX,IZ)*(ALPHA0+BETA0*(CHI(IX,IZ))-CHIINF)*ASCAL(IX,IZ)**4*&
            &(3.*PNEW(IX,IZ)-(ENEW(IX,IZ)+RHONEW(IX,IZ)))*(2.-3.*SINIX**2)*PSI(IX,IZ)**6*R(IZ)**4*SINIX*DR(IZ)*DTT

    END DO
    DO IZ=1,NR
      AA1=-DR(IZ+1)/DR(IZ)/(DR(IZ)+DR(IZ+1))
      AA2= (DR(IZ+1)-DR(IZ))/DR(IZ)/DR(IZ+1)
      AA3= DR(IZ)/DR(IZ+1)/(DR(IZ)+DR(IZ+1))
      BB1=-DTH(IX+1)/DTH(IX)/(DTH(IX)+DTH(IX+1))
      BB2= (DTH(IX+1)-DTH(IX))/DTH(IX)/DTH(IX+1)
      BB3= DTH(IX)/DTH(IX+1)/(DTH(IX)+DTH(IX+1))
      DCHIDR(IX,IZ) = AA1*CHI(IX,IZ-1)+AA3*CHI(IX,IZ+1)+AA2*CHI(IX,IZ)
      DCHIDT(IX,IZ) = BB1*CHI(IX-1,IZ)+BB3*CHI(IX+1,IZ)+BB2*CHI(IX,IZ)
	  IXXS=IXXS-1./8.*(DCHIDR(IX,IZ)**2+DCHIDT(IX,IZ)**2/R(IZ)**2)*R(IZ)**4.*SINIX*(2.-SINIX**2.)*DR(IZ)*DTT 		! scalar field's inertia moment I_xx in the E frame
 	  IZZS=IZZS-1./4.*(DCHIDR(IX,IZ)**2+DCHIDT(IX,IZ)**2/R(IZ)**2)*R(IZ)**4.*SINIX**3*DR(IZ)*DTT         			! scalar field's inertia moment I_zz in the E frame
    END DO
  END DO
  IXX=IXXM+IXXS																															! total inertia moments in the E frame
  IZZ=IZZM+IZZS

  HHTOR=HH
  B2MAXTOR=B2MAx

  END IF


  IF(IPOL.OR.ITWT)THEN

    IF(IPOL .OR. (CUT .LT. 1.)) CUTSURF=WSURF
    DO IX=1,NTH
      SINIX=SIN(TH(IX))
      DTT=DTH(IX)
      IF(NTH .EQ. 1)THEN
        SINIX=2.
        DTT=1.
      END IF

      DO IZ=1,NR
        CALL COVTERM(IX,IZ)
        V2=V3NEW(IX,IZ)**2*ASCAL(IX,IZ)**2*GCOVP
        GLF=1./SQRT(1.-V2)

        VYCOV = V3NEW(IX,IZ)*ASCAL(IX,IZ)**2*GCOVP

        E2 = (EPHI(IX,IZ)**2./GCOVP + EPOLR(IX,IZ)**2./GCOVR + EPOLT(IX,IZ)**2./GCOVT)/ASCAL(IX,IZ)**2
        B2 = ASCAL(IX,IZ)**2*(BPHI(IX,IZ)**2*GCOVP + BPOLR(IX,IZ)**2*GCOVR + BPOLT(IX,IZ)**2*GCOVT)
        B2VEC(IX,IZ)= ASCAL(IX,IZ)**2*(BPHI(IX,IZ)**2*GCOVP + BPOLR(IX,IZ)**2*GCOVR + BPOLT(IX,IZ)**2*GCOVT)
        B2POL = ASCAL(IX,IZ)**2*(BPOLR(IX,IZ)**2*GCOVR + BPOLT(IX,IZ)**2*GCOVT)
        B2TOR = BPHI(IX,IZ)**2*ASCAL(IX,IZ)**2*GCOVP
        B2MAX = MAX(B2,B2MAX)
        B2MAXPOL = MAX(B2POL,B2MAXPOL)
        B2MAXTOR = MAX(B2TOR,B2MAXTOR)
        IF(IX .EQ. 1 .AND. IZ .EQ. 1) BCENT = SQRT(B2)
        IF(IX .EQ. 1 .AND. IZ .EQ. WSURF(1)+1 ) BPOLO = SQRT(B2)
        IF(IX .EQ. 1 .AND. IZ .EQ. WSURF(1)+1 ) EPOLO = SQRT(E2)

        EEBB = ASCAL(IX,IZ)*GP*(EPOLR(IX,IZ)/GCOVR*BPOLT(IX,IZ)-EPOLT(IX,IZ)/GCOVT*BPOLR(IX,IZ))*PSS(IX,IZ)

        BETALOC=PSS(IX,IZ)
        RSTARLOC=(R(IZ)*SINIX*PSI(IX,IZ)**3/PSL(IX,IZ))**2
        CALL OMEGAVALUE(BETALOC,RSTARLOC,OMEGALOC)
        OMGN=OMEGALOC
        IF(IX .EQ. NTH/2+1 .AND. PNEW(IX,IZ) .GT. 1.E-15) THEN
           OMGEQ=OMGN
        END IF



        DET=PSI(IX,IZ)**6*R(IZ)**2*SINIX*DR(IZ)*DTT !!! QUI

        IF(IZ .LE. CUTSURF(IX)) THEN
          MM=MM+2.*PI*ASCAL(IX,IZ)**4*(2.*PNEW(IX,IZ)+(ENEW(IX,IZ)+PNEW(IX,IZ)+RHONEW(IX,IZ))*GLF**2* &
               (1+V2-2*PSI(IX,IZ)**3/PSL(IX,IZ)*R(IZ)*SINIX*SQRT(V2)*PSS(IX,IZ))+B2+E2 &
               -2*PSI(IX,IZ)/PSL(IX,IZ)*EEBB/ASCAL(IX,IZ)**3)*DET*PSL(IX,IZ)/PSI(IX,IZ)
          M0=M0+2*PI*ASCAL(IX,IZ)**3*RHONEW(IX,IZ)/MBARYONFC*GLF*DET
          MP=MP+2*PI*ASCAL(IX,IZ)**4*(ENEW(IX,IZ)+RHONEW(IX,IZ))*GLF*DET
          MPJ=MPJ+2*PI*ASCAL(IX,IZ)**3*(ENEW(IX,IZ)+RHONEW(IX,IZ))*GLF*DET

          TT=TT+  PI*ASCAL(IX,IZ)**4*(ENEW(IX,IZ)+PNEW(IX,IZ)+RHONEW(IX,IZ))*GLF**2*SQRT(V2)*OMGN*R(IZ)*SINIX*DET*PSI(IX,IZ)**2
          IF(OMG .NE. 0) JJ=JJ+2*PI*ASCAL(IX,IZ)**4*((ENEW(IX,IZ)+PNEW(IX,IZ)+RHONEW(IX,IZ))*&
          	   &GLF**2*SQRT(V2)*PSI(IX,IZ)**2*R(IZ)*SINIX + EEBB/PSS(IX,IZ)/ASCAL(IX,IZ)**3)*DET

          FLUX=FLUX+ASCAL(IX,IZ)**2*SQRT(B2)*PSI(IX,IZ)**4.*R(IZ)*DR(IZ)*DTT

          IXXM=IXXM+PI*ASCAL(IX,IZ)**4*(ENEW(IX,IZ)+RHONEW(IX,IZ))*R(IZ)**4.*SINIX*(2.-SINIX**2.)*DR(IZ)*DTT
      	  IZZM=IZZM+2.*PI*ASCAL(IX,IZ)**4*(ENEW(IX,IZ)+RHONEW(IX,IZ))*R(IZ)**4.*SINIX**3.*DR(IZ)*DTT

         AA1=-DR(IZ+1)/DR(IZ)/(DR(IZ)+DR(IZ+1))
         AA2= (DR(IZ+1)-DR(IZ))/DR(IZ)/DR(IZ+1)
         AA3= DR(IZ)/DR(IZ+1)/(DR(IZ)+DR(IZ+1))
         BB1=-DTH(IX+1)/DTH(IX)/(DTH(IX)+DTH(IX+1))
         BB2= (DTH(IX+1)-DTH(IX))/DTH(IX)/DTH(IX+1)
         BB3= DTH(IX)/DTH(IX+1)/(DTH(IX)+DTH(IX+1))
         DCHIDR(IX,IZ) = AA1*CHI(IX,IZ-1)+AA3*CHI(IX,IZ+1)+AA2*CHI(IX,IZ)
         DCHIDT(IX,IZ) = BB1*CHI(IX-1,IZ)+BB3*CHI(IX+1,IZ)+BB2*CHI(IX,IZ)
         IXXS=IXXS-1./8.*(DCHIDR(IX,IZ)**2+DCHIDT(IX,IZ)**2/R(IZ)**2)*R(IZ)**4.*SINIX*(2.-SINIX**2.)*DR(IZ)*DTT
         IZZS=IZZS-1./4.*(DCHIDR(IX,IZ)**2+DCHIDT(IX,IZ)**2/R(IZ)**2)*R(IZ)**4.*SINIX**3*DR(IZ)*DTT

         SCHARGE=SCHARGE+2.*PI*PSL(IX,IZ)/PSI(IX,IZ)*(ALPHA0+BETA0*(CHI(IX,IZ))-CHIINF)*ASCAL(IX,IZ)**4*&
              &(3.*PNEW(IX,IZ)-(ENEW(IX,IZ)+RHONEW(IX,IZ)))*DET
         ! quadrupole of the trace in the E frame (Newtonian)
         TQUAD=TQUAD+2.*PI*(ALPHA0+BETA0*(CHI(IX,IZ))-CHIINF)*ASCAL(IX,IZ)**4*&
              &(3.*PNEW(IX,IZ)-(ENEW(IX,IZ)+RHONEW(IX,IZ)))*(2.-3.*SINIX**2)*R(IZ)**4*SINIX*DR(IZ)*DTT
         TMONO=TMONO+2.*PI*(ALPHA0+BETA0*(CHI(IX,IZ))-CHIINF)*ASCAL(IX,IZ)**4*&
              &(3.*PNEW(IX,IZ)-(ENEW(IX,IZ)+RHONEW(IX,IZ)))*R(IZ)**2*SINIX*DR(IZ)*DTT
         ! quadrupole of the trace in the E frame (Newtonian * alpha * psi**6)
         SQUAD=SQUAD+2.*PI*PSL(IX,IZ)/PSI(IX,IZ)*(ALPHA0+BETA0*(CHI(IX,IZ))-CHIINF)*ASCAL(IX,IZ)**4*&
              &(3.*PNEW(IX,IZ)-(ENEW(IX,IZ)+RHONEW(IX,IZ)))*(2.-3.*SINIX**2)*PSI(IX,IZ)**6*R(IZ)**4*SINIX*DR(IZ)*DTT

      ELSE
         MM=MM+2.*PI*ASCAL(IX,IZ)**4*(B2+E2-2*PSI(IX,IZ)/PSL(IX,IZ)*EEBB/ASCAL(IX,IZ)**3)*DET*PSL(IX,IZ)/PSI(IX,IZ)

         AA1=-DR(IZ+1)/DR(IZ)/(DR(IZ)+DR(IZ+1))
         AA2= (DR(IZ+1)-DR(IZ))/DR(IZ)/DR(IZ+1)
         AA3= DR(IZ)/DR(IZ+1)/(DR(IZ)+DR(IZ+1))
         BB1=-DTH(IX+1)/DTH(IX)/(DTH(IX)+DTH(IX+1))
         BB2= (DTH(IX+1)-DTH(IX))/DTH(IX)/DTH(IX+1)
         BB3= DTH(IX)/DTH(IX+1)/(DTH(IX)+DTH(IX+1))
         DCHIDR(IX,IZ) = AA1*CHI(IX,IZ-1)+AA3*CHI(IX,IZ+1)+AA2*CHI(IX,IZ)
         DCHIDT(IX,IZ) = BB1*CHI(IX-1,IZ)+BB3*CHI(IX+1,IZ)+BB2*CHI(IX,IZ)
         IXXS=IXXS-1./8.*(DCHIDR(IX,IZ)**2+DCHIDT(IX,IZ)**2/R(IZ)**2)*R(IZ)**4.*SINIX*(2.-SINIX**2.)*DR(IZ)*DTT
         IZZS=IZZS-1./4.*(DCHIDR(IX,IZ)**2+DCHIDT(IX,IZ)**2/R(IZ)**2)*R(IZ)**4.*SINIX**3*DR(IZ)*DTT

      END IF
      IXX = IXXM + IXXS
      IZZ = IZZM + IZZS

      IF(IZ .EQ. NR-100 ) ECHRG=ECHRG+ 2.*PI*EPOLR(IX,IZ)*PSI(IX,IZ)**2.*R(IZ)**2.*SINIX*DTT
      JMPHI = JMPHI +2*PI*SQRT(JPHI(IX,IZ)**2*GCOVP)*DET*ASCAL(IX,IZ)**4
      !!!!!! Energy in the J-frame
      HHTOR = HHTOR + PI*ASCAL(IX,IZ)**5*BPHI(IX,IZ)**2*GCOVP*DET + PI*ASCAL(IX,IZ)*EPHI(IX,IZ)**2/GCOVP*DET
      HHPOL = HHPOL + PI*ASCAL(IX,IZ)**5*(BPOLR(IX,IZ)**2*GCOVR + BPOLT(IX,IZ)**2*GCOVT)*DET + &
           PI*ASCAL(IX,IZ)*(EPOLR(IX,IZ)**2/GCOVR + EPOLT(IX,IZ)**2/GCOVT)*DET
      HH = HH + PI*ASCAL(IX,IZ)**5*(BPHI(IX,IZ)**2*GCOVP + BPOLR(IX,IZ)**2*GCOVR + BPOLT(IX,IZ)**2*GCOVT)*DET + &
           PI*ASCAL(IX,IZ)*(EPHI(IX,IZ)**2/GCOVP + EPOLR(IX,IZ)**2/GCOVR + EPOLT(IX,IZ)**2/GCOVT)*DET
      HHE = HHE + PI*ASCAL(IX,IZ)**6*(BPHI(IX,IZ)**2*GCOVP + BPOLR(IX,IZ)**2*GCOVR + BPOLT(IX,IZ)**2*GCOVT)*DET + &
           PI*ASCAL(IX,IZ)**2*(EPHI(IX,IZ)**2/GCOVP + EPOLR(IX,IZ)**2/GCOVR + EPOLT(IX,IZ)**2/GCOVT)*DET

      IF(IZ .LE. WSURF(IX)) THEN
         HHTORIN= HHTORIN + PI*BPHI(IX,IZ)**2*GCOVP*DET*ASCAL(IX,IZ)**5 + PI*EPHI(IX,IZ)**2/GCOVP*DET*ASCAL(IX,IZ)
         JMPHIIN=JMPHIIN +2*PI*SQRT(JPHI(IX,IZ)**2*GCOVP)*DET*ASCAL(IX,IZ)**4
        ENDIF

        IF(IZ .GT. WSURF(IX)) THEN
          HHTOREXT= HHTOREXT + PI*BPHI(IX,IZ)**2*GCOVP*DET*ASCAL(IX,IZ)**5 + PI*EPHI(IX,IZ)**2/GCOVP*DET*ASCAL(IX,IZ)
          JMPHIEXT=JMPHIEXT +2*PI*SQRT((JPHI(IX,IZ)-(ENEW(IX,IZ)+RHONEW(IX,IZ))*KBTT)**2*GCOVP)*DET*ASCAL(IX,IZ)**4
        ENDIF
      END DO
    END DO
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    DIPM = MAXVAL(Aphi(:,NR-10))*4.*R(NR-10)**3./(4.*R(NR-10)+MM)
  END IF

  RPOL=R(WSURF(1))
  REQ=R(WSURF(NTH/2+1))
  REQCHECK=R(WSURF(NTH/2+1)+1)

  IF(IDAT .AND. (.NOT. MPICODE))THEN
    WRITE(6,*)''
	WRITE(6,*)'===================== STELLAR QUANTITIES ====================='
	WRITE(6,*)''
	WRITE(6,*)'Central density (J)','  ',RHOCENT/MBARYONFC  ! Computed on grid = RNEWT (not to be confused with RHOCENT in HYDROEQ)
    WRITE(6,*)'Gravit. mass (E)   ','  ',MM
    WRITE(6,*)'Rest    mass       ','  ',M0
    WRITE(6,*)'Proper  mass (E)   ','  ',MP
    WRITE(6,*)'Proper  mass (J)   ','  ',MPJ
    WRITE(6,*)'Scalar charge (E)  ','  ',SCHARGE!,-scharge/2./mm
    WRITE(6,*)'Rotat.  energy (E) ','  ',TT
    WRITE(6,*)'Angul.  moment. (E)','  ',JJ
    WRITE(6,*)'Magnet. energy (J) ','  ',HH
    WRITE(6,*)'Magnet. energy (E) ','  ',HHE
    WRITE(6,*)'MEnergy  ratio (J/E)','  ',HH/(MP-MM+TT+HHE)
    WRITE(6,*)'MEnergy  ratio (E) ','  ',HHE/(MP-MM+TT+HHE)
    WRITE(6,*)'KEnergy  ratio (E) ','  ',TT/(MP-MM+TT+HHE)

    WRITE(6,*)
    WRITE(6,*)'Equatorial radius (E)    ','  ',REQ,'=',REQ/0.6772030,'KM'
    WRITE(6,*)'Radius ratio (E)         ','  ',RPOL/REQ
    IF(IPOL .AND. OMG .NE. 0) WRITE(6,*) 'N SuperEll.             ','  ', NSPE
    WRITE(6,*)'Circ. radius (E)         ','  ',REQ*PSI(NTH/2+1,WSURF(NTH/2+1))**2.,'='&
    	&,REQ*PSI(NTH/2+1,WSURF(NTH/2+1))**2./0.6772030,'KM'
    WRITE(6,*)'Circ. radius (J)         ','  ',ASCAL(NTH/2+1,WSURF(NTH/2+1))*REQ*PSI(NTH/2+1,WSURF(NTH/2+1))**2.,'='&
    	&,ASCAL(NTH/2+1,WSURF(NTH/2+1))*REQ*PSI(NTH/2+1,WSURF(NTH/2+1))**2./0.6772030,'KM'
    WRITE(6,*)'OMGcen                   ','  ',OMG
    WRITE(6,*)'OMGeq                    ','  ',OMGEQ
    IF(.NOT. IPOL) WRITE(6,*)'Mag. Flux                ','  ',FLUX,'=',SQRT(4.*pi)*FLUX/1.946269d-22,'Wb'
    IF(IPOL.OR.ITWT) WRITE(6,*)'Mag. Dipole (J)          ','  ', DipM
!     WRITE(6,*)'Bmax (J)                 ','  ',SQRT(B2MAX),'=',SQRT(B2MAX*4.*pi)/4.244674d-20,'G'
    WRITE(6,'(A10,A22,ES23.16,A6,ES23.16,A4)')' Bmax (J) ','    ',SQRT(B2MAX),'  =   ',SQRT(B2MAX*4.*pi)/4.244674d-20,'  G'
    WRITE(6,*)'B@pole (J)               ','  ',BPOLO,'=',BPOLO*SQRT(4.*pi)/4.244674d-20,'G'
    WRITE(6,*)'E@pole (J)               ','  ',EPOLO
    !WRITE(6,*)'E.Charge = ', ECHRG, 'VS',ECHRGIN
    WRITE(6,*)'Fluid quadrupole (E)     ','  ',IZZM-IXXM
    WRITE(6,*)'Fluid intertia moment (E)','  ',IZZM
    WRITE(6,*)'Sc. field quadrupole (E) ','  ',IZZS-IXXS
	WRITE(6,*)'Sc. field inertia moment (E) ','  ',IZZS
    WRITE(6,*)'Trace quadrupole (E)     ','  ',TQUAD/(TMONO*REQ**2)
    WRITE(6,*)'Trace monopole (E)       ','  ',TMONO
    WRITE(6,*)'Total quadrupole (E)     ','  ',IZZM+IZZS-IXXM-IXXS
    WRITE(6,*)'Fluid def.rate (E)       ','  ',(IZZM-IXXM)/IZZM
    WRITE(6,*)'Sc. field def. rate (E)  ','  ',(IZZS-IXXS)/IZZS
    WRITE(6,*)'Total def. rate (E)      ','  ',(IZZ-IXX)/IZZM
    WRITE(6,*)' '
    WRITE(6,*)'Poloidal En (J)          ','  ',HHPOL
    WRITE(6,*)'Toroidal En (J)          ','  ',HHTOR
    IF(IMAG) WRITE(6,*) 'Tor./Total (J)          ','  ',HHTOR/HH
    WRITE(6,*)''
	WRITE(6,*)'=============================================================='
	WRITE(6,*)''
  END IF


  ! Write Log file with Run Parameters
  IF(IDAT)THEN
  	CALL DATE_AND_TIME(DATE,TIME,ZONE,VALUES)
    OPEN(12,FILE=trim(adjustl(subdirpath))//'LogFile.dat')

	WRITE(12,*)'Log file for the XNS.f90 program'
	WRITE(12,*)'Ran on date-time-zone: ',DATE,'-',TIME,'-',ZONE
	WRITE(12,*)' '
	IF(STRETCH)THEN
    	WRITE(12,*)'STRETCHED GRID'
    ELSE
    	WRITE(12,*)'UNSTRETCHED GRID'
    END IF
    IF(OMG .GT. 1.E-8)THEN
    	WRITE(12,*)'ROTATING STAR'
		IF(DIFFERENTIAL)THEN
			WRITE(12,*)'DIFFERENTIAL ROTATION'
		ELSE
			WRITE(12,*)'UNIFORM ROTATION'
		END IF
	ELSE
		WRITE(12,*)'NON-ROTATING STAR'
	END IF
	IF(IMAG)THEN
    	WRITE(12,*)'MAGNETIZED STAR'
    	IF(ITOR)WRITE(12,*)'PURELY TOROIDAL MAGNETIC FIELD'
    	IF(IPOL)WRITE(12,*)'PURELY POLOIDAL MAGNETIC FIELD'
    	IF(ITWT)WRITE(12,*)'TWISTED TORUS MAGNETIC FIELD'
    ELSE
    	WRITE(12,*)'UNMAGNETIZED STAR'
    END IF
    IF(GR)THEN
    	WRITE(12,*)'GENERAL RELATIVITY'
    ELSE
    	WRITE(12,*)'SCALAR TENSOR THEORY'
	ENDIF

	WRITE(12,*)'=========================================INPUT========================================================'
	WRITE(12,*)' '
    WRITE(12,*)'>>>>>> GRID SETTINGS : '
    WRITE(12,*)' '
    WRITE(12,*)' NR  = ', NR
    WRITE(12,*)' NTH = ', NTH
    WRITE(12,*) 'NRREG = ', NRREG
    WRITE(12,*) 'STR = ', STRR
    WRITE(12,*)' MLS = ', MLS
    WRITE(12,*)' NGQ = ', NGQ
    WRITE(12,*)' MAXLOOP = ', MAXLOOP
    WRITE(12,*)' MAXSTEP = ', MAXSTEP
    WRITE(12,*)' RELIT = ', RELIT
    WRITE(12,*)' RMIN = ', RMIN
    WRITE(12,*) 'RREG = ', RREG
    WRITE(12,*) 'RMAXSTR = ', RMAXSTR
    WRITE(12,*)' RMAX = ', RMAX
    WRITE(12,*)' RINI = ', RINI
    WRITE(12,*)' '
    WRITE(12,*)'>>>>>> MODEL SETTINGS : '
    WRITE(12,*)' '
    WRITE(12,*)' RHOINI = ', RHOINI
    WRITE(12,*)' RHOINISEQ = ', RHOINISEQ
    WRITE(12,*)' MBARYONFC = ', MBARYONFC
    WRITE(12,*)' GAMMA  = ', GAMMA
    WRITE(12,*)' K1     = ', K1
    WRITE(12,*)' EOSINT     = ', EOSINT
    IF(EOSINT .AND. (RHOCENT .GT. 10**RHOTABMAX))THEN
    	WRITE(12,*)'WARNING: CHOSEN CENTRAL DENSITY'
    	WRITE(12,*)'IS BEYOND THE MAXIMUM TABULATED VALUE OF ',10**RHOTABMAX
    	WRITE(12,*)'POWER-LAW EXTRAPOLATION IS BEING USED'
    ENDIF
    WRITE(12,*)' FILEEOS     = ', FILEEOS
    WRITE(12,*)' ALPHA0    = ', ALPHA0
    WRITE(12,*)' BETA0    = ', BETA0
    WRITE(12,*)' CHIINF    = ', CHIINF
    WRITE(12,*)' OMG = ', OMG
    WRITE(12,*)' PROTDIFF = ', PROTDIFF
    WRITE(12,*)' A2VALUE = ', A2VALUE
    WRITE(12,*)' OMGMAX = ', OMGMAX
    WRITE(12,*)' RMVALUE = ', RMVALUE
    WRITE(12,*)' BCOEF  = ', BCOEF
    WRITE(12,*)' MAGIND = ', MAGIND
	WRITE(12,*)' MLSL = ', MLSL
    WRITE(12,*)' KBPOL = ', KBPOL
    WRITE(12,*)' NPOL = ', NPOL
    WRITE(12,*)' CSI    = ', CSI
	WRITE(12,*)' QNULL = ', QNULL
	WRITE(12,*)' CTP   = ', CTP
    WRITE(12,*)' KBTT = ', KBTT
    WRITE(12,*)' ATWT = ', ATWT
    WRITE(12,*)' ZETA = ', ZETA
    WRITE(12,*)' CUT  = ', CUT
    WRITE(12,*)' '
    WRITE(12,*)'>>>>>> LOGICAL PARAMETERS : '
    WRITE(12,*)' '
    WRITE(12,*)' GR = ', GR
    WRITE(12,*)' STRETCH = ', STRETCH
    WRITE(12,*)' ANALYTIC = ', ANALYTIC
    WRITE(12,*)' MPICODE = ', MPICODE
    WRITE(12,*)' CONVHELP = ', CONVHELP
    WRITE(12,*)' OMGSPACE = ', OMGSPACE
    WRITE(12,*)' JCONSTLAW = ', JCONSTLAW
    WRITE(12,*)' JCMODLAW = ', JCMODLAW
    WRITE(12,*)' URYULAW3 = ', URYULAW3
    WRITE(12,*)' URYULAW4 = ', URYULAW4
    WRITE(12,*)' IMAG = ', IMAG
    WRITE(12,*)' ITOR = ', ITOR
    WRITE(12,*)' IPOL = ', IPOL
    WRITE(12,*)' ITWT = ', ITWT
    WRITE(12,*)' CTP = ', CTP
    WRITE(12,*)' EOSINT = ', EOSINT
    WRITE(12,*)' '
    WRITE(12,*)'>>>>>> RELAXATION/SHOOTING SETTINGS : '
    WRITE(12,*)' '
    WRITE(12,*)' REQMAX =  ',REQMAX
    WRITE(12,*)' QFACTOR = ',QFACTOR
    WRITE(12,*)' QFACTORCHI = ',QFACTORCHI
    WRITE(12,*)' QFACTORMETRIC = ',QFACTORMETRIC
    WRITE(12,*)' QRELAX = ',QRELAX
    WRITE(12,*)' QAPHI =',QAPHI
    WRITE(12,*)' EPS =     ',EPS
    WRITE(12,*)' TOLCONV = ',TOLCONV
	WRITE(12,*)' TOLCHI = ',TOLCHI
	WRITE(12,*)' CONV = ',CONV
    WRITE(12,*)' CONV2 = ',CONV2
    WRITE(12,*)' MMID = ',MMID
    WRITE(12,*)' M = ',M
    WRITE(12,*)' MUIN = ',MUIN
    WRITE(12,*)' '
	WRITE(12,*)'=========================================OUTPUT======================================================='
    WRITE(12,*)' '
    WRITE(12,*)' KOMAR MASS (E) = ', MM
    WRITE(12,*)' REST    MASS = ', M0
    WRITE(12,*)' PROPER  MASS (E) = ', MP
    WRITE(12,*)' PROPER  MASS (J) = ', MPJ
    WRITE(12,*)' SCALAR CHARGE (E) = ', SCHARGE
    WRITE(12,*)' ROTAT.  ENERGY (E)  = ', TT
    WRITE(12,*)' ANGUL.  MOMENT. (E) = ', JJ
    WRITE(12,*)' MAGNET. ENERGY (J)  = ', HH
    WRITE(12,*)' MAGNET. E-ENERGY (E)  = ', HHE
    WRITE(12,*)' MENERGY  RATIO (J/E) = ', HH/(MP-MM+TT+HHE)
    WRITE(12,*)' MENERGY  RATIO (E) = ', HHE/(MP-MM+TT+HHE)
    WRITE(12,*)' KENERGY  RATIO (E) = ', TT/(MP-MM+TT+HHE)
    WRITE(12,*)' '
    WRITE(12,*)'CENTRAL DENSITY (J) = ', RHOCENT/MBARYONFC
    WRITE(12,*)'EQUATORIAL RADIUS (E) = ', REQ,REQ/0.6772030,'KM'
    WRITE(12,*)'RADIUS RATIO (E) = ', RPOL/REQ
    IF(IPOL .AND. OMG .NE. 0) WRITE(12,*) 'NSPE (SUPERELL.) =', NSPE
    WRITE(12,*)'CIRC RADIUS (E) = ', REQ*PSI(NTH/2+1,WSURF(NTH/2+1))**2.,REQ*PSI(NTH/2+1,WSURF(NTH/2+1))**2./0.6772030,'KM'
    WRITE(12,*)'CIRC RADIUS (J) = ', ASCAL(NTH/2+1,WSURF(NTH/2+1))*REQ*PSI(NTH/2+1,WSURF(NTH/2+1))**2.,&
    								&ASCAL(NTH/2+1,WSURF(NTH/2+1))*REQ*PSI(NTH/2+1,WSURF(NTH/2+1))**2./0.6772030,'KM'
    WRITE(12,*)'CENTRAL OMG = ', OMG
    WRITE(12,*)'OMGEQ = ', OMGEQ
    WRITE(12,*)'MAG. FLUX = ', FLUX, SQRT(4.*pi)*FLUX/1.946269d-22,'Wb'
    WRITE(12,*)'MAG. DIPOLE DIPM (J) = ', DipM
    WRITE(12,*)'B@POLE BPOLO (J) = ', BPOLO, BPOLO*SQRT(4.*pi)/4.244674d-20,'G'
    WRITE(12,*)'BMAX (J)  = ', SQRT(B2MAX),SQRT(B2MAX*4.*pi)/4.244674d-20,'G'
    WRITE(12,'(A12,A22,ES23.16,A6,ES23.16,A4)')' Bmax_s (J) ','    ',SQRT(B2MAX),'  =   ',SQRT(B2MAX*4.*pi)/4.244674d-20,'  G'
    WRITE(12,*)'BMAX(POL) (J)  = ', SQRT(B2MAXPOL),SQRT(B2MAXPOL*4.*pi)/4.244674d-20,'G'
    WRITE(12,*)'BMAX(TOR) (J)  = ', SQRT(B2MAXTOR),SQRT(B2MAXTOR*4.*pi)/4.244674d-20,'G'
    WRITE(12,*)'BCENT (J)  = ', BCENT, BCENT*SQRT(4.*pi)/4.244674d-20,'G'
    WRITE(12,*)'E@POLE EPOLO (J)   = ', EPOLO
    WRITE(12,*)'ELEC. CHARGE ECHRG (J) = ', ECHRG, 'VS',ECHRGIN
    WRITE(12,*)'FLUID QUADRUPOLE (E) =      ','  ',IZZM-IXXM
    WRITE(12,*)'FLUID INTERTIA MOMENT (E) = ','  ',IZZM
    WRITE(12,*)'SC. FIELD QUADRUPOLE (E) = ','  ',IZZS-IXXS
	WRITE(12,*)'SC. FIELD INERTIA MOMENT (E) = ','  ',IZZS
	IF(.NOT. GR)WRITE(12,*)'TRACE MONOPOLE (E) =    ','  ', TMONO
    IF(.NOT. GR)WRITE(12,*)'TRACE QUADRUPOLE (E) =    ','  ',TQUAD/(TMONO*REQ**2)
    IF(.NOT. GR)WRITE(12,*)'GRTRACE QUADRUPOLE (E) =    ','  ',SQUAD
    WRITE(12,*)'TOTAL QUADRUPOLE (E) =     ','  ',IZZM+IZZS-IXXM-IXXS
    WRITE(12,*)'FLUID DEF. RATE (E)	=    ','        ',(IZZM-IXXM)/IZZM
    IF(.NOT. GR)WRITE(12,*)'SCL. FIELD DEF. RATE (E) = ','        ',(IZZS-IXXS)/IZZS
    WRITE(12,*)'TOTAL DEF. RATE (E) =    ','        ',(IZZ-IXX)/IZZM
    WRITE(12,*)'HHPOL (J) = ', HHpol
    WRITE(12,*)'HHTOR (J) = ',HHtor
    IF(IMAG) WRITE(12,*)'HHTOR/HH (J) = ',HHtor/HH
    WRITE(12,*)' '
	WRITE(12,*)'======================================================================================================'
    WRITE(12,*)' '
    WRITE(12,*)'NUMBER OF LOOPS DONE = ',ILOOP,'OUT OF ',MAXLOOP

    CLOSE(12)
  END IF

END SUBROUTINE QUANTITIES



! ********************************************************
! ********************************************************

SUBROUTINE SOURCEPOT !!! QUI
!---------------------------------------------------------
! Compute the metric terms and metric derivatives used in
! the source terms of the Maxwell equations/Grad-Shafranov
! equations.
! METERM(:,:)  = Log(Alpha/Psi^2)
! OMTERM(:,:)  = Psi^4*r^2*Sin(th)^2/Alpha^2
! RHOTERM(:,:) = Rho*Enthalpy*Psi^8*A(chi)^4
! VMETTERM(:,:)= Alpha*v^2/(Omega-omega)
! DRMETTERM, DTMETTERM -> derivatives of METERM
! DROMGM, DTOMGM - derivative of the frame-dragging pot.
!--------------------------------------------------------

  USE SYSTEMXNS
  IMPLICIT NONE

  REAL :: A1,A2,A3,B1,B2,B3

  DO IX=1,NTH
    DO IZ=1,NR
      OMGMET(IX,IZ)= -PSS(IX,IZ)
      VLOC(IX,IZ)=MIN((OMG-OMGMET(IX,IZ))*R(IZ)*SIN(TH(IX))/PSL(IX,IZ)*PSI(IX,IZ)**3.,0.9999999)
      GAMLOC(IX,IZ)=1./SQRT(1.-VLOC(IX,IZ)**2.)
      METERM(IX,IZ) = LOG(PSL(IX,IZ)/PSI(IX,IZ)**3)
      OMTERM(IX,IZ) = PSI(IX,IZ)**6.*SIN(TH(IX))**2.*R(IZ)**2/PSL(IX,IZ)**2.
      RHOTERM(IX,IZ) = (RHOSRC(IX,IZ)+GAMMA*PSRC(IX,IZ)/(GAMMA-1.))*PSI(IX,IZ)**8*ASCAL(IX,IZ)**4
! 	  IF((IZ .EQ. 10) .AND. (IX .EQ. 100))WRITE(6,*)'RPOT1',RHOTERM(IX,IZ),RHOSRC(IX,IZ),PSRC(IX,IZ)
	  RHOTERM(IX,IZ) = (ESRC(IX,IZ))*PSI(IX,IZ)**8*ASCAL(IX,IZ)**4
!       IF((IZ .EQ. 10) .AND. (IX .EQ. 100))WRITE(6,*)'RPOT2',RHOTERM(IX,IZ),ESRC(IX,IZ)
      VMETERM(IX,IZ)= ASCAL(IX,IZ)*PSL(IX,IZ)/PSI(IX,IZ)*VLOC(IX,IZ)**2./(OMGMET(IX,IZ)-OMG)
    END DO
  END DO

  !..................
  ! Compute the metric source terms for the vector potential
  !..................
  ! Assume parity on axis
  DO IZ=1,NR
    METERM(0,IZ)     =  METERM(1,IZ)
    METERM(NTH+1,IZ) =  METERM(NTH,IZ)
    GAMLOC(NTH+1,IZ) =  GAMLOC(NTH,IZ)
    GAMLOC(0,IZ)     =  GAMLOC(1,IZ)
    OMGMET(0,IZ)     =  OMGMET(1,IZ)
    OMGMET(NTH+1,IZ) =  OMGMET(NTH,IZ)
  END DO

  ! Assume parity at center
  ! Assume that the metric source term is smooth at the outer boundary.
  DO IX=1,NTH
    METERM(IX,0)    = METERM(NTH-IX+1,1)
    METERM(IX,NR+1) = METERM(IX,NR) + DR(NR)/DR(NR-1)*(METERM(IX,NR)-METERM(IX,NR-1)) !
    GAMLOC(IX,0)    = GAMLOC(NTH-IX+1,1)
    GAMLOC(IX,NR+1) = GAMLOC(IX,NR) + DR(NR)/DR(NR-1)*(GAMLOC(IX,NR)-GAMLOC(IX,NR-1))
    OMGMET(IX,0)    = OMGMET(NTH-IX+1,1)
    OMGMET(IX,NR+1) = OMGMET(IX,NR) + DR(NR)/DR(NR-1)*(OMGMET(IX,NR)-OMGMET(IX,NR-1))
  END DO

  LOGGAM=LOG(GAMLOC**2.)

  DO IX=1,NTH
    DO IZ=1,NR
      A1=-DR(IZ+1)/DR(IZ)/(DR(IZ)+DR(IZ+1))
      A2= (DR(IZ+1)-DR(IZ))/DR(IZ)/DR(IZ+1)
      A3= DR(IZ)/DR(IZ+1)/(DR(IZ)+DR(IZ+1))
      B1=-DTH(IX+1)/DTH(IX)/(DTH(IX)+DTH(IX+1))
      B2= (DTH(IX+1)-DTH(IX))/DTH(IX)/DTH(IX+1)
      B3= DTH(IX)/DTH(IX+1)/(DTH(IX)+DTH(IX+1))

      DRMETTERM(IX,IZ) = A1*METERM(IX,IZ-1)+A3*METERM(IX,IZ+1)+A2*METERM(IX,IZ)
      DTMETTERM(IX,IZ) = B1*METERM(IX-1,IZ)+B3*METERM(IX+1,IZ)+B2*METERM(IX,IZ)

      DRGAML(IX,IZ)=A1*LOGGAM(IX,IZ-1)+A3*LOGGAM(IX,IZ+1)+A2*LOGGAM(IX,IZ)
      DTGAML(IX,IZ)=B1*LOGGAM(IX-1,IZ)+B3*LOGGAM(IX+1,IZ)+B2*LOGGAM(IX,IZ)

      DROMGM(IX,IZ)=A1*OMGMET(IX,IZ-1)+A3*OMGMET(IX,IZ+1)+A2*OMGMET(IX,IZ)
      DTOMGM(IX,IZ)=B1*OMGMET(IX-1,IZ)+B3*OMGMET(IX+1,IZ)+B2*OMGMET(IX,IZ)

    END DO
  END DO


END SUBROUTINE SOURCEPOT

! **************************************************************
! **************************************************************

SUBROUTINE VECPOTPHI !!! QUI

!  ============================================================
!  Purpose : this subroutine solve the Grad-Shafranov eq. for
!    the phi-component of the vector potential in spherical
!    coordinates using coordinate quantities. The procedure and
!    notation are similar to , used to solve the phi-component
!    of the vector poisson equation for the metric.
!  ============================================================
!
!  APHI = azimuthal component of the vector potential
!        (normalized with r*sin(th))
!  RHO  = source term of the Grad-Shafranov Equation
!  APM  = maximum value of APHI on the domain
!  APHIMAX= maximum value of APHI at a distance R(CUT*WSURF)
!  ============================================================

  USE SYSTEMXNS
! #ifdef _OPENMP
!   use omp_lib
! #endif

  IMPLICIT NONE

  INTEGER,PARAMETER :: MAXIT = 1500
  INTEGER :: I,INFO,ITP,ITER
  REAL,DIMENSION(-1:NTH+2,0:NR+1):: RHO,DAPHIDT,DAPHIDR,APHIOLD,APHIOLD2
  REAL :: A1,A2,A3,B1,B2,B3, ERROR
  REAL :: TERM1, TERM2, TERM3, TERM4

  real(kind(1.d0)) :: starttime,endtime

  IF(APM.EQ.0.) APM=1.

  ! Iterate on the potential for relaxation
  DO ITER=1,MAXIT
    !Impose boundary conditions on the potential assuming
    !parity in the axis and at the center. Assume smoothness at the
    !outer boundary.
    DO IZ=1,NR
      APHI(0,IZ) =  -APHI(1,IZ)
      APHI(NTH+1,IZ) =  -APHI(NTH,IZ)
    END DO
    DO IX=1,NTH
      APHI(IX,0) = -APHI(NTH-IX+1,1)
      APHI(IX,NR+1) =  APHI(IX,NR) + DR(NR)/DR(NR-1)*(APHI(IX,NR)-APHI(IX,NR-1))
    END DO

    !Compute derivatives to be used in the source term
    DO IX=1,NTH
      DO IZ=1,NR
        A1=-DR(IZ+1)/DR(IZ)/(DR(IZ)+DR(IZ+1))
        A2= (DR(IZ+1)-DR(IZ))/DR(IZ)/DR(IZ+1)
        A3= DR(IZ)/DR(IZ+1)/(DR(IZ)+DR(IZ+1))
        B1=-DTH(IX+1)/DTH(IX)/(DTH(IX)+DTH(IX+1))
        B2= (DTH(IX+1)-DTH(IX))/DTH(IX)/DTH(IX+1)
        B3= DTH(IX)/DTH(IX+1)/(DTH(IX)+DTH(IX+1))

        DAPHIDR(IX,IZ) = A1*APHI(IX,IZ-1)+A3*APHI(IX,IZ+1)+A2*APHI(IX,IZ)
        DAPHIDT(IX,IZ) = B1*APHI(IX-1,IZ)+B3*APHI(IX+1,IZ)+B2*APHI(IX,IZ)
      END DO
    END DO

    ! Initialize the source term in the domain

    RHO(:,:)=0.
    DO IX=1,NTH
      DO IZ=1,NR
        IF(IPOL)  TERM1 = RHOTERM(IX,IZ)*R(IZ)*SIN(TH(IX))*KBPOL*(1.+Csi*R(IZ)*SIN(TH(IX))*APHI(IX,IZ))
        IF(ITWT)  TERM1 = RHOTERM(IX,IZ)*R(IZ)*SIN(TH(IX))*KBTT
        IF(IZ .GT. WSURF(IX)) TERM1=0.
        TERM2 = DRMETTERM(IX,IZ)*(APHI(IX,IZ)/R(IZ)+DAPHIDR(IX,IZ))
        TERM3 = DTMETTERM(IX,IZ)*(APHI(IX,IZ)*COS(TH(IX))/SIN(TH(IX))+DAPHIDT(IX,IZ))/R(IZ)**2
        IF(IPOL) TERM4=0.
        IF(ITWT) THEN
           TERM4= PSI(IX,IZ)**6./PSL(IX,IZ)**2.*ATWT**2/(ZETA+1.)/(R(IZ)*SIN(TH(IX)))*&
                  MAX( APHI(IX,IZ)*R(IZ)*SIN(TH(IX))-APHIMAX,0.)*&
                 (APHI(IX,IZ)*R(IZ)*SIN(TH(IX))-APHIMAX)**(2.*ZETA)/APM**(2.*ZETA+1.)
        END IF
        RHO(IX,IZ) = -1.*(TERM1+TERM2+TERM3+TERM4)
       END DO
    END DO

    !Solve the equation with a semi-spectral method

    APHIOLD=APHI
    CALL SOLVEAPHI(RHO)

    IF(ITER.GT.1)THEN
       APHI=(QAPHI*APHI+(1.-QAPHI)*APHIOLD)
    END IF

    !Set the surface CUTSURF to control the control the size of the
    !twisted magnetosphere. CUT=1 provides standard TT models
    IF(.NOT. STRETCH) THEN
      CUTSURF=INT(CUT*WSURF)
    ELSE
      DO IX=1,NTH
        IF( CUT*R(WSURF(IX)) .LT. RREG ) THEN
          CUTSURF(IX)=INT(CUT*WSURF(IX))
        ELSE
          CUTSURF(IX)= INT(NRREG+LOG(1.+NRREG*(CUT*R(WSURF(IX))-RREG)*(STRR-1.)/(STRR*RREG))/LOG(STRR))
        END IF
      END DO
    ENDIF

    !Evaluate the maximum value of APHI on CUTSURF
    APHIMAX=0.
    DO IX=1,NTH
       APHIMAX=MAX(APHIMAX,APHI(IX,CUTSURF(IX))*R(CUTSURF(IX))*SIN(TH(IX)))
    END DO

    !Evaluate the maximum of APHI on the entire domain
    APM=0.
    DO IX=1,NTH
      DO IZ=1,NR
         APM=MAX(APM,APHI(IX,IZ)*R(IZ)*SIN(TH(IX)))
      END DO
    END DO

    IF(WCONVA) THEN
       OPEN(12,FILE=trim(adjustl(subdirpath))//'Apconv.dat',access='append')
       WRITE(12,*) ITER, APM
       CLOSE(12)
    END IF

    !Check convergence
    ERROR=MAXVAL(ABS(APHIOLD(1:NTH,1:NR)-APHI(1:NTH,1:NR)))
    IF(VERBOSE .AND. (.NOT. MPICODE))WRITE(6,*) ITER,'MaxAphi & MaxErrAphi =', APM,error
    IF(ERROR .LT. TOLCONV*MAX(KBPOL,KBTT)) EXIT
  END DO

END SUBROUTINE VECPOTPHI

! ********************************************************
! ********************************************************

SUBROUTINE MXWLSOL(ILOOP) !!! QUI

!  =============================================================
!  Purpose: this subroutine solves the Maxwell-Ampere and the
!    Maxwell-Gauss equations. It finally solves for the
!    harmonic function required to impose MHD inside the star.
!    The additional degree of fredoom related to the arbitrary
!    harmonic function correspond to the global net charge of
!    star. If QNULL is set to true, the code searches for
!    globally neutral stars, otherwise it minimizes the electric
!    field at the pole.
!
!  APHI    = azimuthal component of the vector potential
!            (normalized with r*sin(th))
!  ATIM    = global time component of the 4-potential
!  ATIMIN  = time component of the 4-potential inside the star
!  ATIMOUT = time component of the 4-potential outside the star
!  JMXL    = current density
!  RHOEMXL = charge denity
!  RHO     = source term of Maxwell equations
!  APM     = maximum value of APHI on the domain
!  ATMM    = maximum value of ATIM in the domain
!  CC1     = monopolar content of ATIM
!  CONSTANT= arbitrary integral constant for the APHI potential
!  =============================================================

  USE SYSTEMXNS
  IMPLICIT NONE

  INTEGER, PARAMETER :: MAXIT = 100 !1500
  INTEGER :: I, INFO,ILOOP
  REAL, DIMENSION(-1:NTH+2,0:NR+1):: DAPHIDR,DAPHIDT,DATIMDR,DATIMDT
  REAL, DIMENSION(-1:NTH+2,0:NR+1):: RHO,APHIOLD, APHIOLD2, ATIMOLD

  ! Array for the Legendre expansion
  REAL,DIMENSION(NGQ) :: XGQ,WGQ,Y
  REAL,DIMENSION(0:MLS,NGQ) :: PN,PD
  REAL,DIMENSION(0:MLS) :: COEF,PDX,PNX
  REAL,DIMENSION(0:NR+1,0:MLS) :: TAB

  ! Array for the radial solution
  REAL,DIMENSION(NR) :: DC,DC1,DCP,SS1
  REAL,DIMENSION(NR-1) :: DL1,DU1,DLP,DUP
  REAL,DIMENSION(1:NR,0:MLS) :: PHI
  REAL :: A1,A2,A3,A4,A5,A6,B1,B2,B3
  REAL :: DCI,DCF,ERROR,TMP
  REAL :: ATMM,CC1,CONSTANT,CONSTANTTMP
  REAL :: TERM1,TERMa,TERMb,TERM5
  REAL :: TERM2,TERM2a,TERM2b
  REAL :: TERM3,TERM3a,TERM3b
  REAL :: TERM4,TERM4a,TERM4b
  REAL :: ERROR2
  INTEGER :: IL,ITER,ITP

  DO ITER=1, MAXIT
    !Impose boundary conditions on the potentials
    !Assume parity on axis
    DO IZ=1,NR
      APHI(0,IZ) =  -APHI(1,IZ)
      APHI(NTH+1,IZ) =  -APHI(NTH,IZ)
      ATIM(0,IZ) =  ATIM(1,IZ)
      ATIM(NTH+1,IZ) =  ATIM(NTH,IZ)

      ATIMIN(0,IZ) =  ATIMIN(1,IZ)
      ATIMIN(NTH+1,IZ) =  ATIMIN(NTH,IZ)
      ATIMOUT(0,IZ) =  ATIMOUT(1,IZ)
      ATIMOUT(NTH+1,IZ) =  ATIMOUT(NTH,IZ)

    END DO

    !Assume parity at center and smoothness at the outer boundary
    DO IX=1,NTH
      APHI(IX,0) = -APHI(NTH-IX+1,1)
      APHI(IX,NR+1) =  APHI(IX,NR) + DR(NR)/DR(NR-1)*(APHI(IX,NR)-APHI(IX,NR-1))
      ATIM(IX,0) =  ATIM(NTH-IX+1,1)
      ATIM(IX,NR+1) =  ATIM(IX,NR) + DR(NR)/DR(NR-1)*(ATIM(IX,NR)-ATIM(IX,NR-1))

      ATIMIN(IX,0) =  ATIMIN(NTH-IX+1,1)
      ATIMIN(IX,NR+1) =  ATIMIN(IX,NR) + DR(NR)/DR(NR-1)*(ATIMIN(IX,NR)-ATIMIN(IX,NR-1))

      ATIMOUT(IX,0) =  ATIMOUT(NTH-IX+1,1)
      ATIMOUT(IX,NR+1) =  ATIMOUT(IX,NR) + DR(NR)/DR(NR-1)*(ATIMOUT(IX,NR)-ATIMOUT(IX,NR-1))
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

        DAPHIDR(IX,IZ) = A1*APHI(IX,IZ-1)+A3*APHI(IX,IZ+1)+A2*APHI(IX,IZ)
        DAPHIDT(IX,IZ) = B1*APHI(IX-1,IZ)+B3*APHI(IX+1,IZ)+B2*APHI(IX,IZ)

        IF(IZ .LE. WSURF(IX)) THEN
           DATIMDR(IX,IZ) = A1*ATIMIN(IX,IZ-1)+A3*ATIMIN(IX,IZ+1)+A2*ATIMIN(IX,IZ)
           DATIMDT(IX,IZ) = B1*ATIMIN(IX-1,IZ)+B3*ATIMIN(IX+1,IZ)+B2*ATIMIN(IX,IZ)
        ELSE
           DATIMDR(IX,IZ) = A1*ATIMOUT(IX,IZ-1)+A3*ATIMOUT(IX,IZ+1)+A2*ATIMOUT(IX,IZ)
           DATIMDT(IX,IZ) = B1*ATIMOUT(IX-1,IZ)+B3*ATIMOUT(IX+1,IZ)+B2*ATIMOUT(IX,IZ)
        ENDIF
      END DO
    END DO

    !------------------------------------------------------------
    ! Solve Maxwell Ampere equation
    ! -----------------------------------------------------------
    !Initialize the source term for the Maxwell-Ampere equation
    DO IX=1,NTH
      DO IZ=1,NR
        TERM1=RHOTERM(IX,IZ)/PSI(IX,IZ)**8./ASCAL(IX,IZ)**4* GAMLOC(IX,IZ)**2.*KBPOL

        TERM2a=DRGAML(IX,IZ)*(APHI(IX,IZ)/R(IZ)+DAPHIDR(IX,IZ))
        TERM2b=DTGAML(IX,IZ)*(APHI(IX,IZ)*COS(TH(IX))/SIN(TH(IX))+DAPHIDT(IX,IZ))/R(IZ)**2.
        TERM2=(TERM2a+TERM2b)/(ASCAL(IX,IZ)**4*PSI(IX,IZ)**8.*R(IZ)*SIN(TH(IX)))

        TERM3a=DROMGM(IX,IZ)*(APHI(IX,IZ)/R(IZ)+DAPHIDR(IX,IZ))
        TERM3b=DTOMGM(IX,IZ)*(APHI(IX,IZ)*COS(TH(IX))/SIN(TH(IX))+DAPHIDT(IX,IZ))/R(IZ)**2.
        TERM3=(TERM3a+TERM3b)*(OMGMET(IX,IZ)-OMG)*R(IZ)*SIN(TH(IX))/PSL(IX,IZ)**2./PSI(IX,IZ)**2./ASCAL(IX,IZ)**4

        JPHIMXL(IX,IZ)=TERM1-TERM2+TERM3

        !Reset current density outside the star
        IF(IZ .GT. WSURF(IX)) JPHIMXL(IX,IZ)=0.

      END DO
    END DO

    !Evaluate the source term for the Maxwell-Ampere equation
    DO IX=1,NTH
      DO IZ=1,NR
        TERM1=JPHIMXL(IX,IZ)*R(IZ)*SIN(TH(IX))*ASCAL(IX,IZ)**4*PSI(IX,IZ)**8.

        TERMa= APHI(IX,IZ)/R(IZ)+DAPHIDR(IX,IZ)
        TERMb= APHI(IX,IZ)*COS(TH(IX))/SIN(TH(IX))+DAPHIDT(IX,IZ)

        TERM2a=DROMGM(IX,IZ)*(DATIMDR(IX,IZ)+R(IZ)*SIN(TH(IX))*OMGMET(IX,IZ)*TERMa)
        TERM2b=DTOMGM(IX,IZ)*(DATIMDT(IX,IZ)+R(IZ)*SIN(TH(IX))*OMGMET(IX,IZ)*TERMb)/R(IZ)**2.
        TERM2=OMTERM(IX,IX)*(TERM2a+TERM2b)/R(IZ)/SIN(TH(IX))

        TERM3a= DRMETTERM(IX,IZ)*TERMa
        TERM3b= DTMETTERM(IX,IZ)*TERMb/R(IZ)**2

        IF(IZ .GT. WSURF(IX)) TERM1=0. !Redundant

        RHO(IX,IZ)= -TERM1+TERM2-TERM3a-TERM3b
      END DO
    END DO

    !Solve equation by using a semipectral method
    !The solution is stored in APHI declared in the system module

    APHIOLD=APHI
    CALL SOLVEAPHI(RHO)

    !Exploit dumping to stabilize convergence
    IF(ITER.GT.1) APHI= QAPHI*APHI+(1.-QAPHI)*APHIOLD


  !------------------------------------------------------------------------------
  !Solve Maxwell-Gauss equation
  !------------------------------------------------------------------------------

    !Initialize the charge density
    DO IX=1,NTH
      DO IZ=1,NR
        TERM1=VMETERM(IX,IZ)*RHOTERM(IX,IZ)/PSI(IX,IZ)**8/ASCAL(IX,IZ)**4*GAMLOC(IX,IZ)**2.*KBPOL

        TERM2a=DRGAML(IX,IZ)*(APHI(IX,IZ)/R(IZ)+DAPHIDR(IX,IZ))
        TERM2b=DTGAML(IX,IZ)*(APHI(IX,IZ)*COS(TH(IX))/SIN(TH(IX))+DAPHIDT(IX,IZ))/R(IZ)**2.
        TERM2=(TERM2a+TERM2b)/VMETERM(IX,IZ)/PSI(IX,IZ)**4./ASCAL(IX,IZ)**2*R(IZ)*SIN(TH(IX))

        TERM3a=DROMGM(IX,IZ)*(APHI(IX,IZ)/R(IZ)+DAPHIDR(IX,IZ))
        TERM3b=DTOMGM(IX,IZ)*(APHI(IX,IZ)*COS(TH(IX))/SIN(TH(IX))+DAPHIDT(IX,IZ))/R(IZ)**2.
        TERM3=(TERM3a+TERM3b)*R(IZ)*SIN(TH(IX))/PSL(IX,IZ)/PSI(IX,IZ)**3/ASCAL(IX,IZ)**3

        RHOEMXL(IX,IZ)=-TERM1 +TERM2 -TERM3

        !Impose electrovacuum outside the star
        IF(IZ .GT. WSURF(IX)) RHOEMXL(IX,IZ)=0.

      END DO
    END DO

    !Initialize the source of Maxwell-Gauss equation
    DO IX=1,NTH
      DO IZ=1,NR
        TERM1= ASCAL(IX,IZ)**3*PSL(IX,IZ)*PSI(IX,IZ)**3.*RHOEMXL(IX,IZ) + &
           OMGMET(IX,IZ)*ASCAL(IX,IZ)**4*PSI(IX,IZ)**8.*R(IZ)**2.*SIN(TH(IX))**2.*JPHIMXL(IX,IZ)

        TERMa= APHI(IX,IZ)/R(IZ)+DAPHIDR(IX,IZ)
        TERMb= APHI(IX,IZ)*COS(TH(IX))/SIN(TH(IX))+DAPHIDT(IX,IZ)

        TERM2a=DROMGM(IX,IZ)*(DATIMDR(IX,IZ)+OMGMET(IX,IZ)*TERMA*R(IZ)*SIN(TH(IX)))
        TERM2b=DTOMGM(IX,IZ)*(DATIMDT(IX,IZ)+OMGMET(IX,IZ)*TERMB*R(IZ)*SIN(TH(IX)))
        TERM2=OMGMET(IX,IZ)*OMTERM(IX,IZ)*(TERM2a+TERM2b)

        TERM3a=DROMGM(IX,IZ)*TERMA*R(IZ)*SIN(TH(IX))
        TERM3b=DTOMGM(IX,IZ)*TERMB*SIN(TH(IX))/R(IZ)
        TERM3 = TERM3a + TERM3b

        TERM4a=DRMETTERM(IX,IZ)*( DATIMDR(IX,IZ) +2.*OMGMET(IX,IZ)*TERMa*R(IZ)*SIN(TH(IX)) )
        TERM4b=DTMETTERM(IX,IZ)*( DATIMDT(IX,IZ) +2.*OMGMET(IX,IZ)*TERMb*R(IZ)*SIN(TH(IX)) )/R(IZ)**2.
        TERM4 = TERM4a + TERM4b

        TERM5= 2.*OMGMET(IX,IZ)*SIN(TH(IX))*(TERMa + TERMB/R(IZ)/TAN(TH(IX)))

        IF(IZ .GT. WSURF(IX)) TERM1=0. !Redundant

        RHO(IX,IZ)= TERM1 - TERM2 - TERM3 + TERM4 - TERM5

      END DO
    END DO

    !Solve the equation by using a semi-spectral method
    ATIMOLD=ATIM
    CALL SOLVEATIME(RHO,CC1)

    !Find the arbitrary harmonic function depending on the requirement
    !for the net charge of the star and correct the electric potential
    IF(QNULL) THEN
      !Initialize the harmonic function (correct on the surface)
      DO IZ=1,NR
        DO IX=1,NTH
          ATIMARM(IX,IZ)=-OMG*APHI(IX,IZ)*R(IZ)*SIN(TH(IX))-ATIM(IX,IZ) + CONSTANT
        END DO
      END DO

      !Solve Laplacian=0 with boundary condition at the stellar surface
      CALL LAPLACE(CC1,ATIMARM,ATIMIN,ATIMOUT)
      IF(XNSERR.NE.0) THEN
         WRITE(6,*)"ERROR: from HYDROEQ, LAPLACE failure in LUDCMP"
         STOP
      ENDIF
      !Determine the integration constant to minimize the monopolar content of the solution

      ATIMIN=ATIM+ATIMIN
      ATIMOUT=ATIM+ATIMOUT
      ATIM=ATIM +ATIMARM

      CONSTANT=CONSTANT+CC1
    ELSE
      !Initialize the harmonic function (correct on the surface)
      DO IZ=1,NR
        DO IX=1,NTH
          ATIMARM(IX,IZ)=-OMG*APHI(IX,IZ)*R(IZ)*SIN(TH(IX))-ATIM(IX,IZ) + CONSTANT
        END DO
      END DO
      !Solve Laplacian=0 with boundary condition at the stellar surface
      CALL LAPLACE(CC1,ATIMARM,ATIMIN,ATIMOUT)

      ATIMIN=ATIM+ATIMIN
      ATIMOUT=ATIM+ATIMOUT
      ATIM=ATIM +ATIMARM

      !Determine the arbitrary constant to minimize the electric field at the pole
      CONSTANTTMP=MINVAL(ATIM(1,WSURF(1)+1:NR))
      IF(CONSTANTTMP .LE. ATIM(1,WSURF(1)-10)) CONSTANT =CONSTANTTMP
    END IF

    APM=0.
    ATMM=0
    DO IX=1,NTH
      DO IZ=1,NR
        APM=MAX(APM,APHI(IX,IZ)*R(IZ)*SIN(TH(IX)))
        ATMM=MAX(ATMM,ABS(ATIM(IX,IZ)))
      END DO
    END DO

    !Check for convergence
    ERROR=MAXVAL(ABS(APHIOLD(1:NTH,1:NR)-APHI(1:NTH,1:NR)))
    ERROR2=MAXVAL(ABS(ATIMOLD(1:NTH,1:NR)-ATIM(1:NTH,1:NR)))
   ! IF(VERBOSE .AND. (.NOT. MPICODE)) WRITE(6,*) ITER,'MaxAphi,MaxAt & MaxErrAphi =', APM,ATMM,ERROR,ERROR2
    IF((ERROR .LT. TOLCONV*MAX(KBPOL,KBTT)).AND.(ERROR2 .LT. OMG*TOLCONV*MAX(KBPOL,KBTT))) EXIT

  END DO

  IF(WCONVA) THEN
    OPEN(12,FILE=TRIM(ADJUSTL(subdirpath))//'Apconv.dat',access='append')
    WRITE(12,*) ILOOP,APM,ITER,ERROR
    CLOSE(12)
  END IF

  IF(WCONVA) THEN
     OPEN(12,FILE=TRIM(ADJUSTL(subdirpath))//'Atconv.dat',access='append')
     WRITE(12,*) ILOOP,ATMM,ITER,ERROR2
     CLOSE(12)
  END IF

END SUBROUTINE MXWLSOL


!*******************************************************************************
! EQUATION SOLVERS
!*******************************************************************************

SUBROUTINE SOLVEAPHI(RHO) !!! QUI
!  ============================================================
!  Purpose : this subroutine solves for the Grad-Shafranov
!     equation or for the phi-component of the Maxwell-Ampere
!     equation. The procedure and notation are similar to the
!     ones in SHIFTPHI, used to solve the phi-component
!     of the vector poisson equation for the metric.
!
!  RHO = source term
!  ============================================================

  USE SYSTEMXNS
  USE FUNCTIONS
  IMPLICIT NONE

  REAL, DIMENSION(-1:NTH+2,0:NR+1) :: RHO

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

    ! Compute the matrix of Weighted Spherical Harmonics
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

    ! Set the matrix elements (case L=0)
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

    DO IL=1,MLS
      !BC at inner radius (Origin)
      A1= 2./DR(1)/(DR(1)+DR(1+1))
      A4=-DR(1+1)/DR(1)/(DR(1)+DR(1+1))
      DCP(1)=DCI-(A1+2.*A4/R(1))

      !BC at outer radius
      !Various multipole must decay as R^(L+1)
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

    !Compute the matrix of Legendre Polyn in the grid points
    DO IX=1,NTH
      CALL LPN(MLS,XX(IX),PNX(0:MLS),PDX(0:MLS))
      DO I=1,MLS
         PDX(I)=SQRT((2*I+1.)/4./Pi)*PDX(I)
      END DO
      ! Compute the potential on the grid
      DO IZ=1,NR
        APHI(IX,IZ)=-1.*DOT_PRODUCT(PDX(1:MLS)*SIN(TH(IX)),PHI(IZ,1:MLS))
      END DO
    END DO

END SUBROUTINE SOLVEAPHI

! ********************************************************
! ********************************************************

SUBROUTINE SOLVEATIME(RHO,CC1) !!! QUI
  !  ============================================================
  !  Purpose : this subroutine solves for the Maxwell-Gauss
  !     equation. The procedure and notation are similar to the
  !     ones in CONFORMAL subroutine, used to solve for the
  !     conformal factor in XNSMAIN.f90
  !
  !  RHO = Source term
  !  CC1 = monopolar content of the solution
  !  ============================================================

  USE SYSTEMXNS
  USE FUNCTIONS
  IMPLICIT NONE

  REAL, DIMENSION(-1:NTH+2,0:NR+1) :: RHO
  REAL :: CC1

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

     PHI(1:NR,IL)=SS1(1:NR)
  END DO

  ! Compute the matrix of Legendre Polyn in the grid points
  DO IX=1,NTH
    CALL LPN(MLS,XX(IX),PNX(0:MLS),PDX(0:MLS))
    ! Compute the potential on the grid
    DO IZ=1,NR
      ATIM(IX,IZ)=DOT_PRODUCT(PNX(0:MLS),PHI(IZ,0:MLS))
    END DO
  END DO

  !Store the monopolar content at a large radius
  !Used to solve the Laplace equation for the harmonic function
  CC1=PHI(NR-10,0)

END SUBROUTINE SOLVEATIME

! ********************************************************
! ********************************************************

SUBROUTINE LAPLACE(CC1,PHI,PHIIN,PHIOUT)
!  =============================================================
!  Purpose: this subroutine solves for harmonic function PHIARM.
!           It uses LUDCMP and LUBKSB from Numerical Recipes to
!           perform LU decomposition and backward subtitution.
!
!  AELL, BELL = minor and major axis of the surface ellipsoid
!  NSPE       = super-ellipsoid index
!  ELSURF     = fit of the surface
!  PHI(INPUT) = Armonic potential at the stellar surface
!  PHI(OUTPUT)= Armonic potential on the grid
!  PHIIN      = Armonic potential inside the star
!  PHIOUT     = Armonic potential outside the star
!  =============================================================

  USE SYSTEMXNS
  IMPLICIT NONE

  REAL, DIMENSION(-1:NTH+2,0:NR+1) :: PHI, PHIIN, PHIOUT
  REAL                             :: CC1

  !Parameters for the harmonic expansion
  REAL,DIMENSION(0:MLSL,0:MLSL) :: PN,PD
  REAL,DIMENSION(0:MLSL,0:MLSL) :: PNIN,PNOUT,PDIN,PDOUT
  REAL,DIMENSION(0:MLSL) :: WGQ,XGQ,XGQIN,XGQOUT,RXGQIN,RXGQOUT
  REAL,DIMENSION(0:MLSL) :: PNX,PDX

  !Parameters of LU decomposition and BKSB
  REAL, DIMENSION(0:MLSL,0:MLSL) :: MATRIX1, MATRIX2
  REAL, DIMENSION(0:MLSL) :: Y,YOIN,YOOUT
  INTEGER, DIMENSION(0:MLSL/2) :: INDC
  REAL :: DSIGN,YINTRP,DYINT

  !Interpolation and surface
  REAL, DIMENSION(-1:NTH+2) :: PHINTRP
  REAL, DIMENSION (0:NR,0:MLSL) :: PHIINTMP, PHIOUTMP
  REAL :: AELL,BELL,THEL,ATMP,BTMP,RMEAN
  REAL :: NSPE1,NSPE2
  REAL :: CC2,CC3
  !LOGICAL :: CASSINIAN
  INTEGER :: ITHEL,NGQL,NP,RTMP,I,J,CHECK

  !Fit the stellar surface with a superellipse with index NSPE
  AELL=R(WSURF(NTH/2)+1)
  BELL=R(WSURF(1)+1)

  IF(ABS(WSURF(NTH/2)-WSURF(1)) .LE. 3) THEN
    DO IX=1, NTH
      RTMP=WSURF(IX)+1
      ELSURF(IX)=AELL*BELL/SQRT(AELL**2.*COS(TH(IX))**2.+BELL**2.*SIN(TH(IX))**2.)
      !2D interpolation of the potential on the analytical curve
      CALL POLIN2(TH(IX-2:IX+2),R(RTMP-2:RTMP+2),PHI(IX-2:IX+2,RTMP-2:RTMP+2),&
                  5,5,TH(IX),ELSURF(IX),YINTRP,DYINT,CHECK)
      PHINTRP(IX)=YINTRP
    END DO
    NSPE=2.
  ELSE
    THEL=ATAN(AELL/BELL)
    ITHEL=INT(THEL*NTH/PI+0.5)
    IF(THEL-ITHEL .gt. 0.5) ITHEL=ITHEL+1
    NSPE1=LOG(2.)/LOG(BELL/R(WSURF(ITHEL)+1)/COS(TH(ITHEL)))
    NSPE2=LOG(2.)/LOG(AELL/R(WSURF(ITHEL)+1)/SIN(TH(ITHEL)))
    NSPE=(NSPE1+NSPE2)/2.
    DO IX=1,NTH
      RTMP=WSURF(IX)+1
      ATMP=ABS(AELL*COS(TH(IX)))
      BTMP=ABS(BELL*SIN(TH(IX)))
      ELSURF(IX)=AELL*BELL/(ATMP**NSPE+BTMP**NSPE)**(1./NSPE)
      !2D interpolation of the potential on the analytical curve
      CALL POLIN2(TH(IX-2:IX+2),R(RTMP-2:RTMP+2),PHI(IX-2:IX+2,RTMP-2:RTMP+2),&
                  5,5,TH(IX),ELSURF(IX),YINTRP,DYINT,CHECK)
      PHINTRP(IX)=YINTRP
    ENDDO
  ENDIF

  !Alternative fit by using a cassinian oval
  !CASSINIAN=.TRUE.
  !IF(CASSINIAN) THEN
  !   ACAS=sqrt((AELL**2.-BELL**2.)/2.)
  !   BCAS=sqrt((AELL**2.+BELL**2.)/2.)
  !   DO IX=1,NTH
  !      COSC=COS(2*TH(IX))
  !      SENC=SIN(2*TH(IX))
  !      RTMP=WSURF(IX)+1
  !      ELSURF(IX)=SQRT(ACAS**2.*(-COSC+SQRT(BCAS**4./ACAS**4.-SENC**2)))
  !      CALL POLIN2(TH(IX-2:IX+2),R(RTMP-2:RTMP+2),PHI(IX-2:IX+2,RTMP-2:RTMP+2),5,5,TH(IX),ELSURF(IX),YINTRP,DYINT,CHECK)
  !      PHINTRP(IX)=YINTRP
  !   ENDDO
  !ENDIF

  IF(CHECK.EQ.-1) THEN
     WRITE(6,*)"ERROR: from HYDROEQ, LAPLACE failure in POLIN2"
     STOP
  ENDIF

  ! Compute the Gauss-quadrature points XGQ and weights WGQ
  CALL LEGZO(MLSL+1,XGQ(0:MLSL),WGQ(0:MLSL))


  ! Compute the collocation points and the associated Legendre polinomia.
  ! For strongly deformed stars the starting Gauss-quadrature points
  ! are remodulated in an unevenly distributed collection of new collocation
  ! points to avoid aliasing.
  IF(BELL/AELL .GT. 0.9) THEN
     DO I=0,MLSL
        XGQIN(I)=XGQ(I)
        XGQOUT(I)=XGQ(I)
        CALL LPN(MLSL,XGQOUT(I),PNOUT(0:MLSL,I),PDOUT(0:MLSL,I))
        CALL LPN(MLSL,XGQIN(I),PNIN(0:MLSL,I),PDIN(0:MLSL,I))
     ENDDO
  ELSE
     DO I=0,MLSL
        XGQOUT(I)=XGQ(I)/(ABS(XGQ(I))+1.e-9)*SQRT(ABS(XGQ(I)))
        XGQIN(I)=(ABS(XGQ(I)))**1.3*XGQ(I)/ABS(XGQ(I)+1.e-9)
        CALL LPN(MLSL,XGQOUT(I),PNOUT(0:MLSL,I),PDOUT(0:MLSL,I))
        CALL LPN(MLSL,XGQIN(I),PNIN(0:MLSL,I),PDIN(0:MLSL,I))
     END DO
  ENDIF

  !Impose boundary contition for POLINT
  ELSURF(0) = ELSURF(1)      ; ELSURF(-1)= ELSURF(2)
  ELSURF(NTH+1)= ELSURF(NTH) ;ELSURF(NTH+2)= ELSURF(NTH-1)

  PHINTRP(0)= PHINTRP(1)       ; PHINTRP(-1)= PHINTRP(2)
  PHINTRP(NTH+1)= PHINTRP(NTH) ; PHINTRP(NTH+2)= PHINTRP(NTH-1)

  !Interpolate the potential on the new quadrature points
  CALL POLINT(XX(-1:NTH+2),PHINTRP(-1:NTH+2),NTH+4,XGQOUT,YOOUT,MLSL+1)
  CALL POLINT(XX(-1:NTH+2),ELSURF(-1:NTH+2),NTH+4,XGQOUT,RXGQOUT,MLSL+1)

  CALL POLINT(XX(-1:NTH+2),PHINTRP(-1:NTH+2),NTH+4,XGQIN,YOIN,MLSL+1)
  CALL POLINT(XX(-1:NTH+2),ELSURF(-1:NTH+2),NTH+4,XGQIN,RXGQIN,MLSL+1)

  ! Write interpolated values to check interpolation
  IF(CHUP.AND.(WRT.OR.WRTF)) THEN
     OPEN(12,FILE=trim(adjustl(subdirpath))//'CheckInterpl.dat')
     WRITE(12,*)MLSL+1,NR
     DO IX=0,MLSL
        WRITE(12,*) XGQIN(IX),RXGQIN(IX),YOIN(IX),XGQOUT(IX),RXGQOUT(IX),YOOUT(IX)
     END DO
     CLOSE(12)
  ENDIF

  ! Mean Value for the ellipsoidal surface
  RMEAN=SUM(ELSURF(1:NTH))/NTH

  ! Compute the matrix coefficient for the interior solution
  DO I=0,MLSL/2
    DO J=0,MLSL/2
      MATRIX1(I,J)= (RXGQIN(I)/RMEAN)**(2*J)*SQRT((4*J+1.)/4./Pi)*PNIN(2*J,I)
    END DO
  END DO

  ! Compute the matrix coefficient for the exterior solution
  DO I=0,MLSL/2
    DO J=0,MLSL/2
      MATRIX2(I,J)= (RMEAN/RXGQOUT(I))**(2*J+1)*SQRT((4*J+1.)/4./Pi)*PNOUT(2*J,I)
    END DO
  END DO

  ! Find the interior harmonic function
  Y=YOIN
  CALL LUDCMP(MATRIX1(0:MLSL/2,0:MLSL/2),MLSL/2+1,MLSL/2+1,indc,dsign,CHECK)
  CALL LUBKSB(MATRIX1(0:MLSL/2,0:MLSL/2),MLSL/2+1,MLSL/2+1,indc,Y(0:MLSL/2))

  IF(CHECK.EQ.-1) THEN
     WRITE(6,*)"ERROR: from HYDROEQ, LAPLACE failure in LUDCMP"
     STOP
  ENDIF

  DO IZ=1,NR
     PHIINTMP(IZ,0)=Y(0)*SQRT(1./4./Pi)
     DO J=1,MLSL/2
        PHIINTMP(IZ,2*J)=Y(J)*(R(IZ)/RMEAN)**(2*J)*SQRT((4*J+1.)/4./Pi)
        PHIINTMP(IZ,2*J-1)=0.
     END DO
  END DO
  DO IX=1,NTH
     CALL LPN(MLSL,XX(IX),PNX(0:MLSL),PDX(0:MLSL))
     DO IZ=1,NR
        PHIIN(IX,IZ)= DOT_PRODUCT(PNX(0:MLSL),PHIINTMP(IZ,0:MLSL))
     END DO
  END DO

  ! Find the interior exterior harmonic function
  Y=YOOUT
  CALL LUDCMP(MATRIX2(0:MLSL/2,0:MLSL/2),MLSL/2+1,MLSL/2+1,INDC,DSIGN,CHECK)
  CALL LUBKSB(MATRIX2(0:MLSL/2,0:MLSL/2),MLSL/2+1,MLSL/2+1,INDC,Y(0:MLSL/2))

  IF(CHECK.EQ.-1) THEN
     WRITE(6,*)"ERROR: from HYDROEQ, LAPLACE failure in LUDCMP"
     STOP
  ENDIF

  DO  IZ=1,NR
    PHIOUTMP(IZ,0)=Y(0)/(R(IZ)/RMEAN)*SQRT(1./4./Pi)
    DO J=1,MLSL/2
      PHIOUTMP(IZ,2*J)=Y(J)/(R(IZ)/RMEAN)**(2*J+1)*SQRT((4*J+1.)/4./Pi)
      PHIOUTMP(IZ,2*J-1)=0.
    END DO
  END DO

  DO IX=1,NTH
    CALL LPN(MLSL,XX(IX),PNX(0:MLSL),PDX(0:MLSL))
    DO IZ=1,NR
      PHIOUT(IX,IZ) = DOT_PRODUCT(PNX(0:MLSL),PHIOUTMP(IZ,0:MLSL) )
    END DO
  END DO

  ! Save the monopolar component of the harmonic function
  CC2=PHIOUTMP(NR-10,0)

  ! Build the global harmonic potential
  DO IX=1,NTH
    DO IZ=1,NR
      IF(R(IZ) .GE. ELSURF(IX)) THEN
        PHI(IX,IZ)=PHIOUT(IX,IZ)
      ELSE
        PHI(IX,IZ)=PHIIN(IX,IZ)
      ENDIF
    END DO
  END DO

  !Write data to check correct behaviour on the surface
  IF(CHUP.AND.(WRT.OR.WRTF)) THEN
     OPEN(12,FILE=trim(adjustl(subdirpath))//'CheckSurf.dat')
     WRITE(12,*) NTH,AELL,BELL
     DO IX=1,NTH
        WRITE(12,*) R(WSURF(IX)+1), ELSURF(IX),PHI(IX,WSURF(IX)+1),PHINTRP(IX),&
             PHIIN(IX,WSURF(IX)),PHIOUT(IX,WSURF(IX)+1)
     END DO
     CLOSE(12)
  ENDIF

  Y(:)=1.
  CALL LUBKSB(MATRIX2(0:MLSL/2,0:MLSL/2),MLSL/2+1,MLSL/2+1,INDC,Y(0:MLSL/2))

  DO IZ=1,NR
     PHIOUTMP(IZ,0)=Y(0)/(R(IZ)/RMEAN)*SQRT(1./4./Pi)
     DO J=1,MLSL/2
        PHIOUTMP(IZ,2*J)=Y(J)/(R(IZ)/RMEAN)**(2*J+1)*SQRT((4*J+1.)/4./Pi)
        PHIOUTMP(IZ,2*J-1)=0.
     END DO
  END DO
  CC3=PHIOUTMP(NR-10,0)
  CC1=-(CC1+CC2)/CC3

  RETURN

END SUBROUTINE LAPLACE

!*******************************************************************************
!  NR ROUTINES
!*******************************************************************************

SUBROUTINE ludcmp(a,n,np,indx,d,check)

  IMPLICIT NONE

  INTEGER n,np,indx(n),NMAX
  REAL d,a(np,np),TINY
  PARAMETER (NMAX=50,TINY=1.0e-30) !Largest expected n, and a small number.
  INTEGER i,imax,j,k,check
  REAL aamax,dum,sum,vv(NMAX)

  d=1.

  DO i=1,n
     aamax=0.
     DO j=1,n
        IF (ABS(a(i,j)).GT.aamax) aamax=ABS(a(i,j))
     ENDDO
     IF (aamax.EQ.0.) THEN
        PRINT*, 'singular matrix in ludcmp'
        check=-1
        RETURN
     ENDIF
     vv(i)=1./aamax
  ENDDO
  DO j=1,n
     DO i=1,j-1
        sum=a(i,j)
        DO k=1,i-1
           sum=sum-a(i,k)*a(k,j)
        ENDDO
        a(i,j)=sum
     ENDDO
     aamax=0.
     DO i=j,n
        sum=a(i,j)
        DO k=1,j-1
           sum=sum-a(i,k)*a(k,j)
        ENDDO
        a(i,j)=sum
        dum=vv(i)*ABS(sum)
        IF (dum.GE.aamax) THEN
           imax=i
           aamax=dum
        ENDIF
     ENDDO
     IF (j.NE.imax)THEN
        DO k=1,n
           dum=a(imax,k)
           a(imax,k)=a(j,k)
           a(j,k)=dum
        ENDDO
        d=-d
        vv(imax)=vv(j)
     ENDIF
     indx(j)=imax
     IF(a(j,j).EQ.0.) a(j,j)=TINY
     IF(j.NE.n) THEN
        dum=1./a(j,j)
        DO i=j+1,n
           a(i,j)=a(i,j)*dum
        ENDDO
     ENDIF
  ENDDO

  RETURN

END SUBROUTINE ludcmp

SUBROUTINE lubksb(a,n,np,indx,b)
INTEGER n,np,indx(n)
REAL a(np,np),b(n)
INTEGER i,ii,j,ll
REAL sum

  ii=0

  DO i=1,n
    ll=indx(i)
    sum=b(ll)
    b(ll)=b(i)

    IF (ii.NE.0)THEN
      DO j=ii,i-1
        sum=sum-a(i,j)*b(j)
      END DO
    ELSE IF (sum.NE.0.) THEN
       ii=i
    ENDIF
    b(i)=sum
  ENDDO

  DO i=n,1,-1
    sum=b(i)
    DO j=i+1,n
      sum=sum-a(i,j)*b(j)
    ENDDO
    b(i)=sum/a(i,i)
  END DO

  RETURN

END SUBROUTINE lubksb


SUBROUTINE polin2(x1a,x2a,ya,m,n,x1,x2,y,dy,check)

  IMPLICIT NONE

  INTEGER m,n,NMAX,MMAX
  REAL dy,x1,x2,y,x1a(m),x2a(n),ya(m,n)
  PARAMETER (NMAX=20,MMAX=20)
  INTEGER j,k,check
  REAL ymtmp(MMAX),yntmp(NMAX)
  DO j=1,m
    DO k=1,n
      yntmp(k)=ya(j,k)
    ENDDO
    CALL polintNR(x2a,yntmp,n,x2,ymtmp(j),dy,check)
  ENDDO

  CALL polintNR(x1a,ymtmp,m,x1,y,dy,check)

  RETURN

END SUBROUTINE polin2


SUBROUTINE polintNR(xa,ya,n,x,y,dy,check)

  IMPLICIT NONE
  INTEGER n,NMAX
  REAL dy,x,y,xa(n),ya(n)
  PARAMETER (NMAX=10)
  INTEGER i,m,ns,check
  REAL den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)

  ns=1
  dif=ABS(x-xa(1))
  DO i=1,n
    dift=ABS(x-xa(i))
    IF (dift.LT.dif) THEN
      ns=i
      dif=dift
    ENDIF
    c(i)=ya(i)
    d(i)=ya(i)
  ENDDO
  y=ya(ns)
  ns=ns-1
  DO m=1,n-1
    DO i=1,n-m
      ho=xa(i)-x
      hp=xa(i+m)-x
      w=c(i+1)-d(i)
      den=ho-hp
      IF(den.EQ.0.) THEN
         !pause ’failure in polint’
         check=-1
         RETURN
      END IF
      den=w/den
      d(i)=hp*den
      c(i)=ho*den
    ENDDO
    IF (2*ns.LT.n-m)THEN
      dy=c(ns+1)
    ELSE
      dy=d(ns)
      ns=ns-1
    ENDIF
    y=y+dy
  ENDDO
  RETURN

END SUBROUTINE polintNR
