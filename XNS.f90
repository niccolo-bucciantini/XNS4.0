PROGRAM XNS
  !   ============================================================
  !   Purpose :  This progam computes a NS model for a given quantity of interest
  !              (central density, baryonic mass, gravitational mass, etc...)
  !
  !   It invokes XNSMAIN that computes a NS model for a given central density
  !   This is done repeatedly with Newton scheme where the solution having the
  !   desired property is searched for.
  !
  !  Defaults (can be changed by user) - see relevant line.
  !
  !   RHOMIN & RHOMAX are set to
  !         - For Newton-Raphson (1 +/- 0.2 )*RHOINI
  !         - For sequences 0.9E-3 - 2.0E-3
  !
  !   For sequences:
  !         The values of KBMIN,KBMAX, OMGMIN, OMGMAX, RHOMIN, RHOMAX
  !         should be changed by the User - see the relevant lines
  !         Number of models: 
  !         - Minimum number of magnetized sequences = 3 + 1
  !         - Minimum number of density points = 99 + 1
  !         - Minimum number of rotational sequences = 2 + 1
  !   ============================================================
  
  USE SYSTEMXNS
  USE PHYSICS
  USE ROTATION, ONLY: CHECKROTDIFF
  IMPLICIT NONE

  ! #################################################################################
  ! Precompiled for parspace - computation of sequences using multiple CPUs 
#ifdef MPIV
  include 'mpif.h'
  ! use mpi
  ! include '/usr/local/include/mpif.h'
  ! include '/usr/include/openmpi-x86_64/mpif.h'
  ! include '/usr/include/mpich-x86_64/mpif.h'

  INTEGER :: comm = MPI_COMM_WORLD
  INTEGER :: mpierr, mpireal
  INTEGER :: status(MPI_STATUS_SIZE)
  
  REAL, DIMENSION(3) :: wsend, wrecv
  REAL :: omgmaxseq, omgmin, domg, nomg      ! Range of Omega for sequences
  REAL :: kbmax, kbmin, dkb, nkb          ! Range of magnetizations for sequences
  REAL :: nrho1, nrho2               
  REAL :: rhosend, kbpsend, omgsend            
  INTEGER :: ending = 1, ACTPROC
  
  CHARACTER(len=10) :: name1, name2, name3
  CHARACTER(len=200) :: command, textadd
  logical,allocatable,dimension(:) :: vlogical
  
  DOUBLE PRECISION TSTART,TEND
#endif 
  !  Used in the Newton-Raphson scheme
  REAL :: DF, F1, F2, FR, FL
  REAL :: DRHO,DRHOOLD
  REAL :: QUA,QUB,QUR,QUL
  ! Density and models variables
  REAL :: RHOVAR, RHOMAX, RHOMIN
  INTEGER :: NRHOINI
  INTEGER :: NMODELS
  INTEGER :: ERR
  ! Timing variables
  INTEGER :: ICOUNT,ICOUNT0,IRATE,J,ILOOP
  REAL :: SECS
  INTEGER ::  num1, num2
  INTEGER ::  myid, nprocs, id2send

  ! #########################################################
  
  ! Check for consistency of input parameters
  IF((ITOR.AND.IPOL) .OR. (ITOR.AND.ITWT) .OR. (IPOL.AND.ITWT))THEN
    WRITE(6,*)'ERROR: logical parameters for the magnetic structure are conflicting, only one can be set .TRUE.'
    STOP
  END IF
  IF(IPOL.AND.DIFFERENTIAL) THEN
    WRITE(6,*)'ERROR: differential rotation with poloidal magnetic not implemented.'
    STOP 
  ENDIF
  IF(ITWT.AND.OMG.NE.0) THEN
    WRITE(6,*)'ERROR: rotation together with magnetic mixed field not implemented.'
    STOP  
  ENDIF
  IF(CUT.LT.1) THEN
     WRITE(6,*)'WARNING: the parameter -CUT- should not be less than 1'
  ENDIF
  IF(IMAG .AND. ((.NOT. ITOR) .AND. (.NOT. IPOL) .AND. (.NOT. ITWT)))THEN
     WRITE(6,*)'ERROR: the magnetic field configuration has not been specified'
     STOP
  ENDIF
  IF((.NOT. IMAG))THEN
     ITOR=.FALSE.
     IPOL=.FALSE.
     ITWT=.FALSE.
  ENDIF
  IF(GR .AND. (ALPHA0 .NE. 0. .OR. BETA0 .NE. 0.))THEN
     WRITE(6,*)'ERROR: setup for GR but STT parameter are not 0'
     STOP
  END IF
  IF(DIFFERENTIAL) CALL CHECKROTDIFF
  
  ! Start counting system time for log file
  CALL SYSTEM_CLOCK(ICOUNT0)

  ! #################################################################################
  ! Compile if you want to run XNS to converge to a given quantity
#ifdef NWTRPS
  
  IF(QUOC.LT.0 .OR. QUOC .GT. 2) THEN
     WRITE(6,*)' WARNING: please select a quantity for the convergence'
     WRITE(6,*)' of the Newton-Raphson scheme:'
     WRITE(6,*)' QOFC = 0 -> central density'
     WRITE(6,*)' QOFC = 1 -> gravitational mass'
     WRITE(6,*)' QOFC = 2 -> baryonic mass'
  ENDIF
  
  IDAT=.FALSE.
  SUBDIRPATH=''

  ! Inital density for the search
  RHOVAR = RHOINI
  ! Range of density over which to look for solutions (+/- 20%)
  RHOMAX = RHOVAR + 0.2*RHOVAR 
  RHOMIN = RHOVAR - 0.2*RHOVAR
  
  !Check if density range is OK
  WRITE(6,*)'Checking the initial central density range...'
  
  CALL XNSMAIN(RHOMAX,QUR)
  CALL XNSMAIN(RHOMIN,QUL)
  
  FL = QUCONV-QUL
  FR = QUCONV-QUR
  
  IF(FL*FR .GT. 0)THEN
     IF(QUOC .EQ. 0)THEN
        WRITE(6,*)'ERROR: solution is outside the initial  central density range ',RHOMIN,' - ',RHOMAX
        WRITE(6,*)'Selectc a new initial central density or a larger range'
        WRITE(6,*)'At central density = ',RHOMIN,' the quantity is = ',QUL 
        WRITE(6,*)'At central density = ',RHOMAX,' the quantity is = ',QUR 
        STOP
     ELSE
        WRITE(6,*)'ERROR: solution might be a local maxima/minima in the density range ',RHOMIN,' - ',RHOMAX
        WRITE(6,*)'Selectc another variable for convergence'
        WRITE(6,*)'At central density = ',RHOMIN,' the quantity is = ',QUL 
        WRITE(6,*)'At central density = ',RHOMAX,' the quantity is = ',QUR
        STOP
     END IF
  END IF
  
  !Compute the quantity of interest for the initial guess density
  CALL XNSMAIN(RHOVAR,QUA)
  
  F1 = QUCONV-QUA
  J = 1
  WRITE(6,*)'Iter = ',J,', Quantity  =', QUA

  !Compute the numerical derivative with respect to the initial guess density
  CALL XNSMAIN(RHOVAR*1.0001,QUB)
  F2 = QUCONV-QUB
  DF = (F2 - F1)/(RHOVAR*0.0001)

  DRHO = RHOMAX-RHOMIN
  DRHOOLD = DRHO

  DO J = 2, NVALUE
     ! Remain within a bounded search domain
     IF((((RHOVAR-RHOMIN)*DF-F1) * ((RHOVAR-RHOMAX)*DF-F1) .GT. 0) &
        .OR. (ABS(2.*F1) .GT. ABS(DRHOOLD*DF))) THEN
        DRHOOLD = DRHO
        DRHO = 0.5*(RHOMAX - RHOMIN)
        RHOVAR = RHOMIN + 0.5*DRHO
     ELSE 
        DRHOOLD = DRHO
        DRHO= F1/DF
        RHOVAR=RHOVAR-DRHO
     END if

     ! Compute new guess
     CALL XNSMAIN(RHOVAR,QUA)   
     F1 = QUCONV-QUA
     WRITE(6,*)'Iter = ',j,', Quantity =', QUA,QUCONV,F1

     ! Check if conevergent to the desired accuracy
     IF(ABS(F1)/QUCONV .LT. CONVF)THEN
        IDAT=.TRUE.
        CALL XNSMAIN(RHOVAR,QUA)
        WRITE(6,*)''
        WRITE(6,*)'End of Run - Num. iterations = ',J
        WRITE(6,*)'Desired value = ',QUCONV
        WRITE(6,*)'Obtained value = ',QUA
        WRITE(6,*)'Corresponding central density in XNSMAIN = ',RHOVAR
        EXIT
     END IF
     
     ! Compute the numerical derivative 
     CALL XNSMAIN(RHOVAR*1.0001,QUB)

     F2 = QUCONV-QUB
     DF= (F2 - F1)/(RHOVAR*0.0001)
     
     ! Bound the domain of search
     IF(F1 .GT. 0) RHOMIN = RHOVAR
     IF(F1 .LT. 0) RHOMAX = RHOVAR
  END DO

  ! #################################################################################  
  ! Compile for parspace - computation of sequences using multiple CPUs  
#elif MPIV
  
  TSTART = MPI_WTIME()
  IDAT= .TRUE.
  ! Initialize MPI
  CALL MPI_INIT(mpierr)
  CALL MPI_COMM_SIZE(comm,nprocs,mpierr)
  CALL MPI_COMM_RANK(comm,myid,mpierr)
  
  if (kind(1.0)==4) mpireal=mpi_real
  if (kind(1.0)==8) mpireal=mpi_double_precision

  ! Change this values to change the range of sequences
  ! Set the range of sequences for central initial density
  RHOMIN = 0.9E-3
  RHOMAX = 2.0E-3
  ! Set the range of sequences for uniform rotation rate
  OMGMIN = 0.
  OMGMAXSEQ = 0.
  ! Set the range of sequences for magnetization parameter
  IF(IPOL)THEN ! Poloidal
     KBMIN = 0.0
     KBMAX = 0.2
  ELSE IF(ITOR)THEN ! Toroidal
     KBMIN = 0.0
     KBMAX = 0.7
  ELSE IF(ITWT)THEN ! Mixed
     KBMIN = 3.4E-4
     KBMAX = 3.4E-4
  ELSE ! Unmagnetized
     KBMIN = 0.0
     KBMAX = 0.0
  END IF
 
  ! Array to check which CPU has been closed 
  ALLOCATE(VLOGICAL(1:NPROCS-1))
  VLOGICAL = .FALSE.
  
  ! Check the number of magnetic sequences
  IF(KBMIN .EQ. KBMAX)THEN
     NKB = 1.
     NMODELS = 1.
  ELSE ! Minumum number of magnetized sequences = 3 + 1
     NKB= 4.
     NMODELS = NKB+1.
  ENDIF
  
  ! Check the number of density points
  IF(RHOMIN .EQ. RHOMAX)THEN
     NRHO1 = 1.
  ELSE ! Minumum number of density points = 99 + 1
     NRHO1= 19.
     NMODELS = NMODELS*(NRHO1+1.)
  ENDIF
  
  ! Check the number of rotation rates
  IF(OMGMIN .EQ. OMGMAXSEQ)THEN
     NOMG = 1.
  ELSE ! Minumum number of rotation rates = 2 + 1
     NOMG = 2.
     NMODELS = NMODELS*(NOMG + 1.)
  ENDIF
  
  ! Set the density increment
  IF(NRHO1 .GT. 1)THEN
     DRHO=(RHOMAX-RHOMIN)/NRHO1
  ELSE
     DRHO=0.
  ENDIF
  ! Set the magnetic increment
  DKB =(KBMAX-KBMIN)/NKB
  ! Set the rotation increment
  DOMG=(OMGMAXSEQ-OMGMIN)/NOMG
  
  ! Set initial density, magnetization and rotation rate of models to send
  RHOSEND=RHOMIN-DRHO
  IF(NRHO1 .GT. 1)THEN
     KBPSEND=KBMIN
  ELSE
     KBPSEND=KBMIN-DKB
  ENDIF
  OMGSEND=OMGMIN
  
  ! Send the first chunck of initial conditions to the other tasks
  IF(MYID.EQ.0)THEN
     IF(NMODELS .LT. NPROCS-1)THEN
        WRITE(6,*)'ERROR: number of computing processes is greater than the number of models,',nprocs-1,'vs',nmodels
        STOP
     ENDIF
     IF(NMODELS .EQ. 1)THEN
        WRITE(6,*)'ERROR: to compute only one model, please use the serial code'
        STOP
     ENDIF
     ! Cycle over the various CPUs
     DO J=1,NPROCS-1
        ! Set density to send
        RHOSEND=RHOSEND+DRHO
        IF(RHOSEND .GT. (RHOMAX+1.E-6) .OR. (DRHO .EQ. 0))THEN
           RHOSEND=RHOMIN
           ! Set magnetization to send
           KBPSEND=KBPSEND+DKB
           IF((KBPSEND .GT. (KBMAX+1.E-3)) .OR. (DKB .EQ. 0)) THEN
              KBPSEND=KBMIN
              ! Set rotation to send
              OMGSEND=OMGSEND+DOMG
              ! If all range has been sent
              IF(OMGSEND.gt.OMGMAXSEQ .OR. (DOMG .EQ. 0.0))THEN 
                 WRITE(6,*)'All models sent.'
                 RHOSEND=-1.
                 ENDING=0  ! Set parameter for and of run
                 IF(VERBOSE)PRINT*,'ENDING = ',ENDING
              END IF
           END IF
        END IF
        ! Set parameters for each model
        WSEND=[RHOSEND,KBPSEND,OMGSEND]
        
        IF(ENDING .EQ. 1)THEN ! Tell CPU to be active
           IF(VERBOSE)PRINT*,'SEND ENDING TO ',j,ending
           CALL MPI_SEND(ENDING,1,MPI_INTEGER,J,11,COMM,MPIERR)
           WRITE(6,*)'Sending to', j,':', wsend
           ! Send parameters to CPU 
           CALL MPI_SEND(wsend,3,mpireal,j,10,comm,mpierr)! (1) process 0 (root) distributes the initial conditions between the processes  
        END IF
     END DO
     
     ! Iterate till completion
     IF(ENDING .EQ. 1)THEN 
        ACTPROC=NPROCS-1
        DO
           ! Set density to send
           RHOSEND=RHOSEND+DRHO
           IF(RHOSEND.gt.(RHOMAX+1.E-6) .OR. (DRHO .EQ. 0)) THEN
              RHOSEND=RHOMIN
              ! Set magnetization to send
              KBPSEND=KBPSEND+DKB
              IF((KBPSEND.gt. (KBMAX+1.E-3)) .OR. (DKB .EQ. 0)) THEN
                 KBPSEND=KBMIN
                 ! Set rotation to send
                 OMGSEND=OMGSEND+DOMG
                 ! If all range has been sent
                 IF((OMGSEND.gt.OMGMAXSEQ .OR. (DOMG .EQ. 0.0)))then
                    print*,'All models sent.'
                    RHOSEND=-1.
                    ENDING=0  ! Set parameter for and of run
                    if(VERBOSE)print*,'ENDING = ',ENDING
                 end if
              endif
           endif
           ! Set parameters for each model
           wsend=[RHOSEND,KBPSEND,OMGSEND]
           
           ! Receive id of the free task and send new data
           call MPI_RECV(id2send,1,mpi_integer,mpi_any_source,20,comm,status,mpierr)	! (4) root process receives the id of the idling CPU
           if(VERBOSE)print*,'RECV 21: ',id2send,'IS IDLING; WAITING FOR ENDING',ending
           ! Tell CPU to be active or not in case
           CALL MPI_SEND(ENDING,1,MPI_INTEGER,id2send,21,COMM,MPIERR)
           IF(ENDING .EQ. 1)THEN ! Active
              WRITE(6,*), 'More sending to', id2send,':', wsend
              CALL MPI_SEND(wsend,3,mpireal,id2send,30,COMM,MPIERR) ! (5b) root sends new initial data to idling
           END IF
           ! No more model to compute
           IF(ENDING .EQ. 0)THEN
              VLOGICAL(id2send) = .true.
              ACTPROC=ACTPROC-1
              WRITE(6,*)'Number of active processes: ',ACTPROC
              CALL DATE_AND_TIME(DATE,TIME,ZONE,VALUES) 
              WRITE(6,*)'Date-time: ',DATE,'-',TIME
           END IF
           ! If all CPUs are closed then exit  
           IF(ALL(vlogical))EXIT
        END DO
     END IF
     ! If all CPUs are closed then the master can exit 
     if(all(vlogical))THEN
        print*,'Process 0 is finalizing.'
        call mpi_finalize(mpierr)
     endif
     
     ! System clock
     CALL SYSTEM_CLOCK(ICOUNT,IRATE)
     SECS=FLOAT(ICOUNT-ICOUNT0)/FLOAT(IRATE)
     WRITE(6,*)'Elapsed CPU time is ',SECS,'seconds'
     
  ELSE ! if other CPUs
     
     !  Receive first chunk of initial conditions
     CALL MPI_RECV(ENDING,1,MPI_INTEGER,0,11,COMM,STATUS,MPIERR)
     if(VERBOSE)PRINT*,'RECV 11: ',MYID,'ENDING', ending
     
     IF(ENDING .EQ. 1)THEN
        if(VERBOSE)PRINT*,'RECV 10: ',MYID,'WAITING DATA'
        call MPI_RECV(wrecv,3,mpireal,0,10,comm,status,mpierr)	  ! (2) the other processes receive the initial conditions and start 
        if(VERBOSE)PRINT*,'RECV 10: ',MYID,'RECEIVED DATA',wrecv
200     RHOVAR=wrecv(1)      					   ! to compute the solutions
        if(ITOR) BCOEF = wrecv(2)
        if(IPOL) KBPOL = wrecv(2)
        if(ITWT) KBTT  = wrecv(2)
        OMG=wrecv(3)
        
        ! Make folder to write data
        WRITE(name3,'(ES10.4)') wrecv(1)
        WRITE(name2,'(ES10.3)') wrecv(2)
        WRITE(name1,'(ES10.3)') wrecv(3)
        SUBDIRPATH='omg'//trim(adjustl(name1))//'/mag'//trim(adjustl(name2))//'/rho'//trim(adjustl(name3))
        command='mkdir -p '//SUBDIRPATH
        SUBDIRPATH='omg'//trim(adjustl(name1))//'/mag'//trim(adjustl(name2))//'/rho'//trim(adjustl(name3))//'/'
        CALL SYSTEM(command)
        IF(WRECV(2) .EQ. 0.)THEN
           IMAG=.FALSE.
        ELSE
           IMAG=.TRUE.
        ENDIF
        CALL XNSMAIN(RHOVAR,QUA)
        ! Check what is going on when CPU ends
        IF(ENDID)THEN
           WRITE(6,*)WRECV(1),'WARNING:',myid,'is diverging. Lowering qfactors to ',QFACTOR*0.9,'.'
           QFACTOR=QFACTOR*0.9
           QFACTORCHI=QFACTORCHI*0.9
           IF(QFACTOR .LT. 0.3)THEN
              WRITE(6,*)WRECV(1),'WARNING: could not converge. Ending process',myid,'.'
              ENDING=0
           ELSE ! Try again with new QFACTOR
              GOTO 200
           ENDIF
        ENDIF
        ! Gzip datafiles
        command='gzip -f '//trim(SUBDIRPATH)//'Hydroeq.dat'
        CALL SYSTEM(command)
        IF(IMAG)THEN
           command='gzip -f '//trim(SUBDIRPATH)//'Hydroeq_mag.dat'
           CALL SYSTEM(command)
        ENDIF
        
        ! Tell master CPU is free 
        if(VERBOSE)PRINT*,'SEND 20: ',MYID,'FREE TO RECEIVE'
        CALL MPI_SEND(myid,1,mpi_integer,0,20,comm,mpierr)	
     ENDIF
     ! If not more models then CPU can close
     if(ENDING .EQ.0)then
        CALL MPI_SEND(myid,1,mpi_integer,0,20,comm,mpierr)		
        if(VERBOSE)print*,'RECEIVED ENDING SIGNAL FROM PROCESS 0: ',MYID,'IS EXITING'
     endif
     
     do  ! (3) every non-root process sends its id to root
        
        CALL MPI_RECV(ENDING,1,MPI_INTEGER,0,21,COMM,STATUS,MPIERR)
        if(VERBOSE)PRINT*,'RECV 11: ',MYID,'ENDING', ending
        IF(ENDING .EQ. 1)THEN
           if(VERBOSE)PRINT*,'RECV 30: ',MYID,'WAITING DATA'	
           call MPI_RECV(wrecv,3,mpireal,0,30,comm,status,mpierr)! (5a) every non-root receives new initial data
           if(VERBOSE)print*,'RECV 30: ',MYID,'IS RECEVING MORE DATA',wrecv
100        RHOVAR=wrecv(1)
           if(RHOVAR.lt.0.)then
              if(VERBOSE)print*,'RECEIVED SIGNAL FROM PROCESS 0: ',MYID,'IS EXITING'
              exit
           endif
           if(ITOR) BCOEF=wrecv(2)
           if(IPOL) KBPOL=wrecv(2)
           if(ITWT) KBTT=wrecv(2)
           OMG=wrecv(3)
           
           ! Make folder to write data
           WRITE(name3,'(ES10.4)') wrecv(1)
           WRITE(name2,'(ES10.3)') wrecv(2)
           WRITE(name1,'(ES10.3)') wrecv(3)
           SUBDIRPATH='omg'//trim(adjustl(name1))//'/mag'//trim(adjustl(name2))//'/rho'//trim(adjustl(name3))
           command='mkdir -p '//SUBDIRPATH
           SUBDIRPATH='omg'//trim(adjustl(name1))//'/mag'//trim(adjustl(name2))//'/rho'//trim(adjustl(name3))//'/'
           CALL SYSTEM(command)
           IF(WRECV(2) .EQ. 0.)THEN
              IMAG=.FALSE.
           ELSE
              IMAG=.TRUE.
           ENDIF
           CALL XNSMAIN(RHOVAR,QUA)
           ! Check what is going on when CPU ends
           IF(ENDID)THEN
              WRITE(6,*)WRECV(1),'WARNING:',myid,'is diverging. Lowering qfactors to ',QFACTOR*0.9,'.'
              QFACTOR=QFACTOR*0.9
              QFACTORCHI=QFACTORCHI*0.9
              IF(QFACTOR .LT. 0.3)THEN
                 WRITE(6,*)WRECV(1),'WARNING: could not converge. Ending process',myid,'.'
                 ENDING=0
              ELSE ! Try again with new QFACTOR		    
                 GOTO 100
              ENDIF
           ENDIF
           ! Gzip datafiles
           command='gzip -f '//trim(SUBDIRPATH)//'Hydroeq.dat'
           CALL SYSTEM(command)
           IF(IMAG)THEN
              command='gzip -f '//trim(SUBDIRPATH)//'Hydroeq_mag.dat'
              CALL SYSTEM(command)
           ENDIF
           
           ! Tell master CPU is free 
           if(VERBOSE)PRINT*,'SEND 20: ',MYID,'FREE TO RECEIVE'
           CALL MPI_SEND(myid,1,mpi_integer,0,20,comm,mpierr)	
        END IF
        
        ! If no more models then CPU can close
        if(ENDING .EQ.0.)then
           if(VERBOSE)print*,'RECEIVED ENDING SIGNAL FROM PROCESS 0: ',MYID,'IS EXITING'
           exit
        endif
        
     enddo
     
     ! Master closes
     print*,'Process',myid,'is finalizing.'
     call mpi_finalize(mpierr)
     
  endif
  
  TEND = MPI_WTIME()

  ! #################################################################################
  ! Serial run for single model
#else 
  
  IDAT = .TRUE.
  SUBDIRPATH = ''
  RHOVAR = RHOINI
  CALL XNSMAIN(RHOVAR,QUA)
  
  ! System clock
  CALL SYSTEM_CLOCK(ICOUNT,IRATE)
  SECS=FLOAT(ICOUNT-ICOUNT0)/FLOAT(IRATE)
  WRITE(6,*)'Elapsed CPU time is ',SECS,'seconds'
  
#endif

  ! #################################################################################
  
END PROGRAM XNS


