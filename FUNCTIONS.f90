MODULE FUNCTIONS
  !  ============================================================
  !  Purpose : storage module for  generic functions of various use
  !  ============================================================
  
  USE SYSTEMXNS
  USE PHYSICS
  IMPLICIT NONE

CONTAINS
  
! ***************************************************************
! ***************************************************************

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
  
  DO I = 1, N-2
     IF(ABS(D(I)) .GE. ABS(DL(I))) THEN
        ! No row interchange required
        IF(D(I) .NE. 0.) THEN
           FACT = DL(I) / D(I)
           D(I+1) = D(I+1) - FACT*DU(I)
           B(I+1) = B(I+1) - FACT*B(I)
        ELSE
           INFO = I
           RETURN
        END IF
        DL(I) = 0.
     ELSE  ! Interchange rows I and I+1
        FACT = D(I) / DL(I)
        D(I) = DL(I)
        TEMP = D(I+1)
        D(I+1) = DU(I) - FACT*TEMP
        DL(I) = DU(I+1)
        DU(I+1) = -FACT*DL(I)
        DU(I) = TEMP
        TEMP = B(I)
        B(I) = B(I+1)
        B(I+1) = TEMP - FACT*B(I+1)
     END IF
  END DO
  
  IF( N.GT.1 ) THEN
     I = N-1
     IF(ABS(D(I)) .GE. ABS(DL(I))) THEN
        IF(D(I) .NE. 0.) THEN
           FACT = DL(I) / D(I)
           D(I+1) = D(I+1) - FACT*DU(I)
           B(I+1) = B(I+1) - FACT*B(I)
        ELSE
           INFO = I
           RETURN
        END IF
     ELSE
        FACT = D(I) / DL(I)
        D(I) = DL(I)
        TEMP = D(I+1)
        D(I+1) = DU(I) - FACT*TEMP
        DU(I) = TEMP
        TEMP = B(I)
        B(I) = B(I+1)
        B(I+1) = TEMP - FACT*B(I+1)
     END IF
  END IF
  IF(D(N) .EQ. 0.) THEN
     INFO = N
     RETURN
  END IF
  
  ! Back solve with the matrix U from the factorization.
  B(N) = B(N) / D(N)
  IF(N .GT. 1)THEN
     B(N-1) = (B(N-1) - DU(N-1)*B(N)) / D(N-1)
  END IF
  DO I = N-2, 1, -1
     B(I) = (B(I) - DU(I)*B(I+1) - DL(I)* &
          B(I+2)) / D(I)
  END DO
  
  RETURN
  
END SUBROUTINE DGTSV
  
! ***************************************************************
! ***************************************************************
  
SUBROUTINE LUSOLVER(N,A,B,IPIV,JUMP)
  !  =============================================================
  !  Purpose: Solve the linear system A*X = B
  !  where A is an N x N matrix and B is a vector of length N
  !  Use LU decomposition with partial pivoting.
  !  Adapted from LAPACK 3.2 = simplified for this purpose.
  !
  !  N = dimension of A(N,N) and B(N)
  !  A = in input matrix of coefficients, in output LU decomposition
  !  B = in input vector of RHS values, in output vector of solutions X
  !  IPIV = array of pivoting positions
  !  JUMP = Logical argument True => performs only the last LU solution step
  !  SFMIN = Lower floor for singular pivoting
  !
  !  Uses: MYISAMAX, MYSSWAP, SGER, SLASWP, STRSM
  !  =============================================================
  
  IMPLICIT NONE
  LOGICAL :: JUMP
  REAL,PARAMETER :: SFMIN =1.D-20
 ! INTEGER :: MYISAMAX
  INTEGER :: N,J,JP
  INTEGER :: IPIV(1:N)
  REAL :: A(1:N,1:N)
  REAL :: B(1:N)

  IF(.NOT.JUMP)THEN ! Do the LU decomposition
     ! Search for pivoting position
     DO J=1,N
        JP = J-1 + MYISAMAX(N-J+1,A(J:N,J))
        IPIV(J)= JP
        IF (A(JP,J) .NE. 0)THEN
           IF(JP.NE.J ) THEN
              CALL MYSSWAP(N,A(J,1:N),A(JP,1:N))
           END IF
           ! Check for singular points 
           IF(J .LT. N)THEN
              IF(ABS(A(J,J)) .LE. SFMIN) A(J,J) = SFMIN
              A(J+1:N,J) = A(J+1:N,J) / A(J,J) 
           END IF
        ELSE
           WRITE(6,*)'Error in LUSOLVER',A(JP,J),J,JP
        END IF
        
        IF(J.LT.N) THEN
           CALL SGER(N,N-J,A(J+1,J),A(J,J+1),A(J+1,J+1))
        END IF
        
     END DO
  END IF
  ! Once the LU decompostion is available proceed to solution
  ! Perform the exchange over B
  CALL SLASWP(N, B, IPIV)
  CALL STRSM(.False., .False., N, A, B)
  CALL STRSM(.True. , .True. , N, A, B)
  
END SUBROUTINE LUSOLVER

! ***************************************************************
! ***************************************************************

SUBROUTINE MYSSWAP(N,SX,SY)
  !  =============================================================
  !  Purpose: exchange two arrays
  !  Adapted from LAPACK 3.2 = simplified for this purpose.
  !
  !  N = dimension array SX(N) and SY(N)
  !  SX = in input first array, in output second array
  !  SY = in input second array, in output first array
  !  =============================================================
  
  IMPLICIT NONE
  INTEGER :: N
  REAL :: SX(N),SY(N)
  REAL STEMP(N)

  IF (N.LE.0) RETURN
  
  STEMP(1:N) = SX(1:N)
  SX(1:N) = SY(1:N)
  SY(1:N) = STEMP(1:N)

  RETURN
END SUBROUTINE MYSSWAP

! ***************************************************************
! ***************************************************************

INTEGER FUNCTION MYISAMAX(N,SX)
  !  =============================================================
  !  Purpose: return the index of maximum value of an array
  !  Adapted from LAPACK 3.2 = simplified for this purpose.
  !
  !  N = dimension array SX(N)
  !  SX = array to scan for maximum 
  !  =============================================================
  
  IMPLICIT NONE
  INTEGER :: N,I
  REAL :: SMAX
  REAL :: SX(N)
  
  MYISAMAX = 0
  IF (N.LT.1) RETURN
  MYISAMAX = 1
  IF (N.EQ.1) RETURN
  
  SMAX = ABS(SX(1))
  DO I = 2,N
     IF (ABS(SX(I)).GT.SMAX)THEN
        MYISAMAX = I
        SMAX = ABS(SX(I))
     END IF
  END DO
  
  RETURN
END FUNCTION MYISAMAX

! ***************************************************************
! ***************************************************************

SUBROUTINE SGER(NN,N,X,Y,A)
  !  =============================================================
  !  Purpose:  SGER   performs the rank 1 operation A := x*y**T + A,
  !  where  x is an nn element vector, y is an n element vector
  !   and A is an nn by n matrix.
  !
  !  Adapted from LAPACK 3.2 = simplified for this purpose.
  !  =============================================================
  
  IMPLICIT NONE
  INTEGER,PARAMETER:: INCX=1
  INTEGER INCY,N,NN
  INTEGER I,J,JY
  REAL A(NN,*),X(*),Y(*) 
  REAL TEMP

  INCY = NN
  
  IF (INCY.GT.0) THEN
     JY = 1
  ELSE
     JY = 1 - (N-1)*INCY
  END IF
  
  DO  J = 1,N
     IF (Y(JY).NE.0.D0) THEN
        TEMP = -1.0*Y(JY)
        DO I = 1,N
           A(I,J) = A(I,J) + X(I)*TEMP
        END DO
     END IF
     JY = JY + INCY
  END DO
  
  RETURN
END SUBROUTINE SGER

! ***************************************************************
! ***************************************************************

SUBROUTINE SLASWP(N, A, IPIV)
  !  =============================================================
  !  Purpose: exchange elements of an array
  !  Adapted from LAPACK 3.2 = simplified for this purpose.
  !
  !  N = dimension array A(N) and IPIV(N)
  !  A = in input original array, in output array with element echanged
  !  IPIV(N) = element at position N will be moved to position IPIV(N)
  !  =============================================================
  
  IMPLICIT NONE
  INTEGER :: N
  INTEGER :: I,IP
  INTEGER :: IPIV(N)
  REAL :: A(N)
  REAL :: TEMP
  
  DO I = 1,N
     IP = IPIV(I)
     IF(IP.NE.I) THEN
        TEMP = A(I)
        A(I) = A(IP)
        A(IP) = TEMP
     END IF
  END DO
  
  RETURN
END SUBROUTINE SLASWP

! ***************************************************************
! ***************************************************************

SUBROUTINE STRSM(UPPER,NOUNIT,N,A,B)
  !  =============================================================
  !  Purpose: Solve A*X = B where A is a matrix in LU form
  !  Adapted from LAPACK 3.2 = simplified for this purpose.
  !
  !  UPPER = Logical (true = upper diagonal, false = lower diagonal)
  !  NOUNIT = Logical (true unitary matrix)
  !  N = size of matrix A(N,N) and array (B)
  !  A = in input matrix from LU decompostion
  !  B = in input RHS of system, in output solution
  !  =============================================================
  
  IMPLICIT NONE
  INTEGER :: N
  REAL :: A(N,N),B(N)
  LOGICAL :: NOUNIT,UPPER
  INTEGER :: I,J,K,INFO

  ! Solve for upper diagonal matrix
  IF (UPPER) THEN
     DO K = N,1,-1
        IF (B(K).NE.0.D0) THEN
           IF (NOUNIT) B(K) = B(K)/A(K,K)
           DO I = 1,K-1
              B(I) = B(I) - B(K)*A(I,K)
           END DO
        END IF
     END DO
  ELSE
  ! Solve for lower diagonal matrix   
     DO K = 1,N
        IF (B(K).NE.0.D0) THEN
           IF (NOUNIT) B(K) = B(K)/A(K,K)
           DO I = K+1,N
              B(I) = B(I) - B(K)*A(I,K)
           END DO
        END IF
     END DO
  END IF
  
  RETURN
END SUBROUTINE STRSM

! ***************************************************************
! ***************************************************************

SUBROUTINE LAGRANGEINT(N,XP,YP,XINT,YINT,DY)
  !  =============================================================
  !  Purpose: interpolate YP(XP) in XINT using (N-1)th-Polinomial
  !  according to Lagrange Interpolation with Neville algorithm
  !  Adapted from NR = simplified for this purpose.
  !
  !  N = size of array XP(N) and YP(N)
  !  XP = ascissa of input points
  !  YP = ordinata of input point
  !  XINT = value where interpolation is needed
  !  YINT = interpolated value
  !  =============================================================

  IMPLICIT NONE
  INTEGER :: N
  INTEGER :: I,IS,M
  REAL,PARAMETER :: LARGE = 1.D30
  REAL :: XP(N),YP(N)
  REAL :: C(N),D(N)
  REAL :: DELTA
  REAL :: XINT,YINT,DY
  REAL :: CMD

  ! Locate position IS of element of XP closest to XINT
  DELTA = LARGE
  DO I=1,N
     IF (ABS(XINT-XP(I)).LT.DELTA) THEN
        IS = I
        DELTA = ABS(XINT-XP(I))
     ENDIF
     C(I)=YP(I)
     D(I)=YP(I)
  ENDDO
  ! Approximante YINT as XP(IS)
  YINT = YP(IS)

  ! Perform Neville Algorithm
  IS = IS-1
  DO M = 1,N-1
    DO I = 1,N-M
      CMD = C(I+1) - D(I)
      D(I) = (XP(I+M)-XINT)*CMD/(XP(I) - XP(I+M))
      C(I) = (XP(I  )-XINT)*CMD/(XP(I) - XP(I+M))
    ENDDO
    IF (2*IS.LT.N-M)THEN
      DY = C(IS+1)
    ELSE
      DY = D(IS)
      IS = IS-1
    ENDIF
    YINT = YINT + DY
 ENDDO
 
  RETURN
END SUBROUTINE LAGRANGEINT

! ***************************************************************
! ***************************************************************

SUBROUTINE LAGRANGEINT2D(M,N,XP1,XP2,YP,X1,X2,YINT,DY)
  !  =============================================================
  !  Purpose: perform 2D interpolation using (N-1)th-Polinomial
  !  according to Lagrange Interpolation with Neville algorithm,
  !  on structured grid
  !
  !  M = size of 1st coordinates array XP1(N) 
  !  N = size of 2nd coordinates array XP2(N) 
  !  XP1 = ascissa of 1st coordinates points
  !  XP2 = ascissa of 2st coordinates points
  !  YP = values of ordinata  YP(I,J) in point X1P(I) X2P(J)
  !  X1 = value of 1st coordinate where interpolation is needed
  !  X2 = value of 2st coordinate where interpolation is needed
  !  YINT = interpolated value
  !  =============================================================

  IMPLICIT NONE
  INTEGER :: M,N
  REAL :: XP1(M),XP2(N),YP(M,N)
  REAL :: X1,X2,Y,DY,YINT
  INTEGER :: J,K
  REAL :: YNTMP(N),YMTMP(M)
  
  DO J=1,M
    DO K=1,N
      YNTMP(K) = YP(J,K)
    ENDDO
    CALL LAGRANGEINT(N,XP2,YNTMP,X2,YMTMP(J),DY)
    !polintNR(x2a,yntmp,n,x2,ymtmp(j),dy,check)
  ENDDO

  CALL LAGRANGEINT(M,XP1,YMTMP,X1,Y,DY)
  !CALL polintNR(x1a,ymtmp,m,x1,y,dy,check)
  YINT = Y
  
  RETURN

END SUBROUTINE LAGRANGEINT2D

! ***************************************************************
! ***************************************************************

SUBROUTINE LEGZO(N,X,W)

  ! =========================================================
  ! Purpose : Compute the zeros of Legendre polynomial Pn(x)
  !           in the interval [1,-1], and the corresponding
  !           weighting coefficients for Gauss-Legendre
  !           integration
  ! Input :   n    --- Order of the Legendre polynomial
  ! Output:   X(n) --- Zeros of the Legendre polynomial
  !           W(n) --- Corresponding weighting coefficients
  !
  !  From: SPECIAL_FUNCTIONS is a FORTRAN77 library which computes
  !  the value of various special functions, by Shanjie Zhang, Jianming Jin.
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
  !
  !  From: SPECIAL_FUNCTIONS is a FORTRAN77 library which computes
  !  the value of various special functions, by Shanjie Zhang, Jianming Jin.
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
   PF=(2.0D0*K-1.0D0)/K*X*P1-(K-1.0D0)/K*P0	! Bonnet's recursion formula
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

    IF(D2L*D2R .GT. 0)THEN
      Y(I)=YA(JJ(I))+D1L*(X(I)-XA(JJ(I)))+0.5*D2L*(X(I)-XA(JJ(I)))**2.
      Y(I)=Y(I)+YA(JJ(I)+1)+D1R*(X(I)-XA(JJ(I)+1))+0.5*D2R*(X(I)-XA(JJ(I)+1))**2.
      Y(I)=0.5*Y(I)
    ELSE
      Y(I)=YA(JJ(I))+D1C*(X(I)-XA(JJ(I)))
    END IF
  END DO

  RETURN

END SUBROUTINE POLINT

! ***************************************************************************
! ***************************************************************************

SUBROUTINE XNS2ECHO_OUT()

  !  ===============================================
  !  Purpose:   Output file for X-ECHO
  !  ===============================================
  USE SYSTEMXNS
  
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

! ***************************************************************************
! **************************************************************************
  
END MODULE FUNCTIONS
