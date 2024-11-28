!+
MODULE ForsytheMalcolmMoler
! PURPOSE - A collection of Fortran procedures for mathematical computation
!   based on the procedures from the book "Computer Methods for Mathematical
!   Computations", by George E. Forsythe, Michael A. Malcolm, and 
!   Cleve B. Moler. Prentice-Hall, 1977.

!  COLLECTION AUTHORS - George E. Forsythe, Michael A. Malcolm, & 
!    Cleve B. Moler, Stanford University (1977)
!    Ralph L. Carmichael, Public Domain Aeronautical Software

! ORIGINAL AUTHORS - 
!  Decomp,Solve - Forsythe,Malcolm,Moler
!  Spline,Seval - Forsythe,Malcolm,Moler
!  Quanc8 - Forsythe,Malcolm,Moler
!  RKF45 - H.A. Watts AND L.F. Shampine  (Sandia)
!  ZEROIN - Richard Brent
!  FMIN - Van Winjngaarden,Decker,Brent, Et Al
!  SVD - Golub and Reinsch

!  REVISION HISTORY
!   DATE  VERS PERSON STATEMENT OF CHANGES
!   1977   1.0   FMM  Original publication
! 07Jul83  1.1   RLC  Mods for Fortran77
! 09Jul89  1.2   RLC  Added SEVAL3, IMPLICIT NONE
! 17Apr91  1.3   RLC  Fortran 90 style loops
! 10May92  1.4   RLC  Fortran 90 coding
! 23Aug00  1.5   RLC  Made compatible with Elf90
! 26Apr02  1.6   RLC  Added NaturalSpline
! 03May02  1.7   RLC  Charged calling sequence of SVD
! 20May02  1.8   RLC  Added BrentZero - mod of ZeroIn
! 28May02  1.9   RLC  Added errCode to Decomp
! 28Jun02  1.91  RLC  Made cond OPTIONAL in Decomp
! 21Aug02  1.95  RLC  Made single and double precision versions 
! 11Oct02  1.96  RLC  Compute eps each time in RKFS

IMPLICIT NONE
  PUBLIC:: Decomp,Solve
  PUBLIC:: FMMspline,NaturalSpline,Seval,Seval3
  PUBLIC:: Quanc8
  PUBLIC:: Rkf45
  PRIVATE:: RkfsSingle,FehlSingle,RkfsDouble,FehlDouble
  PUBLIC:: Zeroin,BrentZero
  PUBLIC:: Fmin,BrentMin
  PUBLIC:: SVD
  
  INTERFACE Decomp
    MODULE PROCEDURE DecompSingle,DecompDouble
  END INTERFACE
  
  INTERFACE Solve
    MODULE PROCEDURE SolveSingle,SolveDouble
  END INTERFACE
  
  INTERFACE FMMspline
    MODULE PROCEDURE FMMsplineSingle,FMMsplineDouble
  END INTERFACE
  
  INTERFACE NaturalSpline
    MODULE PROCEDURE NaturalSplineSingle,NaturalSplineDouble
  END INTERFACE
  
  INTERFACE Seval
    MODULE PROCEDURE SevalSingle,SevalDouble
  END INTERFACE
  
  INTERFACE Seval3
    MODULE PROCEDURE Seval3Single,Seval3Double
  END INTERFACE
  
  INTERFACE Quanc8
    MODULE PROCEDURE Quanc8Single,Quanc8Double
  END INTERFACE

  INTERFACE RKF45
    MODULE PROCEDURE RKF45Single,RKF45Double
  END INTERFACE
  
  INTERFACE Zeroin
    MODULE PROCEDURE ZeroinSingle,ZeroinDouble
  END INTERFACE
  
  INTERFACE BrentZero
    MODULE PROCEDURE BrentZeroSingle,BrentZeroDouble
  END INTERFACE
  
  INTERFACE Fmin
    MODULE PROCEDURE FminSingle,FminDouble
  END INTERFACE
  
  INTERFACE BrentMin
    MODULE PROCEDURE BrentMinSingle,BrentMinDouble
  END INTERFACE
  
  INTERFACE SVD
    MODULE PROCEDURE SVDSingle,SVDDouble
  END INTERFACE
  


!  INTEGER,PRIVATE:: DBG=3
  CHARACTER(LEN=*),PUBLIC,PARAMETER:: FMM_VERSION = "1.96 (11 October 2002)"
  INTEGER, PRIVATE, PARAMETER :: DP = SELECTED_REAL_KIND(10,50)   ! "double" precision
  INTEGER, PRIVATE, PARAMETER :: SP = SELECTED_REAL_KIND(6,20)   ! "single" precision
  
CONTAINS
  
!+
SUBROUTINE DecompSingle(a, ipvt, errCode, cond)
! ---------------------------------------------------------------------------
! PURPOSE - Matrix triangularization by Gaussian elimination and estimation 
!  of the condition of the matrix. The matrix a is replaced by the LU 
!  decomposition of a rowwise permutation of itself. The array ipvt records
!  the row permutations of the partial pivoting. ipvt(k) is the index of the
!  kth pivot row and ipvt(n) is (-1)**(no. of interchanges). The optional
!  variable cond is an estimate of the conditioning of a. For the linear
!  system Ax=B, changes in A and B may cause changes cond times as large
!  in x. If cond+1.0 == cond, then a is singular to the working precision.
!  The determinant of a can be obtained on output as
!    det = ipvt(n)*a(1,1)*a(2,2)*a(3,3)*...*a(n,n)
! NOTES - If errCode is not zero, the remaining arguments are not to be
!   believed.

  REAL(SP),INTENT(IN OUT),DIMENSION(:,:):: a   ! matrix to be decomposed
  INTEGER,INTENT(OUT),DIMENSION(:):: ipvt   !    index of pivot rows
  INTEGER,INTENT(OUT):: errCode   ! =0 OK; =1 singular
  REAL(SP),INTENT(OUT),OPTIONAL:: cond  ! condition number

  REAL(SP):: anorm
  REAL(SP):: ek
  INTEGER:: j,k,m
  INTEGER,DIMENSION(1) :: mloc  ! receives the result of MAXLOC
  INTEGER:: n
  REAL(SP):: t
  REAL(SP),DIMENSION(SIZE(a,1)) :: work  ! scratch array
  REAL(SP):: ynorm,znorm
  REAL(SP),PARAMETER:: ZERO=0.0, ONE=1.0
!-----------------------------------------------------------------------
  errCode=0
  n=SIZE(a,1)

  IF (n <= 1) THEN
    IF (a(1,1) == ZERO) THEN
      errCode=1
    ELSE
      cond=ONE
    END IF
    RETURN        ! note abnormal RETURN
   END IF

  IF (Present(cond)) THEN
    anorm=ZERO   ! compute 1-norm of a
    DO j=1,n
      anorm=MAX(anorm, SUM(ABS(a(1:n,j))))
    END DO   
  END IF  

  DO k=1,n-1   ! Gaussian elimination with partial pivoting
    mloc= MAXLOC(ABS(a(k:n,k)))   ! pivot row
    m=k-1+mloc(1)
    ipvt(k)=m
    IF (m /= k) ipvt(n)=-ipvt(n)
    t=a(m,k)
    a(m,k)=a(k,k)
    a(k,k)=t
    IF (t /= ZERO) THEN
      t=ONE/t
      a(k+1:n,k)=-t*a(k+1:n,k)
      DO j=k+1,n   ! interchange and eliminate by columns.......
        t=a(m,j)
        a(m,j)=a(k,j)
        a(k,j)=t
        IF (t /= ZERO) a(k+1:n,j)=a(k+1:n,j) + t*a(k+1:n,k)
      END DO
    END IF
  END DO

  DO k=1,n   ! solve (a-transpose)*y=e
    t=ZERO
    IF (k > 1) t=DOT_PRODUCT(a(1:k-1,k), work(1:k-1))
    ek=ONE
    IF (t < ZERO) ek=-ONE
    IF (a(k,k) == ZERO) THEN    ! singular matrix
      errCode=1
      RETURN                            ! note abnormal RETURN
    END IF
    work(k)=-(ek+t)/a(k,k)
  END DO

  DO k=n-1,1,-1 
    t=ZERO
    t=work(k)*SUM(a(k+1:n,k))
    work(k)=t
    m=ipvt(k)
    IF (m /= k) THEN
      t=work(m)
      work(m)=work(k)
      work(k)=t
    END IF
  END DO

  IF (Present(cond)) THEN
    ynorm=SUM(ABS(work(1:n)))
    CALL Solve(a,work,ipvt)   ! solve a*z=y
    znorm=SUM(ABS(work(1:n)))
    cond=anorm*znorm/ynorm   ! estimate condition
    IF (cond < ONE) cond=ONE
  END IF
  
  RETURN
END Subroutine DecompSingle   ! ---------------------------------------------------

!+
SUBROUTINE SolveSingle(a, b, ipvt)
! ---------------------------------------------------------------------------
! PURPOSE - Solve the linear system a*x=b
  REAL(SP),INTENT(IN),DIMENSION(:,:):: a  ! the decomposed matrix from Decomp
  REAL(SP),INTENT(IN OUT),DIMENSION(:):: b ! in:right-hand side; out:solution
  INTEGER,INTENT(IN),DIMENSION(:):: ipvt  ! index of pivot rows (from Decomp)

  INTEGER:: i,k,m
  INTEGER:: n
  REAL(SP):: t
!----------------------------------------------------------------------------
  n=SIZE(a,1)
  DO k=1,n-1   ! forward elimination
    m=ipvt(k)
    t=b(m)
    b(m)=b(k)
    b(k)=t
!    DO i=k+1,n
!      b(i)=b(i)+a(i,k)*t
!    END DO
    b(k+1:n)=b(k+1:n) + t*a(k+1:n,k)
  END DO       

  DO k=n,1,-1      ! back substitution
    b(k)=b(k)/a(k,k)
    t=-b(k)
    DO i=1,k-1
      b(i)=b(i)+a(i,k)*t
    END DO       
  END DO
  RETURN
END Subroutine SolveSingle   ! ----------------------------------------------------

!+
SUBROUTINE FMMsplineSingle(x, y, b, c, d)
! ---------------------------------------------------------------------------
! PURPOSE - Compute the coefficients b,c,d for a cubic interpolating spline
!  so that the interpolated value is given by
!    s(x) = y(k) + b(k)*(x-x(k)) + c(k)*(x-x(k))**2 + d(k)*(x-x(k))**3
!      when x(k) <= x <= x(k+1)
!  The end conditions match the third derivatives of the interpolated curve to
!  the third derivatives of the unique polynomials thru the first four and
!  last four points.
!  Use Seval or Seval3 to evaluate the spline.
   REAL(SP),DIMENSION(:), INTENT(IN)  :: x ! abscissas of knots
   REAL(SP),DIMENSION(:), INTENT(IN)  :: y ! ordinates of knots
   REAL(SP),DIMENSION(:), INTENT(OUT) :: b ! linear coeff
   REAL(SP),DIMENSION(:), INTENT(OUT) :: c ! quadratic coeff.
   REAL(SP),DIMENSION(:), INTENT(OUT) :: d ! cubic coeff.

   INTEGER:: k,n
   REAL(SP):: t
   REAL(SP),PARAMETER:: ZERO=0.0, TWO=2.0, THREE=3.0
!----------------------------------------------------------------------------
  n=SIZE(x)

  IF (n < 3) THEN   ! Straight line - special case for n < 3
    b(1)=ZERO
    IF (n == 2) b(1)=(y(2)-y(1))/(x(2)-x(1))
    c(1)=ZERO
    d(1)=ZERO
    IF (n < 2) RETURN
    b(2)=b(1)
    c(2)=ZERO
    d(2)=ZERO
    RETURN
  END IF

!.....Set up tridiagonal system.........................................
!.    b=diagonal, d=offdiagonal, c=right-hand side
  d(1)=x(2)-x(1)
  c(2)=(y(2)-y(1))/d(1)
  DO k=2,n-1
    d(k)=x(k+1)-x(k)
    b(k)=TWO*(d(k-1)+d(k))
    c(k+1)=(y(k+1)-y(k))/d(k)
    c(k)=c(k+1)-c(k)
  END DO

!.....End conditions.  third derivatives at x(1) and x(n) obtained
!.       from divided differences.......................................
  b(1)=-d(1)
  b(n)=-d(n-1)
  c(1)=ZERO
  c(n)=ZERO
  IF (n > 3) THEN
    c(1)=c(3)/(x(4)-x(2))-c(2)/(x(3)-x(1))
    c(n)=c(n-1)/(x(n)-x(n-2))-c(n-2)/(x(n-1)-x(n-3))
    c(1)=c(1)*d(1)*d(1)/(x(4)-x(1))
    c(n)=-c(n)*d(n-1)*d(n-1)/(x(n)-x(n-3))
  END IF

  DO k=2,n    ! forward elimination
    t=d(k-1)/b(k-1)
    b(k)=b(k)-t*d(k-1)
    c(k)=c(k)-t*c(k-1)
  END DO

  c(n)=c(n)/b(n)   ! back substitution ( makes c the sigma of text)
  DO k=n-1,1,-1
    c(k)=(c(k)-d(k)*c(k+1))/b(k)
  END DO

!.....Compute polynomial coefficients...................................
  b(n)=(y(n)-y(n-1))/d(n-1)+d(n-1)*(c(n-1)+c(n)+c(n))
  DO k=1,n-1
    b(k)=(y(k+1)-y(k))/d(k)-d(k)*(c(k+1)+c(k)+c(k))
    d(k)=(c(k+1)-c(k))/d(k)
    c(k)=THREE*c(k)
  END DO
  c(n)=THREE*c(n)
  d(n)=d(n-1)

  RETURN
END Subroutine FMMsplineSingle   ! ---------------------------------------------------

!+
SUBROUTINE NaturalSplineSingle(x,y,b,c,d)
! ---------------------------------------------------------------------------
! PURPOSE - Construct the natural spline thru a set of points
! NOTES - A natural spline has zero second derivative at both endpoints.

  REAL(SP),INTENT(IN),DIMENSION(:):: x,y   ! coordinates of knots
  REAL(SP),INTENT(OUT),DIMENSION(:):: b,c,d  ! cubic coeff.

  INTEGER:: k,n
  REAL(SP),PARAMETER:: ZERO=0.0, TWO=2.0, THREE=3.0
!-----------------------------------------------------------------------
  n=SIZE(x)
 
  IF (n < 3) THEN   ! Straight line - special case for n < 3
    b(1)=ZERO
    IF (n == 2) b(1)=(y(2)-y(1))/(x(2)-x(1))
    c(1)=ZERO
    d(1)=ZERO
    b(2)=b(1)
    c(2)=ZERO
    d(2)=ZERO
    RETURN
  END IF

  d(1:n-1) = x(2:n)-x(1:n-1)  ! Put the h-array of the text into array d

!.....Set up the upper triangular system in locations 2 thru n-1 of
!        arrays b and c. B holds the diagonal and c the right hand side.
  b(2)=TWO*(d(1)+d(2))
  c(2)=(y(3)-y(2))/d(2)-(y(2)-y(1))/d(1)
  DO  k=3,n-1
    b(k)=TWO*(d(k-1)+d(k))-d(k-1)*d(k-1)/b(k-1)
    c(k)=(y(k+1)-y(k))/d(k)-(y(k)-y(k-1))/d(k-1)-d(k-1)*c(k-1)/b(k-1)
  END DO

  c(n-1)=c(n-1)/b(n-1)   ! Back substitute to get c-array
  DO  k=n-2,2,-1
    c(k)=(c(k)-d(k)*c(k+1))/b(k)
  END DO
  c(1)=ZERO
  c(n)=ZERO   ! c now holds the sigma array of the text 


!.....Compute polynomial coefficients ..................................
  b(n)=(y(n)-y(n-1))/d(n-1)+d(n-1)*(c(n-1)+c(n)+c(n))
  DO  k=1,n-1
    b(k)=(y(k+1)-y(k))/d(k)-d(k)*(c(k+1)+c(k)+c(k))
    d(k)=(c(k+1)-c(k))/d(k)
    c(k)=THREE*c(k)
  END DO
  c(n)=THREE*c(n)
  d(n)=d(n-1)
  RETURN
  
END Subroutine NaturalSplineSingle   ! --------------------------------------------

!+
FUNCTION SevalSingle(u, x,y, b,c,d) RESULT(SevalResult)
! ---------------------------------------------------------------------------
!  PURPOSE - Evaluate the cubic spline function
!     Seval=y(i)+b(i)!(u-x(i))+c(i)*(u-x(i))**2+d(i)*(u-x(i))**3
!           where  x(i) <= u < x(i+1)

!  NOTES- if u<x(1), i=1 is used;if u>x(n), i=n is used

  REAL(SP),INTENT(IN) :: u ! abscissa at which the spline is to be evaluated
  REAL(SP),INTENT(IN),DIMENSION(:):: x ! abscissas of knots
  REAL(SP),INTENT(IN),DIMENSION(:):: y ! ordinates of knots
  REAL(SP),INTENT(IN),DIMENSION(:):: b,c,d ! linear,quadratic,cubic coeff
  REAL(SP):: SevalResult

  INTEGER, SAVE :: i=1
  INTEGER :: j, k, n
  REAL(SP):: dx
!----------------------------------------------------------------------------
  n=SIZE(x)

!.....First check if u is in the same interval found on the
!        last call to Seval.............................................
  IF (  (i<1) .OR. (i >= n) ) i=1
  IF ( (u < x(i))  .OR.  (u >= x(i+1)) ) THEN
    i=1   ! binary search
    j=n+1

    DO
      k=(i+j)/2
      IF (u < x(k)) THEN
        j=k
      ELSE
        i=k
      END IF
      IF (j <= i+1) EXIT
    END DO
  END IF

  dx=u-x(i)   ! evaluate the spline
  SevalResult=y(i)+dx*(b(i)+dx*(c(i)+dx*d(i)))

  RETURN
END Function SevalSingle  ! -------------------------------------------------------

!+
SUBROUTINE Seval3Single(u, x,y, b,c,d, f,fp,fpp,fppp)
! ---------------------------------------------------------------------------
!  PURPOSE - Evaluate the cubic spline function
!     Seval=y(i)+b(i)!(u-x(i))+c(i)*(u-x(i))**2+d(i)*(u-x(i))**3
!           where  x(i) <= u < x(i+1)

!  NOTES- if u<x(1), i=1 is used;if u>x(n), i=n is used

  REAL(SP),INTENT(IN) :: u ! abscissa at which the spline is to be evaluated
  REAL(SP),INTENT(IN),DIMENSION(:):: x ! abscissas of knots
  REAL(SP),INTENT(IN),DIMENSION(:):: y ! ordinates of knots
  REAL(SP),INTENT(IN),DIMENSION(:):: b,c,d ! linear,quadratic,cubic coeff
  REAL(SP),INTENT(OUT),OPTIONAL:: f,fp,fpp,fppp ! function, 1st,2nd,3rd deriv

  INTEGER, SAVE :: i=1
  INTEGER :: j, k, n
  REAL(SP)    :: dx
  REAL(SP),PARAMETER:: TWO=2.0, THREE=3.0, SIX=6.0
!----------------------------------------------------------------------------
  n=SIZE(x)

!.....First check if u is in the same interval found on the
!        last call to Seval.............................................
  IF (  (i<1) .OR. (i >= n) ) i=1
  IF ( (u < x(i))  .OR.  (u >= x(i+1)) ) THEN
    i=1   ! binary search
    j=n+1

    DO
      k=(i+j)/2
      IF (u < x(k)) THEN
        j=k
      ELSE
        i=k
      END IF
      IF (j <= i+1) EXIT
    END DO
  END IF

  dx=u-x(i)   ! evaluate the spline
  IF (Present(f))    f=y(i)+dx*(b(i)+dx*(c(i)+dx*d(i)))
  IF (Present(fp))   fp=b(i)+dx*(TWO*c(i) + dx*THREE*d(i))
  IF (Present(fpp))  fpp=TWO*c(i) + dx*SIX*d(i)
  IF (Present(fppp)) fppp=SIX*d(i)

  RETURN
END Subroutine Seval3Single  ! -------------------------------------------------------

!+
SUBROUTINE Quanc8Single (Fun,a,b,abserr,relerr,result,errest,nofun,flag)
! -----------------------------------------------------------------------------
! PURPOSE - Estimate the integral of FUN(X) from A to B to a user provided 
!  tolerance. This is an automatic adaptive routine based on the 8-panel 
!  Newton-Cotes rule. 

  INTERFACE
    FUNCTION Fun(x) RESULT(f)  ! integrand function subprogram
    IMPLICIT NONE
      INTEGER,PARAMETER:: SP=SELECTED_REAL_KIND(6,20)
      REAL(SP),INTENT(IN)::x
      REAL(SP):: f
     END FUNCTION Fun
   END INTERFACE
   
   REAL(SP),INTENT(IN):: a        !  lower limit of integration
   REAL(SP),INTENT(IN):: b        !  upper limit of integration
   REAL(SP),INTENT(IN):: abserr   !  absolute error tolerance (>0)
   REAL(SP),INTENT(IN):: relerr   !  relative error tolerance (>0)

   REAL(SP),INTENT(OUT):: result   !  approx. value of integral
   REAL(SP),INTENT(OUT):: errest   !  estimate of actual error
   INTEGER,INTENT(OUT):: nofun     !  no. of function evaluations
   REAL(SP),INTENT(OUT):: flag     !  reliability indicator. (=0 is O.K.)

   REAL(SP),PARAMETER:: W0=3956.0_SP/14175.0_SP
   REAL(SP),PARAMETER:: W1=23552.0_SP/14175.0_SP
   REAL(SP),PARAMETER:: W2=-3712.0_SP/14175.0_SP
   REAL(SP),PARAMETER:: W3=41984.0_SP/14175.0_SP
   REAL(SP),PARAMETER:: W4=-18160.0_SP/14175.0_SP

  REAL(SP):: area
  REAL(SP):: cor11  
  REAL(SP):: esterr
  REAL(SP), DIMENSION(0:16):: f,x
  REAL(SP), DIMENSION(8,30):: fsave,xsave
  INTEGER:: i,j
  INTEGER:: lev,nim
  INTEGER:: LEVMAX = 30
  INTEGER,PARAMETER :: LEVMIN = 1
  INTEGER,PARAMETER:: LEVOUT = 6
  INTEGER:: NOFIN     !...Trouble if NOFUN reaches NOFIN
  INTEGER,PARAMETER:: NOMAX = 5000
  REAL(SP):: qprev,qnow,qdiff,qleft
  REAL(SP):: stone,step
  REAL(SP):: tolerr
  REAL(SP),DIMENSION(31):: qright
  REAL(SP),PARAMETER:: ZERO=0.0, HALF=0.5, ONE=1.0, SIXTEEN=16.0
!------------------------------------------------------------------------------

!..... S t a g e   1  (General initialization)..........................
   NOFIN=NOMAX-8*(LEVMAX-LEVOUT+2**(LEVOUT+1))

!...initialize running sums to zero...
  area=ZERO
  cor11=ZERO
  flag=ZERO   ! flag is REAL
  result=ZERO
  errest=ZERO
  nofun=0
  IF (a==b) RETURN
  
!..... S T A G E   2   (initialization for first interval)..............
  lev=0
  nim=1
  x(0)=a
  x(16)=b
  qprev=ZERO
  stone=(b-a)/SIXTEEN
  x(8)= HALF*(x(0)+x(16))
  x(4)= HALF*(x(0)+x(8))
  x(12)=HALF*(x(8)+x(16))
  x(2)= HALF*(x(0)+x(4))
  x(6)= HALF*(x(4)+x(8))
  x(10)=HALF*(x(8)+x(12))
  x(14)=HALF*(x(12)+x(16))

  DO j=0,16,2
    f(j)=Fun(x(j))
  END DO
  nofun=9

!..... S T A G E   3   (central calculation)............................
!...requires qprev,x0,x2,x4,....x(16),f0,f2,f4,....f16.
!...calculates x1,x3,....x15,f1,f3,....f15,qleft,qright,qnow,qdiff,area

30 DO j=1,15,2  ! keeps coming back here until successful
      x(j)=HALF*(x(j-1)+x(j+1))
      f(j)=Fun(x(j))
    END DO
    nofun=nofun+8
    step=(x(16)-x(0))/SIXTEEN
    qleft=(w0*(f(0)+f(8))+w1*(f(1)+f(7))+w2*(f(2)+f(6))+   &
                           w3*(f(3)+f(5))+w4*f(4))*step
    qright(lev+1)=(w0*(f(8)+f(16))+w1*(f(9)+f(15))+w2*(f(10)+f(14))+   &
                    w3*(f(11)+f(13))+w4*f(12))*step
    qnow=qleft+qright(lev+1)
    qdiff=qnow-qprev
    area=area+qdiff

!..... S T A G E   4   (interval convergence test) .....................
    esterr=ABS(qdiff)/1023.0
    tolerr=MAX(abserr,relerr*ABS(area))*(step/stone)
    IF (lev < LEVMIN) GO TO 50
    
    
    IF (lev >= LEVMAX) THEN
      flag=flag+ONE
      GO TO 70
    END IF  
    IF (nofun > NOFIN) THEN      !..... S T A G E   6  (trouble section) 
      nofin=nofin+nofin  ! Number of function values is about to exceed limit
      levmax=levout
      flag=flag+(b-x(0))/(b-a)
      GO TO 70
    END IF
    IF (esterr <= TOLERR) GO TO 70

!..... S T A G E   5   (no convergence).................................
!...locate next interval...
50  nim=nim+nim
    lev=lev+1

    fsave(1:8,lev)=f(9:16)   ! store right hand elements for future use
    xsave(1:8,lev)=x(9:16)
!    DO i=1,8
!      fsave(i,lev)=f(i+8)
!      xsave(i,lev)=x(i+8)
!    END DO

    qprev=qleft
    DO I=1,8   ! assemble left hand elements for immediate use
      f(18-2*i)=f(9-i)
      x(18-2*i)=x(9-i)
    END DO
    GO TO 30


!..... S T A G E   7   (integral converged) ............................
!...add contributions into running sums...
70 result=result+qnow
   errest=errest+esterr
   cor11=cor11+qdiff/1023.0

  DO                                ! 72 loop in text 
    IF (nim == 2*(nim/2)) EXIT   
    nim=nim/2
    lev=lev-1
  END DO
  nim=nim+1
  
  IF (lev <= 0) GO TO 80   ! looks like success!
  qprev=qright(lev)   ! assemble elements required for next interval
  x(0)=x(16)
  f(0)=f(16)
  DO i=1,8
    f(2*i)=fsave(i,lev)
    x(2*i)=xsave(i,lev)
  END DO
  GO TO 30

!..... S T A G E   8  (Finalize and return).............................
80 result=result+cor11
  IF (errest /= ZERO) THEN
    DO  
      IF (ABS(result)+errest > ABS(RESULT) ) EXIT
      errest=errest+errest  ! make sure ERREST is not < roundoff level ....
    END DO
  END IF

  RETURN
END Subroutine QuanC8Single   ! ---------------------------------------------------

!+
SUBROUTINE Rkf45Single (F,y,t,tout,relerr,abserr,iflag,work,iwork)
! ---------------------------------------------------------------------------
! PURPOSE - Integrate a system of neqn first order differential 
!  equations by the Fehlberg fourth-fifth order Runge-Kutta method.
! NOTES - RKF45 is primarily designed to solve non-stiff and mildly stiff
!   differential equations when derivative evaluations are inexpensive.
!   Rkf45 should generally not be used when the user is demanding 
!	high accuracy.

! ABSTRACT - Subroutine Rkf45 integrates a system of neqn first order
!  ordinary differential equations of the form
!	dy(i)/dt = f(t,y(1),y(2),...,y(neqn))
!	where the y(i) are given at t .
! Typically the subroutine is used to integrate from t to tout but it can
! be used as a one-step integrator to advance the solution a single step
! in the direction of tout. On return, the parameters in the call list are
! set for continuing the integration. The user has only to call Rkf45 again
! (and perhaps define a new value for tout). Actually, Rkf45 is an 
! interfacing routine which calls subroutine Rkfs for the solution. 
! Rkfs in turn calls subroutine Fehl which computes an approximate solution
! over one step. 

! Rkf45 uses the Runge-Kutta-Fehlberg(4,5) method described in the reference
! E. Fehlberg:  Low-order Classical Runge-Kutta formulas with Stepsize 
! Control, NASA TR R-315.

! The performance of Rkf45 is illustrated in the reference
! L.F.Shampine,H.A.Watts,S.Davenport: "Solving Non-stiff Ordinary
! Differential Equations - The State of the Art 
!  Sandia Laboratories Report SAND75-0182, to appear in SIAM Review. 

! First call to Rkf45:

! The user must provide storage in his calling program for the arrays in 
! the call list	y(neqn), work(3+6*neqn), iwork(5), supply 
! subroutine F(t,y,yp) and initialize the following parameters-

! y -- vector of initial conditions
! t -- starting point of integration , must be a variable
! tout -- output point at which solution is desired.
!   t=tout is allowed on the first call only, in which case
!   Rkf45 returns with iflag=2 if continuation is possible.
! relerr,abserr -- relative and absolute local error tolerances
!   which must be non-negative. relerr must be a variable while
!   abserr may be a constant. The code should normally not be
!   used with relative error control smaller than about i.e-8 .
!   to avoid limiting precision difficulties the code requires
!   relerr to be larger than an internally computed relative
!   error parameter which is machine dependent. In particular,
!   pure absolute error is not permitted. If a smaller than
!   allowable value of relerr is attempted, rkf45 increases
!   relerr appropriately and returns control to the user before
!   continuing the integration.
! iflag -- +1,-1 indicator to initialize the code for each new problem. 
!   Normal input is +1. the user should set iflag=-1 only when one-step
!   integrator control is essential. In this case, Rkf45 attempts to 
!   advance the solution a single step in the direction of tout each time
!   it is called. Since this mode of operation results in extra computing 
!   overhead, it should be avoided unless needed. 

! Output from Rkf45:

! y -- solution at t
! t - last point reached in integration.
! iflag = 2 -- integration reached tout. Indicates successful return and is 
!              the normal mode for continuing integration.
!	= -2 -- a single successful step in the direction of tout has been 
!               taken. Normal mode for continuing integration one step 
!               at a time.
!	= 3 -- integration was not completed because relative error tolerance 
!              was too small. relerr has been increased	appropriately for 
!              continuing.
!	= 4 -- integration was not completed because more than 3000 derivative
!              evaluations were needed. This is approximately 500 steps.
!	= 5 -- integration was not completed because solution vanished making
!              a pure relative error test impossible. Must use non-zero abserr
!              to continue. Using the one-step integration mode for one step
!              is a good way to proceed.
!	= 6 -- integration was not completed because requested accuracy could 
!              not be achieved using smallest allowable stepsize. User must
!              increase the error tolerance before continued integration can 
!              be attempted.
!	= 7 -- it is likely that Rkf45 is inefficient for solving this problem.
!              Too much output is restricting the natural stepsize choice. 
!              Use the one-step integrator mode.
!	= 8 -- invalid input parameters. This indicator occurs if any of the 
!              following is satisfied -	SIZE(y) <= 0; 
!              t=tout and iflag /= +1 or -1
!	       relerr or abserr < 0.
!              iflag == 0 or < -2 or > 8
! work,iwork -- information which is usually of no interest to the user but 
!    necessary for subsequent calls.
!    work(1),...,work(neqn) contain the first derivatives of the solution 
!    vector y at t. work(neqn+1) contains the stepsize h to be attempted on 
!    the next step. iwork(1) contains the derivative evaluation counter. 

! Subsequent calls to Rkf45 

! Subroutine Rkf45 returns with all information needed to continue the 
! integration. If the integration reached tout, the user need only define
! a new tout and call Rkf45 again. In the one-step integrator mode (iflag=-2)
! the user must keep in mind that each step taken is in the direction of the
! current tout. Upon reaching tout (indicated by changing iflag to 2),the 
! user must then define a new tout and reset iflag to -2 to continue in the 
! one-step integrator mode. 

! If the integration was not completed but the user still wants to continue 
! (iflag=3,4 cases), he just calls Rkf45 again. With iflag=3, the relerr 
! parameter has been adjusted appropriately for continuing the integration.
! In the case of iflag=4 the function counter will be reset to 0 and 
! another 3000 function evaluations are allowed. 

! However,in the case iflag=5, the user must first alter the error criterion
! to use a positive value of abserr before integration can proceed. If he 
! does not, execution is terminated. 

! Also, in the case iflag=6, it is necessary for the user to reset iflag 
! to 2 (or -2 when the one-step integration mode is being used) as well as 
! increasing either abserr,relerr or both before the integration can be 
! continued. If this is not done, execution will be terminated. The 
! occurrence of iflag=6 indicates a trouble spot (solution is changing 
! rapidly,singularity may be present) and it often is inadvisable to continue.

! If (flag=7 is encountered, the user should use the one-step integration mode
! with the stepsize determined by the code or consider switching to the 
! Adams codes DE/STEP,INTRP. If the user insists upon continuing the
! integration with Rkf45, he must reset iflag to 2 before calling Rkf45 again.
! Otherwise, execution will be terminated.

! If iflag=8 is obtained, integration can not be continued unless the invalid
! input parameters are corrected. 

! It should be noted that the arrays work,iwork contain information required
! for subsequent integration. Accordingly, work and iwork should not be 
! altered.

! This interfacing routine merely relieves the user of a long calling list 
! via the splitting apart of two working storage arrays.

INTERFACE
  SUBROUTINE F(arg1,arg2,arg3)
  IMPLICIT NONE
  INTEGER,PARAMETER:: SP = SELECTED_REAL_KIND(6,20)
  REAL(SP),INTENT(IN):: arg1   ! independent variable
  REAL(SP),INTENT(IN),DIMENSION(:):: arg2   ! vector of dependent variables
  REAL(SP),INTENT(OUT),DIMENSION(:):: arg3   ! vector of derivatives
  END Subroutine F
END INTERFACE  


  REAL(SP), INTENT(IN OUT),DIMENSION(:):: y   ! solution vector at t
  REAL(SP), INTENT(IN OUT):: t   ! independent variable
  REAL(SP), INTENT(IN OUT):: tout   ! output point at which solution is desired
  REAL(SP), INTENT(IN OUT):: relerr   ! relative error tolerance
  REAL(SP), INTENT(IN):: abserr   !  absolute error tolerance
  INTEGER, INTENT(IN OUT):: iflag   !   indicator for status of work
  REAL(SP), INTENT(IN OUT),DIMENSION(:):: work   ! work array
  INTEGER, INTENT(IN OUT),DIMENSION(:):: iwork   ! work array

  INTEGER :: k1,k2,k3,k4,k5,k6,k1m
  INTEGER:: neqn
!----------------------------------------------------------------------------
  neqn=SIZE(y)
  k1m=neqn+1   !.....compute indices for the splitting of the work array
  k1=k1m+1
  k2=k1+neqn
  k3=k2+neqn
  k4=k3+neqn
  k5=k4+neqn
  k6=k5+neqn
  
  CALL RkfsSingle(F,neqn,y,t,tout,relerr,abserr,iflag,work(1:neqn),work(neqn+1),  &
    work(k1:k1+neqn-1),work(k2:k2+neqn-1),work(k3:k3+neqn-1),               &
    work(k4:k4+neqn-1), work(k5:k5+neqn-1),work(k6), work(k6+1),            &
    iwork(1),iwork(2),iwork(3),iwork(4),iwork(5))
  RETURN
END Subroutine Rkf45Single   ! ----------------------------------------------------

!+
SUBROUTINE RkfsSingle (F,neqn,y,t,tout,relerr,abserr,iflag,yp,                    &
    h,f1,f2,f3,f4,f5,savre,savae,nfe,kop,init,jflag,kflag)
! ---------------------------------------------------------------------------
! PURPOSE - integrate a system of neqn first order differential equations by 
!   the Fehlberg fourth-fifth order Runge-Kutta method.

! NOTES - Rkfs integrates a system of first order ordinary differential
!  equations as described in the comments for Rkf45. The arrays 
!  yp,f1,f2,f3,f4,and f5 (of dimension at least neqn) and the variables 
!  h,savre,savae,nfe,kop,init,jflag,and kflag are used internally by the 
!  code and appear in the call list to eliminate local retention of variables
!  between calls. Accordingly, they should not be altered. Items of possible 
!  interest are
!   yp - derivative of solution vector at t
!   h - an appropriate stepsize to be used for the next step
!   nfe- counter on the number of derivative function evaluations 

  INTEGER, INTENT(IN):: neqn
  REAL(SP), INTENT(IN OUT),DIMENSION(:):: y   ! solution vector at t
  REAL(SP), INTENT(IN OUT):: t   ! independent variable
  REAL(SP), INTENT(IN):: tout   ! output point at which solution is desired
  REAL(SP), INTENT(IN OUT):: relerr   ! relative error tolerance
  REAL(SP), INTENT(IN):: abserr   !  absolute error tolerance
  INTEGER, INTENT(IN OUT):: iflag  !   indicator for status of work
  REAL(SP), INTENT(IN OUT),DIMENSION(:):: yp   ! derivatives 
  REAL(SP), INTENT(OUT):: h                    ! step size
  REAL(SP), INTENT(IN OUT),DIMENSION(:):: f1
  REAL(SP), INTENT(IN OUT),DIMENSION(:):: f2
  REAL(SP), INTENT(IN OUT),DIMENSION(:):: f3
  REAL(SP), INTENT(IN OUT),DIMENSION(:):: f4
  REAL(SP), INTENT(IN OUT),DIMENSION(:):: f5
  REAL(SP), INTENT(OUT):: savre
  REAL(SP), INTENT(OUT):: savae
  INTEGER, INTENT(OUT):: nfe    ! number of function evaluations
  INTEGER, INTENT(OUT):: kop
  INTEGER, INTENT(IN OUT):: init
  INTEGER, INTENT(IN OUT):: jflag
  INTEGER, INTENT(IN OUT):: kflag
  LOGICAL :: hfaild,output

INTERFACE
  SUBROUTINE F(arg1,arg2,arg3)
  IMPLICIT NONE
  INTEGER,PARAMETER:: SP = SELECTED_REAL_KIND(6,20)
  REAL(SP),INTENT(IN):: arg1
  REAL(SP),INTENT(IN),DIMENSION(:):: arg2
  REAL(SP),INTENT(OUT),DIMENSION(:):: arg3
  END Subroutine F
END INTERFACE  

  REAL(SP):: a,ae,dt,ee,eeoet,esttol,et
  REAL(SP):: eps
  REAL(SP):: hmin
  INTEGER:: k,mflag

  REAL(SP),PARAMETER:: REMIN = 1E-12   ! minimum acceptable value of relerr.  
    ! Attempts to obtain higher accuracy with this routine are usually very
    ! expensive and unsuccessful.

  INTEGER,PARAMETER:: MAXNFE = 3000 ! the expense is controlled by restricting
    ! the number of function evaluations to be approximately MAXNFE. 
    ! As set, this corresponds to about 500 steps.
  REAL(SP):: rer,s 
  REAL(SP):: scale                                   
  REAL(SP):: tol,toln,u26
  REAL(SP):: ypk
  REAL(SP),PARAMETER:: ZERO=0.0, ONE=1.0
!-----------------------------------------------------------------------
  eps=EPSILON(ZERO)
  u26=26.0*eps

!.....check input parameters............................................
  IF (neqn < 1 .OR. relerr < ZERO .OR. abserr < ZERO) THEN
    iflag=8
    RETURN
  END IF  
  mflag=ABS(iflag)
  IF (mflag == 0 .OR. mflag > 8) THEN
    iflag=8
    RETURN
  END IF
  
  IF (mflag == 1) GO TO 50

!.....check continuation posibilities...................................
  IF ((t == tout) .AND. (kflag /= 3)) THEN
    iflag=8   ! invalid flag handler
    RETURN
  END IF  
  IF (mflag /= 2) GO TO 25

!...you get here if iflag=+2 or -2
  IF ((kflag == 3) .OR. (init == 0)) GO TO 45
  IF (kflag == 4) GO TO 40
  IF ((kflag == 5) .AND. (abserr == ZERO)) GO TO 30
  IF ((kflag == 6) .AND. (relerr <= savre) .AND. (abserr <= savae))  GO TO 30
  GO TO 50


!...you get here if iflag=3,4,5,6,7, or 8.....
25 IF (iflag == 3) GO TO 45
   IF (iflag == 4) GO TO 40
   IF ((iflag == 5) .AND. (abserr > ZERO)) GO TO 45


!.....integration cannot be continued since user did not respond to the
!..... instructions pertaining to iflag=5,6,7, or 8
30  STOP   !!! maybe this is too strong


!.....reset function evaluation counter.................................
40 nfe=0
   IF (mflag == 2) GO TO 50


!.....reset flag value from previous call...............................
45 iflag=jflag
   IF (kflag == 3) mflag=ABS(iflag)

50 jflag=iflag   ! save for subsequent calls
   kflag=0       ! set the continuation flag
  savre=relerr   ! save for subsequent calls
  savae=abserr   ! save for subsequent calls


!.....restrict relative error tolerance to be at least as large as
!       2*eps+remin to avoid limiting precision difficulties arising
!       from impossible accuracy requests.
  rer=eps+eps+remin
  IF (relerr >= rer) GO TO 55


!.....relative error tolerance too small................................
  relerr=rer
  iflag=3
  kflag=3
  RETURN

55 dt=tout-t
   IF (mflag == 1) GO TO 60
   IF (init == 0) GO TO 65
   GO TO 80

!.....initialization....................................................
60 init=0
   kop=0
   a=t
   CALL f(a,y,yp)
   nfe=1
   IF (t /= tout) GO TO 65
   iflag=2
   RETURN

65 init=1
   h=ABS(dt)
   toln=ZERO
   DO  k=1,neqn
     tol=relerr*ABS(y(k))+abserr
     IF (tol <= ZERO) CYCLE
     toln=tol
     ypk=ABS(yp(k))
     IF (ypk*h**5 > tol) h=(tol/ypk)**0.2
   END DO

  IF (toln <= ZERO) h=ZERO
  h=MAX(h, u26*MAX(ABS(t), ABS(dt)))
  jflag=SIGN(2,iflag)

!.....set stepsize for integration in the direction from t to tout
80 h=SIGN(h,dt)

!.....test to see if rkf45 is being severely impacted by too many outputs
  IF (ABS(h) >= ABS(dt+dt)) kop=kop+1
  IF (kop /= 100) GO TO 85
  kop=0
  iflag=7
  RETURN

85 IF(ABS(dt) > u26*ABS(t)) GO TO 95

!.....if too close to output point, extrapolate and return
  DO  k=1,neqn
    y(k)=y(k)+dt*yp(k)
  END DO
  a=tout
  CALL F(a,y,yp)
  nfe=nfe+1
  GO TO 300


!.....initialize output point indicator.................................
95 output=.FALSE.


!.....scale the error tolerances to avoid premature underflow in the
!        error tolerance
  scale=2.0/relerr
  ae=scale*abserr


!.....step by step integration..........................................
100 hfaild=.false.
!...set smallest allowable step size...
  hmin=u26*ABS(t)

!.....adjust stepsize if necessary to hit the output point.
!.....look ahead two steps to avoid drastic changes in the stepsize and
!.....thus lessen the impact of output points on the code.
  dt=tout-t
  IF (ABS(dt) >= ABS(h+h)) GO TO 200
  IF (ABS(dt) > ABS(h)) GO TO 150


!.....the next successful step will complete the integration to the
!     output point
  output=.TRUE.
  h=dt
  GO TO 200

150 h=0.5*dt

!.....core integrator for taking a single step..........................

200 IF (nfe <= maxnfe) GO TO 220

!.....too much work...
  iflag=4
  kflag=4
  RETURN


!.....advance an approximate solution over one step of length h
220 CALL FehlSingle(f,y,t,h,yp,f1,f2,f3,f4,f5,f1)
    nfe=nfe+5


!.....compute and test allowable tolerances versus local error estimates
!.....  and remove scaling of tolerances. note that relative error is
!.....  is measured with respect to the average of the magnitudes of the
!.....   solution at the beginning and end of the step.
  eeoet=ZERO
  DO  k=1,neqn
    et=ABS(y(k))+ABS(f1(k))+ae
    IF (et > ZERO) GO TO 240
  
  
!.....Inappropriate error tolerance.....................................
    iflag=5
    RETURN
  
240 ee=ABS((-2090.0*yp(k)+(21970.0*f3(k)-15048.0*f4(k)))+               &
      (22528.0*f2(k)-27360.0*f5(k)))
    eeoet=MAX(eeoet,ee/et)
  END DO
  esttol=ABS(h)*eeoet*scale/752400.0
  IF (esttol <= ONE) GO TO 260


!.....unsuccessful step. reduce the stepsize and start again............
!     (decrease is limited to a factor of 1/10)
  hfaild=.TRUE.
  output=.FALSE.
  s=0.1
  IF (esttol < 59049.0) s=0.9/esttol**0.2
  h=s*h
  IF (ABS(h) > hmin) GO TO 200


!.....requested error unobtainable at smallest allowable stepsize.......
  iflag=6
  kflag=6
  RETURN


!.....successful step.   store solution at t+h and evaluate derivatives
!     there..........................
260 t=t+h
    DO  k=1,neqn
      y(k)=f1(k)
    END DO
    a=t
    CALL F(a,y,yp)
    nfe=nfe+1


! Choose next stepsize.  The increase is limited to a factor of 5.
! If step failure has just occured, next stepsize is not allowed
!  to increase.
  s=5.0
  IF (esttol > 1.889568E-4) s=0.9/esttol**0.2
  IF (hfaild) s=MIN(s,ONE)
  h=SIGN(MAX(s*ABS(h),hmin),h)

!.....E N D   O F   C O R E   I N T E G R A T O R   ....................


!.....Should we take another step? .....................................
  IF (output) GO TO 300
  IF (iflag > 0) GO TO 100


!.....integration successfully completed................................
!...one step mode...
  iflag=-2
  RETURN

!...interval mode...
300 t=tout
    iflag=2
    RETURN
END SUBROUTINE RkfsSingle   ! -----------------------------------------------------

!+
SUBROUTINE FehlSingle (F,y,t,h,yp,f1,f2,f3,f4,f5,s)
! ---------------------------------------------------------------------------
! PURPOSE - Integrate a system of neqn first-order ordinary differential
!  equations of the form dy/dt=F(t,y) [y is a vector], where the initial
!  values of y and the initial values of yp, the derivatives are specified
!  at the starting point t.  Fehl advances the solution over the fixed step h
!  and returns the fifth order (sixth-order locally) solution approximation 
!  at t+h in array s.
INTERFACE
  SUBROUTINE F(arg1,arg2,arg3)
  IMPLICIT NONE
  INTEGER,PARAMETER:: SP = SELECTED_REAL_KIND(6,20)
  REAL(SP),INTENT(IN):: arg1
  REAL(SP),INTENT(IN),DIMENSION(:):: arg2
  REAL(SP),INTENT(OUT),DIMENSION(:):: arg3
  END Subroutine F
END INTERFACE  

!!!  INTEGER, INTENT(IN):: neqn ! number of equations
  REAL(SP), INTENT(IN),DIMENSION(:):: y ! array of length neqn; function at t
  REAL(SP), INTENT(IN):: t   ! starting point
  REAL(SP), INTENT(IN):: h   ! step size
  REAL(SP), INTENT(IN),DIMENSION(:)::  yp ! array of length neqn; derivatives at t
  REAL(SP), INTENT(OUT),DIMENSION(:):: f1 ! array of length neqn for internal use
  REAL(SP), INTENT(OUT),DIMENSION(:):: f2 ! array of length neqn for internal use
  REAL(SP), INTENT(OUT),DIMENSION(:):: f3 ! array of length neqn for internal use
  REAL(SP), INTENT(OUT),DIMENSION(:):: f4 ! array of length neqn for internal use
  REAL(SP), INTENT(OUT),DIMENSION(:):: f5 ! array of length neqn for internal use
  REAL(SP), INTENT(OUT),DIMENSION(:):: s  ! array of length neqn; the results

  REAL(SP),PARAMETER:: C1 = 0.25_SP
  REAL(SP),PARAMETER:: C2 = 3.0_SP/32.0_SP
  REAL(SP),PARAMETER:: C3 = 3.0_SP
  REAL(SP),PARAMETER:: C4 = 3.0_SP/8.0_SP
  REAL(SP),PARAMETER:: C5 = 1.0_SP/2197.0_SP
  REAL(SP),PARAMETER:: C6 = 1932.0_SP
  REAL(SP),PARAMETER:: C7 = 7296.0_SP
  REAL(SP),PARAMETER:: C8 = -7200.0_SP
  REAL(SP),PARAMETER:: C9 = 12.0_SP/13.0_SP
  REAL(SP),PARAMETER:: C10 = 1.0_SP/4104.0_SP
  REAL(SP),PARAMETER:: C11 = 8341.0_SP
  REAL(SP),PARAMETER:: C12 = -845.0_SP
  REAL(SP),PARAMETER:: C13 = 29440.0_SP
  REAL(SP),PARAMETER:: C14 = -32832.0_SP
  REAL(SP),PARAMETER:: C15 = 1.0_SP/20520.0_SP
  REAL(SP),PARAMETER:: C16 = -6080.0_SP
  REAL(SP),PARAMETER:: C17 = 9295.0_SP
  REAL(SP),PARAMETER:: C18 = -5643.0_SP
  REAL(SP),PARAMETER:: C19 = 41040.0_SP
  REAL(SP),PARAMETER:: C20 = -28352.0_SP
  REAL(SP),PARAMETER:: C21 = 0.5_SP
  REAL(SP),PARAMETER:: C22 = 1.0_SP/7618050.0_SP
  REAL(SP),PARAMETER:: C23 = 902880.0_SP
  REAL(SP),PARAMETER:: C24 = 3855735.0_SP
  REAL(SP),PARAMETER:: C25 = -1371249.0_SP
  REAL(SP),PARAMETER:: C26 = 3953664.0_SP
  REAL(SP),PARAMETER:: C27 = 277020.0_SP

  REAL(SP):: ch
!-----------------------------------------------------------------------
  ch=C1*h  
  f5=y+ch*yp  
  CALL F(t+ch, f5,f1)

  ch=C2*h  
  f5=y + ch*(yp+C3*f1)
  CALL F(t+C4*h, f5,f2)  

  ch=C5*h
  f5=y + ch*(C6*yp + (C7*f2+C8*f1) )
  CALL F(t+C9*h, f5,f3)  

  ch=C10*h
  f5=y + ch*( (C11*yp + C12*f3) + (C13*f2+C14*f1) )
  CALL F(t+h, f5,f4)

  ch=C15*h
  f1=y + ch*( C16*yp + (C17*f3+C18*f4) + (C19*f1+C20*f2) )  
  CALL F(t+C21*h, f1, f5)

!.....compute approximate solution at t+h
  ch=C22*h
  s=y + ch*( C23*yp + (C24*f3+C25*f4) + (C26*f2+C27*f5) )  
  RETURN
END SUBROUTINE FehlSingle   ! -----------------------------------------------------


!+
FUNCTION ZeroinSingle (a,b,F,tol) RESULT(z)
! ---------------------------------------------------------------------------
! PURPOSE - Compute a zero of f in the interval (a,b)

  REAL(SP),INTENT(IN):: a,b   ! left and right endpoints of interval
  REAL(SP),INTENT(IN):: tol     ! desired interval of uncertainity 
      
  REAL(SP):: z
INTERFACE
  FUNCTION F(x) RESULT(g)
    IMPLICIT NONE
    INTEGER,PARAMETER:: SP = SELECTED_REAL_KIND(6,20)
    REAL(SP),INTENT(IN):: x
    REAL(SP):: g
  END Function F
END INTERFACE    
  
  INTEGER:: errCode
  INTEGER,PARAMETER:: MAXITER=25
  INTEGER:: neval   ! not used
  REAL(SP):: xZero,fZero
!----------------------------------------------------------------------------
  CALL BrentZeroSingle(a,b,F,tol,MAXITER,neval,errCode,xZero,fZero)
  z=xZero
  RETURN
END Function ZeroinSingle   ! -----------------------------------------------------


!+
SUBROUTINE BrentZeroSingle(ax,bx,F,tol,maxIter,neval,errCode,xZero,fZero)
! ---------------------------------------------------------------------------
! PURPOSE - Compute a zero of F in the interval (ax,bx)

  REAL(SP),INTENT(IN):: ax,bx   ! left and right endpoints of interval
  REAL(SP),INTENT(IN):: tol     ! desired interval of uncertainity 
  INTEGER,INTENT(IN):: maxIter   ! max number of iterations allowed. 25 is good
  INTEGER,INTENT(OUT):: neval
  INTEGER,INTENT(OUT):: errCode   ! =0 is OK; =1 too many iterations
                                  ! =2 if F(ax) and F(bx) have the same sign
  REAL(SP),INTENT(OUT):: xZero,fZero ! the last and best value of the zero                                   
      

INTERFACE
  FUNCTION F(x) RESULT(g)
    IMPLICIT NONE
    INTEGER,PARAMETER:: SP = SELECTED_REAL_KIND(6,20)
    REAL(SP),INTENT(IN):: x
    REAL(SP):: g
  END Function F
END INTERFACE    

  REAL(SP):: a,b,c,d,e,eps
  REAL(SP):: fa,fb,fc,tol1
  INTEGER:: kIter
!  INTEGER:: method   ! =0 bisection; =1 linear; =2 inverse quadratic
  REAL(SP):: xm,p,q,r,s
  REAL(SP),PARAMETER:: ZERO=0.0, ONE=1.0, TWO=2.0, THREE=3.0, HALF=0.5
!----------------------------------------------------------------------------
  eps=EPSILON(ax)
  tol1=ONE+eps

  a=ax   ! initialization
  b=bx
  fa=f(a)
  fb=f(b)
  neval=2
  IF (fa==ZERO) THEN
    xZero=a
    fZero=ZERO
    errCode=0
    RETURN
  END IF  
  IF (fb==ZERO) THEN
    xZero=b
    fZero=ZERO
    errCode=0
    RETURN
  END IF  
  IF (fa*fb > ZERO) THEN
    xZero=HUGE(a)
    fZero=HUGE(a)
    errCode=2
    RETURN
  END IF
! The trivial cases have now been dealt with. On to the real work...

  c=a
  fc=fa
  d=b-a
  e=d
    
  DO kIter=1,maxIter
    IF ( (fb>0 .AND. fc>0) .OR. (fb<0 .AND. fc<0) ) THEN
      c=a  ! we insist that b and c straddle the zero
      fc=fa
      d=b-a
      e=d
    END IF

    IF (ABS(fc) < ABS(fb)) THEN
      a=b    ! we insist that b be the better guess of b and c
      b=c
      c=a
      fa=fb
      fb=fc
      fc=fa
    END IF

    tol1=TWO*eps*ABS(b)+HALF*tol   ! convergence test
    xm=HALF*(c-b)
    IF (ABS(xm) <= tol1 .OR. fb==ZERO) THEN
      xZero=b
      fZero=fb
      errCode=0   ! SUCCESS! The proper way to leave
      RETURN
    END IF

!    WRITE(DBG,*) "BrentZero, start of kIter=",kIter
!    WRITE(DBG,'(A,2ES25.16)') "a,fa=", a,fa
!    WRITE(DBG,'(A,2ES25.16)') "b,fb=", b,fb
!    WRITE(DBG,'(A,2ES25.16)') "c,fc=", c,fc
!    WRITE(DBG,'(A,2ES25.16)') "xm,tol1=", xm,tol1

    IF (ABS(e) < tol1 .OR. ABS(fa) <= ABS(fb) ) THEN
      d=xm   ! bisection
      e=d
!      method=0
    ELSE
      IF (a==c) THEN
        s=fb/fa   ! linear interpolation
        p=TWO*xm*s
        q=ONE-s
!        method=1
      ELSE
        q=fa/fc   ! inverse quadratic interpolation
        r=fb/fc
        s=fb/fa
        p=s*(TWO*xm*q*(q-r)-(b-a)*(r-ONE))
        q=(q-ONE)*(r-ONE)*(s-ONE)
!        method=2
      END IF
      IF (p > ZERO) q=-q   ! adjust signs
      p=ABS(p)
      IF (p+p >= (THREE*xm*q-ABS(tol1*q)) .OR. p+p >= ABS(e*q) ) THEN
        d=xm   ! don't interpolate. Use bisection
        e=d
 !       method=-1
      ELSE
        e=d   ! OK, use interpolation
        d=p/q
      END IF
    END IF  
   
    a=b   ! complete step. a becomes the previous iteration
    fa=fb
   IF (ABS(d) > tol1) THEN
     b=b+d
   ELSE  
     b=b+SIGN(tol1,xm)
   END IF
     
   fb=F(b)   ! the newest and best value (we hope)
   neval=neval+1   ! keep count of the function evaluations
!!   IF ((fb*(fc/ABS(fc))) > ZERO) GO TO 20
!    WRITE(DBG,*) "BrentZero, end of kIter=",kIter, "   method=",method
!    WRITE(DBG,'(A,2ES25.16)') " new b,fb=", b,fb
!    WRITE(DBG,'(A,2ES25.16)') "d,e=", d,e
!    WRITE(DBG,'(A,2ES25.16)') "p,q=", p,q
  END DO
  
! The loop should never terminate. If it does, return the last iteration
!  and set errCode to 1
  xZero=b
  fZero=fb
  errCode=1
  RETURN  

END Subroutine BrentZeroSingle   ! ------------------------------------------------

!+
FUNCTION FminSingle (ax,bx,f,tol) RESULT (xopt)
! ---------------------------------------------------------------------------
! PURPOSE -

  INTERFACE
    FUNCTION F(x) RESULT(ff)  ! integrand function subprogram
    IMPLICIT NONE
      INTEGER,PARAMETER:: SP = SELECTED_REAL_KIND(6,20)
      REAL(SP),INTENT(IN)::x
      REAL(SP):: ff
    END FUNCTION F
  END INTERFACE

  REAL(SP),INTENT(IN):: ax   !  left endpoint of initial interval
  REAL(SP),INTENT(IN):: bx   !  right endpoint of initial interval
  REAL(SP),INTENT(IN):: tol  !  desired length of interval of uncertainity
  REAL(SP) :: xopt ! the result
  
  INTEGER:: errCode
  REAL(SP):: fzero
  INTEGER,PARAMETER:: MAX_ITER = 100
  INTEGER:: neval
!----------------------------------------------------------------------------
  CALL BrentMinSingle(ax,bx,F,tol,MAX_ITER,neval,errCode,xopt,fzero)
  
  RETURN
END Function FminSingle   ! -------------------------------------------------------

!+
SUBROUTINE BrentMinSingle(ax,bx,F,tol,maxIter,neval,errCode,xZero,fZero)
! ---------------------------------------------------------------------------
! PURPOSE -

  INTERFACE
    FUNCTION F(x) RESULT(ff)  ! integrand function subprogram
    IMPLICIT NONE
      INTEGER,PARAMETER:: SP = SELECTED_REAL_KIND(6,20)
      REAL(SP),INTENT(IN)::x
      REAL(SP):: ff
    END FUNCTION F
  END INTERFACE

  REAL(SP),INTENT(IN):: ax   !  left endpoint of initial interval
  REAL(SP),INTENT(IN):: bx   !  right endpoint of initial interval
  REAL(SP),INTENT(IN):: tol  !  desired length of interval of uncertainity
  INTEGER,INTENT(IN):: maxIter   ! max number of iterations allowed. 25 is good
  INTEGER,INTENT(OUT):: neval
  INTEGER,INTENT(OUT):: errCode   ! =0 is OK; =1 too many iterations
  REAL(SP),INTENT(OUT):: xZero,fZero ! the last and best value of the minimum                                  
      
  REAL(SP), PARAMETER :: C=0.381966011 ! (3-Sqrt(5))/2

  REAL(SP) :: a,b,d,e,eps,xm,p,q,r,tol1,tol2,u,v,w
  REAL(SP) :: fu,fv,fw,fx,x
  INTEGER :: iter
  REAL(SP),PARAMETER:: ZERO=0.0, ONE=1.0, TWO=2.0, THREE=3.0, HALF=0.5
!----------------------------------------------------------------------------
  errCode=0
  eps=SQRT(EPSILON(ONE))

  a=ax ! initialization
  b=bx
  v=a+c*(b-a)
  w=v
  x=v
  e=ZERO
  fx=F(x)
  fv=F(v)
  neval=2
  fw=fx

  DO iter=1,maxIter   ! iterate until minimum is found
    xm=HALF*(a+b)
    tol1=eps*ABS(x)+tol/THREE
    tol2=tol1+tol1
    IF (ABS(x-xm)<=(tol2-HALF*(b-a)) ) EXIT    ! check stopping criterion.

    IF (ABS(e) > tol1) THEN   ! is golden section necessary???
      r=(x-w)*(fx-fv)   ! trial parabolic fit
      q=(x-v)*(fx-fw)
      p=(x-v)*q-(x-w)*r
      q=TWO*(q-r)
      IF (q > ZERO) p=-p
      q=ABS(q)
      r=e
      e=d
      IF (ABS(p)>=ABS(HALF*q*r) .OR. p<=q*(a-x) .OR. p>=q*(b-x)) GO TO 40

      d=p/q   ! a parabolic interpolation step
      u=x+d
      IF ((u-a)<tol2 .OR. (b-u)<tol2) D=SIGN(tol1,xm-x)
      GO TO 50
    END IF

!40       CONTINUE
 40 IF (x >= xm) THEN   ! a golden section step
      e=a-x
    ELSE
      e=b-x
    END IF
    d=c*e

50  IF (ABS(d) >= tol1) THEN   ! f must not be evaluated too close to x
      u=x+d
    ELSE
      u=x+SIGN(tol1,d)
    END IF
    fu=F(u)
    neval=neval+1

    IF (fu <= fx) THEN   ! update a,b,v,w, and x
      IF (u >= x) THEN
        a=x
      ELSE
        b=x
      END IF
      v=w
      fv=fw
      w=x
      fw=fx
      x=u
      fx=fu
    ELSE
      IF (u<x) THEN
        a=u
      ELSE
        b=u
      END IF
      IF (fu<=fw .OR. w==x) THEN
        v=w
        fv=fw
        w=u
        fw=fu
      ELSE IF (fu<=fv .OR. v==x .OR. v==w) THEN
        v=u
        fv=fu
      END IF
    END IF
  END DO

   xZero=x   ! end of main loop
   fZero=fx
  RETURN
END Subroutine BrentMinSingle   ! -------------------------------------------------

!+
SUBROUTINE SvdSingle(a,w,matu,u,matv,v,ierr)
! ---------------------------------------------------------------------------
! PURPOSE - Determine the singular value decomposition of a real m by n 
!  rectangular matrix.  Householder bidiagonalization and a variant of the 
!  QR algorithm are used. This subroutine is a translation of the Algol 
!  procedure SVD, Num. Math. 14, 403-420(1970) by Golub and Reinsch.
!  Handbook for Auto. Comp., vol II-Linear algebra, 134-151(1971).
!  This is a modified version of a routine from the Eispack collection by 
!  the Nats project. Modified to eliminate machep.

  REAL(SP), INTENT(IN),DIMENSION(:,:) :: a  ! matrix to be decomposed. 
      ! on output, a is unaltered (unless overwritten by u or v).
  REAL(SP),INTENT(OUT),DIMENSION(:):: w ! w contains the n (non-negative) 
      ! singular values of a (the diagonal elements of s).  They are
      ! unordered. If an error exit is made, the singular values should be 
      ! correct for indices ierr+1,ierr+2,...,n.
  LOGICAL,INTENT(IN):: matu   ! matu should be set to .TRUE. if the u matrix
      ! in the decomposition is desired, and to .FALSE. otherwise.
  REAL(SP),INTENT(OUT),DIMENSION(:,:):: u   ! u contains the matrix u of 
      ! orthogonal column vectors of the decomposition if matu has been set
      ! to .TRUE. Otherwise, u is used as a temporary array.  
      ! u may coincide with a.  If an error exit is made, the columns of u 
      ! corresponding to indices of correct singular values should be correct.
  LOGICAL,INTENT(IN):: matv   ! matv should be set to .TRUE. if the v matrix
      ! in the decomposition is desired, and to .FALSE. otherwise.
  REAL(SP),INTENT(OUT),DIMENSION(:,:):: v   ! v contains the matrix v (orthogonal)
      ! of the decomposition if matv has been set to .TRUE.  Otherwise v is 
      ! not referenced. v may also coincide with a if u is not needed.  If an
      ! error exit is made, the columns of v corresponding to indices of 
      ! correct singular values should be correct.
  INTEGER,INTENT(OUT):: ierr   ! zero for normal return, k if the k-th 
      ! singular value has not been determined after 30 iterations.

  REAL(SP) :: c,f,g,h
  INTEGER :: i,j,k,l, i1,k1
  INTEGER:: its
  INTEGER:: m,n   ! dimensions of a
 
  REAL(SP),DIMENSION(SIZE(a,2)):: rv1    ! rv1 (n) is a temporary storage array.
  REAL(SP):: s,x,y,z,scale,anorm

  REAL(SP),PARAMETER:: ZERO = 0.0, ONE=1.0, TWO=2.0
!----------------------------------------------------------------------------
  m=SIZE(a,1)
  n=SIZE(a,2)
  IF (SIZE(w) < n) THEN
    ierr=-1
    RETURN
  END IF
  
  IF (matu .AND. (SIZE(u,1) < m .OR. SIZE(u,2) < n) ) THEN
    ierr=-2
    RETURN
  END IF
  
  IF (matv .AND. (SIZE(v,1) < n .OR. SIZE(v,2) < n) ) THEN
    ierr=-3
    RETURN
  END IF
    
  ierr = 0

!DO  i = 1, m
!  DO  j = 1, n
!    u(i,j) = a(i,j)
!  END DO
!END DO
  u(1:m,1:n)=a

  g = ZERO
  scale = ZERO
  anorm = ZERO

DO  i = 1, n   ! Householder reduction to bidiagonal form 
!!!  l = i + 1
  rv1(i) = scale * g
  g = ZERO
  s = ZERO
  scale = ZERO
  IF (i > m) GO TO 210
  
  DO  k = i, m
    scale = scale + ABS(u(k,i))
  END DO
  
  IF (scale == ZERO) GO TO 210
  
  DO  k = i, m
    u(k,i) = u(k,i) / scale
    s = s + u(k,i)**2
  END DO
  
  f = u(i,i)
  g = -SIGN(SQRT(s),f)
  h = f * g - s
  u(i,i) = f - g
!!!  IF (i == n) GO TO 190
  
  DO  j = i+1, n
    s = ZERO   
    DO  k = i, m
      s = s + u(k,i) * u(k,j)
    END DO
    
    f = s / h
    
    DO  k = i, m
      u(k,j) = u(k,j) + f * u(k,i)
    END DO
  END DO
  u(i:m,i)=scale*u(i:m,i)
!!!190 DO  k = i, m
!!!      u(k,i) = scale * u(k,i)
!!!    END DO
  
210 w(i) = scale * g
    g = ZERO
    s = ZERO
    scale = ZERO
    IF (i > m .OR. i == n) GO TO 290
  
    DO  k = i+1, n
      scale = scale + ABS(u(i,k))
    END DO
  
    IF (scale == ZERO) GO TO 290
  
    DO  k = i+1, n
      u(i,k) = u(i,k) / scale
      s = s + u(i,k)**2
    END DO
  
    f = u(i,i+1)
    g = -SIGN(SQRT(s),f)
    h = f * g - s
    u(i,i+1) = f - g
  
    DO  k = i+1, n
      rv1(k) = u(i,k) / h
    END DO
  
    IF (i == m) GO TO 270
  
    DO  j = i+1, m
      s = ZERO
      DO  k = i+1, n
        s = s + u(j,k) * u(i,k)
      END DO
      DO  k = i+1, n
        u(j,k) = u(j,k) + s * rv1(k)
      END DO
    END DO
  
270 DO  k = i+1, n
      u(i,k) = scale * u(i,k)
    END DO
  
290 anorm = MAX(anorm,ABS(w(i))+ABS(rv1(i)))
  END DO
!     .......... accumulation of right-hand transformations ..........
  IF (.NOT. matv) GO TO 410
!     .......... for i=n step -1 until 1 do -- ..........
  DO i=n,1,-1
    IF (i == n) GO TO 390
    IF (g == ZERO) GO TO 360
  
    DO  j = i+1, n
      v(j,i) = (u(i,j) / u(i,i+1)) / g   ! double division, YES!
    END DO
  
    DO  j = i+1, n
      s = ZERO
      DO  k = i+1, n
        s = s + u(i,k) * v(k,j)
      END DO
      DO  k = i+1, n
        v(k,j) = v(k,j) + s * v(k,i)
      END DO
    END DO
  
360 v(i,i+1:n)=ZERO
    v(i+1:n,i)=ZERO

!    DO  j = i+1, n
!      v(i,j) = ZERO
!      v(j,i) = ZERO
!    END DO
  
390 v(i,i) = ONE
    g = rv1(i)
!!!    l = i
  END DO
  
!     .......... accumulation of left-hand transformations ..........
410 IF (.NOT. matu) GO TO 510
!     ..........for i=min(m,n) step -1 until 1 do -- ..........
!    mn = n
!    IF (m < n) mn = m
    
    DO i=MIN(m,n),1,-1
!      i = mn + 1 - ii
!      l = i + 1
      g = w(i)
      IF (i == n) GO TO 430
  
      DO  j = i+1, n
        u(i,j) = ZERO
      END DO
  
430 IF (g == ZERO) THEN
      u(i:m,i)=ZERO
    ELSE  
      IF (i == MIN(m,n)) GO TO 460
  
      DO  j = i+1, n
        s = ZERO
        DO  k = i+1, m
          s = s + u(k,i) * u(k,j)
        END DO
        f = (s / u(i,i)) / g   ! double division, YES! Avoid underflow
    
        DO  k = i, m
          u(k,j) = u(k,j) + f * u(k,i)
        END DO
      END DO
  
460 DO  j = i, m
      u(j,i) = u(j,i) / g
    END DO
  
    END IF  

!475 DO  j = i, m
!      u(j,i) = ZERO
!    END DO
  
    u(i,i) = u(i,i) + ONE
  END DO
  
!     .......... diagonalization of the bidiagonal form ..........
!     .......... for k=n step -1 until 1 do -- ..........
510 DO k=n,1,-1
      k1=k-1
!      k1 = n - kk
!      k = k1 + 1
      its = 0
!     .......... test for splitting.
!                for l=k step -1 until 1 do -- ..........
520 DO L=k,1,-1
!      l1 = k - ll
!      l = l1 + 1
      IF (ABS(rv1(L)) + anorm <= anorm) GO TO 565
!     .......... rv1(1) is always zero, so there is no exit
!                through the bottom of the loop ..........
      IF (ABS(w(L-1)) + anorm <= anorm) EXIT
    END DO
!     .......... cancellation of rv1(l) if l greater than 1 ..........
    c = ZERO   ! was statement 540
    s = ONE
  
    DO  i = L, k
      f = s * rv1(i)
      rv1(i) = c * rv1(i)
      IF (ABS(f) + anorm <= anorm) EXIT
      g = w(i)
      h = SQRT(f*f+g*g)
      w(i) = h
      c = g / h
      s = -f / h
      IF (matu) THEN
        DO  j = 1, m
          y = u(j,L-1)
          z = u(j,i)
          u(j,L-1) = y * c + z * s
          u(j,i) = -y * s + z * c
        END DO
      END IF
    END DO
    
!     .......... test for convergence ..........
565 z = w(k)
    IF (L == k) GO TO 650
!     .......... shift from bottom 2 by 2 minor ..........
    IF (its == 30) GO TO 1000
    its = its + 1
    x = w(L)
    y = w(k1)
    g = rv1(k1)
    h = rv1(k)
    f = ((y - z) * (y + z) + (g - h) * (g + h)) / (TWO * h * y)
    g = SQRT(f*f+ONE)
    f = ((x - z) * (x + z) + h * (y / (f + SIGN(g,f)) - h)) / x
!     .......... next qr transformation ..........
    c = ONE
    s = ONE
  
    DO  i1 = L, k1
      i = i1 + 1
      g = rv1(i)
      y = w(i)
      h = s * g
      g = c * g
      z = SQRT(f*f+h*h)
      rv1(i1) = z
      c = f / z
      s = h / z
      f = x * c + g * s
      g = -x * s + g * c
      h = y * s
      y = y * c
      IF (matv) THEN
        DO  j = 1, n
          x = v(j,i1)
          z = v(j,i)
          v(j,i1) = x * c + z * s
          v(j,i) = -x * s + z * c
        END DO
      END IF    
    z = SQRT(f*f+h*h)   ! was statement #575
    w(i1) = z

    IF (z /= ZERO) THEN
      c = f / z   ! rotation can be arbitrary if z is zero
      s = h / z
    END IF  
    f = c * g + s * y
    x = -s * g + c * y
    IF (matu) THEN
      DO  j = 1, m
        y = u(j,i1)
        z = u(j,i)
        u(j,i1) = y * c + z * s
        u(j,i) = -y * s + z * c
      END DO
    END IF
  END DO
  
  rv1(L) = ZERO
  rv1(k) = f
  w(k) = x
  GO TO 520
!     .......... convergence ..........
650 IF (z >= ZERO) CYCLE

    w(k) = -z   !  w(k) is made non-negative 
    IF (matv) v(1:n,k)=-v(1:n,k)
  
!  DO  j = 1, n
!    v(j,k) = -v(j,k)
!  END DO
  
  END DO
  RETURN
  
!  set error -- no convergence to a singular value after 30 iterations ......
1000 ierr = k
  RETURN
END SUBROUTINE SvdSingle

!+
SUBROUTINE DecompDouble(a, ipvt, errCode, cond)
! ---------------------------------------------------------------------------
! PURPOSE - Matrix triangularization by Gaussian elimination and estimation 
!  of the condition of the matrix. The matrix a is replaced by the LU 
!  decomposition of a rowwise permutation of itself. The array ipvt records
!  the row permutations of the partial pivoting. ipvt(k) is the index of the
!  kth pivot row and ipvt(n) is (-1)**(no. of interchanges). The optional
!  variable cond is an estimate of the conditioning of a. For the linear
!  system Ax=B, changes in A and B may cause changes cond times as large
!  in x. If cond+1.0 == cond, then a is singular to the working precision.
!  The determinant of a can be obtained on output as
!    det = ipvt(n)*a(1,1)*a(2,2)*a(3,3)*...*a(n,n)
! NOTES - If errCode is not zero, the remaining arguments are not to be
!   believed.

  REAL(DP),INTENT(IN OUT),DIMENSION(:,:):: a   ! matrix to be decomposed
  INTEGER,INTENT(OUT),DIMENSION(:):: ipvt   !    index of pivot rows
  INTEGER,INTENT(OUT):: errCode   ! =0 OK; =1 singular
  REAL(DP),INTENT(OUT),OPTIONAL:: cond  ! condition number

  REAL(DP):: anorm
  REAL(DP):: ek
  INTEGER:: j,k,m
  INTEGER,DIMENSION(1) :: mloc  ! receives the result of MAXLOC
  INTEGER:: n
  REAL(DP):: t
  REAL(DP),DIMENSION(SIZE(a,1)) :: work  ! scratch array
  REAL(DP):: ynorm,znorm
  REAL(DP),PARAMETER:: ZERO=0.0, ONE=1.0
!-----------------------------------------------------------------------
  errCode=0
  n=SIZE(a,1)

  IF (n <= 1) THEN
    IF (a(1,1) == ZERO) THEN
      errCode=1
    ELSE
      cond=ONE
    END IF
    RETURN        ! note abnormal RETURN
   END IF

  IF (Present(cond)) THEN
    anorm=ZERO   ! compute 1-norm of a
    DO j=1,n
      anorm=MAX(anorm, SUM(ABS(a(1:n,j))))
    END DO   
  END IF  

  DO k=1,n-1   ! Gaussian elimination with partial pivoting
    mloc= MAXLOC(ABS(a(k:n,k)))   ! pivot row
    m=k-1+mloc(1)
    ipvt(k)=m
    IF (m /= k) ipvt(n)=-ipvt(n)
    t=a(m,k)
    a(m,k)=a(k,k)
    a(k,k)=t
    IF (t /= ZERO) THEN
      t=ONE/t
      a(k+1:n,k)=-t*a(k+1:n,k)
      DO j=k+1,n   ! interchange and eliminate by columns.......
        t=a(m,j)
        a(m,j)=a(k,j)
        a(k,j)=t
        IF (t /= ZERO) a(k+1:n,j)=a(k+1:n,j) + t*a(k+1:n,k)
      END DO
    END IF
  END DO

  DO k=1,n   ! solve (a-transpose)*y=e
    t=ZERO
    IF (k > 1) t=DOT_PRODUCT(a(1:k-1,k), work(1:k-1))
    ek=ONE
    IF (t < ZERO) ek=-ONE
    IF (a(k,k) == ZERO) THEN    ! singular matrix
      errCode=1
      RETURN                            ! note abnormal RETURN
    END IF
    work(k)=-(ek+t)/a(k,k)
  END DO

  DO k=n-1,1,-1 
    t=ZERO
    t=work(k)*SUM(a(k+1:n,k))
    work(k)=t
    m=ipvt(k)
    IF (m /= k) THEN
      t=work(m)
      work(m)=work(k)
      work(k)=t
    END IF
  END DO

  IF (Present(cond)) THEN
    ynorm=SUM(ABS(work(1:n)))
    CALL Solve(a,work,ipvt)   ! solve a*z=y
    znorm=SUM(ABS(work(1:n)))
    cond=anorm*znorm/ynorm   ! estimate condition
    IF (cond < ONE) cond=ONE
  END IF
  
  RETURN
END Subroutine DecompDouble   ! ---------------------------------------------------

!+
SUBROUTINE SolveDouble(a, b, ipvt)
! ---------------------------------------------------------------------------
! PURPOSE - Solve the linear system a*x=b
  REAL(DP),INTENT(IN),DIMENSION(:,:):: a  ! the decomposed matrix from Decomp
  REAL(DP),INTENT(IN OUT),DIMENSION(:):: b ! in:right-hand side; out:solution
  INTEGER,INTENT(IN),DIMENSION(:):: ipvt  ! index of pivot rows (from Decomp)

  INTEGER:: i,k,m
  INTEGER:: n
  REAL(DP):: t
!----------------------------------------------------------------------------
  n=SIZE(a,1)
  DO k=1,n-1   ! forward elimination
    m=ipvt(k)
    t=b(m)
    b(m)=b(k)
    b(k)=t
!    DO i=k+1,n
!      b(i)=b(i)+a(i,k)*t
!    END DO
    b(k+1:n)=b(k+1:n) + t*a(k+1:n,k)
  END DO       

  DO k=n,1,-1      ! back substitution
    b(k)=b(k)/a(k,k)
    t=-b(k)
    DO i=1,k-1
      b(i)=b(i)+a(i,k)*t
    END DO       
  END DO
  RETURN
END Subroutine SolveDouble   ! ----------------------------------------------------

!+
SUBROUTINE FMMsplineDouble(x, y, b, c, d)
! ---------------------------------------------------------------------------
! PURPOSE - Compute the coefficients b,c,d for a cubic interpolating spline
!  so that the interpolated value is given by
!    s(x) = y(k) + b(k)*(x-x(k)) + c(k)*(x-x(k))**2 + d(k)*(x-x(k))**3
!      when x(k) <= x <= x(k+1)
!  The end conditions match the third derivatives of the interpolated curve to
!  the third derivatives of the unique polynomials thru the first four and
!  last four points.
!  Use Seval or Seval3 to evaluate the spline.
   REAL(DP),DIMENSION(:), INTENT(IN)  :: x ! abscissas of knots
   REAL(DP),DIMENSION(:), INTENT(IN)  :: y ! ordinates of knots
   REAL(DP),DIMENSION(:), INTENT(OUT) :: b ! linear coeff
   REAL(DP),DIMENSION(:), INTENT(OUT) :: c ! quadratic coeff.
   REAL(DP),DIMENSION(:), INTENT(OUT) :: d ! cubic coeff.

   INTEGER:: k,n
   REAL(DP):: t
   REAL(DP),PARAMETER:: ZERO=0.0, TWO=2.0, THREE=3.0
!----------------------------------------------------------------------------
  n=SIZE(x)

  IF (n < 3) THEN   ! Straight line - special case for n < 3
    b(1)=ZERO
    IF (n == 2) b(1)=(y(2)-y(1))/(x(2)-x(1))
    c(1)=ZERO
    d(1)=ZERO
    IF (n < 2) RETURN
    b(2)=b(1)
    c(2)=ZERO
    d(2)=ZERO
    RETURN
  END IF

!.....Set up tridiagonal system.........................................
!.    b=diagonal, d=offdiagonal, c=right-hand side
  d(1)=x(2)-x(1)
  c(2)=(y(2)-y(1))/d(1)
  DO k=2,n-1
    d(k)=x(k+1)-x(k)
    b(k)=TWO*(d(k-1)+d(k))
    c(k+1)=(y(k+1)-y(k))/d(k)
    c(k)=c(k+1)-c(k)
  END DO

!.....End conditions.  third derivatives at x(1) and x(n) obtained
!.       from divided differences.......................................
  b(1)=-d(1)
  b(n)=-d(n-1)
  c(1)=ZERO
  c(n)=ZERO
  IF (n > 3) THEN
    c(1)=c(3)/(x(4)-x(2))-c(2)/(x(3)-x(1))
    c(n)=c(n-1)/(x(n)-x(n-2))-c(n-2)/(x(n-1)-x(n-3))
    c(1)=c(1)*d(1)*d(1)/(x(4)-x(1))
    c(n)=-c(n)*d(n-1)*d(n-1)/(x(n)-x(n-3))
  END IF

  DO k=2,n    ! forward elimination
    t=d(k-1)/b(k-1)
    b(k)=b(k)-t*d(k-1)
    c(k)=c(k)-t*c(k-1)
  END DO

  c(n)=c(n)/b(n)   ! back substitution ( makes c the sigma of text)
  DO k=n-1,1,-1
    c(k)=(c(k)-d(k)*c(k+1))/b(k)
  END DO

!.....Compute polynomial coefficients...................................
  b(n)=(y(n)-y(n-1))/d(n-1)+d(n-1)*(c(n-1)+c(n)+c(n))
  DO k=1,n-1
    b(k)=(y(k+1)-y(k))/d(k)-d(k)*(c(k+1)+c(k)+c(k))
    d(k)=(c(k+1)-c(k))/d(k)
    c(k)=THREE*c(k)
  END DO
  c(n)=THREE*c(n)
  d(n)=d(n-1)

  RETURN
END Subroutine FMMsplineDouble   ! ---------------------------------------------------

!+
SUBROUTINE NaturalSplineDouble(x,y,b,c,d)
! ---------------------------------------------------------------------------
! PURPOSE - Construct the natural spline thru a set of points
! NOTES - A natural spline has zero second derivative at both endpoints.

  REAL(DP),INTENT(IN),DIMENSION(:):: x,y   ! coordinates of knots
  REAL(DP),INTENT(OUT),DIMENSION(:):: b,c,d  ! cubic coeff.

  INTEGER:: k,n
  REAL(DP),PARAMETER:: ZERO=0.0, TWO=2.0, THREE=3.0
!-----------------------------------------------------------------------
  n=SIZE(x)
 
  IF (n < 3) THEN   ! Straight line - special case for n < 3
    b(1)=ZERO
    IF (n == 2) b(1)=(y(2)-y(1))/(x(2)-x(1))
    c(1)=ZERO
    d(1)=ZERO
    b(2)=b(1)
    c(2)=ZERO
    d(2)=ZERO
    RETURN
  END IF

  d(1:n-1) = x(2:n)-x(1:n-1)  ! Put the h-array of the text into array d

!.....Set up the upper triangular system in locations 2 thru n-1 of
!        arrays b and c. B holds the diagonal and c the right hand side.
  b(2)=TWO*(d(1)+d(2))
  c(2)=(y(3)-y(2))/d(2)-(y(2)-y(1))/d(1)
  DO  k=3,n-1
    b(k)=TWO*(d(k-1)+d(k))-d(k-1)*d(k-1)/b(k-1)
    c(k)=(y(k+1)-y(k))/d(k)-(y(k)-y(k-1))/d(k-1)-d(k-1)*c(k-1)/b(k-1)
  END DO

  c(n-1)=c(n-1)/b(n-1)   ! Back substitute to get c-array
  DO  k=n-2,2,-1
    c(k)=(c(k)-d(k)*c(k+1))/b(k)
  END DO
  c(1)=ZERO
  c(n)=ZERO   ! c now holds the sigma array of the text 


!.....Compute polynomial coefficients ..................................
  b(n)=(y(n)-y(n-1))/d(n-1)+d(n-1)*(c(n-1)+c(n)+c(n))
  DO  k=1,n-1
    b(k)=(y(k+1)-y(k))/d(k)-d(k)*(c(k+1)+c(k)+c(k))
    d(k)=(c(k+1)-c(k))/d(k)
    c(k)=THREE*c(k)
  END DO
  c(n)=THREE*c(n)
  d(n)=d(n-1)
  RETURN
  
END Subroutine NaturalSplineDouble   ! --------------------------------------------

!+
FUNCTION SevalDouble(u, x,y, b,c,d) RESULT(SevalResult)
! ---------------------------------------------------------------------------
!  PURPOSE - Evaluate the cubic spline function
!     Seval=y(i)+b(i)!(u-x(i))+c(i)*(u-x(i))**2+d(i)*(u-x(i))**3
!           where  x(i) <= u < x(i+1)

!  NOTES- if u<x(1), i=1 is used;if u>x(n), i=n is used

  REAL(DP),INTENT(IN) :: u ! abscissa at which the spline is to be evaluated
  REAL(DP),INTENT(IN),DIMENSION(:):: x ! abscissas of knots
  REAL(DP),INTENT(IN),DIMENSION(:):: y ! ordinates of knots
  REAL(DP),INTENT(IN),DIMENSION(:):: b,c,d ! linear,quadratic,cubic coeff
  REAL(DP):: SevalResult

  INTEGER, SAVE :: i=1
  INTEGER :: j, k, n
  REAL(DP):: dx
!----------------------------------------------------------------------------
  n=SIZE(x)

!.....First check if u is in the same interval found on the
!        last call to Seval.............................................
  IF (  (i<1) .OR. (i >= n) ) i=1
  IF ( (u < x(i))  .OR.  (u >= x(i+1)) ) THEN
    i=1   ! binary search
    j=n+1

    DO
      k=(i+j)/2
      IF (u < x(k)) THEN
        j=k
      ELSE
        i=k
      END IF
      IF (j <= i+1) EXIT
    END DO
  END IF

  dx=u-x(i)   ! evaluate the spline
  SevalResult=y(i)+dx*(b(i)+dx*(c(i)+dx*d(i)))

  RETURN
END Function SevalDouble  ! -------------------------------------------------------

!+
SUBROUTINE Seval3Double(u, x,y, b,c,d, f,fp,fpp,fppp)
! ---------------------------------------------------------------------------
!  PURPOSE - Evaluate the cubic spline function
!     Seval=y(i)+b(i)!(u-x(i))+c(i)*(u-x(i))**2+d(i)*(u-x(i))**3
!           where  x(i) <= u < x(i+1)

!  NOTES- if u<x(1), i=1 is used;if u>x(n), i=n is used

  REAL(DP),INTENT(IN) :: u ! abscissa at which the spline is to be evaluated
  REAL(DP),INTENT(IN),DIMENSION(:):: x ! abscissas of knots
  REAL(DP),INTENT(IN),DIMENSION(:):: y ! ordinates of knots
  REAL(DP),INTENT(IN),DIMENSION(:):: b,c,d ! linear,quadratic,cubic coeff
  REAL(DP),INTENT(OUT),OPTIONAL:: f,fp,fpp,fppp ! function, 1st,2nd,3rd deriv

  INTEGER, SAVE :: i=1
  INTEGER :: j, k, n
  REAL(DP)    :: dx
  REAL(DP),PARAMETER:: TWO=2.0, THREE=3.0, SIX=6.0
!----------------------------------------------------------------------------
  n=SIZE(x)

!.....First check if u is in the same interval found on the
!        last call to Seval.............................................
  IF (  (i<1) .OR. (i >= n) ) i=1
  IF ( (u < x(i))  .OR.  (u >= x(i+1)) ) THEN
    i=1   ! binary search
    j=n+1

    DO
      k=(i+j)/2
      IF (u < x(k)) THEN
        j=k
      ELSE
        i=k
      END IF
      IF (j <= i+1) EXIT
    END DO
  END IF

  dx=u-x(i)   ! evaluate the spline
  IF (Present(f))    f=y(i)+dx*(b(i)+dx*(c(i)+dx*d(i)))
  IF (Present(fp))   fp=b(i)+dx*(TWO*c(i) + dx*THREE*d(i))
  IF (Present(fpp))  fpp=TWO*c(i) + dx*SIX*d(i)
  IF (Present(fppp)) fppp=SIX*d(i)

  RETURN
END Subroutine Seval3Double  ! -------------------------------------------------------

!+
SUBROUTINE Quanc8Double (Fun,a,b,abserr,relerr,result,errest,nofun,flag)
! -----------------------------------------------------------------------------
! PURPOSE - Estimate the integral of FUN(X) from A to B to a user provided 
!  tolerance. This is an automatic adaptive routine based on the 8-panel 
!  Newton-Cotes rule. 

  INTERFACE
    FUNCTION Fun(x) RESULT(f)  ! integrand function subprogram
    IMPLICIT NONE
      INTEGER,PARAMETER:: DP=SELECTED_REAL_KIND(10,50)
      REAL(DP),INTENT(IN)::x
      REAL(DP):: f
     END FUNCTION Fun
   END INTERFACE
   
   REAL(DP),INTENT(IN):: a        !  lower limit of integration
   REAL(DP),INTENT(IN):: b        !  upper limit of integration
   REAL(DP),INTENT(IN):: abserr   !  absolute error tolerance (>0)
   REAL(DP),INTENT(IN):: relerr   !  relative error tolerance (>0)

   REAL(DP),INTENT(OUT):: result   !  approx. value of integral
   REAL(DP),INTENT(OUT):: errest   !  estimate of actual error
   INTEGER,INTENT(OUT):: nofun     !  no. of function evaluations
   REAL(DP),INTENT(OUT):: flag     !  reliability indicator. (=0 is O.K.)

   REAL(DP),PARAMETER:: W0=3956.0_DP/14175.0_DP
   REAL(DP),PARAMETER:: W1=23552.0_DP/14175.0_DP
   REAL(DP),PARAMETER:: W2=-3712.0_DP/14175.0_DP
   REAL(DP),PARAMETER:: W3=41984.0_DP/14175.0_DP
   REAL(DP),PARAMETER:: W4=-18160.0_DP/14175.0_DP

  REAL(DP):: area
  REAL(DP):: cor11  
  REAL(DP):: esterr
  REAL(DP), DIMENSION(0:16):: f,x
  REAL(DP), DIMENSION(8,30):: fsave,xsave
  INTEGER:: i,j
  INTEGER:: lev,nim
  INTEGER:: LEVMAX = 30
  INTEGER,PARAMETER :: LEVMIN = 1
  INTEGER,PARAMETER:: LEVOUT = 6
  INTEGER:: NOFIN     !...Trouble if NOFUN reaches NOFIN
  INTEGER,PARAMETER:: NOMAX = 5000
  REAL(DP):: qprev,qnow,qdiff,qleft
  REAL(DP):: stone,step
  REAL(DP):: tolerr
  REAL(DP),DIMENSION(31):: qright
  REAL(DP),PARAMETER:: ZERO=0.0, HALF=0.5, ONE=1.0, SIXTEEN=16.0
!------------------------------------------------------------------------------

!..... S t a g e   1  (General initialization)..........................
   NOFIN=NOMAX-8*(LEVMAX-LEVOUT+2**(LEVOUT+1))

!...initialize running sums to zero...
  area=ZERO
  cor11=ZERO
  flag=ZERO   ! flag is REAL
  result=ZERO
  errest=ZERO
  nofun=0
  IF (a==b) RETURN
  
!..... S T A G E   2   (initialization for first interval)..............
  lev=0
  nim=1
  x(0)=a
  x(16)=b
  qprev=ZERO
  stone=(b-a)/SIXTEEN
  x(8)= HALF*(x(0)+x(16))
  x(4)= HALF*(x(0)+x(8))
  x(12)=HALF*(x(8)+x(16))
  x(2)= HALF*(x(0)+x(4))
  x(6)= HALF*(x(4)+x(8))
  x(10)=HALF*(x(8)+x(12))
  x(14)=HALF*(x(12)+x(16))

  DO j=0,16,2
    f(j)=Fun(x(j))
  END DO
  nofun=9

!..... S T A G E   3   (central calculation)............................
!...requires qprev,x0,x2,x4,....x(16),f0,f2,f4,....f16.
!...calculates x1,x3,....x15,f1,f3,....f15,qleft,qright,qnow,qdiff,area

30 DO j=1,15,2  ! keeps coming back here until successful
      x(j)=HALF*(x(j-1)+x(j+1))
      f(j)=Fun(x(j))
    END DO
    nofun=nofun+8
    step=(x(16)-x(0))/SIXTEEN
    qleft=(w0*(f(0)+f(8))+w1*(f(1)+f(7))+w2*(f(2)+f(6))+   &
                           w3*(f(3)+f(5))+w4*f(4))*step
    qright(lev+1)=(w0*(f(8)+f(16))+w1*(f(9)+f(15))+w2*(f(10)+f(14))+   &
                    w3*(f(11)+f(13))+w4*f(12))*step
    qnow=qleft+qright(lev+1)
    qdiff=qnow-qprev
    area=area+qdiff

!..... S T A G E   4   (interval convergence test) .....................
    esterr=ABS(qdiff)/1023.0
    tolerr=MAX(abserr,relerr*ABS(area))*(step/stone)
    IF (lev < LEVMIN) GO TO 50
    
    
    IF (lev >= LEVMAX) THEN
      flag=flag+ONE
      GO TO 70
    END IF  
    IF (nofun > NOFIN) THEN      !..... S T A G E   6  (trouble section) 
      nofin=nofin+nofin  ! Number of function values is about to exceed limit
      levmax=levout
      flag=flag+(b-x(0))/(b-a)
      GO TO 70
    END IF
    IF (esterr <= TOLERR) GO TO 70

!..... S T A G E   5   (no convergence).................................
!...locate next interval...
50  nim=nim+nim
    lev=lev+1

    fsave(1:8,lev)=f(9:16)   ! store right hand elements for future use
    xsave(1:8,lev)=x(9:16)
!    DO i=1,8
!      fsave(i,lev)=f(i+8)
!      xsave(i,lev)=x(i+8)
!    END DO

    qprev=qleft
    DO I=1,8   ! assemble left hand elements for immediate use
      f(18-2*i)=f(9-i)
      x(18-2*i)=x(9-i)
    END DO
    GO TO 30


!..... S T A G E   7   (integral converged) ............................
!...add contributions into running sums...
70 result=result+qnow
   errest=errest+esterr
   cor11=cor11+qdiff/1023.0

  DO                                ! 72 loop in text 
    IF (nim == 2*(nim/2)) EXIT   
    nim=nim/2
    lev=lev-1
  END DO
  nim=nim+1
  
  IF (lev <= 0) GO TO 80   ! looks like success!
  qprev=qright(lev)   ! assemble elements required for next interval
  x(0)=x(16)
  f(0)=f(16)
  DO i=1,8
    f(2*i)=fsave(i,lev)
    x(2*i)=xsave(i,lev)
  END DO
  GO TO 30

!..... S T A G E   8  (Finalize and return).............................
80 result=result+cor11
  IF (errest /= ZERO) THEN
    DO  
      IF (ABS(result)+errest > ABS(RESULT) ) EXIT
      errest=errest+errest  ! make sure ERREST is not < roundoff level ....
    END DO
  END IF

  RETURN
END Subroutine QuanC8Double   ! ---------------------------------------------------

!+
SUBROUTINE Rkf45Double (F,y,t,tout,relerr,abserr,iflag,work,iwork)
! ---------------------------------------------------------------------------
! PURPOSE - Integrate a system of neqn first order differential 
!  equations by the Fehlberg fourth-fifth order Runge-Kutta method.
! NOTES - RKF45 is primarily designed to solve non-stiff and mildly stiff
!   differential equations when derivative evaluations are inexpensive.
!   Rkf45 should generally not be used when the user is demanding 
!	high accuracy.

! ABSTRACT - Subroutine Rkf45 integrates a system of neqn first order
!  ordinary differential equations of the form
!	dy(i)/dt = f(t,y(1),y(2),...,y(neqn))
!	where the y(i) are given at t .
! Typically the subroutine is used to integrate from t to tout but it can
! be used as a one-step integrator to advance the solution a Double step
! in the direction of tout. On return, the parameters in the call list are
! set for continuing the integration. The user has only to call Rkf45 again
! (and perhaps define a new value for tout). Actually, Rkf45 is an 
! interfacing routine which calls subroutine Rkfs for the solution. 
! Rkfs in turn calls subroutine Fehl which computes an approximate solution
! over one step. 

! Rkf45 uses the Runge-Kutta-Fehlberg(4,5) method described in the reference
! E. Fehlberg:  Low-order Classical Runge-Kutta formulas with Stepsize 
! Control, NASA TR R-315.

! The performance of Rkf45 is illustrated in the reference
! L.F.Shampine,H.A.Watts,S.Davenport: "Solving Non-stiff Ordinary
! Differential Equations - The State of the Art 
!  Sandia Laboratories Report SAND75-0182, to appear in SIAM Review. 

! First call to Rkf45:

! The user must provide storage in his calling program for the arrays in 
! the call list	y(neqn), work(3+6*neqn), iwork(5), supply 
! subroutine F(t,y,yp) and initialize the following parameters-

! y -- vector of initial conditions
! t -- starting point of integration , must be a variable
! tout -- output point at which solution is desired.
!   t=tout is allowed on the first call only, in which case
!   Rkf45 returns with iflag=2 if continuation is possible.
! relerr,abserr -- relative and absolute local error tolerances
!   which must be non-negative. relerr must be a variable while
!   abserr may be a constant. The code should normally not be
!   used with relative error control smaller than about i.e-8 .
!   to avoid limiting precision difficulties the code requires
!   relerr to be larger than an internally computed relative
!   error parameter which is machine dependent. In particular,
!   pure absolute error is not permitted. If a smaller than
!   allowable value of relerr is attempted, rkf45 increases
!   relerr appropriately and returns control to the user before
!   continuing the integration.
! iflag -- +1,-1 indicator to initialize the code for each new problem. 
!   Normal input is +1. the user should set iflag=-1 only when one-step
!   integrator control is essential. In this case, Rkf45 attempts to 
!   advance the solution a Double step in the direction of tout each time
!   it is called. Since this mode of operation results in extra computing 
!   overhead, it should be avoided unless needed. 

! Output from Rkf45:

! y -- solution at t
! t - last point reached in integration.
! iflag = 2 -- integration reached tout. Indicates successful return and is 
!              the normal mode for continuing integration.
!	= -2 -- a Double successful step in the direction of tout has been 
!               taken. Normal mode for continuing integration one step 
!               at a time.
!	= 3 -- integration was not completed because relative error tolerance 
!              was too small. relerr has been increased	appropriately for 
!              continuing.
!	= 4 -- integration was not completed because more than 3000 derivative
!              evaluations were needed. This is approximately 500 steps.
!	= 5 -- integration was not completed because solution vanished making
!              a pure relative error test impossible. Must use non-zero abserr
!              to continue. Using the one-step integration mode for one step
!              is a good way to proceed.
!	= 6 -- integration was not completed because requested accuracy could 
!              not be achieved using smallest allowable stepsize. User must
!              increase the error tolerance before continued integration can 
!              be attempted.
!	= 7 -- it is likely that Rkf45 is inefficient for solving this problem.
!              Too much output is restricting the natural stepsize choice. 
!              Use the one-step integrator mode.
!	= 8 -- invalid input parameters. This indicator occurs if any of the 
!              following is satisfied -	SIZE(y) <= 0; 
!              t=tout and iflag /= +1 or -1
!	       relerr or abserr < 0.
!              iflag == 0 or < -2 or > 8
! work,iwork -- information which is usually of no interest to the user but 
!    necessary for subsequent calls.
!    work(1),...,work(neqn) contain the first derivatives of the solution 
!    vector y at t. work(neqn+1) contains the stepsize h to be attempted on 
!    the next step. iwork(1) contains the derivative evaluation counter. 

! Subsequent calls to Rkf45 

! Subroutine Rkf45 returns with all information needed to continue the 
! integration. If the integration reached tout, the user need only define
! a new tout and call Rkf45 again. In the one-step integrator mode (iflag=-2)
! the user must keep in mind that each step taken is in the direction of the
! current tout. Upon reaching tout (indicated by changing iflag to 2),the 
! user must then define a new tout and reset iflag to -2 to continue in the 
! one-step integrator mode. 

! If the integration was not completed but the user still wants to continue 
! (iflag=3,4 cases), he just calls Rkf45 again. With iflag=3, the relerr 
! parameter has been adjusted appropriately for continuing the integration.
! In the case of iflag=4 the function counter will be reset to 0 and 
! another 3000 function evaluations are allowed. 

! However,in the case iflag=5, the user must first alter the error criterion
! to use a positive value of abserr before integration can proceed. If he 
! does not, execution is terminated. 

! Also, in the case iflag=6, it is necessary for the user to reset iflag 
! to 2 (or -2 when the one-step integration mode is being used) as well as 
! increasing either abserr,relerr or both before the integration can be 
! continued. If this is not done, execution will be terminated. The 
! occurrence of iflag=6 indicates a trouble spot (solution is changing 
! rapidly,singularity may be present) and it often is inadvisable to continue.

! If (flag=7 is encountered, the user should use the one-step integration mode
! with the stepsize determined by the code or consider switching to the 
! Adams codes DE/STEP,INTRP. If the user insists upon continuing the
! integration with Rkf45, he must reset iflag to 2 before calling Rkf45 again.
! Otherwise, execution will be terminated.

! If iflag=8 is obtained, integration can not be continued unless the invalid
! input parameters are corrected. 

! It should be noted that the arrays work,iwork contain information required
! for subsequent integration. Accordingly, work and iwork should not be 
! altered.

! This interfacing routine merely relieves the user of a long calling list 
! via the splitting apart of two working storage arrays.

INTERFACE
  SUBROUTINE F(arg1,arg2,arg3)
  IMPLICIT NONE
  INTEGER,PARAMETER:: DP = SELECTED_REAL_KIND(10,50)
  REAL(DP),INTENT(IN):: arg1   ! independent variable
  REAL(DP),INTENT(IN),DIMENSION(:):: arg2   ! vector of dependent variables
  REAL(DP),INTENT(OUT),DIMENSION(:):: arg3   ! vector of derivatives
  END Subroutine F
END INTERFACE  


  REAL(DP), INTENT(IN OUT),DIMENSION(:):: y   ! solution vector at t
  REAL(DP), INTENT(IN OUT):: t   ! independent variable
  REAL(DP), INTENT(IN OUT):: tout   ! output point at which solution is desired
  REAL(DP), INTENT(IN OUT):: relerr   ! relative error tolerance
  REAL(DP), INTENT(IN):: abserr   !  absolute error tolerance
  INTEGER, INTENT(IN OUT):: iflag   !   indicator for status of work
  REAL(DP), INTENT(IN OUT),DIMENSION(:):: work   ! work array
  INTEGER, INTENT(IN OUT),DIMENSION(:):: iwork   ! work array

  INTEGER :: k1,k2,k3,k4,k5,k6,k1m
  INTEGER:: neqn
!----------------------------------------------------------------------------
  neqn=SIZE(y)
  k1m=neqn+1   !.....compute indices for the splitting of the work array
  k1=k1m+1
  k2=k1+neqn
  k3=k2+neqn
  k4=k3+neqn
  k5=k4+neqn
  k6=k5+neqn
  
  CALL RkfsDouble(F,neqn,y,t,tout,relerr,abserr,iflag,work(1:neqn),work(neqn+1),  &
    work(k1:k1+neqn-1),work(k2:k2+neqn-1),work(k3:k3+neqn-1),               &
    work(k4:k4+neqn-1), work(k5:k5+neqn-1),work(k6), work(k6+1),            &
    iwork(1),iwork(2),iwork(3),iwork(4),iwork(5))
  RETURN
END Subroutine Rkf45Double   ! ----------------------------------------------------

!+
SUBROUTINE RkfsDouble (F,neqn,y,t,tout,relerr,abserr,iflag,yp,                    &
    h,f1,f2,f3,f4,f5,savre,savae,nfe,kop,init,jflag,kflag)
! ---------------------------------------------------------------------------
! PURPOSE - integrate a system of neqn first order differential equations by 
!   the Fehlberg fourth-fifth order Runge-Kutta method.

! NOTES - Rkfs integrates a system of first order ordinary differential
!  equations as described in the comments for Rkf45. The arrays 
!  yp,f1,f2,f3,f4,and f5 (of dimension at least neqn) and the variables 
!  h,savre,savae,nfe,kop,init,jflag,and kflag are used internally by the 
!  code and appear in the call list to eliminate local retention of variables
!  between calls. Accordingly, they should not be altered. Items of possible 
!  interest are
!   yp - derivative of solution vector at t
!   h - an appropriate stepsize to be used for the next step
!   nfe- counter on the number of derivative function evaluations 

  INTEGER, INTENT(IN):: neqn
  REAL(DP), INTENT(IN OUT),DIMENSION(:):: y   ! solution vector at t
  REAL(DP), INTENT(IN OUT):: t   ! independent variable
  REAL(DP), INTENT(IN):: tout   ! output point at which solution is desired
  REAL(DP), INTENT(IN OUT):: relerr   ! relative error tolerance
  REAL(DP), INTENT(IN):: abserr   !  absolute error tolerance
  INTEGER, INTENT(IN OUT):: iflag  !   indicator for status of work
  REAL(DP), INTENT(IN OUT),DIMENSION(:):: yp   ! derivatives 
  REAL(DP), INTENT(OUT):: h                    ! step size
  REAL(DP), INTENT(IN OUT),DIMENSION(:):: f1
  REAL(DP), INTENT(IN OUT),DIMENSION(:):: f2
  REAL(DP), INTENT(IN OUT),DIMENSION(:):: f3
  REAL(DP), INTENT(IN OUT),DIMENSION(:):: f4
  REAL(DP), INTENT(IN OUT),DIMENSION(:):: f5
  REAL(DP), INTENT(OUT):: savre
  REAL(DP), INTENT(OUT):: savae
  INTEGER, INTENT(OUT):: nfe    ! number of function evaluations
  INTEGER, INTENT(OUT):: kop
  INTEGER, INTENT(IN OUT):: init
  INTEGER, INTENT(IN OUT):: jflag
  INTEGER, INTENT(IN OUT):: kflag
  LOGICAL :: hfaild,output

INTERFACE
  SUBROUTINE F(arg1,arg2,arg3)
  IMPLICIT NONE
  INTEGER,PARAMETER:: DP = SELECTED_REAL_KIND(10,50)
  REAL(DP),INTENT(IN):: arg1
  REAL(DP),INTENT(IN),DIMENSION(:):: arg2
  REAL(DP),INTENT(OUT),DIMENSION(:):: arg3
  END Subroutine F
END INTERFACE  

  REAL(DP):: a,ae,dt,ee,eeoet,esttol,et
  REAL(DP):: eps
  REAL(DP):: hmin
  INTEGER:: k,mflag

  REAL(DP),PARAMETER:: REMIN = 1E-12   ! minimum acceptable value of relerr.  
    ! Attempts to obtain higher accuracy with this routine are usually very
    ! expensive and unsuccessful.

  INTEGER,PARAMETER:: MAXNFE = 3000 ! the expense is controlled by restricting
    ! the number of function evaluations to be approximately MAXNFE. 
    ! As set, this corresponds to about 500 steps.
  REAL(DP):: rer,s 
  REAL(DP):: scale                                   
  REAL(DP):: tol,toln,u26
  REAL(DP):: ypk
  REAL(DP),PARAMETER:: ZERO=0.0, ONE=1.0
!-----------------------------------------------------------------------
  eps=EPSILON(ZERO)
  u26=26.0*eps

!.....check input parameters............................................
  IF (neqn < 1 .OR. relerr < ZERO .OR. abserr < ZERO) THEN
    iflag=8
    RETURN
  END IF  
  mflag=ABS(iflag)
  IF (mflag == 0 .OR. mflag > 8) THEN
    iflag=8
    RETURN
  END IF
  
  IF (mflag == 1) GO TO 50

!.....check continuation posibilities...................................
  IF ((t == tout) .AND. (kflag /= 3)) THEN
    iflag=8   ! invalid flag handler
    RETURN
  END IF  
  IF (mflag /= 2) GO TO 25

!...you get here if iflag=+2 or -2
  IF ((kflag == 3) .OR. (init == 0)) GO TO 45
  IF (kflag == 4) GO TO 40
  IF ((kflag == 5) .AND. (abserr == ZERO)) GO TO 30
  IF ((kflag == 6) .AND. (relerr <= savre) .AND. (abserr <= savae))  GO TO 30
  GO TO 50


!...you get here if iflag=3,4,5,6,7, or 8.....
25 IF (iflag == 3) GO TO 45
   IF (iflag == 4) GO TO 40
   IF ((iflag == 5) .AND. (abserr > ZERO)) GO TO 45


!.....integration cannot be continued since user did not respond to the
!..... instructions pertaining to iflag=5,6,7, or 8
30  STOP   !!! maybe this is too strong


!.....reset function evaluation counter.................................
40 nfe=0
   IF (mflag == 2) GO TO 50


!.....reset flag value from previous call...............................
45 iflag=jflag
   IF (kflag == 3) mflag=ABS(iflag)

50 jflag=iflag   ! save for subsequent calls
   kflag=0       ! set the continuation flag
  savre=relerr   ! save for subsequent calls
  savae=abserr   ! save for subsequent calls


!.....restrict relative error tolerance to be at least as large as
!       2*eps+remin to avoid limiting precision difficulties arising
!       from impossible accuracy requests.
  rer=eps+eps+remin
  IF (relerr >= rer) GO TO 55


!.....relative error tolerance too small................................
  relerr=rer
  iflag=3
  kflag=3
  RETURN

55 dt=tout-t
   IF (mflag == 1) GO TO 60
   IF (init == 0) GO TO 65
   GO TO 80

!.....initialization....................................................
60 init=0
   kop=0
   a=t
   CALL f(a,y,yp)
   nfe=1
   IF (t /= tout) GO TO 65
   iflag=2
   RETURN

65 init=1
   h=ABS(dt)
   toln=ZERO
   DO  k=1,neqn
     tol=relerr*ABS(y(k))+abserr
     IF (tol <= ZERO) CYCLE
     toln=tol
     ypk=ABS(yp(k))
     IF (ypk*h**5 > tol) h=(tol/ypk)**0.2
   END DO

  IF (toln <= ZERO) h=ZERO
  h=MAX(h, u26*MAX(ABS(t), ABS(dt)))
  jflag=SIGN(2,iflag)

!.....set stepsize for integration in the direction from t to tout
80 h=SIGN(h,dt)

!.....test to see if rkf45 is being severely impacted by too many outputs
  IF (ABS(h) >= ABS(dt+dt)) kop=kop+1
  IF (kop /= 100) GO TO 85
  kop=0
  iflag=7
  RETURN

85 IF(ABS(dt) > u26*ABS(t)) GO TO 95

!.....if too close to output point, extrapolate and return
  DO  k=1,neqn
    y(k)=y(k)+dt*yp(k)
  END DO
  a=tout
  CALL F(a,y,yp)
  nfe=nfe+1
  GO TO 300


!.....initialize output point indicator.................................
95 output=.FALSE.


!.....scale the error tolerances to avoid premature underflow in the
!        error tolerance
  scale=2.0/relerr
  ae=scale*abserr


!.....step by step integration..........................................
100 hfaild=.false.
!...set smallest allowable step size...
  hmin=u26*ABS(t)

!.....adjust stepsize if necessary to hit the output point.
!.....look ahead two steps to avoid drastic changes in the stepsize and
!.....thus lessen the impact of output points on the code.
  dt=tout-t
  IF (ABS(dt) >= ABS(h+h)) GO TO 200
  IF (ABS(dt) > ABS(h)) GO TO 150


!.....the next successful step will complete the integration to the
!     output point
  output=.TRUE.
  h=dt
  GO TO 200

150 h=0.5*dt

!.....core integrator for taking a Double step..........................

200 IF (nfe <= maxnfe) GO TO 220

!.....too much work...
  iflag=4
  kflag=4
  RETURN


!.....advance an approximate solution over one step of length h
220 CALL FehlDouble(f,y,t,h,yp,f1,f2,f3,f4,f5,f1)
    nfe=nfe+5


!.....compute and test allowable tolerances versus local error estimates
!.....  and remove scaling of tolerances. note that relative error is
!.....  is measured with respect to the average of the magnitudes of the
!.....   solution at the beginning and end of the step.
  eeoet=ZERO
  DO  k=1,neqn
    et=ABS(y(k))+ABS(f1(k))+ae
    IF (et > ZERO) GO TO 240
  
  
!.....Inappropriate error tolerance.....................................
    iflag=5
    RETURN
  
240 ee=ABS((-2090.0*yp(k)+(21970.0*f3(k)-15048.0*f4(k)))+               &
      (22528.0*f2(k)-27360.0*f5(k)))
    eeoet=MAX(eeoet,ee/et)
  END DO
  esttol=ABS(h)*eeoet*scale/752400.0
  IF (esttol <= ONE) GO TO 260


!.....unsuccessful step. reduce the stepsize and start again............
!     (decrease is limited to a factor of 1/10)
  hfaild=.TRUE.
  output=.FALSE.
  s=0.1
  IF (esttol < 59049.0) s=0.9/esttol**0.2
  h=s*h
  IF (ABS(h) > hmin) GO TO 200


!.....requested error unobtainable at smallest allowable stepsize.......
  iflag=6
  kflag=6
  RETURN


!.....successful step.   store solution at t+h and evaluate derivatives
!     there..........................
260 t=t+h
    DO  k=1,neqn
      y(k)=f1(k)
    END DO
    a=t
    CALL F(a,y,yp)
    nfe=nfe+1


! Choose next stepsize.  The increase is limited to a factor of 5.
! If step failure has just occured, next stepsize is not allowed
!  to increase.
  s=5.0
  IF (esttol > 1.889568E-4) s=0.9/esttol**0.2
  IF (hfaild) s=MIN(s,ONE)
  h=SIGN(MAX(s*ABS(h),hmin),h)

!.....E N D   O F   C O R E   I N T E G R A T O R   ....................


!.....Should we take another step? .....................................
  IF (output) GO TO 300
  IF (iflag > 0) GO TO 100


!.....integration successfully completed................................
!...one step mode...
  iflag=-2
  RETURN

!...interval mode...
300 t=tout
    iflag=2
    RETURN
END SUBROUTINE RkfsDouble   ! -----------------------------------------------------

!+
SUBROUTINE FehlDouble (F,y,t,h,yp,f1,f2,f3,f4,f5,s)
! ---------------------------------------------------------------------------
! PURPOSE - Integrate a system of neqn first-order ordinary differential
!  equations of the form dy/dt=F(t,y) [y is a vector], where the initial
!  values of y and the initial values of yp, the derivatives are specified
!  at the starting point t.  Fehl advances the solution over the fixed step h
!  and returns the fifth order (sixth-order locally) solution approximation 
!  at t+h in array s.
INTERFACE
  SUBROUTINE F(arg1,arg2,arg3)
  IMPLICIT NONE
  INTEGER,PARAMETER:: DP = SELECTED_REAL_KIND(10,50)
  REAL(DP),INTENT(IN):: arg1
  REAL(DP),INTENT(IN),DIMENSION(:):: arg2
  REAL(DP),INTENT(OUT),DIMENSION(:):: arg3
  END Subroutine F
END INTERFACE  

!!!  INTEGER, INTENT(IN):: neqn ! number of equations
  REAL(DP), INTENT(IN),DIMENSION(:):: y ! array of length neqn; function at t
  REAL(DP), INTENT(IN):: t   ! starting point
  REAL(DP), INTENT(IN):: h   ! step size
  REAL(DP), INTENT(IN),DIMENSION(:)::  yp ! array of length neqn; derivatives at t
  REAL(DP), INTENT(OUT),DIMENSION(:):: f1 ! array of length neqn for internal use
  REAL(DP), INTENT(OUT),DIMENSION(:):: f2 ! array of length neqn for internal use
  REAL(DP), INTENT(OUT),DIMENSION(:):: f3 ! array of length neqn for internal use
  REAL(DP), INTENT(OUT),DIMENSION(:):: f4 ! array of length neqn for internal use
  REAL(DP), INTENT(OUT),DIMENSION(:):: f5 ! array of length neqn for internal use
  REAL(DP), INTENT(OUT),DIMENSION(:):: s  ! array of length neqn; the results

  REAL(DP),PARAMETER:: C1 = 0.25_DP
  REAL(DP),PARAMETER:: C2 = 3.0_DP/32.0_DP
  REAL(DP),PARAMETER:: C3 = 3.0_DP
  REAL(DP),PARAMETER:: C4 = 3.0_DP/8.0_DP
  REAL(DP),PARAMETER:: C5 = 1.0_DP/2197.0_DP
  REAL(DP),PARAMETER:: C6 = 1932.0_DP
  REAL(DP),PARAMETER:: C7 = 7296.0_DP
  REAL(DP),PARAMETER:: C8 = -7200.0_DP
  REAL(DP),PARAMETER:: C9 = 12.0_DP/13.0_DP
  REAL(DP),PARAMETER:: C10 = 1.0_DP/4104.0_DP
  REAL(DP),PARAMETER:: C11 = 8341.0_DP
  REAL(DP),PARAMETER:: C12 = -845.0_DP
  REAL(DP),PARAMETER:: C13 = 29440.0_DP
  REAL(DP),PARAMETER:: C14 = -32832.0_DP
  REAL(DP),PARAMETER:: C15 = 1.0_DP/20520.0_DP
  REAL(DP),PARAMETER:: C16 = -6080.0_DP
  REAL(DP),PARAMETER:: C17 = 9295.0_DP
  REAL(DP),PARAMETER:: C18 = -5643.0_DP
  REAL(DP),PARAMETER:: C19 = 41040.0_DP
  REAL(DP),PARAMETER:: C20 = -28352.0_DP
  REAL(DP),PARAMETER:: C21 = 0.5_DP
  REAL(DP),PARAMETER:: C22 = 1.0_DP/7618050.0_DP
  REAL(DP),PARAMETER:: C23 = 902880.0_DP
  REAL(DP),PARAMETER:: C24 = 3855735.0_DP
  REAL(DP),PARAMETER:: C25 = -1371249.0_DP
  REAL(DP),PARAMETER:: C26 = 3953664.0_DP
  REAL(DP),PARAMETER:: C27 = 277020.0_DP

  REAL(DP):: ch
!-----------------------------------------------------------------------
  ch=C1*h  
  f5=y+ch*yp  
  CALL F(t+ch, f5,f1)

  ch=C2*h  
  f5=y + ch*(yp+C3*f1)
  CALL F(t+C4*h, f5,f2)  

  ch=C5*h
  f5=y + ch*(C6*yp + (C7*f2+C8*f1) )
  CALL F(t+C9*h, f5,f3)  

  ch=C10*h
  f5=y + ch*( (C11*yp + C12*f3) + (C13*f2+C14*f1) )
  CALL F(t+h, f5,f4)

  ch=C15*h
  f1=y + ch*( C16*yp + (C17*f3+C18*f4) + (C19*f1+C20*f2) )  
  CALL F(t+C21*h, f1, f5)

!.....compute approximate solution at t+h
  ch=C22*h
  s=y + ch*( C23*yp + (C24*f3+C25*f4) + (C26*f2+C27*f5) )  
  RETURN
END SUBROUTINE FehlDouble   ! -----------------------------------------------------

!+
FUNCTION ZeroinDouble (a,b,F,tol) RESULT(z)
! ---------------------------------------------------------------------------
! PURPOSE - Compute a zero of f in the interval (a,b)

  REAL(DP),INTENT(IN):: a,b   ! left and right endpoints of interval
  REAL(DP),INTENT(IN):: tol     ! desired interval of uncertainity 
      
  REAL(DP):: z
INTERFACE
  FUNCTION F(x) RESULT(g)
    IMPLICIT NONE
    INTEGER,PARAMETER:: DP = SELECTED_REAL_KIND(10,50)
    REAL(DP),INTENT(IN):: x
    REAL(DP):: g
  END Function F
END INTERFACE    
  
  INTEGER:: errCode
  INTEGER,PARAMETER:: MAXITER=25
  INTEGER:: neval   ! not used
  REAL(DP):: xZero,fZero
!----------------------------------------------------------------------------
  CALL BrentZeroDouble(a,b,F,tol,MAXITER,neval,errCode,xZero,fZero)
  z=xZero
  RETURN
END Function ZeroinDouble   ! -----------------------------------------------------

!+
SUBROUTINE BrentZeroDouble(ax,bx,F,tol,maxIter,neval,errCode,xZero,fZero)
! ---------------------------------------------------------------------------
! PURPOSE - Compute a zero of F in the interval (ax,bx)

  REAL(DP),INTENT(IN):: ax,bx   ! left and right endpoints of interval
  REAL(DP),INTENT(IN):: tol     ! desired interval of uncertainity 
  INTEGER,INTENT(IN):: maxIter   ! max number of iterations allowed. 25 is good
  INTEGER,INTENT(OUT):: neval
  INTEGER,INTENT(OUT):: errCode   ! =0 is OK; =1 too many iterations
                                  ! =2 if F(ax) and F(bx) have the same sign
  REAL(DP),INTENT(OUT):: xZero,fZero ! the last and best value of the zero                                   
      
INTERFACE
  FUNCTION F(x) RESULT(g)
    IMPLICIT NONE
    INTEGER,PARAMETER:: DP = SELECTED_REAL_KIND(10,50)
    REAL(DP),INTENT(IN):: x
    REAL(DP):: g
  END Function F
END INTERFACE    

  REAL(DP):: a,b,c,d,e,eps
  REAL(DP):: fa,fb,fc,tol1
  INTEGER:: kIter
!  INTEGER:: method   ! =0 bisection; =1 linear; =2 inverse quadratic
  REAL(DP):: xm,p,q,r,s
  REAL(DP),PARAMETER:: ZERO=0.0, ONE=1.0, TWO=2.0, THREE=3.0, HALF=0.5
!----------------------------------------------------------------------------
  eps=EPSILON(ax)
  tol1=ONE+eps

  a=ax   ! initialization
  b=bx
  fa=f(a)
  fb=f(b)
  neval=2
  IF (fa==ZERO) THEN
    xZero=a
    fZero=ZERO
    errCode=0
    RETURN
  END IF  
  IF (fb==ZERO) THEN
    xZero=b
    fZero=ZERO
    errCode=0
    RETURN
  END IF  
  IF (fa*fb > ZERO) THEN
    xZero=HUGE(a)
    fZero=HUGE(a)
    errCode=2
    RETURN
  END IF
! The trivial cases have now been dealt with. On to the real work...

  c=a
  fc=fa
  d=b-a
  e=d
    
  DO kIter=1,maxIter
    IF ( (fb>0 .AND. fc>0) .OR. (fb<0 .AND. fc<0) ) THEN
      c=a  ! we insist that b and c straddle the zero
      fc=fa
      d=b-a
      e=d
    END IF

    IF (ABS(fc) < ABS(fb)) THEN
      a=b    ! we insist that b be the better guess of b and c
      b=c
      c=a
      fa=fb
      fb=fc
      fc=fa
    END IF

    tol1=TWO*eps*ABS(b)+HALF*tol   ! convergence test
    xm=HALF*(c-b)
    IF (ABS(xm) <= tol1 .OR. fb==ZERO) THEN
      xZero=b
      fZero=fb
      errCode=0   ! SUCCESS! The proper way to leave
      RETURN
    END IF

!    WRITE(DBG,*) "BrentZero, start of kIter=",kIter
!    WRITE(DBG,'(A,2ES25.16)') "a,fa=", a,fa
!    WRITE(DBG,'(A,2ES25.16)') "b,fb=", b,fb
!    WRITE(DBG,'(A,2ES25.16)') "c,fc=", c,fc
!    WRITE(DBG,'(A,2ES25.16)') "xm,tol1=", xm,tol1

    IF (ABS(e) < tol1 .OR. ABS(fa) <= ABS(fb) ) THEN
      d=xm   ! bisection
      e=d
!      method=0
    ELSE
      IF (a==c) THEN
        s=fb/fa   ! linear interpolation
        p=TWO*xm*s
        q=ONE-s
!        method=1
      ELSE
        q=fa/fc   ! inverse quadratic interpolation
        r=fb/fc
        s=fb/fa
        p=s*(TWO*xm*q*(q-r)-(b-a)*(r-ONE))
        q=(q-ONE)*(r-ONE)*(s-ONE)
!        method=2
      END IF
      IF (p > ZERO) q=-q   ! adjust signs
      p=ABS(p)
      IF (p+p >= (THREE*xm*q-ABS(tol1*q)) .OR. p+p >= ABS(e*q) ) THEN
        d=xm   ! don't interpolate. Use bisection
        e=d
 !       method=-1
      ELSE
        e=d   ! OK, use interpolation
        d=p/q
      END IF
    END IF  
   
    a=b   ! complete step. a becomes the previous iteration
    fa=fb
   IF (ABS(d) > tol1) THEN
     b=b+d
   ELSE  
     b=b+SIGN(tol1,xm)
   END IF
     
   fb=F(b)   ! the newest and best value (we hope)
   neval=neval+1   ! keep count of the function evaluations
!!   IF ((fb*(fc/ABS(fc))) > ZERO) GO TO 20
!    WRITE(DBG,*) "BrentZero, end of kIter=",kIter, "   method=",method
!    WRITE(DBG,'(A,2ES25.16)') " new b,fb=", b,fb
!    WRITE(DBG,'(A,2ES25.16)') "d,e=", d,e
!    WRITE(DBG,'(A,2ES25.16)') "p,q=", p,q
  END DO
  
! The loop should never terminate. If it does, return the last iteration
!  and set errCode to 1
  xZero=b
  fZero=fb
  errCode=1
  RETURN  

END Subroutine BrentZeroDouble   ! ------------------------------------------------

!+
FUNCTION FminDouble (ax,bx,f,tol) RESULT (xopt)
! ---------------------------------------------------------------------------
! PURPOSE -

  INTERFACE
    FUNCTION F(x) RESULT(ff)  ! integrand function subprogram
    IMPLICIT NONE
      INTEGER,PARAMETER:: DP = SELECTED_REAL_KIND(10,50)
      REAL(DP),INTENT(IN)::x
      REAL(DP):: ff
    END FUNCTION F
  END INTERFACE

  REAL(DP),INTENT(IN):: ax   !  left endpoint of initial interval
  REAL(DP),INTENT(IN):: bx   !  right endpoint of initial interval
  REAL(DP),INTENT(IN):: tol  !  desired length of interval of uncertainity
  REAL(DP) :: xopt ! the result
  
  INTEGER:: errCode
  REAL(DP):: fzero
  INTEGER,PARAMETER:: MAX_ITER = 100
  INTEGER:: neval
!----------------------------------------------------------------------------
  CALL BrentMinDouble(ax,bx,F,tol,MAX_ITER,neval,errCode,xopt,fzero)
  
  RETURN
END Function FminDouble   ! -------------------------------------------------------

!+
SUBROUTINE BrentMinDouble(ax,bx,F,tol,maxIter,neval,errCode,xZero,fZero)
! ---------------------------------------------------------------------------
! PURPOSE -

  INTERFACE
    FUNCTION F(x) RESULT(ff)  ! integrand function subprogram
    IMPLICIT NONE
      INTEGER,PARAMETER:: DP = SELECTED_REAL_KIND(10,50)
      REAL(DP),INTENT(IN)::x
      REAL(DP):: ff
    END FUNCTION F
  END INTERFACE

  REAL(DP),INTENT(IN):: ax   !  left endpoint of initial interval
  REAL(DP),INTENT(IN):: bx   !  right endpoint of initial interval
  REAL(DP),INTENT(IN):: tol  !  desired length of interval of uncertainity
  INTEGER,INTENT(IN):: maxIter   ! max number of iterations allowed. 25 is good
  INTEGER,INTENT(OUT):: neval
  INTEGER,INTENT(OUT):: errCode   ! =0 is OK; =1 too many iterations
  REAL(DP),INTENT(OUT):: xZero,fZero ! the last and best value of the minimum                                  
      
  REAL(DP), PARAMETER :: C=0.381966011 ! (3-Sqrt(5))/2

  REAL(DP) :: a,b,d,e,eps,xm,p,q,r,tol1,tol2,u,v,w
  REAL(DP) :: fu,fv,fw,fx,x
  INTEGER :: iter
  REAL(DP),PARAMETER:: ZERO=0.0, ONE=1.0, TWO=2.0, THREE=3.0, HALF=0.5
!----------------------------------------------------------------------------
  errCode=0
  eps=SQRT(EPSILON(ONE))

  a=ax ! initialization
  b=bx
  v=a+c*(b-a)
  w=v
  x=v
  e=ZERO
  fx=F(x)
  fv=F(v)
  neval=2
  fw=fx

  DO iter=1,maxIter   ! iterate until minimum is found
    xm=HALF*(a+b)
    tol1=eps*ABS(x)+tol/THREE
    tol2=tol1+tol1
    IF (ABS(x-xm)<=(tol2-HALF*(b-a)) ) EXIT    ! check stopping criterion.

    IF (ABS(e) > tol1) THEN   ! is golden section necessary???
      r=(x-w)*(fx-fv)   ! trial parabolic fit
      q=(x-v)*(fx-fw)
      p=(x-v)*q-(x-w)*r
      q=TWO*(q-r)
      IF (q > ZERO) p=-p
      q=ABS(q)
      r=e
      e=d
      IF (ABS(p)>=ABS(HALF*q*r) .OR. p<=q*(a-x) .OR. p>=q*(b-x)) GO TO 40

      d=p/q   ! a parabolic interpolation step
      u=x+d
      IF ((u-a)<tol2 .OR. (b-u)<tol2) D=SIGN(tol1,xm-x)
      GO TO 50
    END IF

!40       CONTINUE
 40 IF (x >= xm) THEN   ! a golden section step
      e=a-x
    ELSE
      e=b-x
    END IF
    d=c*e

50  IF (ABS(d) >= tol1) THEN   ! f must not be evaluated too close to x
      u=x+d
    ELSE
      u=x+SIGN(tol1,d)
    END IF
    fu=F(u)
    neval=neval+1

    IF (fu <= fx) THEN   ! update a,b,v,w, and x
      IF (u >= x) THEN
        a=x
      ELSE
        b=x
      END IF
      v=w
      fv=fw
      w=x
      fw=fx
      x=u
      fx=fu
    ELSE
      IF (u<x) THEN
        a=u
      ELSE
        b=u
      END IF
      IF (fu<=fw .OR. w==x) THEN
        v=w
        fv=fw
        w=u
        fw=fu
      ELSE IF (fu<=fv .OR. v==x .OR. v==w) THEN
        v=u
        fv=fu
      END IF
    END IF
  END DO

   xZero=x   ! end of main loop
   fZero=fx
  RETURN
END Subroutine BrentMinDouble   ! -------------------------------------------------

!+
SUBROUTINE SvdDouble(a,w,matu,u,matv,v,ierr)
! ---------------------------------------------------------------------------
! PURPOSE - Determine the singular value decomposition of a real m by n 
!  rectangular matrix.  Householder bidiagonalization and a variant of the 
!  QR algorithm are used. This subroutine is a translation of the Algol 
!  procedure SVD, Num. Math. 14, 403-420(1970) by Golub and Reinsch.
!  Handbook for Auto. Comp., vol II-Linear algebra, 134-151(1971).
!  This is a modified version of a routine from the Eispack collection by 
!  the Nats project. Modified to eliminate machep.

  REAL(DP), INTENT(IN),DIMENSION(:,:) :: a  ! matrix to be decomposed. 
      ! on output, a is unaltered (unless overwritten by u or v).
  REAL(DP),INTENT(OUT),DIMENSION(:):: w ! w contains the n (non-negative) 
      ! singular values of a (the diagonal elements of s).  They are
      ! unordered. If an error exit is made, the singular values should be 
      ! correct for indices ierr+1,ierr+2,...,n.
  LOGICAL,INTENT(IN):: matu   ! matu should be set to .TRUE. if the u matrix
      ! in the decomposition is desired, and to .FALSE. otherwise.
  REAL(DP),INTENT(OUT),DIMENSION(:,:):: u   ! u contains the matrix u of 
      ! orthogonal column vectors of the decomposition if matu has been set
      ! to .TRUE. Otherwise, u is used as a temporary array.  
      ! u may coincide with a.  If an error exit is made, the columns of u 
      ! corresponding to indices of correct singular values should be correct.
  LOGICAL,INTENT(IN):: matv   ! matv should be set to .TRUE. if the v matrix
      ! in the decomposition is desired, and to .FALSE. otherwise.
  REAL(DP),INTENT(OUT),DIMENSION(:,:):: v   ! v contains the matrix v (orthogonal)
      ! of the decomposition if matv has been set to .TRUE.  Otherwise v is 
      ! not referenced. v may also coincide with a if u is not needed.  If an
      ! error exit is made, the columns of v corresponding to indices of 
      ! correct singular values should be correct.
  INTEGER,INTENT(OUT):: ierr   ! zero for normal return, k if the k-th 
      ! singular value has not been determined after 30 iterations.

  REAL(DP) :: c,f,g,h
  INTEGER :: i,j,k,l, i1,k1
  INTEGER:: its
  INTEGER:: m,n   ! dimensions of a
 
  REAL(DP),DIMENSION(SIZE(a,2)):: rv1    ! rv1 (n) is a temporary storage array.
  REAL(DP):: s,x,y,z,scale,anorm

  REAL(DP),PARAMETER:: ZERO = 0.0, ONE=1.0, TWO=2.0
!----------------------------------------------------------------------------
  m=SIZE(a,1)
  n=SIZE(a,2)
  IF (SIZE(w) < n) THEN
    ierr=-1
    RETURN
  END IF
  
  IF (matu .AND. (SIZE(u,1) < m .OR. SIZE(u,2) < n) ) THEN
    ierr=-2
    RETURN
  END IF
  
  IF (matv .AND. (SIZE(v,1) < n .OR. SIZE(v,2) < n) ) THEN
    ierr=-3
    RETURN
  END IF
    
  ierr = 0

!DO  i = 1, m
!  DO  j = 1, n
!    u(i,j) = a(i,j)
!  END DO
!END DO
  u(1:m,1:n)=a

  g = ZERO
  scale = ZERO
  anorm = ZERO

DO  i = 1, n   ! Householder reduction to bidiagonal form 
!!!  l = i + 1
  rv1(i) = scale * g
  g = ZERO
  s = ZERO
  scale = ZERO
  IF (i > m) GO TO 210
  
  DO  k = i, m
    scale = scale + ABS(u(k,i))
  END DO
  
  IF (scale == ZERO) GO TO 210
  
  DO  k = i, m
    u(k,i) = u(k,i) / scale
    s = s + u(k,i)**2
  END DO
  
  f = u(i,i)
  g = -SIGN(SQRT(s),f)
  h = f * g - s
  u(i,i) = f - g
!!!  IF (i == n) GO TO 190
  
  DO  j = i+1, n
    s = ZERO   
    DO  k = i, m
      s = s + u(k,i) * u(k,j)
    END DO
    
    f = s / h
    
    DO  k = i, m
      u(k,j) = u(k,j) + f * u(k,i)
    END DO
  END DO
  u(i:m,i)=scale*u(i:m,i)
!!!190 DO  k = i, m
!!!      u(k,i) = scale * u(k,i)
!!!    END DO
  
210 w(i) = scale * g
    g = ZERO
    s = ZERO
    scale = ZERO
    IF (i > m .OR. i == n) GO TO 290
  
    DO  k = i+1, n
      scale = scale + ABS(u(i,k))
    END DO
  
    IF (scale == ZERO) GO TO 290
  
    DO  k = i+1, n
      u(i,k) = u(i,k) / scale
      s = s + u(i,k)**2
    END DO
  
    f = u(i,i+1)
    g = -SIGN(SQRT(s),f)
    h = f * g - s
    u(i,i+1) = f - g
  
    DO  k = i+1, n
      rv1(k) = u(i,k) / h
    END DO
  
    IF (i == m) GO TO 270
  
    DO  j = i+1, m
      s = ZERO
      DO  k = i+1, n
        s = s + u(j,k) * u(i,k)
      END DO
      DO  k = i+1, n
        u(j,k) = u(j,k) + s * rv1(k)
      END DO
    END DO
  
270 DO  k = i+1, n
      u(i,k) = scale * u(i,k)
    END DO
  
290 anorm = MAX(anorm,ABS(w(i))+ABS(rv1(i)))
  END DO
!     .......... accumulation of right-hand transformations ..........
  IF (.NOT. matv) GO TO 410
!     .......... for i=n step -1 until 1 do -- ..........
  DO i=n,1,-1
    IF (i == n) GO TO 390
    IF (g == ZERO) GO TO 360
  
    DO  j = i+1, n
      v(j,i) = (u(i,j) / u(i,i+1)) / g   ! double division, YES!
    END DO
  
    DO  j = i+1, n
      s = ZERO
      DO  k = i+1, n
        s = s + u(i,k) * v(k,j)
      END DO
      DO  k = i+1, n
        v(k,j) = v(k,j) + s * v(k,i)
      END DO
    END DO
  
360 v(i,i+1:n)=ZERO
    v(i+1:n,i)=ZERO

!    DO  j = i+1, n
!      v(i,j) = ZERO
!      v(j,i) = ZERO
!    END DO
  
390 v(i,i) = ONE
    g = rv1(i)
!!!    l = i
  END DO
  
!     .......... accumulation of left-hand transformations ..........
410 IF (.NOT. matu) GO TO 510
!     ..........for i=min(m,n) step -1 until 1 do -- ..........
!    mn = n
!    IF (m < n) mn = m
    
    DO i=MIN(m,n),1,-1
!      i = mn + 1 - ii
!      l = i + 1
      g = w(i)
      IF (i == n) GO TO 430
  
      DO  j = i+1, n
        u(i,j) = ZERO
      END DO
  
430 IF (g == ZERO) THEN
      u(i:m,i)=ZERO
    ELSE  
      IF (i == MIN(m,n)) GO TO 460
  
      DO  j = i+1, n
        s = ZERO
        DO  k = i+1, m
          s = s + u(k,i) * u(k,j)
        END DO
        f = (s / u(i,i)) / g   ! double division, YES! Avoid underflow
    
        DO  k = i, m
          u(k,j) = u(k,j) + f * u(k,i)
        END DO
      END DO
  
460 DO  j = i, m
      u(j,i) = u(j,i) / g
    END DO
  
    END IF  

!475 DO  j = i, m
!      u(j,i) = ZERO
!    END DO
  
    u(i,i) = u(i,i) + ONE
  END DO
  
!     .......... diagonalization of the bidiagonal form ..........
!     .......... for k=n step -1 until 1 do -- ..........
510 DO k=n,1,-1
      k1=k-1
!      k1 = n - kk
!      k = k1 + 1
      its = 0
!     .......... test for splitting.
!                for l=k step -1 until 1 do -- ..........
520 DO L=k,1,-1
!      l1 = k - ll
!      l = l1 + 1
      IF (ABS(rv1(L)) + anorm <= anorm) GO TO 565
!     .......... rv1(1) is always zero, so there is no exit
!                through the bottom of the loop ..........
      IF (ABS(w(L-1)) + anorm <= anorm) EXIT
    END DO
!     .......... cancellation of rv1(l) if l greater than 1 ..........
    c = ZERO   ! was statement 540
    s = ONE
  
    DO  i = L, k
      f = s * rv1(i)
      rv1(i) = c * rv1(i)
      IF (ABS(f) + anorm <= anorm) EXIT
      g = w(i)
      h = SQRT(f*f+g*g)
      w(i) = h
      c = g / h
      s = -f / h
      IF (matu) THEN
        DO  j = 1, m
          y = u(j,L-1)
          z = u(j,i)
          u(j,L-1) = y * c + z * s
          u(j,i) = -y * s + z * c
        END DO
      END IF
    END DO
    
!     .......... test for convergence ..........
565 z = w(k)
    IF (L == k) GO TO 650
!     .......... shift from bottom 2 by 2 minor ..........
    IF (its == 30) GO TO 1000
    its = its + 1
    x = w(L)
    y = w(k1)
    g = rv1(k1)
    h = rv1(k)
    f = ((y - z) * (y + z) + (g - h) * (g + h)) / (TWO * h * y)
    g = SQRT(f*f+ONE)
    f = ((x - z) * (x + z) + h * (y / (f + SIGN(g,f)) - h)) / x
!     .......... next qr transformation ..........
    c = ONE
    s = ONE
  
    DO  i1 = L, k1
      i = i1 + 1
      g = rv1(i)
      y = w(i)
      h = s * g
      g = c * g
      z = SQRT(f*f+h*h)
      rv1(i1) = z
      c = f / z
      s = h / z
      f = x * c + g * s
      g = -x * s + g * c
      h = y * s
      y = y * c
      IF (matv) THEN
        DO  j = 1, n
          x = v(j,i1)
          z = v(j,i)
          v(j,i1) = x * c + z * s
          v(j,i) = -x * s + z * c
        END DO
      END IF    
    z = SQRT(f*f+h*h)   ! was statement #575
    w(i1) = z

    IF (z /= ZERO) THEN
      c = f / z   ! rotation can be arbitrary if z is zero
      s = h / z
    END IF  
    f = c * g + s * y
    x = -s * g + c * y
    IF (matu) THEN
      DO  j = 1, m
        y = u(j,i1)
        z = u(j,i)
        u(j,i1) = y * c + z * s
        u(j,i) = -y * s + z * c
      END DO
    END IF
  END DO
  
  rv1(L) = ZERO
  rv1(k) = f
  w(k) = x
  GO TO 520
!     .......... convergence ..........
650 IF (z >= ZERO) CYCLE

    w(k) = -z   !  w(k) is made non-negative 
    IF (matv) v(1:n,k)=-v(1:n,k)
  
!  DO  j = 1, n
!    v(j,k) = -v(j,k)
!  END DO
  
  END DO
  RETURN
  
!  set error -- no convergence to a singular value after 30 iterations ......
1000 ierr = k
  RETURN
END SUBROUTINE SvdDouble


END Module ForsytheMalcolmMoler   ! =========================================
