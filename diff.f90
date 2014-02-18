!DECK DIFF

! Code converted using TO_F90 by Alan Miller
! Date: 2014-02-17  Time: 20:46:07

SUBROUTINE diff(iord,x0,xmin,xmax,f,eps,acc,deriv,error,ifail)

!             NUMERICAL DIFFERENTIATION OF USER DEFINED FUNCTION

!                         DAVID KAHANER, NBS (GAITHERSBURG)

!  THE PROCEDURE DIFFERENTIATE CALCULATES THE FIRST, SECOND OR
!   THIRD ORDER DERIVATIVE OF A FUNCTION BY USING NEVILLE'S PROCESS TO
!   EXTRAPOLATE FROM A SEQUENCE OF SIMPLE POLYNOMIAL APPROXIMATIONS BASED ON
!   INTERPOLATING POINTS DISTRIBUTED SYMMETRICALLY ABOUT X0 (OR LYING ONLY ON
!   ONE SIDE OF X0 SHOULD THIS BE NECESSARY).  IF THE SPECIFIED TOLERANCE IS
!   NON-ZERO THEN THE PROCEDURE ATTEMPTS TO SATISFY THIS ABSOLUTE OR RELATIVE
!   ACCURACY REQUIREMENT, WHILE IF IT IS UNSUCCESSFUL OR IF THE TOLERANCE IS
!   SET TO ZERO THEN THE RESULT HAVING THE MINIMUM ACHIEVABLE ESTIMATED ERROR
!   IS RETURNED INSTEAD.

! INPUT PARAMETERS:
! IORD = 1, 2 OR 3 SPECIFIES THAT THE FIRST, SECOND OR THIRD ORDER
!   DERIVATIVE,RESPECTIVELY, IS REQUIRED.
! X0 IS THE POINT AT WHICH THE DERIVATIVE OF THE FUNCTION IS TO BE CALCULATED.
! XMIN, XMAX RESTRICT THE INTERPOLATING POINTS TO LIE IN [XMIN, XMAX], WHICH
!   SHOULD BE THE LARGEST INTERVAL INCLUDING X0 IN WHICH THE FUNCTION IS
!   CALCULABLE AND CONTINUOUS.
! F, A REAL PROCEDURE SUPPLIED BY THE USER, MUST YIELD THE VALUE OF THE
!   FUNCTION AT X FOR ANY X IN [XMIN, XMAX] WHEN CALLED BY F(X).
! EPS DENOTES THE TOLERANCE, EITHER ABSOLUTE OR RELATIVE.  EPS=0 SPECIFIES THAT
!   THE ERROR IS TO BE MINIMISED, WHILE EPS>0 OR EPS<0 SPECIFIES THAT THE
!   ABSOLUTE OR RELATIVE ERROR, RESPECTIVELY, MUST NOT EXCEED ABS(EPS) IF
!   POSSIBLE.  THE ACCURACY REQUIREMENT SHOULD NOT BE MADE STRICTER THAN
!   NECESSARY, SINCE THE AMOUNT OF COMPUTATION TENDS TO INCREASE AS
!   THE MAGNITUDE OF EPS DECREASES, AND IS PARTICULARLY HIGH WHEN EPS=0.
! ACC DENOTES THAT THE ABSOLUTE (ACC>0) OR RELATIVE (ACC<0) ERRORS IN THE
!   COMPUTED VALUES OF THE FUNCTION ARE MOST UNLIKELY TO EXCEED ABS(ACC), WHICH
!   SHOULD BE AS SMALL AS POSSIBLE.  IF THE USER CANNOT ESTIMATE ACC WITH
!   COMPLETE CONFIDENCE, THEN IT SHOULD BE SET TO ZERO.

! OUTPUT PARAMETERS:
! DERIV IS THE CALCULATED VALUE OF THE DERIVATIVE.
! ERROR IS AN ESTIMATED UPPER BOUND ON THE MAGNITUDE OF THE ABSOLUTE ERROR IN
!   THE CALCULATED RESULT.  IT SHOULD ALWAYS BE EXAMINED, SINCE IN EXTREME CASE
!   MAY INDICATE THAT THERE ARE NO CORRECT SIGNIFICANT DIGITS IN THE VALUE
!   RETURNED FOR DERIVATIVE.
! IFAIL WILL HAVE ONE OF THE FOLLOWING VALUES ON EXIT:
!   0   THE PROCEDURE WAS SUCCESSFUL.
!   1   THE ESTIMATED ERROR IN THE RESULT EXCEEDS THE (NON-ZERO) REQUESTED
!          ERROR, BUT THE MOST ACCURATE RESULT POSSIBLE HAS BEEN RETURNED.
!   2   INPUT DATA INCORRECT (DERIVATIVE AND ERROR WILL BE UNDEFINED).
!   3   THE INTERVAL [XMIN, XMAX] IS TOO SMALL (DERIVATIVE AND ERROR WILL BE
!          UNDEFINED);


INTEGER, INTENT(IN)                      :: iord
REAL, INTENT(IN)                         :: x0
REAL, INTENT(IN)                         :: xmin
REAL, INTENT(IN)                         :: xmax
REAL, INTENT(IN)                         :: f
REAL, INTENT(IN)                         :: eps
REAL, INTENT(IN OUT)                     :: acc
REAL, INTENT(OUT)                        :: deriv
REAL, INTENT(OUT)                        :: error
INTEGER, INTENT(OUT)                     :: ifail
EXTERNAL f
REAL :: beta,beta4,h,h0,h1,h2,  &
    newh1,newh2,heval,hprev,baseh,hacc1,hacc2,nhacc1,  &
    nhacc2,minh,maxh,maxh1,maxh2,tderiv,f0,twof0,f1,f2,f3,f4,fmax,  &
    maxfun,pmaxf,df1,deltaf,pdelta,z,zpower,c0f0,c1,c2,c3,dnew,dprev,  &
    re,te,newerr,temerr,newacc,pacc1,pacc2,facc1,facc2,acc0,  &
    acc1,acc2,relacc,twoinf,twosup,s,  &
    d(10),denom(10),e(10),minerr(10),maxf(0:10),SAVE(0:13), storef(-45:45),factor

INTEGER :: eta,inf,sup,i,j,k,n,nmax,method,signh,fcount, init
LOGICAL :: ignore(10),contin,saved


! ETA IS THE MINIMUM NUMBER OF SIGNIFICANT BINARY DIGITS (APART FROM THE
! SIGN DIGIT) USED TO REPRESENT THE MANTISSA OF REAL NUMBERS. IT SHOULD
! BE DEVREASED BY ONE IF THE COMPUTER TRUNCATES RATHER THAN ROUNDS.
! INF, SUP ARE THE LARGEST POSSIBLE POSITIVE INTEGERS SUBJECT TO
! 2**(-INF), -2**(-INF), 2**SUP, AND -2**SUP ALL BEING REPRESENTABLE REAL
! NUMBERS.
eta=i1mach(11) - 1
inf=-i1mach(12) - 2
sup=i1mach(13)-1
IF(iord < 1 .OR. iord > 3 .OR. xmax <= xmin .OR.  &
      x0 > xmax .OR. x0 < xmin) THEN
  ifail = 2
  RETURN
END IF

twoinf = 2.**(-inf)
twosup = 2.**sup
factor = 2**(FLOAT((inf+sup))/30.)
IF(factor < 256.)factor=256.
maxh1 = xmax - x0
signh = 1
IF(x0-xmin <= maxh1)THEN
  maxh2 = x0 - xmin
ELSE
  maxh2 = maxh1
  maxh1 = x0 - xmin
  signh = -1
END IF
relacc = 2.**(1-eta)
maxh1 = (1.-relacc)*maxh1
maxh2 = (1.-relacc)*maxh2
s=128.*twoinf
IF(ABS(x0) > 128.*twoinf*2.**eta) s = ABS(x0)*2.**(-eta)
IF(maxh1 < s)THEN
!         INTERVAL TOO SMALL
  ifail =3
  RETURN
END IF
IF(acc < 0.) THEN
  IF(-acc > relacc)relacc = -acc
  acc = 0.
END IF

!     DETERMINE THE SMALLEST SPACING AT WHICH THE CALCULATED
!     FUNCTION VALUES ARE UNEQUAL NEAR X0.

f0 = f(x0)
twof0 = f0 + f0
IF(ABS(x0) > twoinf*2.**eta) THEN
  h = ABS(x0)*2.**(-eta)
  z = 2.
ELSE
  h = twoinf
  z = 64.
END IF
df1 = f(x0+signh*h) - f0
80 IF(df1 /= 0. .OR. z*h > maxh1) GO TO 100
h = z*h
df1 = f(x0+signh*h) - f0
IF(z /= 2.) THEN
  IF(df1 /= 0.) THEN
    h = h/z
    z = 2.
    df1 = 0.
  ELSE
    IF(z*h > maxh1) z = 2.
  END IF
END IF
GO TO 80
100 CONTINUE

IF(df1 == 0.) THEN
!         CONSTANT FUNCTION
  deriv = 0.
  error = 0.
  ifail = 0
  RETURN
END IF
IF(h > maxh1/128.) THEN
!         MINIMUM H TOO LARGE
  ifail = 3
  RETURN
END IF

h = 8.*h
h1 = signh*h
h0 = h1
h2 = -h1
minh = 2.**(-MIN(inf,sup)/iord)
IF(minh < h) minh = h
IF(iord == 1) s = 8.
IF(iord == 2) s = 9.*SQRT(3.)
IF(iord == 3) s = 27.
IF(minh > maxh1/s) THEN
  ifail = 3
  RETURN
END IF
IF(minh > maxh2/s .OR. maxh2 < 128.*twoinf) THEN
  method = 1
ELSE
  method = 2
END IF

!     METHOD 1 USES 1-SIDED FORMULAE, AND METHOD 2 SYMMETRIC.
!         NOW ESTIMATE ACCURACY OF CALCULATED FUNCTION VALUES.

IF(method /= 2 .OR. iord == 2) THEN
  IF(x0 /= 0.) THEN
    CALL faccur(0.,-h1,acc0,x0,f,twoinf,f0,f1)
  ELSE
    acc0 = 0.
  END IF
END IF

IF(ABS(h1) > twosup/128.) THEN
  hacc1 = twosup
ELSE
  hacc1 = 128.*h1
END IF

IF(ABS(hacc1)/4. < minh) THEN
  hacc1 = 4.*signh*minh
ELSE IF(ABS(hacc1) > maxh1) THEN
  hacc1 = signh*maxh1
END IF
f1 = f(x0+hacc1)
CALL faccur(hacc1,h1,acc1,x0,f,twoinf,f0,f1)
IF(method == 2) THEN
  hacc2 = -hacc1
  IF(ABS(hacc2) > maxh2) hacc2 = -signh * maxh2
  f1 = f(x0 + hacc2)
  CALL faccur(hacc2,h2,acc2,x0,f,twoinf,f0,f1)
END IF
nmax = 8
IF(eta > 36) nmax = 10
n = -1
fcount = 0
deriv = 0.
error = twosup
init = 3
contin = .true.

130 n = n+1
IF(.NOT. contin) GO TO 800

IF(init == 3) THEN
!         CALCULATE COEFFICIENTS FOR DIFFERENTIATION FORMULAE
!             AND NEVILLE EXTRAPOLATION ALGORITHM
  IF(iord == 1) THEN
    beta=2.
  ELSE IF(method == 2)THEN
    beta = SQRT(2.)
  ELSE
    beta = SQRT(3.)
  END IF
  beta4 = beta**4.
  z = beta
  IF(method == 2) z = z**2
  zpower = 1.
  DO  k = 1,nmax
    zpower = z*zpower
    denom(k) = zpower-1
  END DO
  IF(method == 2 .AND. iord == 1) THEN
    e(1) = 5.
    e(2) = 6.3
    DO  i = 3,nmax
      e(i) = 6.81
    END DO
  ELSE IF((method /= 2.AND.iord == 1) .OR. (method == 2.AND.  &
        iord == 2)) THEN
    e(1) = 10.
    e(2) = 16.
    e(3) = 20.36
    e(4) = 23.
    e(5) = 24.46
    DO  i = 6,nmax
      e(i) = 26.
    END DO
    IF(method == 2.AND.iord == 2) THEN
      DO  i = 1,nmax
        e(i)=2*e(i)
      END DO
    END IF
  ELSE IF(method /= 2.AND.iord == 2) THEN
    e(1) = 17.78
    e(2) = 30.06
    e(3) = 39.66
    e(4) = 46.16
    e(5) = 50.26
    DO  i = 6,nmax
      e(i) = 55.
    END DO
  ELSE IF(method == 2.AND.iord == 3) THEN
    e(1) = 25.97
    e(2) = 41.22
    e(3) = 50.95
    e(4) = 56.4
    e(5) = 59.3
    DO  i = 6,nmax
      e(i) = 62.
    END DO
  ELSE
    e(1) = 24.5
    e(2) = 40.4
    e(3) = 52.78
    e(4) = 61.2
    e(5) = 66.55
    DO  i = 6,nmax
      e(i) = 73.
    END DO
    c0f0 = -twof0/(3.*beta)
    c1 = 3./(3.*beta-1.)
    c2 = -1./(3.*(beta-1.))
    c3 = 1./(3.*beta*(5.-2.*beta))
  END IF
END IF


IF(init >= 2) THEN
!         INITIALIZATION OF STEPLENGTHS, ACCURACY AND OTHER
!             PARAMETERS
  
  heval = signh*minh
  h = heval
  baseh = heval
  maxh = maxh2
  IF(method == 1)maxh = maxh1
  DO  k = 1,nmax
    minerr(k) = twosup
    ignore(k) = .false.
  END DO
  IF(method == 1) newacc = acc1
  IF(method == -1) newacc = acc2
  IF(method == 2) newacc = (acc1+acc2)/2.
  IF(newacc < acc) newacc = acc
  IF((method /= 2 .OR. iord == 2) .AND. newacc < acc0) newacc = acc0
  IF(method /= -1) THEN
    facc1 = acc1
    nhacc1 = hacc1
    newh1 = h1
  END IF
  IF(method /= 1) THEN
    facc2 = acc2
    nhacc2 = hacc2
    newh2 = h2
  ELSE
    facc2 = 0.
    nhacc2 = 0.
  END IF
  init = 1
  j = 0
  saved = .false.
END IF

!     CALCULATE NEW OR INITIAL FUNCTION VALUES

IF(init == 1 .AND. (n == 0 .OR. iord == 1) .AND.  &
      .NOT.(method == 2 .AND. fcount >= 45)) THEN
  IF(method == 2) THEN
    fcount = fcount + 1
    f1 = f(x0+heval)
    storef(fcount) = f1
    f2 = f(x0-heval)
    storef(-fcount) = f2
  ELSE
    j = j+1
    IF(j <= fcount) THEN
      f1 = storef(j*method)
    ELSE
      f1 = f(x0+heval)
    END IF
  END IF
ELSE
  f1 = f(x0+heval)
  IF(method == 2) f2 = f(x0-heval)
END IF
IF(n == 0) THEN
  IF(method == 2 .AND. iord == 3) THEN
    pdelta = f1-f2
    pmaxf = (ABS(f1)+ABS(f2))/2.
    heval = beta*heval
    f1 = f(x0+heval)
    f2 = f(x0-heval)
    deltaf = f1-f2
    maxfun = (ABS(f1)+ABS(f2))/2.
    heval = beta*heval
    f1 = f(x0+heval)
    f2 = f(x0-heval)
  ELSE IF(method /= 2 .AND. iord >= 2) THEN
    IF(iord == 2) THEN
      f3 = f1
    ELSE
      f4 = f1
      heval = beta*heval
      f3 = f(x0+heval)
    END IF
    heval = beta*heval
    f2 = f(x0+heval)
    heval = beta*heval
    f1 = f(x0+heval)
  END IF
END IF

!     EVALUATE A NEW APPROXIMATION DNEW TO THE DERIVATIVE

IF(n > nmax) THEN
  n = nmax
  DO  i = 1,n
    maxf(i-1) = maxf(i)
  END DO
END IF
IF(method == 2) THEN
  maxf(n) = (ABS(f1)+ABS(f2))/2.
  IF(iord == 1) THEN
    dnew = (f1-f2)/2.
  ELSE IF(iord == 2) THEN
    dnew = f1+f2-twof0
  ELSE
    dnew = -pdelta
    pdelta = deltaf
    deltaf = f1-f2
    dnew = dnew + .5*deltaf
    IF(maxf(n) < pmaxf) maxf(n) = pmaxf
    pmaxf = maxfun
    maxfun = (ABS(f1)+ABS(f2))/2.
  END IF
ELSE
  maxf(n) = ABS(f1)
  IF(iord == 1) THEN
    dnew = f1-f0
  ELSE IF(iord == 2) THEN
    dnew = (twof0-3*f3+f1)/3.
    IF(maxf(n) < ABS(f3)) maxf(n) = ABS(f3)
    f3 = f2
    f2 = f1
  ELSE
    dnew = c3*f1+c2*f2+c1*f4+c0f0
    IF(maxf(n) < ABS(f2)) maxf(n) = ABS(f2)
    IF(maxf(n) < ABS(f4)) maxf(n) = ABS(f4)
    f4 = f3
    f3 = f2
    f2 = f1
  END IF
END IF
IF(ABS(h) > 1) THEN
  dnew = dnew/h**iord
ELSE
  IF(128.*ABS(dnew) > twosup*ABS(h)**iord) THEN
    dnew = twosup/128.
  ELSE
    dnew = dnew/h**iord
  END IF
END IF

IF(init == 0) THEN
!         UPDATE ESTIMATED ACCURACY OF FUNCTION VALUES
  newacc = acc
  IF((method /= 2 .OR. iord == 2) .AND. newacc < acc0) newacc = acc0
  IF(method /= -1 .AND. ABS(nhacc1) <= 1.125*ABS(heval)/beta4) THEN
    nhacc1 = heval
    pacc1 = facc1
    CALL faccur(nhacc1,newh1,facc1,x0,f,twoinf,f0,f1)
    IF(facc1 < pacc1) facc1=(3*facc1+pacc1)/4.
  END IF
  IF(method /= 1 .AND. ABS(nhacc2) <= 1.125*ABS(heval)/beta4) THEN
    IF(method == 2) THEN
      f1 = f2
      nhacc2 = -heval
    ELSE
      nhacc2 = heval
    END IF
    pacc2 = facc2
    CALL faccur(nhacc2,newh2,facc2,x0,f,twoinf,f0,f1)
    IF(facc2 < pacc2) facc2 = (3*facc2+pacc2)/4.
  END IF
  IF(method == 1 .AND. newacc < facc1) newacc = facc1
  IF(method == -1 .AND. newacc < facc2) newacc = facc2
  IF(method == 2 .AND. newacc < (facc1+facc2)/2.) newacc = (facc1+facc2)/2.
END IF

!     EVALUATE SUCCESSIVE ELEMENTS OF THE CURRENT ROW IN THE NEVILLE
!     ARRAY, ESTIMATING AND EXAMINING THE TRUNCATION AND ROUNDING
!     ERRORS IN EACH

contin = n < nmax
hprev = ABS(h)
fmax = maxf(n)
IF((method /= 2 .OR. iord == 2) .AND. fmax < ABS(f0)) fmax = ABS(f0)

DO  k = 1,n
  dprev = d(k)
  d(k) = dnew
  dnew = dprev+(dprev-dnew)/denom(k)
  te = ABS(dnew-d(k))
  IF(fmax < maxf(n-k)) fmax = maxf(n-k)
  hprev = hprev/beta
  IF(newacc >= relacc*fmax) THEN
    re = newacc*e(k)
  ELSE
    re = relacc*fmax*e(k)
  END IF
  IF(re /= 0.) THEN
    IF(hprev > 1) THEN
      re = re/hprev**iord
    ELSE IF(2*re > twosup*hprev**iord) THEN
      re = twosup/2.
    ELSE
      re = re/hprev**iord
    END IF
  END IF
  newerr = te+re
  IF(te > re) newerr = 1.25*newerr
  IF(.NOT. ignore(k)) THEN
    IF((init == 0 .OR. (k == 2 .AND. .NOT.ignore(1)))  &
          .AND. newerr < error) THEN
      deriv = d(k)
      error = newerr
    END IF
    IF(init == 1 .AND. n == 1) THEN
      tderiv = d(1)
      temerr = newerr
    END IF
    IF(minerr(k) < twosup/4) THEN
      s = 4*minerr(k)
    ELSE
      s = twosup
    END IF
    IF(te > re .OR. newerr > s) THEN
      ignore(k) = .true.
    ELSE
      contin = .true.
    END IF
    IF(newerr < minerr(k)) minerr(k) = newerr
    IF(init == 1 .AND. n == 2 .AND. k == 1 .AND. .NOT.ignore(1)) THEN
      IF(newerr < temerr) THEN
        tderiv = d(1)
        temerr = newerr
      END IF
      IF(temerr < error) THEN
        deriv = tderiv
        error = temerr
      END IF
    END IF
  END IF
END DO

IF(n < nmax) d(n+1) = dnew
IF(eps < 0.) THEN
  s = ABS(eps*deriv)
ELSE
  s = eps
END IF
IF(error <= s) THEN
  contin = .false.
ELSE IF(init == 1 .AND. (n == 2 .OR. ignore(1))) THEN
  IF((ignore(1) .OR. ignore(2)) .AND. saved) THEN
    saved = .false.
    n = 2
    h = beta * SAVE(0)
    heval = beta*SAVE(1)
    maxf(0) = SAVE(2)
    maxf(1) = SAVE(3)
    maxf(2) = SAVE(4)
    d(1) = SAVE(5)
    d(2) = SAVE(6)
    d(3) = SAVE(7)
    minerr(1) = SAVE(8)
    minerr(2) = SAVE(9)
    IF(method == 2 .AND. iord == 3) THEN
      pdelta = SAVE(10)
      deltaf = SAVE(11)
      pmaxf = SAVE(12)
      maxfun = SAVE(13)
    ELSE IF(method /= 2 .AND. iord >= 2) THEN
      f2 = SAVE(10)
      f3 = SAVE(11)
      IF(iord == 3) f4 = SAVE(12)
    END IF
    init = 0
    ignore(1) = .false.
    ignore(2) = .false.
  ELSE IF(.NOT. (ignore(1) .OR. ignore(2)) .AND. n == 2  &
        .AND. beta4*factor*ABS(heval) <= maxh) THEN
!             SAVE ALL CURRENT VALUES IN CASE OF RETURN TO
!                 CURRENT POINT
    saved = .true.
    SAVE(0) = h
    SAVE(1) = heval
    SAVE(2) = maxf(0)
    SAVE(3) = maxf(1)
    SAVE(4) = maxf(2)
    SAVE(5) = d(1)
    SAVE(6) = d(2)
    SAVE(7) = d(3)
    SAVE(8) = minerr(1)
    SAVE(9) = minerr (2)
    IF(method == 2 .AND. iord == 3) THEN
      SAVE(10) = pdelta
      SAVE(11) = deltaf
      SAVE(12) = pmaxf
      SAVE(13) = maxfun
    ELSE IF(method /= 2 .AND. iord >= 2) THEN
      SAVE(10) = f2
      SAVE(11) = f3
      IF(iord == 3) SAVE(12) = f4
    END IF
    h = factor*baseh
    heval = h
    baseh = h
    n = -1
  ELSE
    init = 0
    h = beta*h
    heval = beta*heval
  END IF
ELSE IF(contin .AND. beta*ABS(heval) <= maxh) THEN
  h = beta*h
  heval = beta*heval
ELSE IF(method /= 1) THEN
  contin = .true.
  IF(method == 2) THEN
    init = 3
    method = -1
    IF(iord /= 2) THEN
      IF(x0 /= 0.) THEN
        CALL faccur(0.,-h0,acc0,x0,f,twoinf,f0,f1)
      ELSE
        acc0 = 0.
      END IF
    END IF
  ELSE
    init = 2
    method = 1
  END IF
  n = -1
  signh = -signh
ELSE
  contin = .false.
END IF
GO TO 130
800 IF(eps < 0.) THEN
  s = ABS(eps*deriv)
ELSE
  s = eps
END IF
ifail = 0
IF(eps /= 0. .AND. error > s) ifail = 1
RETURN
END SUBROUTINE diff
!DECK FACCUR

SUBROUTINE faccur(h0,h1,facc,x0,f,twoinf,f0,f1)

REAL, INTENT(IN OUT)                     :: h0
REAL, INTENT(OUT)                        :: h1
REAL, INTENT(OUT)                        :: facc
REAL, INTENT(IN OUT)                     :: x0
REAL, INTENT(IN)                         :: f
REAL, INTENT(IN)                         :: twoinf
REAL, INTENT(IN)                         :: f0
REAL, INTENT(IN)                         :: f1
REAL :: a0,a1,f00,f2,deltaf,t0,t1, df(5), f1
INTEGER :: j
EXTERNAL f

t0 = 0.
t1 = 0.
IF(h0 /= 0.) THEN
  IF(x0+h0 /= 0.) THEN
    f00 = f1
  ELSE
    h0 = 0.875*h0
    f00 = f(x0+h0)
  END IF
  IF(ABS(h1) >= 32.*twoinf) h1 = h1/8.
  IF(16.*ABS(h1) > ABS(h0)) h1 = SIGN(h1,1.)*ABS(h0)/16.
  IF(f(x0+h0-h1) == f00) THEN
    IF(256.*ABS(h1) <= ABS(h0)) THEN
      h1 = 2.*h1
      10             IF(f(x0+h0-h1) /= f00 .OR. 256.*ABS(h1) > ABS(h0))  &
          GO TO 20
      h1 = 2.*h1
      GO TO 10
      20             h1 = 8.*h1
      
    ELSE
      h1 = SIGN(h1,1.)*ABS(h0)/16.
    END IF
  ELSE
    IF(256.*twoinf <= ABS(h0)) THEN
      30             IF(f(x0+h0-h1/2.) == f00 .OR. ABS(h1) < 4.*twoinf)  &
          GO TO 40
      h1 = h1/2.
      GO TO 30
      40             CONTINUE
      h1 = 8.*h1
      IF(16.*ABS(h1) > ABS(h0)) h1 = SIGN(h1,1.) *ABS(h0)/16.
    ELSE
      h1 = SIGN(h1,1.)*ABS(h0)/16.
    END IF
  END IF
ELSE
  f00 = f0
END IF

DO  j = 1,5
  f2 = f(x0+h0-FLOAT(2*j-1)*h1)
  df(j) = f2 - f00
  t0 = t0+df(j)
  t1 = t1+FLOAT(2*j-1)*df(j)
END DO
a0 = (33.*t0-5.*t1)/73.
a1 = (-5.*t0+1.2*t1)/73.
facc = ABS(a0)
DO  j = 1,5
  deltaf = ABS(df(j)-(a0+FLOAT(2*j-1)*a1))
  IF(facc < deltaf) facc = deltaf
END DO
facc = 2.*facc
RETURN
END SUBROUTINE faccur
