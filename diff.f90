!deck diff

	module diff_module
	
	implicit none
	
	
	integer, parameter :: sp = selected_real_kind(6, 37)	!single
	integer, parameter :: dp = selected_real_kind(15, 307)	!double
	integer, parameter :: qp = selected_real_kind(33, 4931)	!quad
	
	!using double precision by default:
	integer,parameter,public :: wp = sp
	
	abstract interface
		function func(x) result(fx)
			import :: wp
			implicit none
			real(wp),intent(in) :: x
			real(wp) :: fx	
		end function func	
	end interface
	
	
	contains
	

! code converted using to_f90 by alan miller
! date: 2014-02-17  time: 20:23:56

subroutine diff(iord,x0,xmin,xmax,f,eps,acc,deriv,error,ifail)

!             numerical differentiation of user defined function
!
!                         david kahaner, nbs (gaithersburg)
!
!  the procedure differentiate calculates the first, second or
!   third order derivative of a function by using neville's process to
!   extrapolate from a sequence of simple polynomial approximations based on
!   interpolating points distributed symmetrically about x0 (or lying only on
!   one side of x0 should this be necessary).  if the specified tolerance is
!   non-zero then the procedure attempts to satisfy this absolute or relative
!   accuracy requirement, while if it is unsuccessful or if the tolerance is
!   set to zero then the result having the minimum achievable estimated error
!   is returned instead.
!
! input parameters:
! iord = 1, 2 or 3 specifies that the first, second or third order
!   derivative,respectively, is required.
! x0 is the point at which the derivative of the function is to be calculated.
! xmin, xmax restrict the interpolating points to lie in [xmin, xmax], which
!   should be the largest interval including x0 in which the function is
!   calculable and continuous.
! f, a real(wp) procedure supplied by the user, must yield the value of the
!   function at x for any x in [xmin, xmax] when called by f(x).
! eps denotes the tolerance, either absolute or relative.  eps=0 specifies that
!   the error is to be minimised, while eps>0 or eps<0 specifies that the
!   absolute or relative error, respectively, must not exceed abs(eps) if
!   possible.  the accuracy requirement should not be made stricter than
!   necessary, since the amount of computation tends to increase as
!   the magnitude of eps decreases, and is particularly high when eps=0.
! acc denotes that the absolute (acc>0) or relative (acc<0) errors in the
!   computed values of the function are most unlikely to exceed abs(acc), which
!   should be as small as possible.  if the user cannot estimate acc with
!   complete confidence, then it should be set to zero.
!
! output parameters:
! deriv is the calculated value of the derivative.
! error is an estimated upper bound on the magnitude of the absolute error in
!   the calculated result.  it should always be examined, since in extreme case
!   may indicate that there are no correct significant digits in the value
!   returned for derivative.
! ifail will have one of the following values on exit:
!   0   the procedure was successful.
!   1   the estimated error in the result exceeds the (non-zero) requested
!          error, but the most accurate result possible has been returned.
!   2   input data incorrect (derivative and error will be undefined).
!   3   the interval [xmin, xmax] is too small (derivative and error will be
!          undefined);

integer, intent(in)          :: iord
real(wp), intent(in)         :: x0
real(wp), intent(in)         :: xmin
real(wp), intent(in)         :: xmax
real(wp), intent(in)         :: eps
real(wp), intent(in out)     :: acc
real(wp), intent(out)        :: deriv
real(wp), intent(out)        :: error
integer, intent(out)         :: ifail
procedure(func) :: f

real(wp) :: beta,beta4,h,h0,h1,h2,  &
    newh1,newh2,heval,hprev,baseh,hacc1,hacc2,nhacc1,  &
    nhacc2,minh,maxh,maxh1,maxh2,tderiv,f0,twof0,f1,f2,f3,f4,fmax,  &
    maxfun,pmaxf,df1,deltaf,pdelta,z,zpower,c0f0,c1,c2,c3,dnew,dprev,  &
    re,te,newerr,temerr,newacc,pacc1,pacc2,facc1,facc2,acc0,  &
    acc1,acc2,relacc,twoinf,twosup,s,  &
    d(10),denom(10),e(10),minerr(10),maxf(0:10),save(0:13), storef(-45:45),factor
integer :: eta,inf,sup,i,j,k,n,nmax,method,signh,fcount, init
logical :: ignore(10),contin,saved
real(wp) :: dummy1,dummy2

! eta is the minimum number of significant binary digits (apart from the
! sign digit) used to represent the mantissa of real(wp) numbers. it should
! be devreased by one if the computer truncates rather than rounds.
! inf, sup are the largest possible positive integers subject to
! 2**(-inf), -2**(-inf), 2**sup, and -2**sup all being representable real(wp)
! numbers.
eta=i1mach(11) - 1
inf=-i1mach(12) - 2
sup=i1mach(13) - 1
if(iord < 1 .or. iord > 3 .or. xmax <= xmin .or.  &
      x0 > xmax .or. x0 < xmin) then
  ifail = 2
  return
end if

twoinf = 2.**(-inf)
twosup = 2.**sup
factor = 2**(float((inf+sup))/30.)
if(factor < 256.)factor=256.
maxh1 = xmax - x0
signh = 1
if(x0-xmin <= maxh1)then
  maxh2 = x0 - xmin
else
  maxh2 = maxh1
  maxh1 = x0 - xmin
  signh = -1
end if
relacc = 2.**(1-eta)
maxh1 = (1.-relacc)*maxh1
maxh2 = (1.-relacc)*maxh2
s=128.*twoinf
if(abs(x0) > 128.*twoinf*2.**eta) s = abs(x0)*2.**(-eta)
if(maxh1 < s)then
!         interval too small
  ifail =3
  return
end if
if(acc < 0.) then
  if(-acc > relacc)relacc = -acc
  acc = 0.
end if

!     determine the smallest spacing at which the calculated
!     function values are unequal near x0.

f0 = f(x0)
twof0 = f0 + f0
if(abs(x0) > twoinf*2.**eta) then
  h = abs(x0)*2.**(-eta)
  z = 2.
else
  h = twoinf
  z = 64.
end if
df1 = f(x0+signh*h) - f0
80 if(df1 /= 0. .or. z*h > maxh1) go to 100
h = z*h
df1 = f(x0+signh*h) - f0
if(z /= 2.) then
  if(df1 /= 0.) then
    h = h/z
    z = 2.
    df1 = 0.
  else
    if(z*h > maxh1) z = 2.
  end if
end if
go to 80
100 continue

if(df1 == 0.) then
!         constant function
  deriv = 0.
  error = 0.
  ifail = 0
  return
end if
if(h > maxh1/128.) then
!         minimum h too large
  ifail = 3
  return
end if

h = 8.*h
h1 = signh*h
h0 = h1
h2 = -h1
minh = 2.**(-min(inf,sup)/iord)
if(minh < h) minh = h
if(iord == 1) s = 8.
if(iord == 2) s = 9.*sqrt(3.)
if(iord == 3) s = 27.
if(minh > maxh1/s) then
  ifail = 3
  return
end if
if(minh > maxh2/s .or. maxh2 < 128.*twoinf) then
  method = 1
else
  method = 2
end if

!     method 1 uses 1-sided formulae, and method 2 symmetric.
!         now estimate accuracy of calculated function values.

if(method /= 2 .or. iord == 2) then
  if(x0 /= 0.) then
    dummy1 = 0.
    dummy2 = -h1
    call faccur(dummy1,dummy2,acc0,x0,f,twoinf,f0,f1)
  else
    acc0 = 0.
  end if
end if

if(abs(h1) > twosup/128.) then
  hacc1 = twosup
else
  hacc1 = 128.*h1
end if

if(abs(hacc1)/4. < minh) then
  hacc1 = 4.*signh*minh
else if(abs(hacc1) > maxh1) then
  hacc1 = signh*maxh1
end if
f1 = f(x0+hacc1)
call faccur(hacc1,h1,acc1,x0,f,twoinf,f0,f1)
if(method == 2) then
  hacc2 = -hacc1
  if(abs(hacc2) > maxh2) hacc2 = -signh * maxh2
  f1 = f(x0 + hacc2)
  call faccur(hacc2,h2,acc2,x0,f,twoinf,f0,f1)
end if
nmax = 8
if(eta > 36) nmax = 10
n = -1
fcount = 0
deriv = 0.
error = twosup
init = 3
contin = .true.

130 n = n+1
if(.not. contin) go to 800

if(init == 3) then
!         calculate coefficients for differentiation formulae
!             and neville extrapolation algorithm
  if(iord == 1) then
    beta=2.
  else if(method == 2)then
    beta = sqrt(2.)
  else
    beta = sqrt(3.)
  end if
  beta4 = beta**4.
  z = beta
  if(method == 2) z = z**2
  zpower = 1.
  do  k = 1,nmax
    zpower = z*zpower
    denom(k) = zpower-1
  end do
  if(method == 2 .and. iord == 1) then
    e(1) = 5.
    e(2) = 6.3
    do  i = 3,nmax
      e(i) = 6.81
    end do
  else if((method /= 2.and.iord == 1) .or. (method == 2.and.  &
        iord == 2)) then
    e(1) = 10.
    e(2) = 16.
    e(3) = 20.36
    e(4) = 23.
    e(5) = 24.46
    do  i = 6,nmax
      e(i) = 26.
    end do
    if(method == 2.and.iord == 2) then
      do  i = 1,nmax
        e(i)=2*e(i)
      end do
    end if
  else if(method /= 2.and.iord == 2) then
    e(1) = 17.78
    e(2) = 30.06
    e(3) = 39.66
    e(4) = 46.16
    e(5) = 50.26
    do  i = 6,nmax
      e(i) = 55.
    end do
  else if(method == 2.and.iord == 3) then
    e(1) = 25.97
    e(2) = 41.22
    e(3) = 50.95
    e(4) = 56.4
    e(5) = 59.3
    do  i = 6,nmax
      e(i) = 62.
    end do
  else
    e(1) = 24.5
    e(2) = 40.4
    e(3) = 52.78
    e(4) = 61.2
    e(5) = 66.55
    do  i = 6,nmax
      e(i) = 73.
    end do
    c0f0 = -twof0/(3.*beta)
    c1 = 3./(3.*beta-1.)
    c2 = -1./(3.*(beta-1.))
    c3 = 1./(3.*beta*(5.-2.*beta))
  end if
end if


if(init >= 2) then
!         initialization of steplengths, accuracy and other
!             parameters
  
  heval = signh*minh
  h = heval
  baseh = heval
  maxh = maxh2
  if(method == 1)maxh = maxh1
  do  k = 1,nmax
    minerr(k) = twosup
    ignore(k) = .false.
  end do
  if(method == 1) newacc = acc1
  if(method == -1) newacc = acc2
  if(method == 2) newacc = (acc1+acc2)/2.
  if(newacc < acc) newacc = acc
  if((method /= 2 .or. iord == 2) .and. newacc < acc0) newacc = acc0
  if(method /= -1) then
    facc1 = acc1
    nhacc1 = hacc1
    newh1 = h1
  end if
  if(method /= 1) then
    facc2 = acc2
    nhacc2 = hacc2
    newh2 = h2
  else
    facc2 = 0.
    nhacc2 = 0.
  end if
  init = 1
  j = 0
  saved = .false.
end if

!     calculate new or initial function values

if(init == 1 .and. (n == 0 .or. iord == 1) .and.  &
      .not.(method == 2 .and. fcount >= 45)) then
  if(method == 2) then
    fcount = fcount + 1
    f1 = f(x0+heval)
    storef(fcount) = f1
    f2 = f(x0-heval)
    storef(-fcount) = f2
  else
    j = j+1
    if(j <= fcount) then
      f1 = storef(j*method)
    else
      f1 = f(x0+heval)
    end if
  end if
else
  f1 = f(x0+heval)
  if(method == 2) f2 = f(x0-heval)
end if
if(n == 0) then
  if(method == 2 .and. iord == 3) then
    pdelta = f1-f2
    pmaxf = (abs(f1)+abs(f2))/2.
    heval = beta*heval
    f1 = f(x0+heval)
    f2 = f(x0-heval)
    deltaf = f1-f2
    maxfun = (abs(f1)+abs(f2))/2.
    heval = beta*heval
    f1 = f(x0+heval)
    f2 = f(x0-heval)
  else if(method /= 2 .and. iord >= 2) then
    if(iord == 2) then
      f3 = f1
    else
      f4 = f1
      heval = beta*heval
      f3 = f(x0+heval)
    end if
    heval = beta*heval
    f2 = f(x0+heval)
    heval = beta*heval
    f1 = f(x0+heval)
  end if
end if

!     evaluate a new approximation dnew to the derivative

if(n > nmax) then
  n = nmax
  do  i = 1,n
    maxf(i-1) = maxf(i)
  end do
end if
if(method == 2) then
  maxf(n) = (abs(f1)+abs(f2))/2.
  if(iord == 1) then
    dnew = (f1-f2)/2.
  else if(iord == 2) then
    dnew = f1+f2-twof0
  else
    dnew = -pdelta
    pdelta = deltaf
    deltaf = f1-f2
    dnew = dnew + .5*deltaf
    if(maxf(n) < pmaxf) maxf(n) = pmaxf
    pmaxf = maxfun
    maxfun = (abs(f1)+abs(f2))/2.
  end if
else
  maxf(n) = abs(f1)
  if(iord == 1) then
    dnew = f1-f0
  else if(iord == 2) then
    dnew = (twof0-3*f3+f1)/3.
    if(maxf(n) < abs(f3)) maxf(n) = abs(f3)
    f3 = f2
    f2 = f1
  else
    dnew = c3*f1+c2*f2+c1*f4+c0f0
    if(maxf(n) < abs(f2)) maxf(n) = abs(f2)
    if(maxf(n) < abs(f4)) maxf(n) = abs(f4)
    f4 = f3
    f3 = f2
    f2 = f1
  end if
end if
if(abs(h) > 1) then
  dnew = dnew/h**iord
else
  if(128.*abs(dnew) > twosup*abs(h)**iord) then
    dnew = twosup/128.
  else
    dnew = dnew/h**iord
  end if
end if

if(init == 0) then
!         update estimated accuracy of function values
  newacc = acc
  if((method /= 2 .or. iord == 2) .and. newacc < acc0) newacc = acc0
  if(method /= -1 .and. abs(nhacc1) <= 1.125*abs(heval)/beta4) then
    nhacc1 = heval
    pacc1 = facc1
    call faccur(nhacc1,newh1,facc1,x0,f,twoinf,f0,f1)
    if(facc1 < pacc1) facc1=(3*facc1+pacc1)/4.
  end if
  if(method /= 1 .and. abs(nhacc2) <= 1.125*abs(heval)/beta4) then
    if(method == 2) then
      f1 = f2
      nhacc2 = -heval
    else
      nhacc2 = heval
    end if
    pacc2 = facc2
    call faccur(nhacc2,newh2,facc2,x0,f,twoinf,f0,f1)
    if(facc2 < pacc2) facc2 = (3*facc2+pacc2)/4.
  end if
  if(method == 1 .and. newacc < facc1) newacc = facc1
  if(method == -1 .and. newacc < facc2) newacc = facc2
  if(method == 2 .and. newacc < (facc1+facc2)/2.) newacc = (facc1+facc2)/2.
end if

!     evaluate successive elements of the current row in the neville
!     array, estimating and examining the truncation and rounding
!     errors in each

contin = n < nmax
hprev = abs(h)
fmax = maxf(n)
if((method /= 2 .or. iord == 2) .and. fmax < abs(f0)) fmax = abs(f0)

do  k = 1,n
  dprev = d(k)
  d(k) = dnew
  dnew = dprev+(dprev-dnew)/denom(k)
  te = abs(dnew-d(k))
  if(fmax < maxf(n-k)) fmax = maxf(n-k)
  hprev = hprev/beta
  if(newacc >= relacc*fmax) then
    re = newacc*e(k)
  else
    re = relacc*fmax*e(k)
  end if
  if(re /= 0.) then
    if(hprev > 1) then
      re = re/hprev**iord
    else if(2*re > twosup*hprev**iord) then
      re = twosup/2.
    else
      re = re/hprev**iord
    end if
  end if
  newerr = te+re
  if(te > re) newerr = 1.25*newerr
  if(.not. ignore(k)) then
    if((init == 0 .or. (k == 2 .and. .not.ignore(1)))  &
          .and. newerr < error) then
      deriv = d(k)
      error = newerr
    end if
    if(init == 1 .and. n == 1) then
      tderiv = d(1)
      temerr = newerr
    end if
    if(minerr(k) < twosup/4) then
      s = 4*minerr(k)
    else
      s = twosup
    end if
    if(te > re .or. newerr > s) then
      ignore(k) = .true.
    else
      contin = .true.
    end if
    if(newerr < minerr(k)) minerr(k) = newerr
    if(init == 1 .and. n == 2 .and. k == 1 .and. .not.ignore(1)) then
      if(newerr < temerr) then
        tderiv = d(1)
        temerr = newerr
      end if
      if(temerr < error) then
        deriv = tderiv
        error = temerr
      end if
    end if
  end if
end do

if(n < nmax) d(n+1) = dnew
if(eps < 0.) then
  s = abs(eps*deriv)
else
  s = eps
end if
if(error <= s) then
  contin = .false.
else if(init == 1 .and. (n == 2 .or. ignore(1))) then
  if((ignore(1) .or. ignore(2)) .and. saved) then
    saved = .false.
    n = 2
    h = beta * save(0)
    heval = beta*save(1)
    maxf(0) = save(2)
    maxf(1) = save(3)
    maxf(2) = save(4)
    d(1) = save(5)
    d(2) = save(6)
    d(3) = save(7)
    minerr(1) = save(8)
    minerr(2) = save(9)
    if(method == 2 .and. iord == 3) then
      pdelta = save(10)
      deltaf = save(11)
      pmaxf = save(12)
      maxfun = save(13)
    else if(method /= 2 .and. iord >= 2) then
      f2 = save(10)
      f3 = save(11)
      if(iord == 3) f4 = save(12)
    end if
    init = 0
    ignore(1) = .false.
    ignore(2) = .false.
  else if(.not. (ignore(1) .or. ignore(2)) .and. n == 2  &
        .and. beta4*factor*abs(heval) <= maxh) then
!             save all current values in case of return to
!                 current point
    saved = .true.
    save(0) = h
    save(1) = heval
    save(2) = maxf(0)
    save(3) = maxf(1)
    save(4) = maxf(2)
    save(5) = d(1)
    save(6) = d(2)
    save(7) = d(3)
    save(8) = minerr(1)
    save(9) = minerr (2)
    if(method == 2 .and. iord == 3) then
      save(10) = pdelta
      save(11) = deltaf
      save(12) = pmaxf
      save(13) = maxfun
    else if(method /= 2 .and. iord >= 2) then
      save(10) = f2
      save(11) = f3
      if(iord == 3) save(12) = f4
    end if
    h = factor*baseh
    heval = h
    baseh = h
    n = -1
  else
    init = 0
    h = beta*h
    heval = beta*heval
  end if
else if(contin .and. beta*abs(heval) <= maxh) then
  h = beta*h
  heval = beta*heval
else if(method /= 1) then
  contin = .true.
  if(method == 2) then
    init = 3
    method = -1
    if(iord /= 2) then
      if(x0 /= 0.) then
        dummy1 = 0.
        dummy2 = -h0       
        call faccur(dummy1,dummy2,acc0,x0,f,twoinf,f0,f1)
      else
        acc0 = 0.
      end if
    end if
  else
    init = 2
    method = 1
  end if
  n = -1
  signh = -signh
else
  contin = .false.
end if
go to 130
800 if(eps < 0.) then
  s = abs(eps*deriv)
else
  s = eps
end if
ifail = 0
if(eps /= 0. .and. error > s) ifail = 1
return
end subroutine diff
!deck faccur

subroutine faccur(h0,h1,facc,x0,f,twoinf,f0,f1)

real(wp), intent(in out)                     :: h0
real(wp), intent(out)                        :: h1
real(wp), intent(out)                        :: facc
real(wp), intent(in)                         :: x0
real(wp), intent(in)                         :: twoinf
real(wp), intent(in)                         :: f0
real(wp), intent(in)                         :: f1
real(wp) :: a0,a1,f00,f2,deltaf,t0,t1, df(5)
integer :: j
procedure(func) :: f

t0 = 0.
t1 = 0.
if(h0 /= 0.) then
  if(x0+h0 /= 0.) then
    f00 = f1
  else
    h0 = 0.875*h0
    f00 = f(x0+h0)
  end if
  if(abs(h1) >= 32.*twoinf) h1 = h1/8.
  if(16.*abs(h1) > abs(h0)) h1 = sign(h1,1.)*abs(h0)/16.
  if(f(x0+h0-h1) == f00) then
    if(256.*abs(h1) <= abs(h0)) then
      h1 = 2.*h1
      10             if(f(x0+h0-h1) /= f00 .or. 256.*abs(h1) > abs(h0))  &
          go to 20
      h1 = 2.*h1
      go to 10
      20             h1 = 8.*h1
      
    else
      h1 = sign(h1,1.)*abs(h0)/16.
    end if
  else
    if(256.*twoinf <= abs(h0)) then
      30             if(f(x0+h0-h1/2.) == f00 .or. abs(h1) < 4.*twoinf)  &
          go to 40
      h1 = h1/2.
      go to 30
      40             continue
      h1 = 8.*h1
      if(16.*abs(h1) > abs(h0)) h1 = sign(h1,1.) *abs(h0)/16.
    else
      h1 = sign(h1,1.)*abs(h0)/16.
    end if
  end if
else
  f00 = f0
end if

do  j = 1,5
  f2 = f(x0+h0-float(2*j-1)*h1)
  df(j) = f2 - f00
  t0 = t0+df(j)
  t1 = t1+float(2*j-1)*df(j)
end do
a0 = (33.*t0-5.*t1)/73.
a1 = (-5.*t0+1.2*t1)/73.
facc = abs(a0)
do  j = 1,5
  deltaf = abs(df(j)-(a0+float(2*j-1)*a1))
  if(facc < deltaf) facc = deltaf
end do
facc = 2.*facc
return
end subroutine faccur

	integer function i1mach (i)
	implicit none
	integer :: i
	real(wp),parameter :: x = 1.0_wp

	select case(i)
	case (11)
		i1mach = digits(x)		!t, the number of base-b digits.
	case (12)
		i1mach = minexponent(x)	!emin, the smallest exponent e.
	case (13)
		i1mach = maxexponent(x)	!emax, the largest exponent e.
	end select
	
	end function i1mach


	subroutine test_case()
	implicit none
	
	integer,parameter  :: iord = 1
	real(wp),parameter :: x0    = 0.12345_wp
	real(wp),parameter :: xmin  = 0.0_wp
	real(wp),parameter :: xmax  = 1.0_wp
	real(wp),parameter :: eps   = 1.0e-5_wp

	real(wp) :: acc
	real(wp) :: deriv
	real(wp) :: error
	integer  :: ifail

	acc = 0.0_wp
	
	call diff(iord,x0,xmin,xmax,test_function,eps,acc,deriv,error,ifail)
	
	write(*,*) ''
	write(*,*) 'solution          :',deriv
	write(*,*) 'estimated error   :',error
	write(*,*) 'ifail             :',ifail
	write(*,*) ''
	write(*,*) 'actual derivative :', cos(x0)
	write(*,*) 'actual error      :', cos(x0) - deriv
	write(*,*) ''
	
	contains
	
		function test_function(x) result(fx)
		implicit none
		real(wp),intent(in) :: x
		real(wp) :: fx
		
		fx = sin(x)
		
		end function test_function
	
	end subroutine test_case

	end module diff_module

	
	program test
	
	use diff_module
	
	call test_case()
	
	end program test
