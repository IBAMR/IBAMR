module quadglobal
use qdmodule
implicit none
integer ndebug, ndigits, nerror, nquadl
end module

! program tquadgsq
subroutine f_main

!   David H. Bailey      2004-12-16
!   This is Quad-Double Fortran-90 version.

!   This work was supported by the Director, Office of Science, Division
!   of Mathematical, Information, and Computational Sciences of the
!   U.S. Department of Energy under contract number DE-AC03-76SF00098.

!   This program demonstrates the quadrature routine 'quadgsq', which employs
!   Gaussian quadrature.  The function quadgsq is suitable to integrate
!   a function that is continuous, infinitely differentiable and integrable on a
!   finite open interval.  It can also be used for certain integrals on infinite
!   intervals, by making a suitable change of variable -- see below.  While
!   this routine is usually more efficient than quaderf or quadts for functions
!   that are regular on a closed interval, it is not very effective for 
!   functions with a singularity at one or both of the endpoints.

!   The function(s) to be integrated is(are) defined in external function
!   subprogram(s) -- see the sample function subprograms below.  The name(s) of
!   the function subprogram(s) must be included in appropriate type and external
!   statements in the main program.

!   Note that an integral of a function on an infinite interval can be
!   converted to an integral on a finite interval by means of a suitable
!   change of variable.  Example (here the notation "inf" means infinity):

!   Int_0^inf f(t) dt  =  Int_0^1 f(t) dt + Int_1^inf f(t) dt
!                      =  Int_0^1 f(t) dt + Int_0^1 f(1/t)/t^2 dt

!   See examples below.

!   Inputs set in parameter statement below:
!   kdebug Debug level setting.  Default = 2.
!   ndp    Digits of precision.  May not exceed mpipl in file mpmod90.f.
!            In some cases, ndp must be significantly greater than the desired
!            tolerance in the result-- see the examples below.
!   neps   Log10 of the desired tolerance in the result (negative integer).
!   nq1    Max number of phases in quadrature routine; adding 1 increases
!            (possibly doubles) the number of accurate digits in the result,
!            but also roughly doubles the run time and memory.  nq1 > 2.
!   nq2    Space parameter for wk and xk arrays in the calling program.  By
!            default it is set to 8 * 2^nq1.  Increase nq2 if directed by a 
!            message produced in initqerf.  Note that the dimension of the
!            wk and xk arrays starts with -1, so the length of these arrays is
!            (nq2+2) * 4 eight-byte words.

use qdmodule
use quadglobal
implicit none
integer i, kdebug, ndp, neps, nq1, nq2, n1
parameter (kdebug = 2, ndp = 64, neps = -64, nq1 = 8, nq2 = 8 * 2 ** nq1 +100)
double precision dplog10q, d1, d2, second, tm0, tm1
type (qd_real) err, quadgsq, fun01, fun02, fun03, fun04, fun05, fun06, fun07, &
  fun08, fun09, fun10, fun11, fun12, fun13, fun14, fun15a, fun15b, &
  t1, t2, t3, t4, wk(-1:nq2), xk(-1:nq2), x1, x2
external quadgsq, fun01, fun02, fun03, fun04, fun05, fun06, fun07, fun08, &
  fun09, fun10, fun11, fun12, fun13, fun14, fun15a, fun15b, second
integer*4 old_cw

call f_fpu_fix_start (old_cw)

ndebug = kdebug
ndigits = ndp
nerror = 0
nquadl = nq1
write (6, 1) ndigits, neps, nquadl
1 format ('Quadgsq test'/'Digits =',i6,'  Epsilon =',i6,'   Quadlevel =',i6)

!   Initialize quadrature tables wk and xk (weights and abscissas).

tm0 = second ()
call initqgsq (nq1, nq2, wk, xk)
tm1 = second ()
if (nerror > 0) stop
write (6, 2) tm1 - tm0
2 format ('Quadrature initialization completed: cpu time =',f12.6)

!   Begin quadrature tests.

write (6, 11)
11 format (/'Continuous functions on finite itervals:'//&
  'Problem 1: Int_0^1 t*log(1+t) dt = 1/4')
x1 = 0.d0
x2 = 1.d0
tm0 = second ()
t1 = quadgsq (fun01, x1, x2, nq1, nq2, wk, xk)
tm1 = second ()
write (6, 3) tm1 - tm0
3 format ('Quadrature completed: CPU time =',f12.6/'Result =')
call qdwrite (6, t1)
t2 = 0.25d0
call decmdq (t2 - t1, d1, n1)
write (6, 4) d1, n1
4 format ('Actual error =',f10.6,'x10^',i5)

write (6, 12)
12 format (/'Problem 2: Int_0^1 t^2*arctan(t) dt = (pi - 2 + 2*log(2))/12')
x1 = 0.d0
x2 = 1.d0
tm0 = second ()
t1 = quadgsq (fun02, x1, x2, nq1, nq2, wk, xk)
tm1 = second ()
write (6, 3) tm1 - tm0
call qdwrite (6, t1)
t2 = (qdpi() - 2.d0 + 2.d0 * log (qdreal (2.d0))) / 12.d0
call decmdq (t2 - t1, d1, n1)
write (6, 4) d1, n1

write (6, 13)
13 format (/'Problem 3: Int_0^(pi/2) e^t*cos(t) dt = 1/2*(e^(pi/2) - 1)')
x1 = 0.d0
x2 = 0.5d0 * qdpi()
tm0 = second ()
t1 = quadgsq (fun03, x1, x2, nq1, nq2, wk, xk)
tm1 = second ()
write (6, 3) tm1 - tm0
call qdwrite (6, t1)
t2 = 0.5d0 * (exp (0.5d0 * qdpi()) - 1.d0)
call decmdq (t2 - t1, d1, n1)
write (6, 4) d1, n1

write (6, 14)
14 format (/ &
  'Problem 4: Int_0^1 arctan(sqrt(2+t^2))/((1+t^2)sqrt(2+t^2)) dt = 5*Pi^2/96')
x1 = 0.d0
x2 = 1.d0
tm0 = second ()
t1 = quadgsq (fun04, x1, x2, nq1, nq2, wk, xk)
tm1 = second ()
write (6, 3) tm1 - tm0
call qdwrite (6, t1)
t2 = 5.d0 * qdpi()**2 / 96.d0
call decmdq (t2 - t1, d1, n1)
write (6, 4) d1, n1

write (6, 15)
15 format (/&
  'Continuous functions on finite itervals, but non-diff at an endpoint'// &
  'Problem 5: Int_0^1 sqrt(t)*log(t) dt = -4/9')
x1 = 0.d0
x2 = 1.d0
tm0 = second ()
t1 = quadgsq (fun05, x1, x2, nq1, nq2, wk, xk)
tm1 = second ()
write (6, 3) tm1 - tm0
call qdwrite (6, t1)
t2 = qdreal (-4.d0) / 9.d0
call decmdq (t2 - t1, d1, n1)
write (6, 4) d1, n1

write (6, 16)
16 format (/'Problem 6: Int_0^1 sqrt(1-t^2) dt = pi/4')
x1 = 0.d0
x2 = 1.d0
tm0 = second ()
t1 = quadgsq (fun06, x1, x2, nq1, nq2, wk, xk)
tm1 = second ()
write (6, 3) tm1 - tm0
call qdwrite (6, t1)
t2 = 0.25d0 * qdpi()
call decmdq (t2 - t1, d1, n1)
write (6, 4) d1, n1

write (6, 17)
17 format (/&
  'Functions on finite intervals with integrable singularity at an endpoint.'//&
  'Problem 7: Int_0^1 t/sqrt(1-t^2) dt = 1')
x1 = 0.d0
x2 = 1.d0
tm0 = second ()
t1 = quadgsq (fun07, x1, x2, nq1, nq2, wk, xk)
tm1 = second ()
write (6, 3) tm1 - tm0
call qdwrite (6, t1)
t2 = 1.d0
call decmdq (t2 - t1, d1, n1)
write (6, 4) d1, n1

write (6, 18)
18 format (/'Problem 8: Int_0^1 log(t)^2 dt = 2')
x1 = 0.d0
x2 = 1.d0
tm0 = second ()
t1 = quadgsq (fun08, x1, x2, nq1, nq2, wk, xk)
tm1 = second ()
write (6, 3) tm1 - tm0
call qdwrite (6, t1)
t2 = 2.d0
call decmdq (t2 - t1, d1, n1)
write (6, 4) d1, n1

write (6, 19)
19 format (/'Problem 9: Int_0^(pi/2) log(cos(t)) dt = -pi*log(2)/2')
x1 = 0.d0
x2 = 0.5d0 * qdpi()
tm0 = second ()
t1 = quadgsq (fun09, x1, x2, nq1, nq2, wk, xk)
tm1 = second ()
write (6, 3) tm1 - tm0
call qdwrite (6, t1)
t2 = -0.5d0 * qdpi() * log (qdreal (2.d0))
call decmdq (t2 - t1, d1, n1)
write (6, 4) d1, n1

write (6, 20)
20 format (/'Problem 10: Int_0^(pi/2) sqrt(tan(t)) dt = pi*sqrt(2)/2')
x1 = 0.d0
x2 = 0.5d0 * qdpi()
tm0 = second ()
t1 = quadgsq (fun10, x1, x2, nq1, nq2, wk, xk)
tm1 = second ()
write (6, 3) tm1 - tm0
call qdwrite (6, t1)
t2 = 0.5d0 * qdpi() * sqrt (qdreal (2.d0))
call decmdq (t2 - t1, d1, n1)
write (6, 4) d1, n1

write (6, 21)
21 format (/&
  'Functions on an infinite interval (requiring a two-step solution'//&
  'Problem 11: Int_0^inf 1/(1+t^2) dt = pi/2')
x1 = 0.d0
x2 = 1.d0
tm0 = second ()
t1 = quadgsq (fun11, x1, x2, nq1, nq2, wk, xk)
tm1 = second ()
write (6, 3) tm1 - tm0
call qdwrite (6, t1)
t2 = 0.5d0 * qdpi()
call decmdq (t2 - t1, d1, n1)
write (6, 4) d1, n1

write (6, 22)
22 format (/'Problem 12: Int_0^inf e^(-t)/sqrt(t) dt = sqrt(pi)')
x1 = 0.d0
x2 = 1.d0
tm0 = second ()
t1 = quadgsq (fun12, x1, x2, nq1, nq2, wk, xk)
tm1 = second ()
write (6, 3) tm1 - tm0
call qdwrite (6, t1)
t2 = sqrt (qdpi())
call decmdq (t2 - t1, d1, n1)
write (6, 4) d1, n1

write (6, 23)
23 format (/'Problem 13: Int_0^inf e^(-t^2/2) dt = sqrt(pi/2)')
x1 = 0.d0
x2 = 1.d0
tm0 = second ()
t1 = quadgsq (fun13, x1, x2, nq1, nq2, wk, xk)
tm1 = second ()
write (6, 3) tm1 - tm0
call qdwrite (6, t1)
t2 = sqrt (0.5d0 * qdpi())
call decmdq (t2 - t1, d1, n1)
write (6, 4) d1, n1

write (6, 24)
24 format (/&
  'Oscillatory functions on an infinite interval.'//&
  'Problem 14: Int_0^inf e^(-t)*cos(t) dt = 1/2')
x1 = 0.d0
x2 = 1.d0
tm0 = second ()
t1 = quadgsq (fun14, x1, x2, nq1, nq2, wk, xk)
tm1 = second ()
write (6, 3) tm1 - tm0
call qdwrite (6, t1)
t2 = 0.5d0
call decmdq (t2 - t1, d1, n1)
write (6, 4) d1, n1

write (6, 25)
25 format (/'Problem 15: Int_0^inf sin(t)/t = pi/2')
x1 = 0.d0
x2 = qdpi()
tm0 = second ()
t1 = quadgsq (fun15a, x1, x2, nq1, nq2, wk, xk)
x2 = 1.d0 / qdpi()
t2 = quadgsq (fun15b, x1, x2, nq1, nq2, wk, xk)
t3 = t1 + 40320.d0 * t2 - 1.d0 / qdpi() + 2.d0 / qdpi() ** 3 &
  - 24.d0 / qdpi() ** 5 + 720.d0 / qdpi() ** 7
tm1 = second ()
write (6, 3) tm1 - tm0
t4 = 0.5d0 * qdpi()
call decmdq (t4 - t3, d1, n1)
write (6, 4) d1, n1
write (6, 26)
26 format ('Prob 15 error may be 40,000 X higher than estimated error.')

call f_fpu_fix_end (old_cw)
stop
end

function fun01 (t)

!   fun01(t) = t * log(1+t)

use qdmodule
implicit none
type (qd_real) fun01, t

fun01 = t * log (1.d0 + t)
return
end

function fun02 (t)

!   fun02(t) = t^2 * arctan(t)

use qdmodule
implicit none
type (qd_real) fun02, t

fun02 = t ** 2 * atan (t)
return
end

function fun03 (t)

!   fun03(t) = e^t * cos(t)

use qdmodule
implicit none
type (qd_real) fun03, t

fun03 = exp(t) * cos(t)
return
end

function fun04 (t)

!   fun04(t) = arctan(sqrt(2+t^2))/((1+t^2)sqrt(2+t^2))

use qdmodule
implicit none
type (qd_real) fun04, t, t1

t1 = sqrt (2.d0 + t**2)
fun04 = atan(t1)/((1.d0 + t**2)*t1)
return
end

function fun05 (t)

!    fun05(t) = sqrt(t)*log(t)

use qdmodule
implicit none
type (qd_real) fun05, t

fun05 = sqrt (t) * log (t)
return
end

function fun06 (t)

!    fun06(t) = sqrt(1-t^2)

use qdmodule
implicit none
type (qd_real) fun06, t

fun06 = sqrt (1.d0 - t**2)
return
end

function fun07 (t)

!   fun07(t) = t / sqrt(1-t^2)

use qdmodule
implicit none
type (qd_real) fun07, t

fun07 = t / sqrt (1.d0 - t**2)
return
end

function fun08 (t)

!   fun08(t) = log(t)^2

use qdmodule
implicit none
type (qd_real) fun08, t

fun08 = log (t) ** 2
return
end

function fun09 (t)

!   fun09(t) = log (cos (t))

use qdmodule
implicit none
type (qd_real) fun09, t

fun09 = log (cos (t))
return
end

function fun10 (t)

!   fun10(t) = sqrt(tan(t))

use qdmodule
implicit none
type (qd_real) fun10, t

fun10 = sqrt (tan (t))
return
end

function fun11 (t)

!   fun11(t) = 1/(u^2(1+(1/u-1)^2)) = 1/(1 - 2*u + u^2)

use qdmodule
implicit none
type (qd_real) fun11, t

fun11 = 1.d0 / (1.d0 - 2.d0 * t + 2.d0 * t ** 2)
return
end

function fun12 (t)

!   fun12(t) = e^(-(1/t-1)) / sqrt(t^3 - t^4)

use qdmodule
implicit none
type (qd_real) fun12, t, t1

t1 = 1.d0 / t - 1.d0
fun12 = exp (-t1) / sqrt (t ** 3 - t ** 4)
return
end

function fun13 (t)

!   fun13(t) = e^(-(1/t-1)^2/2) / t^2

use qdmodule
implicit none
type (qd_real) fun13, t, t1

t1 = 1.d0 / t - 1.d0
fun13 = exp (-0.5d0 * t1 ** 2) / t ** 2
return
end

function fun14 (t)

!   fun14(t) = e^(-(1/t-1)) * cos (1/t-1) / t^2

use qdmodule
implicit none
type (qd_real) fun14, t, t1

t1 = 1.d0 / t - 1.d0
fun14 = exp (-t1) * cos (t1) / t ** 2
return
end

function fun15a (t)

!   fun15a(t) = sin(t)/t

use qdmodule
use quadglobal
implicit none
type (qd_real) fun15a, t

fun15a = sin (t) / t
return
end

function fun15b (t)

!   fun15b(t) = t^7 * sin(1/t)

use qdmodule
use quadglobal
implicit none
type (qd_real) fun15b, t

if (abs (t) > 1.d-10) then
  fun15b = t**7 * sin (1.d0 / t)
else
  fun15b = 0.d0
endif
return
end

subroutine initqgsq (nq1, nq2, wk, xk)

!   This subroutine initializes the quadrature arays xk and wk for Gaussian
!   quadrature.  It employs a Newton iteration scheme with a dynamic precision
!   level.  The argument nq2, which is the space allocated for wk and xk in
!   the calling program, should be at least 8 * 2^nq1 + 100, although a higher
!   value may be required, depending on precision level.  Monitor the space
!   figure given in the message below during initialization to be certain.
!   David H Bailey    2002-11-04

use qdmodule
use quadglobal
implicit none
integer i, ierror, ik0, is, j, j1, k, n, nq1, nq2, nwp, nws
double precision pi
parameter (pi = 3.141592653589793238d0)
type (qd_real) eps, r, t1, t2, t3, t4, t5, wk(-1:nq2), xk(-1:nq2)
parameter (ik0 = 100)

if (ndebug >= 1) then
  write (6, 1)
1 format ('initqgsq: Gaussian quadrature initialization')
endif

eps = 1.d-64
wk(-1) = dble (nq1)
wk(0) = 0.d0
xk(0) = 0.d0
wk(1) = dble (nq1)
xk(1) = dble (ik0)
i = ik0

do j = 2, ik0
  wk(j) = 0.d0
  xk(j) = 0.d0
enddo

do k = 1, nq1
  if (ndebug >= 2) write (6, *) k, i, nq2
  n = 3 * 2 ** (k + 1)

  do j = 1, n / 2

!   Compute a double precision estimate of the root.

    is = 0
    r = cos ((pi * (j - 0.25d0)) / (n + 0.5d0))

!   Compute the j-th root of the n-degree Legendre polynomial using Newton's
!   iteration.

100 continue

    t1 = 1.d0
    t2 = 0.d0

    do j1 = 1, n
      t3 = t2
      t2 = t1
      t1 = (dble (2 * j1 - 1) * r * t2 - dble (j1 - 1) * t3) / dble (j1)
    enddo

    t4 = dble (n) * (r * t1 - t2) / (r ** 2  - 1.d0)
    t5 = r
    r = r - t1 / t4

!   Once convergence is achieved at nwp = 3, then start doubling (almost) the
!   working precision level at each iteration until full precision is reached.

    if (abs (r - t5) > eps) goto 100

    i = i + 1
    if (i > nq2) goto 110
    xk(i) = r
    t4 = dble (n) * (r * t1 - t2) / (r ** 2  - 1.d0)
    wk(i) = 2.d0 / ((1.d0 - r ** 2) * t4 ** 2)
  enddo

  xk(k+1) = dble (i)
enddo

xk(-1) = dble (i)
if (ndebug >= 2) then
  write (6, 2) i
2 format ('initqerf: Table spaced used =',i8)
endif
goto 130

110 continue

write (6, 3) nq2
3 format ('initqgsq: Table space parameter is too small; value =',i8)
nerror = 92
goto 130

120 continue

nerror = ierror + 100
write (6, 4) nerror
4 format ('initqgsq: Error in quadrature initialization; code =',i5)

130 continue

return
end

function quadgsq (fun, x1, x2, nq1, nq2, wk, xk)

!   This routine computes the integral of the function fun on the interval
!   [0, 1], with up to nq1 iterations, with a target tolerance of eps.
!   wk and xk are precomputed tables of weights and abscissas, each of size
!   nq2 and of type qd_real.  The function fun is not evaluated at x1 or x2.
!   David H. Bailey    2002-11-04

use qdmodule
use quadglobal
implicit none
integer i, ierror, ik0, k, j, n, nds, nq1, nq2
double precision d1, d2, d3, d4, dplog10q
type (qd_real) a, b, c10, quadgsq, eps, eps1, eps2, err, fmx, fun, h, sum, &
  s1, s2, s3, t1, t2, t3, wk(-1:nq2), xk(-1:nq2), x1, x2, xx1, xx2
external fun, dplog10q
parameter (ik0 = 100)

a = 0.5d0 * (x2 - x1)
b = 0.5d0 * (x2 + x1)
s1 = 0.d0
s2 = 0.d0
c10 = 10.d0
eps = 1.d-64
eps1 = 1.d-61
  if (wk(-1) < dble (nq1)) then
  write (6, 1) nq1
1 format ('quadgsq: quadrature arrays have not been initialized; nq1 =',i6)
  nerror = 70
  goto 130
endif

do k = 1, nq1
  n = 3 * 2 ** (k + 1)
  s3 = s2
  s2 = s1
  fmx = 0.d0
  sum = 0.d0
  i = dble (xk(k))

  do j = 1, n / 2
    i = i + 1
    xx1 = - a * xk(i) + b
    xx2 = a * xk(i) + b
    if (xx1 > x1) then
      t1 = fun (xx1)
    else
      t1 = 0.d0
    endif
    if (xx2 < x2 .and. j + k > 2) then
      t2 = fun (xx2)
    else
      t2 = 0.d0
    endif
    sum = sum + wk(i) * (t1 + t2)
    fmx = max (fmx, abs (t1), abs (t2))
  enddo

  s1 =  a * sum
  eps2 = fmx * eps
  d1 = dplog10q (abs (s1 - s2))
  d2 = dplog10q (abs (s1 - s3))
  d3 = dplog10q (eps2) - 1.d0

  if (k <= 2) then
    err = 1.d0
  elseif (d1 .eq. -9999.d0) then
    err = 0.d0
  else
    d4 = min (0.d0, max (d1 ** 2 / d2, 2.d0 * d1, d3))
    err = c10 ** nint (d4)
  endif

!   Output current integral approximation and error estimate, to 56 dp.

  if (ndebug >= 2) then
    write (6, 2) k, nq1, nint (dplog10q (abs (err)))
2   format ('quadgsq: Iteration',i3,' of',i3,'; est error = 10^',i5, &
      '; approx value =')
    call qdwrite (6, s1)
  endif
  if (k >= 3 .and. err < eps1) goto 130
  if (k >= 3 .and. err < eps2) goto 110
enddo

write (6, 3) nint (dplog10q (abs (err))), nquadl
3 format ('quadgsq: Estimated error = 10^',i5/&
  'Increase Quadlevel for greater accuracy. Current Quadlevel =',i4)
if (err > 1.d-20) then
  write (6, 4)
4 format ('quadgsq: Poor results may be due to singularities at endpoints.'/&
  'If so, try the erf or tanh-sinh quadrature routines (Quadtype = 2 or 3).')
endif
goto 130

110 continue

write (6, 5) nint (dplog10q (abs (err))), ndigits
5 format ('quadgsq: Estimated error = 10^',i5/&
  'Increase working prec (Digits) for greater accuracy. Current Digits =',i4)
goto 130

120 continue

nerror = ierror + 100
write (6, 6) nerror
6 format ('quadgsq: Error in quadrature calculation; code =',i5)
s1 = 0.d0

130 continue

quadgsq = s1
return
end

function dplog10q (a)

!   For input MP value a, this routine returns a DP approximation to log10 (a).

use qdmodule
implicit none
integer ia
double precision da, dplog10q, t1
type (qd_real) a

da = a
if (da .eq. 0.d0) then
  dplog10q = -9999.d0
else
  dplog10q = log10 (abs (da))
endif

100 continue
return
end

subroutine decmdq (a, b, ib)

!   For input MP value a, this routine returns DP b and integer ib such that 
!   a = b * 10^ib, with 1 <= abs (b) < 10 for nonzero a.

use qdmodule
implicit none
integer ia, ib
double precision da, b, t1, xlt
parameter (xlt = 0.3010299956639812d0)
type (qd_real) a

da = a
if (da .ne. 0.d0) then
  t1 = log10 (abs (da))
  ib = t1
  if (t1 .lt. 0.d0) ib = ib - 1
  b = sign (10.d0 ** (t1 - ib), da)
else
  b = 0.d0
  ib = 0
endif

return
end
