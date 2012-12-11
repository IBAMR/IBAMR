module quadglobal

!   This work was supported by the Director, Office of Science, Division
!   of Mathematical, Information, and Computational Sciences of the
!   U.S. Department of Energy under contract number DE-AC02-05CH11231.

!   This program demonstrates the 2D quadrature routine 'quadtsq2d'.

!   David H. Bailey      2007-03-23

!   The function quadtsq2d is suitable to integrate any function that is 
!   continuous, infinitely differentiable and integrable on a 2D finite open 
!   interval.  Singularities or vertical derivatives are permitted on the 
!   boundaries.  This can also be used for certain integrals on infinite
!   intervals, by making a suitable change of variable.

!   The function(s) to be integrated is(are) defined in external function
!   subprogram(s) -- see the sample function subprograms below.  The name(s) of
!   the function subprogram(s) must be included in appropriate type and external
!   statements in the main program.

!   Inputs set in parameter statement below:
!   ndebug Debug level setting.  Default = 2.
!   ndigits1  Primary working precision.  With QD precision, set to 64.
!   nepsilson1  Log10 of the desired tolerance.  Normally set to - ndigits1.
!   nquadl  Max number of phases in quadrature routine; adding 1 increases
!           (possibly doubles) the number of accurate digits in the result,
!           but also roughly *quadruples* the run time.
!   nq2     Space parameter for wk and xk arrays in the calling program.  By
!           default it is set to 8 * 2^nq1.  Increase nq2 if directed by a 
!           message in subroutine initqts.

use qdmodule
implicit none
integer ndebug, nquadl, ndigits1, nepsilon1, nwords1, nq1, nq2, nqmx
parameter (ndebug = 2, nquadl = 6, ndigits1 = 64, nepsilon1 = -64, &
  nwords1 = 2, nq1 = nquadl, nq2 = 8 * 2**nq1)
type (qd_real) xk(-nq2:nq2), wk(-nq2:nq2)
end module

! program tquadtsq2d
subroutine f_main
use qdmodule
use quadglobal
implicit none
integer i, n, n1
double precision dplog10q, d1, d2, second, tm0, tm1
type (qd_real) cat, catalan, err, quadtsq2d, fun01, fun02, fun03, fun04, &
  t1, t2, t3, t4, x1, x2, y1, y2
external quadtsq2d, catalan, fun01, fun02, fun03, fun04, second
integer*4 old_cw

!  This line must be present in DD and QD main programs.

call f_fpu_fix_start (old_cw)

write (6, 1) ndigits1, nepsilon1, nquadl
1 format ('Quadtsq2d test'/'Digits =',i6,'  Epsilon =',i6,'   Quadlevel =',i6)

!   Initialize quadrature tables wk and xk (weights and abscissas).

tm0 = second ()
call initqts
tm1 = second ()
write (6, 2) tm1 - tm0
2 format ('Quadrature initialization completed: cpu time =',f12.6)
cat = catalan ()

!   Begin quadrature tests.

write (6, 11)
11 format (/ &
  'Problem 1: Int_-1^1 Int_-1^1 1/(1+x^2+y^2) dx dy = 4*log(2+sqrt(3))-2*pi/3')
x1 = -1.d0
x2 = 1.d0
y1 = -1.d0
y2 = 1.d0
tm0 = second ()
t1 = quadtsq2d (fun01, x1, x2, y1, y2)
tm1 = second ()
write (6, 3) tm1 - tm0
3 format ('Quadrature completed: CPU time =',f12.6/'Result =')
call qdwrite (6, t1)
t2 = 4.d0 * log (2.d0 + sqrt (qdreal (3.d0))) - 2.d0 * qdpi () / 3.d0
call decmdq (t2 - t1, d1, n1)
write (6, 4) d1, n1
4 format ('Actual error =',f10.6,'x10^',i5)

write (6, 12)
12 format (/&
  'Problem 2: Int_0^pi Int_0^pi log (2-cos(s)-cos(t)) = 4*pi*cat- pi^2*log(2)')
x1 = 0.d0
x2 = qdpi()
y1 = 0.d0
y2 = qdpi()
tm0 = second ()
t1 = quadtsq2d (fun02, x1, x2, y1, y2)
tm1 = second ()
t2 = 4.d0 * qdpi() * cat - qdpi()**2 * log (qdreal (2.d0))
write (6, 3) tm1 - tm0
call qdwrite (6, t1)
call decmdq (t2 - t1, d1, n1)
write (6, 4) d1, n1

write (6, 13)
13 format (/&
  'Problem 3: Int_0^inf Int_0^inf sqrt(x^2+xy+y^2) * exp(-x-y) = 1 + 3/4*log(3)')
x1 = 0.d0
x2 = 1.d0
y1 = 0.d0
y2 = 1.d0
tm0 = second ()
t1 = quadtsq2d (fun03, x1, x2, y1, y2)
tm1 = second ()
t2 = 1.d0 + 0.75d0 * log (qdreal (3.d0))
write (6, 3) tm1 - tm0
call qdwrite (6, t1)
call decmdq (t2 - t1, d1, n1)
write (6, 4) d1, n1

write (6, 14)
14 format (/&
  'Problem 4: Int_0^1 Int_0^1 1/(sqrt((1-x)*(1-y))*(x+y)) dx dy = 4*cat')
x1 = 0.d0
x2 = 1.d0
y1 = 0.d0
y2 = 1.d0
tm0 = second ()
t1 = quadtsq2d (fun04, x1, x2, y1, y2)
tm1 = second ()
t2 = 4.d0 * cat
write (6, 3) tm1 - tm0
call qdwrite (6, t1)
call decmdq (t2 - t1, d1, n1)
write (6, 4) d1, n1

call f_fpu_fix_end (old_cw)
stop
end

function fun01 (s, t)

!   fun01 (s,t) = 1/sqrt[1+s^2+t^2]

use qdmodule
implicit none
type (qd_real) fun01, s, t

fun01 = 1.d0 / sqrt (1.d0 + s**2 + t**2)
return
end

function fun02 (s, t)

!   fun02 (s,t) = log (2 - cos(s) - cos(t))

use qdmodule
implicit none
type (qd_real) fun02, s, t, t1

t1 = 2.d0 - cos (s) - cos (t)
if (t1 > 0.d0) then
  fun02 = log (2.d0 - cos (s) - cos (t))
else
  fun02 = 0.d0
endif
return
end

function fun03 (s, t)

!   fun03 (s,t) = ((1/s-1)^2 + (1/s-1)*(1/t-1) + (1/t-1)^2) 
!                / (s^2 * t^2 * exp(1/s + 1/t - 2)

use qdmodule
implicit none
type (qd_real) fun03, s, t, s1, t1, sq
external dplog10q

if (s > 3.d-3 .and. t > 3.d-3) then
  s1 = 1.d0 / s - 1.d0
  t1 = 1.d0 / t - 1.d0
  sq = sqrt (s1**2 + s1 * t1 + t1**2)
  fun03 = sq / (s**2 * t**2) * exp (-s1 - t1)
else
  fun03 = 0.d0
endif

return
end

function fun04 (s, t)

!   fun04 (s,t) = 1/(sqrt((1-s)*(1-t)) * (s+t))

use qdmodule
implicit none
type (qd_real) fun04, s, t

fun04 = 1.d0 / (sqrt ((1.d0 - s) * (1.d0 - t)) * (s + t))
return
end

subroutine initqts

!   This subroutine initializes the quadrature arays xk and wk for
!   tanh-sinh quadrature.  The argument nq2 is set to 8 * 2^nq1 in module
!   quadmod above.  Increase nq2 if directed by a message produced below.
!   One option here is to set

!   x(k) = 1.d0 / (t3 * t4) 

!   instead of the line below.  This change stores 1 - absicssas in x(k),
!   which can be used to improve quadrature result accuracy in some cases.

use qdmodule
use quadglobal
implicit none
integer i, ierror, iprint, j, k, k1
double precision eps, h
parameter (iprint = 1000)
type (qd_real) pi, p2, t1, t2, t3, t4, t5, u1, u2

pi = qdpi()
eps = 10.d0 ** nepsilon1
if (ndebug >= 1) then
  write (6, 1)
1 format ('initqts: Tanh-sinh quadrature initialization')
endif

p2 = 0.5d0 * pi
h = 0.5d0 ** nq1

do k = 0, nq2
  if (ndebug >= 2 .and. mod (k, iprint) == 0) write (6, *) k, nq2
  t1 = dble (k) * h

!   xk(k) = 1 - tanh (u1) = 1 /(e^u1 * cosh (u1))
!   wk(k) = u2 / cosh (u1)^2
!   where u1 = pi/2 * cosh (t1), u2 = pi/2 * sinh (t1)

  t2 = exp (t1)
  u1 = 0.5d0 * p2 * (t2 + 1.d0 / t2)
  u2 = 0.5d0 * p2 * (t2 - 1.d0 / t2)
  t3 = exp (u2)
  t4 = 0.5d0 * (t3 + 1.d0 / t3)
  xk(k) = 1.d0 - 1.d0 / (t3 * t4)
  wk(k) = u1 / t4 ** 2
  if (k > 0) then
    xk(-k) = - xk(k)
    wk(-k) = wk(k)
  endif
  if (wk(k) < eps) goto 100
enddo

write (6, 2) nq2
2 format ('initqts: Table space parameter is too small; value =',i8)
stop

100 continue

nqmx = k
if (ndebug >= 2) then
  write (6, 3) k
3 format ('initqts: Table spaced used =',i8)
  write (6, *) 'final xk, wk ='
  call qdwrite (6, xk(k), wk(k))
endif

return
end

function quadtsq2d (fun, x1, x2, y1, y2)

!   This performs 2-D tanh-sinh quadrature.  No attempt is made in this code
!   to estimate errors.
!   David H. Bailey     2007-03-23

use qdmodule
use quadglobal
implicit none
integer i, j, k, n
double precision h
type (qd_real) ax, bx, ay, by, quadtsq2d, fun, s1, t1, t2, &
  x1, x2, xx1, y1, y2, yy1
external fun

ax = 0.5d0 * (x2 - x1)
bx = 0.5d0 * (x2 + x1)
ay = 0.5d0 * (y2 - y1)
by = 0.5d0 * (y2 + y1)

if (nqmx == 0) then
  write (6, 1)
1 format ('quadtsq2d: quadrature arrays have not been initialized')
  stop
endif
h = 0.5d0 ** nq1
s1 = 0.d0

do k = -nqmx, nqmx
!   write (6, *) k, nqmx
  yy1 = ay * xk(k) + by

  do j = -nqmx, nqmx
    xx1 = ax * xk(j) + bx
    t1 = fun (xx1, yy1)
    s1 = s1 + wk(j) * wk(k) * t1
  enddo
enddo

quadtsq2d = ax * ay * h**2 * s1
return
end

function catalan ()
use qdmodule
implicit none
integer k
real*8 dk, eps
type (qd_real) catalan, c1, c2, c4, c8, r16, t1, t2, t3
type (qd_real) x1, x2, x3, x4, x5, x6

c1 = 1.d0
c2 = 2.d0
c4 = 4.d0
c8 = 8.d0
r16 = 1.d0 / 16.d0
t1 = 0.d0
t2 = 1.d0
eps = 1.d-64

do k = 0, 10000000
  dk = k
  t3 = t2 * (c8 / (8.d0 * dk + 1.d0) ** 2 + c8 / (8.d0 * dk + 2.d0) ** 2 &
       + c4 / (8.d0 * dk + 3.d0) ** 2 - c2 / (8.d0 * dk + 5.d0) ** 2 &
       - c2 / (8.d0 * dk + 6.d0) ** 2 - c1 / (8.d0 * dk + 7.d0) ** 2)
  t1 = t1 + t3
  t2 = r16 * t2
  if (t3 < 1.d-5 * eps) goto 100
enddo

write (6, *) 'catalan: error - contact author'

100 continue

catalan = 1.d0 / 8.d0 * qdpi() * log (c2) + 1.d0 / 16.d0 * t1
return
end

function dplog10q (a)

!   For input MP value a, this routine returns a DP approximation to log10 (a).

use qdmodule
implicit none
integer ia
double precision da, dplog10q, t1
type (qd_real) a

! call mpmdc (a%mpr, da, ia)
da = a
ia = 0
if (da .eq. 0.d0) then
  dplog10q = -9999.d0
else
  dplog10q = log10 (abs (da)) + ia * log10 (2.d0)
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

! call mpmdc (a%mpr, da, ia)
da = a
ia = 0
if (da .ne. 0.d0) then
  t1 = xlt * ia + log10 (abs (da))
  ib = t1
  if (t1 .lt. 0.d0) ib = ib - 1
  b = sign (10.d0 ** (t1 - ib), da)
else
  b = 0.d0
  ib = 0
endif

return
end
