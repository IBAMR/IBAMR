c
c     Copyright (c) 2002-2010 Boyce Griffith
c
c     Permission is hereby granted, free of charge, to any person
c     obtaining a copy of this software and associated documentation
c     files (the "Software"), to deal in the Software without
c     restriction, including without limitation the rights to use, copy,
c     modify, merge, publish, distribute, sublicense, and/or sell copies
c     of the Software, and to permit persons to whom the Software is
c     furnished to do so, subject to the following conditions:
c
c     The above copyright notice and this permission notice shall be
c     included in all copies or substantial portions of the Software.
c
c     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
c     EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
c     MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
c     NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
c     HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
c     WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
c     OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
c     DEALINGS IN THE SOFTWARE.
c
dnl Process this file with m4 to produce FORTRAN source code
define(REAL,`double precision')dnl
define(INTEGER,`integer')dnl
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     The MUSCL limited 4th order centered difference approximation to
c     dQ/dx.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      REAL function muscldiff(Q)
c
      implicit none
include(TOP_SRCDIR/src/fortran/const.i)dnl
c
c     Input.
c
      REAL Q(-2:2)
c
c     Local variables.
c
      REAL dQ_lim,dQf_left,dQf_rght
c
c     The MUSCL limited 4th order centered difference approximation to
c     dQ/dx.
c
      if ( ((Q(1)-Q(0))*(Q(0)-Q(-1))).gt.0.d0 ) then
         if ( ((Q(0)-Q(-1))*(Q(-1)-Q(-2))).gt.0.d0 ) then
            dQ_lim = 2.d0*dmin1(dabs(Q(0)-Q(-1)),dabs(Q(-1)-Q(-2)))
            dQf_left = dmin1(0.5d0*dabs(Q(0)-Q(-2)),dQ_lim)

            if ( (Q(0)-Q(-2)).lt.0.d0 ) dQf_left = -dQf_left
         else
            dQf_left = 0.d0
         endif

         if ( ((Q(2)-Q(1))*(Q(1)-Q(0))).gt.0.d0 ) then
            dQ_lim = 2.d0*dmin1(dabs(Q(2)-Q(1)),dabs(Q(1)-Q(0)))
            dQf_rght = dmin1(0.5d0*dabs(Q(2)-Q(0)),dQ_lim)

            if ( (Q(2)-Q(0)).lt.0.d0 ) dQf_rght = -dQf_rght
         else
            dQf_rght = 0.d0
         endif

         dQ_lim = 2.d0*dmin1(dabs(Q(1)-Q(0)),dabs(Q(0)-Q(-1)))

         muscldiff = dmin1(twothird*dabs( Q( 1)-0.25d0*dQf_rght
     &                                   -Q(-1)-0.25d0*dQf_left),dQ_lim)

         if ( (Q(1)-Q(-1)).lt.0.d0 ) then
            muscldiff = -muscldiff
         endif
      else
         muscldiff = 0.d0
      endif
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     The minmod function of two arguments.
c
c     The minmod function is a function of two or more arguments that
c     takes the value of the argument with the smallest modulus if all
c     arguments have the same sign.  Otherwise it takes the value zero.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      REAL function minmod2(a,b)
c
      implicit none
include(TOP_SRCDIR/src/fortran/const.i)dnl
c
c     Input.
c
      REAL a,b
c
c     minmod(a,b)
c
      if     ( (a.ge.0.0d0).and.(b.ge.0.0d0) ) then
         minmod2 = dmin1(a,b)
      elseif ( (a.le.0.0d0).and.(b.le.0.0d0) ) then
         minmod2 = dmax1(a,b)
      else
         minmod2 = 0.d0
      endif
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     The minmod function of three arguments.
c
c     The minmod function is a function of two or more arguments that
c     takes the value of the argument with the smallest modulus if all
c     arguments have the same sign.  Otherwise it takes the value zero.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      REAL function minmod3(a,b,c)
c
      implicit none
include(TOP_SRCDIR/src/fortran/const.i)dnl
c
c     Input.
c
      REAL a,b,c
c
c     minmod(a,b,c)
c
      if     ( (a.ge.0.0d0).and.(b.ge.0.0d0).and.(c.ge.0.d0) ) then
         minmod3 = dmin1(a,b,c)
      elseif ( (a.le.0.0d0).and.(b.le.0.0d0).and.(c.le.0.d0) ) then
         minmod3 = dmax1(a,b,c)
      else
         minmod3 = 0.d0
      endif
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
