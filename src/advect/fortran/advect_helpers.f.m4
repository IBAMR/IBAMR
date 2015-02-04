c
c     Copyright (c) 2002-2014, Boyce Griffith
c     All rights reserved.
c
c     Redistribution and use in source and binary forms, with or without
c     modification, are permitted provided that the following conditions
c     are met:
c
c        * Redistributions of source code must retain the above
c          copyright notice, this list of conditions and the following
c          disclaimer.
c
c        * Redistributions in binary form must reproduce the above
c          copyright notice, this list of conditions and the following
c          disclaimer in the documentation and/or other materials
c          provided with the distribution.
c
c        * Neither the name of The University of North Carolina nor the
c          names of its contributors may be used to endorse or promote
c          products derived from this software without specific prior
c          written permission.
c
c     THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
c     CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
c     INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
c     MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
c     DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS
c     BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
c     EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
c     TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
c     DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
c     ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
c     TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
c     THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
c     SUCH DAMAGE.
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
      REAL twothird
      parameter (twothird=0.66666666666667d0)
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
c
c     Input.
c
      REAL a,b,c
c
c     minmod(a,b,c)
c
      if     ( (a.ge.0.0d0).and.(b.ge.0.0d0).and.(c.ge.0.0d0) ) then
         minmod3 = dmin1(a,b,c)
      elseif ( (a.le.0.0d0).and.(b.le.0.0d0).and.(c.le.0.0d0) ) then
         minmod3 = dmax1(a,b,c)
      else
         minmod3 = 0.d0
      endif
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     The maxmod function of two arguments.
c
c     The maxmod function is a function of two or more arguments that
c     takes the value of the argument with the largest modulus if all
c     arguments have the same sign.  Otherwise it takes the value zero.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      REAL function maxmod2(a,b)
c
      implicit none
c
c     Input.
c
      REAL a,b
c
c     maxmod(a,b)
c
      if     ( (a.ge.0.0d0).and.(b.ge.0.0d0) ) then
         maxmod2 = dmax1(a,b)
      elseif ( (a.le.0.0d0).and.(b.le.0.0d0) ) then
         maxmod2 = dmin1(a,b)
      else
         maxmod2 = 0.d0
      endif
c
      return
      end
c
