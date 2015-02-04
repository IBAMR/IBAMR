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
define(REAL,`double precision')dnl
define(INTEGER,`integer')dnl
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Returns the greatest integer less than or equal to x.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      function lagrangian_floor(x)
c
      implicit none
      integer lagrangian_floor
      double precision x
c
      lagrangian_floor = int(x)
      if ( x.lt.0.d0 ) lagrangian_floor = lagrangian_floor - 1
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Returns the interpolation weight phi(r) for the piecewise linear
c     "hat" function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      function lagrangian_piecewise_linear_delta(r)
c
      implicit none
      double precision lagrangian_piecewise_linear_delta,r
c
      if ( r.lt.0.d0 ) then
         r = -r
      endif

      if ( r.lt.1.d0 ) then
         lagrangian_piecewise_linear_delta = 1.d0-r
      else
         lagrangian_piecewise_linear_delta = 0.d0
      endif
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Returns the interpolation weight phi(r) for the broadened 4-point
c     piecewise linear "hat" function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      function lagrangian_wide4_piecewise_linear_delta(r)
c
      implicit none
      external lagrangian_piecewise_linear_delta
      double precision lagrangian_piecewise_linear_delta
      double precision lagrangian_wide4_piecewise_linear_delta,r
c
      lagrangian_wide4_piecewise_linear_delta =
     &     0.5d0*lagrangian_piecewise_linear_delta(0.5d0*r)
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Returns the interpolation weight phi(r) for the piecewise cubic
c     function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      function lagrangian_piecewise_cubic_delta(r)
c
      implicit none
      double precision lagrangian_piecewise_cubic_delta,r
c
      if ( r.lt.0.d0 ) then
         r = -r
      endif

      if ( r.lt.1.d0 ) then
         lagrangian_piecewise_cubic_delta =
     &        1.d0-0.5d0*r-r*r+0.5d0*r*r*r
      else if ( r .lt. 2.d0 ) then
         lagrangian_piecewise_cubic_delta =
     &        1.d0-(11.d0/6.d0)*r+r*r-(1.d0/6.d0)*r*r*r
      else
         lagrangian_piecewise_cubic_delta = 0.d0
      endif
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Returns the interpolation weight phi(r) for the broadened 8-point
c     piecewise cubic function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      function lagrangian_wide8_piecewise_cubic_delta(r)
c
      implicit none
      external lagrangian_piecewise_cubic_delta
      double precision lagrangian_piecewise_cubic_delta
      double precision lagrangian_wide8_piecewise_cubic_delta,r
c
      lagrangian_wide8_piecewise_cubic_delta =
     &     0.5d0*lagrangian_piecewise_cubic_delta(0.5d0*r)
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Returns the interpolation weight phi(r) for the IB 3-point delta
c     function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      function lagrangian_ib_3_delta(r)
c
      implicit none
      double precision lagrangian_ib_3_delta,r
      double precision sixth, third
      parameter (sixth=0.16666666666667d0)
      parameter (third=0.333333333333333d0)
c
      if ( r.lt.0.d0 ) then
         r = -r
      endif

      if ( r.lt.0.5d0 ) then
         lagrangian_ib_3_delta =
     &        third*(1.d0+sqrt(1.d0-3.d0*r*r))
      else if ( r.lt.1.5d0 ) then
         lagrangian_ib_3_delta =
     &        sixth*(5.d0-3.d0*r-sqrt(1.d0-3.d0*(1.d0-r)*(1.d0-r)))
      else
         lagrangian_ib_3_delta = 0.d0
      endif
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Returns the interpolation weight phi(r) for the broadened 6-point
c     version of the IB 3-point delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      function lagrangian_wide6_ib_3_delta(r)
c
      implicit none
      external lagrangian_ib_3_delta
      double precision lagrangian_ib_3_delta
      double precision lagrangian_wide6_ib_3_delta,r
c
      lagrangian_wide6_ib_3_delta = 0.5d0*lagrangian_ib_3_delta(0.5d0*r)
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Returns the interpolation weight phi(r) for the IB 4-point delta
c     function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      function lagrangian_ib_4_delta(r)
c
      implicit none
      double precision lagrangian_ib_4_delta,r,t2,t6
c
      if ( r.lt.0.d0 ) then
         r = -r
      endif

      if ( r.lt.1.d0 ) then
         t2 = r * r
         t6 = sqrt(-0.4D1 * t2 + 0.4D1 * r + 0.1D1)
         lagrangian_ib_4_delta = -r / 0.4D1 + 0.3D1 / 0.8D1 + t6 / 0.8D1
      else if ( r.lt.2.d0 ) then
         t2 = r * r
         t6 = sqrt(0.12D2 * r - 0.7D1 - 0.4D1 * t2)
         lagrangian_ib_4_delta = -r / 0.4D1 + 0.5D1 / 0.8D1 - t6 / 0.8D1
      else
         lagrangian_ib_4_delta = 0.d0
      endif
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Returns the interpolation weight phi(r) for the broadened 8-point
c     version of the IB 4-point delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      function lagrangian_wide8_ib_4_delta(r)
c
      implicit none
      external lagrangian_ib_4_delta
      double precision lagrangian_ib_4_delta
      double precision lagrangian_wide8_ib_4_delta,r
c
      lagrangian_wide8_ib_4_delta = 0.5d0*lagrangian_ib_4_delta(0.5d0*r)
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Returns the interpolation weight phi(r) for the broadened 16-point
c     version of the IB 4-point delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      function lagrangian_wide16_ib_4_delta(r)
c
      implicit none
      external lagrangian_ib_4_delta
      double precision lagrangian_ib_4_delta
      double precision lagrangian_wide16_ib_4_delta,r
c
      lagrangian_wide16_ib_4_delta =
     &     0.25d0*lagrangian_ib_4_delta(0.25d0*r)
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Initializes the lookup table for the interpolation weight phi(r)
c     for the IB 4-point delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_ib_4_table_init(ib_4_table,NTABLE)
c
      implicit none
      double precision lagrangian_ib_4_delta,ib_4_table(0:NTABLE),x
      integer NTABLE,k
c
      do k = 0,NTABLE
         x = 2.d0*dble(k)/dble(NTABLE)
         ib_4_table(k) = lagrangian_ib_4_delta(x)
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Returns the interpolation weights f(x,1/2), f(x,3/2), f(x,5/2),
c     and f(x,7/2) for the one-sided IB four point delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_one_sided_ib_4_delta(f,x)
c
      implicit none
      double precision f(0:3),N,x
c
      f(0) = 0.75d0 - 0.25d0*x
      N = -4.d0*x*x+16.d0*x-14.d0
      if ( N.gt.0.d0 ) then
         f(0) = f(0) - 0.125d0*sqrt(N)
      endif
      f(1) =  1.5d0 - 0.5d0*x - f(0)
      f(2) =  0.5d0           - f(0)
      f(3) = -1.0d0 + 0.5d0*x + f(0)
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Returns the interpolation weights f(x,3/2), f(x,5/2), f(x,7/2),
c     and f(x,9/2) for the one-sided IB four point delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_alt_one_sided_ib_4_delta(f,x)
c
      implicit none
      double precision f(0:3),N,x
c
      f(0) = 1.d0 - 0.25d0*x
      N = -4.d0*x*x+24.d0*x-34.d0
      if ( N.gt.0.d0 ) then
         f(0) = f(0) - 0.125d0*sqrt(N)
      endif
      f(1) =  2.0d0 - 0.5d0*x - f(0)
      f(2) =  0.5d0           - f(0)
      f(3) = -1.5d0 + 0.5d0*x + f(0)
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Returns the interpolation weight phi(r) for the IB 6-point delta
c     function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      function lagrangian_ib_6_delta(r)
c
      implicit none
      double precision lagrangian_ib_6_delta,r,t2,t4,t9,t16
c
      if ( r.lt.0.d0 ) then
         r = -r
      endif

      if ( r.lt.1.d0 ) then
         t2 = r * r
         t4 = t2 * r
         t9 = t2 * t2
         t16 = sqrt(0.729D3 + 0.4752D4 * r - 0.2244D4 * t2 - 0.4680D4 *
     &        t4 + 0.1500D4 * t9 + 0.1008D4 * t9 * r - 0.336D3 * t9 *
     &        t2)
         lagrangian_ib_6_delta =
     &        0.61D2 / 0.112D3 - 0.11D2 / 0.42D2 * r -
     &        0.11D2 / 0.56D2 * t2 + t4 / 0.12D2 + t16 / 0.336D3
      elseif ( r.lt.2.d0 ) then
         t2 = r * r
         t4 = t2 * r
         t9 = t2 * t2
         t16 = sqrt(-0.1431D4 - 0.3744D4 * r + 0.5676D4 * t2 + 0.6120D4
     &        * t4 + 0.3024D4 * t9 * r - 0.8580D4 * t9 - 0.336D3 * t9 *
     &        t2)
         lagrangian_ib_6_delta =
     &        r / 0.84D2 + 0.117D3 / 0.224D3 - 0.23D2
     &        / 0.112D3 * t2 + t4 / 0.24D2 - t16 / 0.224D3
      elseif ( r.lt.3.d0 ) then
         t2 = r * r
         t4 = t2 * r
         t9 = t2 * t2
         t16 = sqrt(-0.10071D5 + 0.54720D5 * r - 0.99444D5 * t2 +
     &        0.77400D5 * t4 + 0.5040D4 * t9 * r - 0.28740D5 * t9 -
     &        0.336D3 * t9 * t2)
         lagrangian_ib_6_delta =
     &        -0.97D2 / 0.84D2 * r + 0.209D3 / 0.224D3
     &        + 0.45D2 / 0.112D3 * t2 - t4 / 0.24D2 + t16 / 0.672D3
      else
         lagrangian_ib_6_delta = 0.d0
      endif
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Initializes the lookup table for the interpolation weight phi(r)
c     for the IB 6-point delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_ib_6_table_init(ib_6_table,NTABLE)
c
      implicit none
      double precision lagrangian_ib_6_delta,ib_6_table(0:NTABLE),x
      integer NTABLE,k
c
      do k = 0,NTABLE
         x = 3.d0*dble(k)/dble(NTABLE)
         ib_6_table(k) = lagrangian_ib_6_delta(x)
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Returns the interpolation weights f(x,1/2), f(x,3/2), f(x,5/2),
c     f(x,7/2), f(x,9/2), and f(x,11/2) for the one-sided IB six point
c     delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_one_sided_ib_6_delta(f,x)
c
      implicit none
      double precision f(0:5),N,t3,t6,x
c
      t3 = x * x
      f(0) = 0.181D3 / 0.112D3 - 0.1067D4 / 0.672D3 * x + 0.13D2 /
     &     0.28D2 * t3 - t3 * x / 0.24D2
      t6 = t3 * t3
      N = -0.73926D5 + 0.228222D6 * x - 0.265269D6 * t3 + 0.148320D6 *
     &     t3 * x - 0.42600D5 * t6 + 0.6048D4 * t6 * x - 0.336D3 * t6 *
     &     t3
      if ( N.gt.0.d0 ) then
         f(0) = f(0) + sqrt(N) / 0.672D3
      endif
      f(1) = 0.19D2 / 0.4D1 - 0.185D3 / 0.48D2 * x - 0.3D1 * f(0) + t3 -
     &     t3 * x / 0.12D2
      f(2) = 0.131D3 / 0.24D2 * x - 0.73D2 / 0.16D2 + 0.2D1 * f(0) -
     &     0.7D1 / 0.4D1 * t3 + t3 * x / 0.6D1
      f(3) = 0.2D1 * f(0) + 0.7D1 / 0.4D1 * x - t3 / 0.4D1 - 0.39D2 /
     &     0.16D2
      f(4) = -0.3D1 * f(0) - 0.131D3 / 0.24D2 * x + 0.81D2 / 0.16D2 +
     &     0.7D1 / 0.4D1 * t3 - t3 * x / 0.6D1
      f(5) = -0.29D2 / 0.16D2 + 0.101D3 / 0.48D2 * x + f(0) - 0.3D1 /
     &     0.4D1 * t3 + t3 * x / 0.12D2
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
