c ---------------------------------------------------------------------
c
c Copyright (c) 2011 - 2020 by the IBAMR developers
c All rights reserved.
c
c This file is part of IBAMR.
c
c IBAMR is free software and is distributed under the 3-clause BSD
c license. The full text of the license can be found in the file
c COPYRIGHT at the top level directory of IBAMR.
c
c ---------------------------------------------------------------------

dnl This file only exists to be included into
dnl lagrangian_interaction2d.f and lagrangian_interaction3d.f: i.e., it is
dnl never explicitly compiled on its own so that the delta function kernels
dnl can be inlined. This is necessary since Fortran compilers will (unless
dnl we use something fancy like LTO) only do inlining of functions and
dnl subroutines defined in the same translation unit.

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
c     Returns the interpolation weight phi(r) for the 3-point B-spline.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      function lagrangian_bspline_3_delta(x)
c
      implicit none
      double precision lagrangian_bspline_3_delta,modx,r,r2,x
c
      modx = dabs(x)
      r = modx + 1.5d0
      r2 = r*r
      if (modx .le. 0.5d0) then
         lagrangian_bspline_3_delta = 0.5d0*(-2.d0*r2 + 6.d0*r - 3.d0)
      elseif (modx .le. 1.5d0) then
         lagrangian_bspline_3_delta = 0.5d0*(r2 - 6.d0*r + 9.d0)
      else
         lagrangian_bspline_3_delta = 0.d0
      endif
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Returns the interpolation weight phi(r) for the 4-point B-spline.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      function lagrangian_bspline_4_delta(x)
c
      implicit none
      double precision lagrangian_bspline_4_delta,modx,r,r2,r3,x
c
      modx = dabs(x)
      r = modx + 2.d0
      r2 = r*r
      r3 = r2*r
      if (modx .le. 1.d0) then
         lagrangian_bspline_4_delta = (1.d0/6.d0)*(3.d0*r3 - 24.d0*r2 +
     &        60.d0*r - 44.d0)
      elseif (modx .le. 2.d0) then
         lagrangian_bspline_4_delta = (1.d0/6.d0)*(-r3 + 12.d0*r2 -
     &        48.d0*r + 64.d0)
      else
         lagrangian_bspline_4_delta = 0.d0
      endif
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Returns the interpolation weight phi(r) for the 5-point B-spline.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      function lagrangian_bspline_5_delta(x)
c
      implicit none
      double precision lagrangian_bspline_5_delta,modx,r,r2,r3,r4,x
c
      modx = dabs(x)
      r = modx + 2.5d0;
      r2 = r*r
      r3 = r2*r
      r4 = r3*r
      if (modx .le. 0.5d0) then
         lagrangian_bspline_5_delta = (1.d0/24.d0)*(6.d0*r4 - 60.d0*r3 +
     &        210.d0*r2 - 300.d0*r + 155.d0)
      elseif (modx .le. 1.5d0) then
         lagrangian_bspline_5_delta = (1.d0/24.d0)*(-4.d0*r4 + 60.d0*r3
     &        - 330.d0*r2 + 780.d0*r - 655.d0)
      elseif (modx .le. 2.5d0) then
         lagrangian_bspline_5_delta = (1.d0/24.d0)*(r4 - 20.d0*r3 +
     &        150.d0*r2 - 500.d0*r + 625.d0)
      else
         lagrangian_bspline_5_delta = 0.d0
      endif
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Returns the interpolation weight phi(r) for the 6-point B-spline.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      function lagrangian_bspline_6_delta(x)
c
      implicit none
      double precision lagrangian_bspline_6_delta,modx,r,r2,r3,r4,r5,x
c
      modx = dabs(x)
      r = modx + 3.d0;
      r2 = r*r
      r3 = r2*r
      r4 = r3*r
      r5 = r4*r
      if (modx .le. 1.d0) then
         lagrangian_bspline_6_delta = (1.d0/60.d0)*(2193.d0 - 3465.d0*r
     &        + 2130.d0*r2 - 630.d0*r3 + 90.d0*r4 - 5.d0*r5)
      elseif (modx .le. 2.d0) then
         lagrangian_bspline_6_delta = (1.d0/120.d0)*(-10974.d0 +
     &        12270.d0*r - 5340.d0*r2 + 1140.d0*r3 - 120.d0*r4 + 5.d0
     &        *r5)
      elseif (modx .le. 3.d0) then
         lagrangian_bspline_6_delta = (1.d0/120.d0)*(7776.d0 - 6480.d0*r
     $        + 2160.d0*r2 - 360.d0*r3 + 30.d0*r4 - r5)
      else
         lagrangian_bspline_6_delta = 0.d0
      endif
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
      integer NTABLE,k
      double precision lagrangian_ib_4_delta,ib_4_table(0:NTABLE),x
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
c     Returns the interpolation weight phi(r) for the 5-point IB delta
c     with three continuous derivatives
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      function lagrangian_ib_5_delta(t)
c
      implicit none
      double precision lagrangian_ib_5_delta,r,x,t,K,phi
c
      x = dabs(t)
      K = (38.0d0 - sqrt(69.0d0))/60.0d0
c
      if (x .le. 0.5d0) then
c
          r = x
c
          phi = (136.0d0 - 40.0d0*K - 40.0d0*r**2 + sqrt(2.0d0)
     &          *sqrt(3123.0d0 - 6840.0d0*K + 3600.0d0*(K**2)
     &          - 12440.0d0*(r**2) + 25680.0d0*K*(r**2)
     &          - 12600.0d0*(K**2)*(r**2) + 8080.0d0*(r**4)
     &          - 8400.0d0*K*(r**4) - 1400.0d0*(r**6)))
     &          /280.0d0
c
         lagrangian_ib_5_delta = phi
c
      elseif (x .le. 1.5d0) then
c
          r = x - 1.0d0
c
          phi = (136.0d0 - 40.0d0*K - 40.0d0*r**2 + sqrt(2.0d0)
     &          *sqrt(3123.0d0 - 6840.0d0*K + 3600.0d0*(K**2)
     &          - 12440.0d0*(r**2) + 25680.0d0*K*(r**2)
     &          - 12600.0d0*(K**2)*(r**2) + 8080.0d0*(r**4)
     &          - 8400.0d0*K*(r**4) - 1400.0d0*(r**6)))
     &          /280.0d0
c
         lagrangian_ib_5_delta = (4.0d0 - 4.0d0*phi - K
     &          - 4.0d0*r + 3.0d0*K*r - r**2 + r**3)/6.0d0
c
      elseif (x .le. 2.5d0) then
c
          r = x - 2.0d0
c
          phi = (136.0d0 - 40.0d0*K - 40.0d0*r**2 + sqrt(2.0d0)
     &          *sqrt(3123.0d0 - 6840.0d0*K + 3600.0d0*(K**2)
     &          - 12440.0d0*(r**2) + 25680.0d0*K*(r**2)
     &          - 12600.0d0*(K**2)*(r**2) + 8080.0d0*(r**4)
     &          - 8400.0d0*K*(r**4) - 1400.0d0*(r**6)))
     &          /280.0d0
c
         lagrangian_ib_5_delta = (-2.0d0 + 2.0d0*phi + 2*K +
     &          r - 3.0d0*K*r + 2.0d0*r**2 - r**3)/12.0d0
c
      else
         lagrangian_ib_5_delta = 0.d0
      endif
c
      return
      end
