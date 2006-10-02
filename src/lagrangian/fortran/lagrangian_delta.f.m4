ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Returns the greatest integer less than or equal to x.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      function lagrangian_floor(x)
c
      implicit none
c
      double precision x
      integer lagrangian_floor
c
      lagrangian_floor = int(x)
      if ( x.lt.0.d0 ) lagrangian_floor = lagrangian_floor - 1
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Returns the interpolation weight phi(|r|) for the IB four point
c     delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      function lagrangian_ib4_delta(r)
c
      implicit none
c
      double precision lagrangian_ib4_delta
c
      double precision r,r1,r2
      double precision eighth
      parameter (eighth=0.125d0)
c
      r1 = dabs(r)
      r2 = r1*r1

      if ( r1.lt.1.d0 ) then
         lagrangian_ib4_delta =
     &        eighth*(3.d0 - 2.d0*r1 + dsqrt( 1.d0+ 4.d0*r1-4.d0*r2))
      else if ( r1.lt.2.d0 ) then
         lagrangian_ib4_delta =
     &        eighth*(5.d0 - 2.d0*r1 - dsqrt(-7.d0+12.d0*r1-4.d0*r2))
      else
         lagrangian_ib4_delta = 0.d0
      endif
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Returns the interpolation weight phi(|r|) for the IB four point
c     delta function, spread out over 8 meshwidths.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      function lagrangian_wib4_delta(r)
c
      implicit none
c
      double precision lagrangian_wib4_delta
c
      double precision r,r1,r2
      double precision sixteenth
      parameter (sixteenth=0.0625d0)
c
      r1 = dabs(0.5d0*r)
      r2 = r1*r1

      if ( r1.lt.1.d0 ) then
         lagrangian_wib4_delta =
     &        sixteenth*(3.d0 - 2.d0*r1 + dsqrt( 1.d0+ 4.d0*r1-4.d0*r2))
      else if ( r1.lt.2.d0 ) then
         lagrangian_wib4_delta =
     &        sixteenth*(5.d0 - 2.d0*r1 - dsqrt(-7.d0+12.d0*r1-4.d0*r2))
      else
         lagrangian_wib4_delta = 0.d0
      endif
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Returns the interpolation weight phi(|r|) for the IB six point
c     delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      function lagrangian_ib6_delta(r)
c
      implicit none
c
      double precision lagrangian_ib6_delta
c
      double precision r,r1,r2,r3,r4,r5,r6
      double precision fac00,fac01,fac02,fac03,fac04
      double precision fac10,fac11,fac12,fac13,fac14
      double precision fac20,fac21,fac22,fac23,fac24
      parameter (fac00= 5.446428571428571d-01)
      parameter (fac01=-2.619047619047619d-01)
      parameter (fac02=-1.964285714285714d-01)
      parameter (fac03= 8.333333333333333d-02)
      parameter (fac04= 5.154913117764516d-03)
      parameter (fac10= 5.223214285714286d-01)
      parameter (fac11= 1.190476190476190d-02)
      parameter (fac12=-2.053571428571428d-01)
      parameter (fac13= 4.166666666666666d-02)
      parameter (fac14=-7.732369676646774d-03)
      parameter (fac20= 9.330357142857143d-01)
      parameter (fac21=-1.154761904761905d+00)
      parameter (fac22= 4.017857142857143d-01)
      parameter (fac23=-4.166666666666666d-02)
      parameter (fac24= 2.577456558882258d-03)
c
      r1 = dabs(r)
      r2 = r1*r1
      r3 = r1*r2
      r4 = r1*r3
      r5 = r1*r4
      r6 = r1*r5

      if ( r1.lt.1.d0 ) then
         lagrangian_ib6_delta =
     &        fac00 + fac01*r1 + fac02*r2 + fac03*r3 +
     &        fac04*dsqrt(243.d0 + 1584.d0*r1 - 748.d0*r2
     &        - 1560.d0*r3 + 500.d0*r4 + 336.d0*r5 - 112.d0*r6)
      else if ( r1.lt.2.d0 ) then
         lagrangian_ib6_delta =
     &        fac10 + fac11*r1 + fac12*r2 + fac13*r3 +
     &        fac14*dsqrt(-477.d0 - 1248.d0*r1 + 1892.d0*r2
     &        + 2040.d0*r3 - 2860.d0*r4 + 1008.d0*r5 - 112.d0*r6)
      else if ( r1.lt.3.d0 ) then
         lagrangian_ib6_delta =
     &        fac20 + fac21*r1 + fac22*r2 + fac23*r3 +
     &        fac24*dsqrt(-3357.d0 + 18240.d0*r1 - 33148.d0*r2
     &        + 25800.d0*r3 - 9580.d0*r4 + 1680.d0*r5 - 112.d0*r6)
      else
         lagrangian_ib6_delta = 0.d0
      endif
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Returns the interpolation weight phi(|r|) for piecewise cubic
c     interpolation.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      function lagrangian_pwcubic_delta(r)
c
      implicit none
c
      double precision lagrangian_pwcubic_delta
c
      double precision r,r1,r2,r3
      double precision sixth,elevensixth
      parameter (sixth=0.16666666666667d0)
      parameter (elevensixth=1.833333333333333d0)
c
      r1 = dabs(r)
      r2 = r1*r1
      r3 = r2*r1

      if ( r1.lt.1.d0 ) then
         lagrangian_pwcubic_delta =
     &        1.d0 - 0.5d0*r1 - r2 + 0.5d0*r3
      else if ( r1.lt.2.d0 ) then
         lagrangian_pwcubic_delta =
     &        1.d0 - elevensixth*r1 + r2 - sixth*r3
      else
         lagrangian_pwcubic_delta = 0.d0
      endif
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
