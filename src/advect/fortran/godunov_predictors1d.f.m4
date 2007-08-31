dnl Process this file with m4 to produce FORTRAN source code
define(NDIM,1)dnl
define(REAL,`double precision')dnl
define(INTEGER,`integer')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim1d.i)dnl
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     A Godunov predictor used to predict face and time centered values
c     from cell centered values using a Taylor expansion about each cell
c     center.
c
c     The predictor assumes that Q satisfies an equation of the form
c
c          dQ/dt + u * grad Q = 0
c
c     i.e. Q satisfies an advection equation that is not in conservation
c     form.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine godunov_predict1d(
     &     dx,dt,
     &     limiter,
     &     ifirst0,ilast0,
     &     nQgc0,
     &     Q,
     &     nugc0,
     &     nqhalfgc0,
     &     u0,
     &     qhalf0)
c
      implicit none
include(TOP_SRCDIR/src/fortran/const.i)dnl
include(TOP_SRCDIR/src/advect/fortran/limitertypes.i)dnl
c
c     Functions.
c
      REAL muscldiff,minmod3
c
c     Input.
c
      INTEGER ifirst0,ilast0

      INTEGER nQgc0

      INTEGER nugc0
      INTEGER nqhalfgc0

      INTEGER limiter

      REAL dx(0:NDIM-1),dt

      REAL Q(CELL1dVECG(ifirst,ilast,nQgc))

      REAL u0(FACE1dVECG(ifirst,ilast,nugc))
c
c     Input/Output.
c
      REAL qhalf0(FACE1dVECG(ifirst,ilast,nqhalfgc))
c
c     Local variables.
c
      INTEGER ic0
      REAL Qx,qL,qR
      REAL unorm
c
c     Predict face centered values using a Taylor expansion about each
c     cell center.
c
      if     ( limiter.eq.second_order ) then
c     Employ second order slopes (no limiting).
         Qx = half*(Q(ifirst0-1+1)-Q(ifirst0-1-1))
      elseif ( limiter.eq.fourth_order ) then
         Qx = twothird*(Q(ifirst0-1+1)-Q(ifirst0-1-1))
     &        - sixth*half*(Q(ifirst0-1+2)-Q(ifirst0-1-2))
      elseif ( limiter.eq.mc_limited ) then
c     Employ van Leer's MC limiter.
         Qx = minmod3(
     &        0.5d0*(Q(ifirst0-1+1)-Q(ifirst0-1-1)),
     &        2.0d0*(Q(ifirst0-1  )-Q(ifirst0-1-1)),
     &        2.0d0*(Q(ifirst0-1+1)-Q(ifirst0-1  )))
      elseif ( limiter.eq.muscl_limited ) then
c     Employ Colella's MUSCL limiter.
         Qx = muscldiff(Q(ifirst0-1-2))
      endif

      unorm = 0.5d0*(u0(ifirst0-1  )+u0(ifirst0-1+1))
!     unorm = fourth*fourth*
!    &     ( 9.d0*(u0(ifirst0-1  )+u0(ifirst0-1+1))
!    &     - 1.d0*(u0(ifirst0-1-1)+u0(ifirst0-1+2)) )

      do ic0 = ifirst0-1,ilast0
         qL = Q(ic0  )
     &        + 0.5d0*(+1.d0-unorm*dt/dx(0))*Qx

         if     ( limiter.eq.second_order ) then
            Qx = 0.5d0*(Q(ic0+1+1)-Q(ic0+1-1))
         elseif ( limiter.eq.fourth_order ) then
            Qx = twothird*(Q(ic0+1+1)-Q(ic0+1-1))
     &           - sixth*half*(Q(ic0+1+2)-Q(ic0+1-2))
         elseif ( limiter.eq.mc_limited ) then
            Qx = minmod3(
     &           0.5d0*(Q(ic0+1+1)-Q(ic0+1-1)),
     &           2.0d0*(Q(ic0+1  )-Q(ic0+1-1)),
     &           2.0d0*(Q(ic0+1+1)-Q(ic0+1  )))
         elseif ( limiter.eq.muscl_limited ) then
            Qx = muscldiff(Q(ic0+1-2))
         endif

         unorm = 0.5d0*(u0(ic0+1)+u0(ic0+2))
!        unorm = fourth*fourth*
!    &        ( 9.d0*(u0(ic0+1)+u0(ic0+2))
!    &        - 1.d0*(u0(ic0  )+u0(ic0+3)) )

         qR = Q(ic0+1)
     &        + 0.5d0*(-1.d0-unorm*dt/dx(0))*Qx

         qhalf0(ic0+1) =
     &        0.5d0*(qL+qR)+sign(1.d0,u0(ic0+1))*0.5d0*(qL-qR)
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     A Godunov predictor used to predict face and time centered values
c     from cell centered values using a Taylor expansion about each cell
c     center.
c
c     The predictor assumes that Q satisfies an equation of the form
c
c          dQ/dt + u * grad Q = F
c
c     i.e. Q satisfies an advection equation that is not in conservation
c     form.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine godunov_predictwithsource1d(
     &     dx,dt,
     &     limiter,
     &     ifirst0,ilast0,
     &     nQgc0,
     &     nFgc0,
     &     Q,
     &     F,
     &     nugc0,
     &     nqhalfgc0,
     &     u0,
     &     qhalf0)
c
      implicit none
include(TOP_SRCDIR/src/fortran/const.i)dnl
include(TOP_SRCDIR/src/advect/fortran/limitertypes.i)dnl
c
c     Functions.
c
      REAL muscldiff,minmod3
c
c     Input.
c
      INTEGER ifirst0,ilast0

      INTEGER nQgc0
      INTEGER nFgc0

      INTEGER nugc0
      INTEGER nqhalfgc0

      INTEGER limiter

      REAL dx(0:NDIM-1),dt

      REAL Q(CELL1dVECG(ifirst,ilast,nQgc))
      REAL F(CELL1dVECG(ifirst,ilast,nFgc))

      REAL u0(FACE1dVECG(ifirst,ilast,nugc))
c
c     Input/Output.
c
      REAL qhalf0(FACE1dVECG(ifirst,ilast,nqhalfgc))
c
c     Local variables.
c
      INTEGER ic0
      REAL Qx,qL,qR
      REAL unorm
c
c     Predict face centered values using a Taylor expansion about each
c     cell center.
c
      if     ( limiter.eq.second_order ) then
c     Employ second order slopes (no limiting).
         Qx = half*(Q(ifirst0-1+1)-Q(ifirst0-1-1))
      elseif ( limiter.eq.fourth_order ) then
         Qx = twothird*(Q(ifirst0-1+1)-Q(ifirst0-1-1))
     &        - sixth*half*(Q(ifirst0-1+2)-Q(ifirst0-1-2))
      elseif ( limiter.eq.mc_limited ) then
c     Employ van Leer's MC limiter.
         Qx = minmod3(
     &        0.5d0*(Q(ifirst0-1+1)-Q(ifirst0-1-1)),
     &        2.0d0*(Q(ifirst0-1  )-Q(ifirst0-1-1)),
     &        2.0d0*(Q(ifirst0-1+1)-Q(ifirst0-1  )))
      elseif ( limiter.eq.muscl_limited ) then
c     Employ Colella's MUSCL limiter.
         Qx = muscldiff(Q(ifirst0-1-2))
      endif

      unorm = 0.5d0*(u0(ifirst0-1  )+u0(ifirst0-1+1))
!     unorm = fourth*fourth*
!    &     ( 9.d0*(u0(ifirst0-1  )+u0(ifirst0-1+1))
!    &     - 1.d0*(u0(ifirst0-1-1)+u0(ifirst0-1+2)) )

      do ic0 = ifirst0-1,ilast0
         qL = Q(ic0  )
     &        + 0.5d0*(+1.d0-unorm*dt/dx(0))*Qx
     &        + 0.5d0*dt*F(ic0  )

         if     ( limiter.eq.second_order ) then
            Qx = 0.5d0*(Q(ic0+1+1)-Q(ic0+1-1))
         elseif ( limiter.eq.fourth_order ) then
            Qx = twothird*(Q(ic0+1+1)-Q(ic0+1-1))
     &           - sixth*half*(Q(ic0+1+2)-Q(ic0+1-2))
         elseif ( limiter.eq.mc_limited ) then
            Qx = minmod3(
     &           0.5d0*(Q(ic0+1+1)-Q(ic0+1-1)),
     &           2.0d0*(Q(ic0+1  )-Q(ic0+1-1)),
     &           2.0d0*(Q(ic0+1+1)-Q(ic0+1  )))
         elseif ( limiter.eq.muscl_limited ) then
            Qx = muscldiff(Q(ic0+1-2))
         endif

         unorm = 0.5d0*(u0(ic0+1)+u0(ic0+2))
!        unorm = fourth*fourth*
!    &        ( 9.d0*(u0(ic0+1)+u0(ic0+2))
!    &        - 1.d0*(u0(ic0  )+u0(ic0+3)) )

         qR = Q(ic0+1)
     &        + 0.5d0*(-1.d0-unorm*dt/dx(0))*Qx
     &        + 0.5d0*dt*F(ic0+1)

         qhalf0(ic0+1) =
     &        0.5d0*(qL+qR)+sign(1.d0,u0(ic0+1))*0.5d0*(qL-qR)
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
