dnl Process this file with m4 to produce FORTRAN source code
define(NDIM,2)dnl
define(REAL,`double precision')dnl
define(INTEGER,`integer')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim2d.i)dnl
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
      subroutine godunov_predict2d(
     &     dx,dt,
     &     limiter,
     &     ifirst0,ilast0,ifirst1,ilast1,
     &     nQgc0,nQgc1,
     &     Q,Qscratch1,
     &     nugc0,nugc1,
     &     nqhalfgc0,nqhalfgc1,
     &     u0,u1,
     &     qtemp0,qtemp1,
     &     qhalf0,qhalf1)
c
      implicit none
include(TOP_SRCDIR/src/fortran/const.i)dnl
include(TOP_SRCDIR/src/advect/fortran/limitertypes.i)dnl
c
c     Input.
c
      INTEGER ifirst0,ilast0,ifirst1,ilast1

      INTEGER nQgc0,nQgc1

      INTEGER nugc0,nugc1
      INTEGER nqhalfgc0,nqhalfgc1

      INTEGER limiter

      REAL dx(0:NDIM-1),dt

      REAL Q(CELL2dVECG(ifirst,ilast,nQgc))
      REAL Qscratch1(ifirst1-nQgc1:ilast1+nQgc1,
     &               ifirst0-nQgc0:ilast0+nQgc0)

      REAL u0(FACE2d0VECG(ifirst,ilast,nugc))
      REAL u1(FACE2d1VECG(ifirst,ilast,nugc))

      REAL qtemp0(FACE2d0VECG(ifirst,ilast,nqhalfgc))
      REAL qtemp1(FACE2d1VECG(ifirst,ilast,nqhalfgc))
c
c     Input/Output.
c
      REAL qhalf0(FACE2d0VECG(ifirst,ilast,nqhalfgc))
      REAL qhalf1(FACE2d1VECG(ifirst,ilast,nqhalfgc))
c
c     Local variables.
c
      INTEGER ic0,ic1
c
c     For ease of implementation, we make a copy of Q with permuted
c     indices.
c
      do ic1 = ifirst1-nQgc1,ilast1+nQgc1
         do ic0 = ifirst0-nQgc0,ilast0+nQgc0
            Qscratch1(ic1,ic0) = Q(ic0,ic1)
         enddo
      enddo
c
c     Compute temporary predicted values on cell faces.
c
c     In this computation, normal derivatives are approximated by
c     (limited) centered differences.  Transverse derivatives are not
c     included.
c
      call godunov_predictnormal2d( ! predict values on the x-faces
     &     dx(0),dt,
     &     limiter,
     &     ifirst0,ilast0,ifirst1,ilast1,
     &     nQgc0,nQgc1,
     &     Q,
     &     nugc0,nugc1,
     &     nqhalfgc0,nqhalfgc1,
     &     u0,
     &     qtemp0)

      call godunov_predictnormal2d( ! predict values on the y-faces
     &     dx(1),dt,
     &     limiter,
     &     ifirst1,ilast1,ifirst0,ilast0,
     &     nQgc1,nQgc0,
     &     Qscratch1,
     &     nugc1,nugc0,
     &     nqhalfgc1,nqhalfgc0,
     &     u1,
     &     qtemp1)
c
c     Compute final predicted values on cell faces.
c
c     This computation approximates transverse derivatives by centered
c     differences of the "temporary" predicted values.
c
      call godunov_transversefix2d( ! update values on the x-faces
     &     dx(1),dt,
     &     ifirst0,ilast0,ifirst1,ilast1,
     &     nugc0,nugc1,
     &     nqhalfgc0,nqhalfgc1,
     &     u0,u1,
     &     qtemp0,qtemp1,
     &     qhalf0)

      call godunov_transversefix2d( ! update values on the y-faces
     &     dx(0),dt,
     &     ifirst1,ilast1,ifirst0,ilast0,
     &     nugc1,nugc0,
     &     nqhalfgc1,nqhalfgc0,
     &     u1,u0,
     &     qtemp1,qtemp0,
     &     qhalf1)
c
      return
      end
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
      subroutine godunov_predictwithsource2d(
     &     dx,dt,
     &     limiter,
     &     ifirst0,ilast0,ifirst1,ilast1,
     &     nQgc0,nQgc1,
     &     nFgc0,nFgc1,
     &     Q,Qscratch1,
     &     F,Fscratch1,
     &     nugc0,nugc1,
     &     nqhalfgc0,nqhalfgc1,
     &     u0,u1,
     &     qtemp0,qtemp1,
     &     qhalf0,qhalf1)
c
      implicit none
include(TOP_SRCDIR/src/fortran/const.i)dnl
include(TOP_SRCDIR/src/advect/fortran/limitertypes.i)dnl
c
c     Input.
c
      INTEGER ifirst0,ilast0,ifirst1,ilast1

      INTEGER nQgc0,nQgc1
      INTEGER nFgc0,nFgc1

      INTEGER nugc0,nugc1
      INTEGER nqhalfgc0,nqhalfgc1

      INTEGER limiter

      REAL dx(0:NDIM-1),dt

      REAL Q(CELL2dVECG(ifirst,ilast,nQgc))
      REAL Qscratch1(ifirst1-nQgc1:ilast1+nQgc1,
     &               ifirst0-nQgc0:ilast0+nQgc0)

      REAL F(CELL2dVECG(ifirst,ilast,nFgc))
      REAL Fscratch1(ifirst1-nFgc1:ilast1+nFgc1,
     &               ifirst0-nFgc0:ilast0+nFgc0)

      REAL u0(FACE2d0VECG(ifirst,ilast,nugc))
      REAL u1(FACE2d1VECG(ifirst,ilast,nugc))

      REAL qtemp0(FACE2d0VECG(ifirst,ilast,nqhalfgc))
      REAL qtemp1(FACE2d1VECG(ifirst,ilast,nqhalfgc))
c
c     Input/Output.
c
      REAL qhalf0(FACE2d0VECG(ifirst,ilast,nqhalfgc))
      REAL qhalf1(FACE2d1VECG(ifirst,ilast,nqhalfgc))
c
c     Local variables.
c
      INTEGER ic0,ic1
c
c     For ease of implementation, we make copies of Q and F with
c     permuted indices.
c
      do ic1 = ifirst1-nQgc1,ilast1+nQgc1
         do ic0 = ifirst0-nQgc0,ilast0+nQgc0
            Qscratch1(ic1,ic0) = Q(ic0,ic1)
         enddo
      enddo

      do ic1 = ifirst1-nFgc1,ilast1+nFgc1
         do ic0 = ifirst0-nFgc0,ilast0+nFgc0
            Fscratch1(ic1,ic0) = F(ic0,ic1)
         enddo
      enddo
c
c     Compute temporary predicted values on cell faces.
c
c     In this computation, normal derivatives are approximated by
c     (limited) centered differences.  Transverse derivatives are not
c     included.
c
      call godunov_predictnormalwithsource2d( ! predict values on the x-faces
     &     dx(0),dt,
     &     limiter,
     &     ifirst0,ilast0,ifirst1,ilast1,
     &     nQgc0,nQgc1,
     &     nFgc0,nFgc1,
     &     Q,
     &     F,
     &     nugc0,nugc1,
     &     nqhalfgc0,nqhalfgc1,
     &     u0,
     &     qtemp0)

      call godunov_predictnormalwithsource2d( ! predict values on the y-faces
     &     dx(1),dt,
     &     limiter,
     &     ifirst1,ilast1,ifirst0,ilast0,
     &     nQgc1,nQgc0,
     &     nFgc1,nFgc0,
     &     Qscratch1,
     &     Fscratch1,
     &     nugc1,nugc0,
     &     nqhalfgc1,nqhalfgc0,
     &     u1,
     &     qtemp1)
c
c     Compute final predicted values on cell faces.
c
c     This computation approximates transverse derivatives by centered
c     differences of the "temporary" predicted values.
c
      call godunov_transversefix2d( ! update values on the x-faces
     &     dx(1),dt,
     &     ifirst0,ilast0,ifirst1,ilast1,
     &     nugc0,nugc1,
     &     nqhalfgc0,nqhalfgc1,
     &     u0,u1,
     &     qtemp0,qtemp1,
     &     qhalf0)

      call godunov_transversefix2d( ! update values on the y-faces
     &     dx(0),dt,
     &     ifirst1,ilast1,ifirst0,ilast0,
     &     nugc1,nugc0,
     &     nqhalfgc1,nqhalfgc0,
     &     u1,u0,
     &     qtemp1,qtemp0,
     &     qhalf1)
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Subtract off the gradient of a scalar to enforce
c     incompressibility for a predicted value.
c
c     NOTE: The gradtype variable specifies which component of the
c     gradient we should be using to enforce incompressibility.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine godunov_incompressibilityfix2d(
     &     gradtype,
     &     ifirst0,ilast0,ifirst1,ilast1,
     &     ngradgc0,ngradgc1,
     &     nqhalfgc0,nqhalfgc1,
     &     grad0,grad1,
     &     qhalf0,qhalf1)
c
      implicit none
include(TOP_SRCDIR/src/fortran/const.i)dnl
c
c     Input.
c
      INTEGER ifirst0,ilast0,ifirst1,ilast1

      INTEGER ngradgc0,ngradgc1
      INTEGER nqhalfgc0,nqhalfgc1

      INTEGER gradtype

      REAL grad0(FACE2d0VECG(ifirst,ilast,ngradgc))
      REAL grad1(FACE2d1VECG(ifirst,ilast,ngradgc))
c
c     Input/Output.
c
      REAL qhalf0(FACE2d0VECG(ifirst,ilast,nqhalfgc))
      REAL qhalf1(FACE2d1VECG(ifirst,ilast,nqhalfgc))
c
c     Local variables.
c
      INTEGER gradtype0,gradtype1
c
      if ( gradtype.eq.0 ) then
         gradtype0 = 0
         gradtype1 = 1
      else
         gradtype0 = 1
         gradtype1 = 0
      endif

      call godunov_gradfix2d(
     &     gradtype0,
     &     ifirst0,ilast0,ifirst1,ilast1,
     &     ngradgc0,ngradgc1,
     &     nqhalfgc0,nqhalfgc1,
     &     grad0,grad1,
     &     qhalf0)

      call godunov_gradfix2d(
     &     gradtype1,
     &     ifirst1,ilast1,ifirst0,ilast0,
     &     ngradgc1,ngradgc0,
     &     nqhalfgc1,nqhalfgc0,
     &     grad1,grad0,
     &     qhalf1)
c
      return
      end
c
      subroutine godunov_predictnormal2d(
     &     dx0,dt,
     &     limiter,
     &     ifirst0,ilast0,ifirst1,ilast1,
     &     nQgc0,nQgc1,
     &     Q,
     &     nugc0,nugc1,
     &     nqhalfgc0,nqhalfgc1,
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
      INTEGER ifirst0,ilast0,ifirst1,ilast1

      INTEGER nQgc0,nQgc1

      INTEGER nugc0,nugc1
      INTEGER nqhalfgc0,nqhalfgc1

      INTEGER limiter

      REAL dx0,dt

      REAL Q(CELL2dVECG(ifirst,ilast,nQgc))

      REAL u0(FACE2d0VECG(ifirst,ilast,nugc))
c
c     Input/Output.
c
      REAL qhalf0(FACE2d0VECG(ifirst,ilast,nqhalfgc))
c
c     Local variables.
c
      INTEGER ic0,ic1
      REAL Qx,qL,qR
      REAL unorm
c
c     Predict face centered values using a Taylor expansion about each
c     cell center.
c
c     (Limited) centered differences are used to approximate normal
c     derivatives.  Transverse derivatives are NOT included.
c
      Qx = 0.d0

      do ic1 = ifirst1-1,ilast1+1

         if     ( limiter.eq.second_order ) then
c     Employ second order slopes (no limiting).
            Qx = half*(Q(ifirst0-1+1,ic1)-Q(ifirst0-1-1,ic1))
         elseif ( limiter.eq.fourth_order ) then
            Qx = twothird*(Q(ifirst0-1+1,ic1)-Q(ifirst0-1-1,ic1))
     &           - sixth*half*(Q(ifirst0-1+2,ic1)-Q(ifirst0-1-2,ic1))
         elseif ( limiter.eq.mc_limited ) then
c     Employ van Leer's MC limiter.
            Qx = minmod3(
     &           0.5d0*(Q(ifirst0-1+1,ic1)-Q(ifirst0-1-1,ic1)),
     &           2.0d0*(Q(ifirst0-1  ,ic1)-Q(ifirst0-1-1,ic1)),
     &           2.0d0*(Q(ifirst0-1+1,ic1)-Q(ifirst0-1  ,ic1)))
         elseif ( limiter.eq.muscl_limited ) then
c     Employ Colella's MUSCL limiter.
            Qx = muscldiff(Q(ifirst0-1-2,ic1))
         endif

!     unorm = 0.5d0*(u0(ifirst0-1  ,ic1)+u0(ifirst0-1+1,ic1))
         unorm = fourth*fourth*
     &        ( 9.d0*(u0(ifirst0-1  ,ic1)+u0(ifirst0-1+1,ic1))
     &        - 1.d0*(u0(ifirst0-1-1,ic1)+u0(ifirst0-1+2,ic1)) )

         do ic0 = ifirst0-1,ilast0
            qL = Q(ic0  ,ic1)
     &           + 0.5d0*(1.d0-unorm*dt/dx0)*Qx

            if     ( limiter.eq.second_order ) then
               Qx = 0.5d0*(Q(ic0+1+1,ic1)-Q(ic0+1-1,ic1))
            elseif ( limiter.eq.fourth_order ) then
               Qx = twothird*(Q(ic0+1+1,ic1)-Q(ic0+1-1,ic1))
     &              - sixth*half*(Q(ic0+1+2,ic1)-Q(ic0+1-2,ic1))
            elseif ( limiter.eq.mc_limited ) then
               Qx = minmod3(
     &              0.5d0*(Q(ic0+1+1,ic1)-Q(ic0+1-1,ic1)),
     &              2.0d0*(Q(ic0+1  ,ic1)-Q(ic0+1-1,ic1)),
     &              2.0d0*(Q(ic0+1+1,ic1)-Q(ic0+1  ,ic1)))
            elseif ( limiter.eq.muscl_limited ) then
               Qx = muscldiff(Q(ic0+1-2,ic1))
            endif

!     unorm = 0.5d0*(u0(ic0+1,ic1)+u0(ic0+2,ic1))
            unorm = fourth*fourth*
     &           ( 9.d0*(u0(ic0+1,ic1)+u0(ic0+2,ic1))
     &           - 1.d0*(u0(ic0  ,ic1)+u0(ic0+3,ic1)) )

            qR = Q(ic0+1,ic1)
     &           - 0.5d0*(1.d0+unorm*dt/dx0)*Qx

            qhalf0(ic0+1,ic1) =
     &           0.5d0*(qL+qR)+sign(1.d0,u0(ic0+1,ic1))*0.5d0*(qL-qR)
         enddo

      enddo
c
      return
      end
c
      subroutine godunov_predictnormalwithsource2d(
     &     dx0,dt,
     &     limiter,
     &     ifirst0,ilast0,ifirst1,ilast1,
     &     nQgc0,nQgc1,
     &     nFgc0,nFgc1,
     &     Q,
     &     F,
     &     nugc0,nugc1,
     &     nqhalfgc0,nqhalfgc1,
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
      INTEGER ifirst0,ilast0,ifirst1,ilast1

      INTEGER nQgc0,nQgc1
      INTEGER nFgc0,nFgc1

      INTEGER nugc0,nugc1
      INTEGER nqhalfgc0,nqhalfgc1

      INTEGER limiter

      REAL dx0,dt

      REAL Q(CELL2dVECG(ifirst,ilast,nQgc))
      REAL F(CELL2dVECG(ifirst,ilast,nFgc))

      REAL u0(FACE2d0VECG(ifirst,ilast,nugc))
c
c     Input/Output.
c
      REAL qhalf0(FACE2d0VECG(ifirst,ilast,nqhalfgc))
c
c     Local variables.
c
      INTEGER ic0,ic1
      REAL Qx,qL,qR
      REAL unorm
c
c     Predict face centered values using a Taylor expansion about each
c     cell center.
c
c     (Limited) centered differences are used to approximate normal
c     derivatives.  Transverse derivatives are NOT included.
c
      Qx = 0.d0

      do ic1 = ifirst1-1,ilast1+1

         if     ( limiter.eq.second_order ) then
c     Employ second order slopes (no limiting).
            Qx = half*(Q(ifirst0-1+1,ic1)-Q(ifirst0-1-1,ic1))
         elseif ( limiter.eq.fourth_order ) then
            Qx = twothird*(Q(ifirst0-1+1,ic1)-Q(ifirst0-1-1,ic1))
     &           - sixth*half*(Q(ifirst0-1+2,ic1)-Q(ifirst0-1-2,ic1))
         elseif ( limiter.eq.mc_limited ) then
c     Employ van Leer's MC limiter.
            Qx = minmod3(
     &           0.5d0*(Q(ifirst0-1+1,ic1)-Q(ifirst0-1-1,ic1)),
     &           2.0d0*(Q(ifirst0-1  ,ic1)-Q(ifirst0-1-1,ic1)),
     &           2.0d0*(Q(ifirst0-1+1,ic1)-Q(ifirst0-1  ,ic1)))
         elseif ( limiter.eq.muscl_limited ) then
c     Employ Colella's MUSCL limiter.
            Qx = muscldiff(Q(ifirst0-1-2,ic1))
         endif

!     unorm = 0.5d0*(u0(ifirst0-1,ic1)+u0(ifirst0-1+1,ic1))
         unorm = fourth*fourth*
     &        ( 9.d0*(u0(ifirst0-1  ,ic1)+u0(ifirst0-1+1,ic1))
     &        - 1.d0*(u0(ifirst0-1-1,ic1)+u0(ifirst0-1+2,ic1)) )

         do ic0 = ifirst0-1,ilast0
            qL = Q(ic0  ,ic1)
     &           + 0.5d0*(1.d0-unorm*dt/dx0)*Qx
     &           + 0.5d0*dt*F(ic0  ,ic1)

            if     ( limiter.eq.second_order ) then
               Qx = 0.5d0*(Q(ic0+1+1,ic1)-Q(ic0+1-1,ic1))
            elseif ( limiter.eq.fourth_order ) then
               Qx = twothird*(Q(ic0+1+1,ic1)-Q(ic0+1-1,ic1))
     &              - sixth*half*(Q(ic0+1+2,ic1)-Q(ic0+1-2,ic1))
            elseif ( limiter.eq.mc_limited ) then
               Qx = minmod3(
     &              0.5d0*(Q(ic0+1+1,ic1)-Q(ic0+1-1,ic1)),
     &              2.0d0*(Q(ic0+1  ,ic1)-Q(ic0+1-1,ic1)),
     &              2.0d0*(Q(ic0+1+1,ic1)-Q(ic0+1  ,ic1)))
            elseif ( limiter.eq.muscl_limited ) then
               Qx = muscldiff(Q(ic0+1-2,ic1))
            endif

!     unorm = 0.5d0*(u0(ic0+1,ic1)+u0(ic0+2,ic1))
            unorm = fourth*fourth*
     &           ( 9.d0*(u0(ic0+1,ic1)+u0(ic0+2,ic1))
     &           - 1.d0*(u0(ic0  ,ic1)+u0(ic0+3,ic1)) )

            qR = Q(ic0+1,ic1)
     &           - 0.5d0*(1.d0+unorm*dt/dx0)*Qx
     &           + 0.5d0*dt*F(ic0+1,ic1)

            qhalf0(ic0+1,ic1) =
     &           0.5d0*(qL+qR)+sign(1.d0,u0(ic0+1,ic1))*0.5d0*(qL-qR)
         enddo

      enddo
c
      return
      end
c
      subroutine godunov_transversefix2d(
     &     dx1,dt,
     &     ifirst0,ilast0,ifirst1,ilast1,
     &     nugc0,nugc1,
     &     nqhalfgc0,nqhalfgc1,
     &     u0,u1,
     &     qtemp0,qtemp1,
     &     qhalf0)
c
      implicit none
include(TOP_SRCDIR/src/fortran/const.i)dnl
c
c     Input.
c
      INTEGER ifirst0,ilast0,ifirst1,ilast1

      INTEGER nugc0,nugc1
      INTEGER nqhalfgc0,nqhalfgc1

      REAL dx1,dt

      REAL u0(FACE2d0VECG(ifirst,ilast,nugc))
      REAL u1(FACE2d1VECG(ifirst,ilast,nugc))

      REAL qtemp0(FACE2d0VECG(ifirst,ilast,nqhalfgc))
      REAL qtemp1(FACE2d1VECG(ifirst,ilast,nqhalfgc))
c
c     Input/Output.
c
      REAL qhalf0(FACE2d0VECG(ifirst,ilast,nqhalfgc))
c
c     Local variables.
c
      INTEGER ic0,ic1
      REAL Qy,qL_diff,qR_diff
      REAL vtan
c
c     Add transverse derivitives by taking centered differences of
c     temporary predicted values.
c
      do ic1 = ifirst1,ilast1

!     vtan = 0.5d0*(u1(ic1,ifirst0-1)+u1(ic1+1,ifirst0-1))
         vtan = fourth*fourth*
     &        ( 9.d0*(u1(ic1  ,ifirst0-1)+u1(ic1+1,ifirst0-1))
     &        - 1.d0*(u1(ic1-1,ifirst0-1)+u1(ic1+2,ifirst0-1)) )

         do ic0 = ifirst0-1,ilast0
            Qy = qtemp1(ic1+1,ic0)-qtemp1(ic1,ic0)

            qL_diff =
     &           - 0.5d0*dt*vtan*Qy/dx1

!     vtan = 0.5d0*(u1(ic1,ic0+1)+u1(ic1+1,ic0+1))
            vtan = fourth*fourth*
     &           ( 9.d0*(u1(ic1  ,ic0+1)+u1(ic1+1,ic0+1))
     &           - 1.d0*(u1(ic1-1,ic0+1)+u1(ic1+2,ic0+1)) )

            Qy = qtemp1(ic1+1,ic0+1)-qtemp1(ic1,ic0+1)

            qR_diff =
     &           - 0.5d0*dt*vtan*Qy/dx1

            qhalf0(ic0+1,ic1) = qtemp0(ic0+1,ic1) +
     &           0.5d0*(qL_diff+qR_diff)+
     &           sign(1.d0,u0(ic0+1,ic1))*0.5d0*(qL_diff-qR_diff)
         enddo

      enddo
c
      return
      end
c
      subroutine godunov_gradfix2d(
     &     gradtype,
     &     ifirst0,ilast0,ifirst1,ilast1,
     &     ngradgc0,ngradgc1,
     &     nqhalfgc0,nqhalfgc1,
     &     grad0,grad1,
     &     qhalf0)
c
      implicit none
include(TOP_SRCDIR/src/fortran/const.i)dnl
c
c     Input.
c
      INTEGER ifirst0,ilast0,ifirst1,ilast1

      INTEGER ngradgc0,ngradgc1
      INTEGER nqhalfgc0,nqhalfgc1

      INTEGER gradtype

      REAL grad0(FACE2d0VECG(ifirst,ilast,ngradgc))
      REAL grad1(FACE2d1VECG(ifirst,ilast,ngradgc))
c
c     Input/Output.
c
      REAL qhalf0(FACE2d0VECG(ifirst,ilast,nqhalfgc))
c
c     Local variables.
c
      INTEGER ic0,ic1
      REAL g
c
c     Fix predicted values to account for the inclusion of the gradient
c     of a scalar to enforce incompressibility.
c
      do ic1 = ifirst1,ilast1
         do ic0 = ifirst0-1,ilast0
            if ( gradtype.eq.0 ) then
               g = grad0(ic0+1,ic1)
            else
               g = 0.25d0*(
     &              grad1(ic1,ic0  ) + grad1(ic1+1,ic0  ) +
     &              grad1(ic1,ic0+1) + grad1(ic1+1,ic0+1) )
            endif

            qhalf0(ic0+1,ic1) = qhalf0(ic0+1,ic1)-g
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
