dnl Process this file with m4 to produce FORTRAN source code
define(NDIM,3)dnl
define(REAL,`double precision')dnl
define(INTEGER,`integer')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim3d.i)dnl
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
      subroutine godunov_predict3d(
     &     dx,dt,
     &     usefullctu,limiter,
     &     ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &     nQgc0,nQgc1,nQgc2,
     &     Q,Qscratch1,Qscratch2,
     &     nugc0,nugc1,nugc2,
     &     nqhalfgc0,nqhalfgc1,nqhalfgc2,
     &     u0,u1,u2,
     &     qtemp0,qtemp1,qtemp2,
     &     qhalf0,qhalf1,qhalf2)
c
      implicit none
include(TOP_SRCDIR/src/fortran/const.i)dnl
include(TOP_SRCDIR/src/advect/fortran/limitertypes.i)dnl
c
c     Input.
c
      INTEGER ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2

      INTEGER nQgc0,nQgc1,nQgc2

      INTEGER nugc0,nugc1,nugc2
      INTEGER nqhalfgc0,nqhalfgc1,nqhalfgc2

      LOGICAL usefullctu

      INTEGER limiter

      REAL dx(0:NDIM-1), dt

      REAL Q(CELL3dVECG(ifirst,ilast,nQgc))
      REAL Qscratch1(ifirst1-nQgc1:ilast1+nQgc1,
     &               ifirst2-nQgc2:ilast2+nQgc2,
     &               ifirst0-nQgc0:ilast0+nQgc0)
      REAL Qscratch2(ifirst2-nQgc2:ilast2+nQgc2,
     &               ifirst0-nQgc0:ilast0+nQgc0,
     &               ifirst1-nQgc1:ilast1+nQgc1)

      REAL u0(FACE3d0VECG(ifirst,ilast,nugc))
      REAL u1(FACE3d1VECG(ifirst,ilast,nugc))
      REAL u2(FACE3d2VECG(ifirst,ilast,nugc))

      REAL qtemp0(FACE3d0VECG(ifirst,ilast,nqhalfgc))
      REAL qtemp1(FACE3d1VECG(ifirst,ilast,nqhalfgc))
      REAL qtemp2(FACE3d2VECG(ifirst,ilast,nqhalfgc))
c
c     Input/Output.
c
      REAL qhalf0(FACE3d0VECG(ifirst,ilast,nqhalfgc))
      REAL qhalf1(FACE3d1VECG(ifirst,ilast,nqhalfgc))
      REAL qhalf2(FACE3d2VECG(ifirst,ilast,nqhalfgc))
c
c     Local variables.
c
      INTEGER ic0,ic1,ic2
c
c     For ease of implementation, we make copies of Q with permuted
c     indices.
c
      do ic2 = ifirst2-nQgc2,ilast2+nQgc2
         do ic1 = ifirst1-nQgc1,ilast1+nQgc1
            do ic0 = ifirst0-nQgc0,ilast0+nQgc0
               Qscratch1(ic1,ic2,ic0) = Q(ic0,ic1,ic2)
               Qscratch2(ic2,ic0,ic1) = Q(ic0,ic1,ic2)
            enddo
         enddo
      enddo
c
c     Compute temporary predicted values on cell faces.
c
c     In this computation, normal derivatives are approximated by
c     (limited) centered differences.  Transverse derivatives are not
c     included.
c
      call godunov_predict_normal3d( ! predict values on the x-faces
     &     dx(0),dt,
     &     limiter,
     &     ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &     nQgc0,nQgc1,nQgc2,
     &     Q,
     &     nugc0,nugc1,nugc2,
     &     nqhalfgc0,nqhalfgc1,nqhalfgc2,
     &     u0,
     &     qtemp0)

      call godunov_predict_normal3d( ! predict values on the y-faces
     &     dx(1),dt,
     &     limiter,
     &     ifirst1,ilast1,ifirst2,ilast2,ifirst0,ilast0,
     &     nQgc1,nQgc2,nQgc0,
     &     Qscratch1,
     &     nugc1,nugc2,nugc0,
     &     nqhalfgc1,nqhalfgc2,nqhalfgc0,
     &     u1,
     &     qtemp1)

      call godunov_predict_normal3d( ! predict values on the z-faces
     &     dx(2),dt,
     &     limiter,
     &     ifirst2,ilast2,ifirst0,ilast0,ifirst1,ilast1,
     &     nQgc2,nQgc0,nQgc1,
     &     Qscratch2,
     &     nugc2,nugc0,nugc1,
     &     nqhalfgc2,nqhalfgc0,nqhalfgc1,
     &     u2,
     &     qtemp2)
c
c     Compute final predicted values on cell faces.
c
c     This computation approximates transverse derivatives by centered
c     differences of the "temporary" predicted values.  Full corner
c     transport upwinding (increasing the largest stable timestep but
c     not order of accuracy) is somewhat expensive and optional.
c
      if ( usefullctu ) then
c
c     Include full corner transport upwinding.
c
         call godunov_transverse_ctu_fix3d( ! update values on the x-faces
     &        dx(1),dx(2),dt,
     &        ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &        nugc0,nugc1,nugc2,
     &        nqhalfgc0,nqhalfgc1,nqhalfgc2,
     &        u0,u1,u2,
     &        qtemp0,qtemp1,qtemp2,
     &        qhalf0)

         call godunov_transverse_ctu_fix3d( ! update values on the y-faces
     &        dx(2),dx(0),dt,
     &        ifirst1,ilast1,ifirst2,ilast2,ifirst0,ilast0,
     &        nugc1,nugc2,nugc0,
     &        nqhalfgc1,nqhalfgc2,nqhalfgc0,
     &        u1,u2,u0,
     &        qtemp1,qtemp2,qtemp0,
     &        qhalf1)

         call godunov_transverse_ctu_fix3d( ! update values on the y-faces
     &        dx(0),dx(1),dt,
     &        ifirst2,ilast2,ifirst0,ilast0,ifirst1,ilast1,
     &        nugc2,nugc0,nugc1,
     &        nqhalfgc2,nqhalfgc0,nqhalfgc1,
     &        u2,u0,u1,
     &        qtemp2,qtemp0,qtemp1,
     &        qhalf2)
      else
c
c     Do not include full corner transport upwinding.
c
         call godunov_transverse_fix3d( ! update values on the x-faces
     &        dx(1),dx(2),dt,
     &        ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &        nugc0,nugc1,nugc2,
     &        nqhalfgc0,nqhalfgc1,nqhalfgc2,
     &        u0,u1,u2,
     &        qtemp0,qtemp1,qtemp2,
     &        qhalf0)

         call godunov_transverse_fix3d( ! update values on the y-faces
     &        dx(2),dx(0),dt,
     &        ifirst1,ilast1,ifirst2,ilast2,ifirst0,ilast0,
     &        nugc1,nugc2,nugc0,
     &        nqhalfgc1,nqhalfgc2,nqhalfgc0,
     &        u1,u2,u0,
     &        qtemp1,qtemp2,qtemp0,
     &        qhalf1)

         call godunov_transverse_fix3d( ! update values on the y-faces
     &        dx(0),dx(1),dt,
     &        ifirst2,ilast2,ifirst0,ilast0,ifirst1,ilast1,
     &        nugc2,nugc0,nugc1,
     &        nqhalfgc2,nqhalfgc0,nqhalfgc1,
     &        u2,u0,u1,
     &        qtemp2,qtemp0,qtemp1,
     &        qhalf2)
      endif
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     A Godunov predictor used to predict face and time centered valuess
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
      subroutine godunov_predict_with_source3d(
     &     dx,dt,
     &     usefullctu,limiter,
     &     ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &     nQgc0,nQgc1,nQgc2,
     &     nFgc0,nFgc1,nFgc2,
     &     Q,Qscratch1,Qscratch2,
     &     F,Fscratch1,Fscratch2,
     &     nugc0,nugc1,nugc2,
     &     nqhalfgc0,nqhalfgc1,nqhalfgc2,
     &     u0,u1,u2,
     &     qtemp0,qtemp1,qtemp2,
     &     qhalf0,qhalf1,qhalf2)
c
      implicit none
include(TOP_SRCDIR/src/fortran/const.i)dnl
include(TOP_SRCDIR/src/advect/fortran/limitertypes.i)dnl
c
c     Input.
c
      INTEGER ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2

      INTEGER nQgc0,nQgc1,nQgc2
      INTEGER nFgc0,nFgc1,nFgc2

      INTEGER nugc0,nugc1,nugc2
      INTEGER nqhalfgc0,nqhalfgc1,nqhalfgc2

      LOGICAL usefullctu

      INTEGER limiter

      REAL dx(0:NDIM-1), dt

      REAL Q(CELL3dVECG(ifirst,ilast,nQgc))
      REAL Qscratch1(ifirst1-nQgc1:ilast1+nQgc1,
     &               ifirst2-nQgc2:ilast2+nQgc2,
     &               ifirst0-nQgc0:ilast0+nQgc0)
      REAL Qscratch2(ifirst2-nQgc2:ilast2+nQgc2,
     &               ifirst0-nQgc0:ilast0+nQgc0,
     &               ifirst1-nQgc1:ilast1+nQgc1)

      REAL F(CELL3dVECG(ifirst,ilast,nFgc))
      REAL Fscratch1(ifirst1-nFgc1:ilast1+nFgc1,
     &               ifirst2-nFgc2:ilast2+nFgc2,
     &               ifirst0-nFgc0:ilast0+nFgc0)
      REAL Fscratch2(ifirst2-nFgc2:ilast2+nFgc2,
     &               ifirst0-nFgc0:ilast0+nFgc0,
     &               ifirst1-nFgc1:ilast1+nFgc1)

      REAL u0(FACE3d0VECG(ifirst,ilast,nugc))
      REAL u1(FACE3d1VECG(ifirst,ilast,nugc))
      REAL u2(FACE3d2VECG(ifirst,ilast,nugc))

      REAL qtemp0(FACE3d0VECG(ifirst,ilast,nqhalfgc))
      REAL qtemp1(FACE3d1VECG(ifirst,ilast,nqhalfgc))
      REAL qtemp2(FACE3d2VECG(ifirst,ilast,nqhalfgc))
c
c     Input/Output.
c
      REAL qhalf0(FACE3d0VECG(ifirst,ilast,nqhalfgc))
      REAL qhalf1(FACE3d1VECG(ifirst,ilast,nqhalfgc))
      REAL qhalf2(FACE3d2VECG(ifirst,ilast,nqhalfgc))
c
c     Local variables.
c
      INTEGER ic0,ic1,ic2
c
c     For ease of implementation, we make copies of Q and F with
c     permuted indices.
c
      do ic2 = ifirst2-nQgc2,ilast2+nQgc2
         do ic1 = ifirst1-nQgc1,ilast1+nQgc1
            do ic0 = ifirst0-nQgc0,ilast0+nQgc0
               Qscratch1(ic1,ic2,ic0) = Q(ic0,ic1,ic2)
               Qscratch2(ic2,ic0,ic1) = Q(ic0,ic1,ic2)
            enddo
         enddo
      enddo

      do ic2 = ifirst2-nFgc2,ilast2+nFgc2
         do ic1 = ifirst1-nFgc1,ilast1+nFgc1
            do ic0 = ifirst0-nFgc0,ilast0+nFgc0
               Fscratch1(ic1,ic2,ic0) = F(ic0,ic1,ic2)
               Fscratch2(ic2,ic0,ic1) = F(ic0,ic1,ic2)
            enddo
         enddo
      enddo
c
c     Compute temporary predicted values on cell faces.
c
c     In this computation, normal derivatives are approximated by
c     (limited) centered differences.  Transverse derivatives are not
c     included.
c
      call godunov_predict_normal_with_source3d( ! predict values on the x-faces
     &     dx(0),dt,
     &     limiter,
     &     ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &     nQgc0,nQgc1,nQgc2,
     &     nFgc0,nFgc1,nFgc2,
     &     Q,
     &     F,
     &     nugc0,nugc1,nugc2,
     &     nqhalfgc0,nqhalfgc1,nqhalfgc2,
     &     u0,
     &     qtemp0)

      call godunov_predict_normal_with_source3d( ! predict values on the y-faces
     &     dx(1),dt,
     &     limiter,
     &     ifirst1,ilast1,ifirst2,ilast2,ifirst0,ilast0,
     &     nQgc1,nQgc2,nQgc0,
     &     nFgc1,nFgc2,nFgc0,
     &     Qscratch1,
     &     Fscratch1,
     &     nugc1,nugc2,nugc0,
     &     nqhalfgc1,nqhalfgc2,nqhalfgc0,
     &     u1,
     &     qtemp1)

      call godunov_predict_normal_with_source3d( ! predict values on the z-faces
     &     dx(2),dt,
     &     limiter,
     &     ifirst2,ilast2,ifirst0,ilast0,ifirst1,ilast1,
     &     nQgc2,nQgc0,nQgc1,
     &     nFgc2,nFgc0,nFgc1,
     &     Qscratch2,
     &     Fscratch2,
     &     nugc2,nugc0,nugc1,
     &     nqhalfgc2,nqhalfgc0,nqhalfgc1,
     &     u2,
     &     qtemp2)
c
c     Compute final predicted values on cell faces.
c
c     This computation approximates transverse derivatives by centered
c     differences of the "temporary" predicted values.  Full corner
c     transport upwinding (increasing the largest stable timestep but
c     not order of accuracy) is somewhat expensive and optional.
c
      if ( usefullctu ) then
c
c     Include full corner transport upwinding.
c
         call godunov_transverse_ctu_fix3d( ! update values on the x-faces
     &        dx(1),dx(2),dt,
     &        ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &        nugc0,nugc1,nugc2,
     &        nqhalfgc0,nqhalfgc1,nqhalfgc2,
     &        u0,u1,u2,
     &        qtemp0,qtemp1,qtemp2,
     &        qhalf0)

         call godunov_transverse_ctu_fix3d( ! update values on the y-faces
     &        dx(2),dx(0),dt,
     &        ifirst1,ilast1,ifirst2,ilast2,ifirst0,ilast0,
     &        nugc1,nugc2,nugc0,
     &        nqhalfgc1,nqhalfgc2,nqhalfgc0,
     &        u1,u2,u0,
     &        qtemp1,qtemp2,qtemp0,
     &        qhalf1)

         call godunov_transverse_ctu_fix3d( ! update values on the y-faces
     &        dx(0),dx(1),dt,
     &        ifirst2,ilast2,ifirst0,ilast0,ifirst1,ilast1,
     &        nugc2,nugc0,nugc1,
     &        nqhalfgc2,nqhalfgc0,nqhalfgc1,
     &        u2,u0,u1,
     &        qtemp2,qtemp0,qtemp1,
     &        qhalf2)
      else
c
c     Do not include full corner transport upwinding.
c
         call godunov_transverse_fix3d( ! update values on the x-faces
     &        dx(1),dx(2),dt,
     &        ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &        nugc0,nugc1,nugc2,
     &        nqhalfgc0,nqhalfgc1,nqhalfgc2,
     &        u0,u1,u2,
     &        qtemp0,qtemp1,qtemp2,
     &        qhalf0)

         call godunov_transverse_fix3d( ! update values on the y-faces
     &        dx(2),dx(0),dt,
     &        ifirst1,ilast1,ifirst2,ilast2,ifirst0,ilast0,
     &        nugc1,nugc2,nugc0,
     &        nqhalfgc1,nqhalfgc2,nqhalfgc0,
     &        u1,u2,u0,
     &        qtemp1,qtemp2,qtemp0,
     &        qhalf1)

         call godunov_transverse_fix3d( ! update values on the y-faces
     &        dx(0),dx(1),dt,
     &        ifirst2,ilast2,ifirst0,ilast0,ifirst1,ilast1,
     &        nugc2,nugc0,nugc1,
     &        nqhalfgc2,nqhalfgc0,nqhalfgc1,
     &        u2,u0,u1,
     &        qtemp2,qtemp0,qtemp1,
     &        qhalf2)
      endif
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
      subroutine godunov_incompressibility_fix3d(
     &     gradtype,
     &     ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &     ngradgc0,ngradgc1,ngradgc2,
     &     nqhalfgc0,nqhalfgc1,nqhalfgc2,
     &     grad0,grad1,grad2,
     &     qhalf0,qhalf1,qhalf2)
c
      implicit none
include(TOP_SRCDIR/src/fortran/const.i)dnl
c
c     Input.
c
      INTEGER ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2

      INTEGER ngradgc0,ngradgc1,ngradgc2
      INTEGER nqhalfgc0,nqhalfgc1,nqhalfgc2

      INTEGER gradtype

      REAL grad0(FACE3d0VECG(ifirst,ilast,ngradgc))
      REAL grad1(FACE3d1VECG(ifirst,ilast,ngradgc))
      REAL grad2(FACE3d2VECG(ifirst,ilast,ngradgc))
c
c     Input/Output.
c
      REAL qhalf0(FACE3d0VECG(ifirst,ilast,nqhalfgc))
      REAL qhalf1(FACE3d1VECG(ifirst,ilast,nqhalfgc))
      REAL qhalf2(FACE3d2VECG(ifirst,ilast,nqhalfgc))
c
c     Local variables.
c
      INTEGER gradtype0,gradtype1,gradtype2
c
      if ( gradtype.eq.0 ) then
         gradtype0 = 0
         gradtype1 = 1
         gradtype2 = 2
      elseif ( gradtype.eq.1 ) then
         gradtype1 = 0
         gradtype2 = 1
         gradtype0 = 2
      else
         gradtype2 = 0
         gradtype0 = 1
         gradtype1 = 2
      endif

      call godunov_grad_fix3d(
     &     gradtype0,
     &     ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &     ngradgc0,ngradgc1,ngradgc2,
     &     nqhalfgc0,nqhalfgc1,nqhalfgc2,
     &     grad0,grad1,grad2,
     &     qhalf0)

      call godunov_grad_fix3d(
     &     gradtype1,
     &     ifirst1,ilast1,ifirst2,ilast2,ifirst0,ilast0,
     &     ngradgc1,ngradgc2,ngradgc0,
     &     nqhalfgc1,nqhalfgc2,nqhalfgc0,
     &     grad1,grad2,grad0,
     &     qhalf1)

      call godunov_grad_fix3d(
     &     gradtype2,
     &     ifirst2,ilast2,ifirst0,ilast0,ifirst1,ilast1,
     &     ngradgc2,ngradgc0,ngradgc1,
     &     nqhalfgc2,nqhalfgc0,nqhalfgc1,
     &     grad2,grad0,grad1,
     &     qhalf2)
c
      return
      end
c
      subroutine godunov_predict_normal3d(
     &     dx0,dt,
     &     limiter,
     &     ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &     nQgc0,nQgc1,nQgc2,
     &     Q,
     &     nugc0,nugc1,nugc2,
     &     nqhalfgc0,nqhalfgc1,nqhalfgc2,
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
      INTEGER ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2

      INTEGER nQgc0,nQgc1,nQgc2

      INTEGER nugc0,nugc1,nugc2
      INTEGER nqhalfgc0,nqhalfgc1,nqhalfgc2

      INTEGER limiter

      REAL dx0,dt

      REAL Q(CELL3dVECG(ifirst,ilast,nQgc))

      REAL u0(FACE3d0VECG(ifirst,ilast,nugc))
c
c     Input/Output.
c
      REAL qhalf0(FACE3d0VECG(ifirst,ilast,nqhalfgc))
c
c     Local variables.
c
      INTEGER ic0,ic1,ic2
      REAL Qx,qL,qR
      REAL unorm
c
c     Predict face centered values using a Taylor expansion about each
c     cell center.
c
c     (Limited) centered differences are used to approximate normal
c     derivatives.  Transverse derivatives are NOT included.
c
      do ic2 = ifirst2-1,ilast2+1
         do ic1 = ifirst1-1,ilast1+1

      if     ( limiter.eq.second_order ) then
c     Employ second order slopes (no limiting).
         Qx = half*(Q(ifirst0-1+1,ic1,ic2)-Q(ifirst0-1-1,ic1,ic2))
      elseif ( limiter.eq.fourth_order ) then
         Qx = twothird*(Q(ifirst0-1+1,ic1,ic2)-Q(ifirst0-1-1,ic1,ic2))
     &      - sixth*half*(Q(ifirst0-1+2,ic1,ic2)-Q(ifirst0-1-2,ic1,ic2))
      elseif ( limiter.eq.mc_limited ) then
c     Employ van Leer's MC limiter.
         Qx = minmod3(
     &        0.5d0*(Q(ifirst0-1+1,ic1,ic2)-Q(ifirst0-1-1,ic1,ic2)),
     &        2.0d0*(Q(ifirst0-1  ,ic1,ic2)-Q(ifirst0-1-1,ic1,ic2)),
     &        2.0d0*(Q(ifirst0-1+1,ic1,ic2)-Q(ifirst0-1  ,ic1,ic2)))
      elseif ( limiter.eq.muscl_limited ) then
c     Employ Colella's MUSCL limiter.
         Qx = muscldiff(Q(ifirst0-1-2,ic1,ic2))
      else
         Qx = 0.d0
      endif

      unorm = 0.5d0*(u0(ifirst0-1,ic1,ic2)+u0(ifirst0-1+1,ic1,ic2))
!     unorm = fourth*fourth*
!    &     ( 9.d0*(u0(ifirst0-1  ,ic1,ic2)+u0(ifirst0-1+1,ic1,ic2))
!    &     - 1.d0*(u0(ifirst0-1-1,ic1,ic2)+u0(ifirst0-1+2,ic1,ic2)) )

            do ic0 = ifirst0-1,ilast0
               qL = Q(ic0  ,ic1,ic2)
     &              + 0.5d0*(1.d0-unorm*dt/dx0)*Qx

               if     ( limiter.eq.second_order ) then
                  Qx = 0.5d0*(Q(ic0+1+1,ic1,ic2)-Q(ic0+1-1,ic1,ic2))
               elseif ( limiter.eq.fourth_order ) then
                  Qx = twothird*(Q(ic0+1+1,ic1,ic2)-Q(ic0+1-1,ic1,ic2))
     &              - sixth*half*(Q(ic0+1+2,ic1,ic2)-Q(ic0+1-2,ic1,ic2))
               elseif ( limiter.eq.mc_limited ) then
             Qx = minmod3(
     &                 0.5d0*(Q(ic0+1+1,ic1,ic2)-Q(ic0+1-1,ic1,ic2)),
     &                 2.0d0*(Q(ic0+1  ,ic1,ic2)-Q(ic0+1-1,ic1,ic2)),
     &                 2.0d0*(Q(ic0+1+1,ic1,ic2)-Q(ic0+1  ,ic1,ic2)))
               elseif ( limiter.eq.muscl_limited ) then
                  Qx = muscldiff(Q(ic0+1-2,ic1,ic2))
               endif

               unorm = 0.5d0*(u0(ic0+1,ic1,ic2)+u0(ic0+2,ic1,ic2))
!              unorm = fourth*fourth*
!    &              ( 9.d0*(u0(ic0+1,ic1,ic2)+u0(ic0+2,ic1,ic2))
!    &              - 1.d0*(u0(ic0  ,ic1,ic2)+u0(ic0+3,ic1,ic2)) )

               qR = Q(ic0+1,ic1,ic2)
     &              - 0.5d0*(1.d0+unorm*dt/dx0)*Qx

               qhalf0(ic0+1,ic1,ic2) =
     &              0.5d0*(qL+qR) +
     &              sign(1.d0,u0(ic0+1,ic1,ic2))*0.5d0*(qL-qR)
            enddo

         enddo
      enddo
c
      return
      end
c
      subroutine godunov_predict_normal_with_source3d(
     &     dx0,dt,
     &     limiter,
     &     ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &     nQgc0,nQgc1,nQgc2,
     &     nFgc0,nFgc1,nFgc2,
     &     Q,
     &     F,
     &     nugc0,nugc1,nugc2,
     &     nqhalfgc0,nqhalfgc1,nqhalfgc2,
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
      INTEGER ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2

      INTEGER nQgc0,nQgc1,nQgc2
      INTEGER nFgc0,nFgc1,nFgc2

      INTEGER nugc0,nugc1,nugc2
      INTEGER nqhalfgc0,nqhalfgc1,nqhalfgc2

      INTEGER limiter

      REAL dx0,dt

      REAL Q(CELL3dVECG(ifirst,ilast,nQgc))
      REAL F(CELL3dVECG(ifirst,ilast,nFgc))

      REAL u0(FACE3d0VECG(ifirst,ilast,nugc))
c
c     Input/Output.
c
      REAL qhalf0(FACE3d0VECG(ifirst,ilast,nqhalfgc))
c
c     Local variables.
c
      INTEGER ic0,ic1,ic2
      REAL Qx,qL,qR
      REAL unorm
c
c     Predict face centered values using a Taylor expansion about each
c     cell center.
c
c     (Limited) centered differences are used to approximate normal
c     derivatives.  Transverse derivatives are NOT included.
c
      do ic2 = ifirst2-1,ilast2+1
         do ic1 = ifirst1-1,ilast1+1

      if     ( limiter.eq.second_order ) then
c     Employ second order slopes (no limiting).
         Qx = half*(Q(ifirst0-1+1,ic1,ic2)-Q(ifirst0-1-1,ic1,ic2))
      elseif ( limiter.eq.fourth_order ) then
         Qx = twothird*(Q(ifirst0-1+1,ic1,ic2)-Q(ifirst0-1-1,ic1,ic2))
     &      - sixth*half*(Q(ifirst0-1+2,ic1,ic2)-Q(ifirst0-1-2,ic1,ic2))
      elseif ( limiter.eq.mc_limited ) then
c     Employ van Leer's MC limiter.
         Qx = minmod3(
     &        0.5d0*(Q(ifirst0-1+1,ic1,ic2)-Q(ifirst0-1-1,ic1,ic2)),
     &        2.0d0*(Q(ifirst0-1  ,ic1,ic2)-Q(ifirst0-1-1,ic1,ic2)),
     &        2.0d0*(Q(ifirst0-1+1,ic1,ic2)-Q(ifirst0-1  ,ic1,ic2)))
      elseif ( limiter.eq.muscl_limited ) then
c     Employ Colella's MUSCL limiter.
         Qx = muscldiff(Q(ifirst0-1-2,ic1,ic2))
      else
         Qx = 0.d0
      endif

      unorm = 0.5d0*(u0(ifirst0-1,ic1,ic2)+u0(ifirst0-1+1,ic1,ic2))
!     unorm = fourth*fourth*
!    &        ( 9.d0*(u0(ifirst0-1  ,ic1,ic2)+u0(ifirst0-1+1,ic1,ic2))
!    &        - 1.d0*(u0(ifirst0-1-1,ic1,ic2)+u0(ifirst0-1+2,ic1,ic2)) )

            do ic0 = ifirst0-1,ilast0
               qL = Q(ic0  ,ic1,ic2)
     &              + 0.5d0*(1.d0-unorm*dt/dx0)*Qx
     &              + 0.5d0*dt*F(ic0  ,ic1,ic2)

               if     ( limiter.eq.second_order ) then
                  Qx = 0.5d0*(Q(ic0+1+1,ic1,ic2)-Q(ic0+1-1,ic1,ic2))
               elseif ( limiter.eq.fourth_order ) then
                  Qx = twothird*(Q(ic0+1+1,ic1,ic2)-Q(ic0+1-1,ic1,ic2))
     &              - sixth*half*(Q(ic0+1+2,ic1,ic2)-Q(ic0+1-2,ic1,ic2))
               elseif ( limiter.eq.mc_limited ) then
             Qx = minmod3(
     &                 0.5d0*(Q(ic0+1+1,ic1,ic2)-Q(ic0+1-1,ic1,ic2)),
     &                 2.0d0*(Q(ic0+1  ,ic1,ic2)-Q(ic0+1-1,ic1,ic2)),
     &                 2.0d0*(Q(ic0+1+1,ic1,ic2)-Q(ic0+1  ,ic1,ic2)))
               elseif ( limiter.eq.muscl_limited ) then
                  Qx = muscldiff(Q(ic0+1-2,ic1,ic2))
               endif

               unorm = 0.5d0*(u0(ic0+1,ic1,ic2)+u0(ic0+2,ic1,ic2))
!              unorm = fourth*fourth*
!    &              ( 9.d0*(u0(ic0+1,ic1,ic2)+u0(ic0+2,ic1,ic2))
!    &              - 1.d0*(u0(ic0  ,ic1,ic2)+u0(ic0+3,ic1,ic2)) )

               qR = Q(ic0+1,ic1,ic2)
     &              - 0.5d0*(1.d0+unorm*dt/dx0)*Qx
     &              + 0.5d0*dt*F(ic0+1,ic1,ic2)

               qhalf0(ic0+1,ic1,ic2) =
     &              0.5d0*(qL+qR) +
     &              sign(1.d0,u0(ic0+1,ic1,ic2))*0.5d0*(qL-qR)
            enddo

         enddo
      enddo
c
      return
      end
c
      subroutine godunov_transverse_fix3d(
     &     dx1,dx2,dt,
     &     ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &     nugc0,nugc1,nugc2,
     &     nqhalfgc0,nqhalfgc1,nqhalfgc2,
     &     u0,u1,u2,
     &     qtemp0,qtemp1,qtemp2,
     &     qhalf0)
c
      implicit none
include(TOP_SRCDIR/src/fortran/const.i)dnl
c
c     Input.
c
      INTEGER ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2

      INTEGER nugc0,nugc1,nugc2
      INTEGER nqhalfgc0,nqhalfgc1,nqhalfgc2

      REAL dx1,dx2,dt

      REAL u0(FACE3d0VECG(ifirst,ilast,nugc))
      REAL u1(FACE3d1VECG(ifirst,ilast,nugc))
      REAL u2(FACE3d2VECG(ifirst,ilast,nugc))

      REAL qtemp0(FACE3d0VECG(ifirst,ilast,nqhalfgc))
      REAL qtemp1(FACE3d1VECG(ifirst,ilast,nqhalfgc))
      REAL qtemp2(FACE3d2VECG(ifirst,ilast,nqhalfgc))
c
c     Input/Output.
c
      REAL qhalf0(FACE3d0VECG(ifirst,ilast,nqhalfgc))
c
c     Local variables.
c
      INTEGER ic0,ic1,ic2
      REAL Qy,Qz,qL_diff,qR_diff
      REAL vtan,wtan
c
c     Add transverse derivitives by taking centered differences of
c     temporary predicted values.
c
c     This computation DOES NOT include full corner transport upwinding.
c
      do ic2 = ifirst2,ilast2
         do ic1 = ifirst1,ilast1

            vtan = 0.5d0*(u1(ic1,ic2,ifirst0-1)+u1(ic1+1,ic2,ifirst0-1))
            wtan = 0.5d0*(u2(ic2,ifirst0-1,ic1)+u2(ic2+1,ifirst0-1,ic1))

!           vtan = fourth*fourth*
!    &        ( 9.d0*(u1(ic1  ,ic2,ifirst0-1)+u1(ic1+1,ic2,ifirst0-1))
!    &        - 1.d0*(u1(ic1-1,ic2,ifirst0-1)+u1(ic1+2,ic2,ifirst0-1)) )
!           wtan = fourth*fourth*
!    &        ( 9.d0*(u2(ic2  ,ifirst0-1,ic1)+u2(ic2+1,ifirst0-1,ic1))
!    &        - 1.d0*(u2(ic2-1,ifirst0-1,ic1)+u2(ic2+2,ifirst0-1,ic1)) )

            do ic0 = ifirst0-1,ilast0
               Qy = qtemp1(ic1+1,ic2,ic0)-qtemp1(ic1,ic2,ic0)
               Qz = qtemp2(ic2+1,ic0,ic1)-qtemp2(ic2,ic0,ic1)

               qL_diff =
     &              - 0.5d0*dt*vtan*Qy/dx1
     &              - 0.5d0*dt*wtan*Qz/dx2

               vtan = 0.5d0*(u1(ic1,ic2,ic0+1)+u1(ic1+1,ic2,ic0+1))
               wtan = 0.5d0*(u2(ic2,ic0+1,ic1)+u2(ic2+1,ic0+1,ic1))

!              vtan = fourth*fourth*
!    &              ( 9.d0*(u1(ic1  ,ic2,ic0+1)+u1(ic1+1,ic2,ic0+1))
!    &              - 1.d0*(u1(ic1-1,ic2,ic0+1)+u1(ic1+2,ic2,ic0+1)) )
!              wtan = fourth*fourth*
!    &              ( 9.d0*(u2(ic2  ,ic0+1,ic1)+u2(ic2+1,ic0+1,ic1))
!    &              - 1.d0*(u2(ic2-1,ic0+1,ic1)+u2(ic2+2,ic0+1,ic1)) )

               Qy = qtemp1(ic1+1,ic2,ic0+1)-qtemp1(ic1,ic2,ic0+1)
               Qz = qtemp2(ic2+1,ic0+1,ic1)-qtemp2(ic2,ic0+1,ic1)

               qR_diff =
     &              - 0.5d0*dt*vtan*Qy/dx1
     &              - 0.5d0*dt*wtan*Qz/dx2

               qhalf0(ic0+1,ic1,ic2) = qtemp0(ic0+1,ic1,ic2) +
     &              0.5d0*(qL_diff+qR_diff)+
     &              sign(1.d0,u0(ic0+1,ic1,ic2))*0.5d0*(qL_diff-qR_diff)
            enddo

         enddo
      enddo
c
      return
      end
c
      subroutine godunov_transverse_ctu_fix3d(
     &     dx1,dx2,dt,
     &     ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &     nugc0,nugc1,nugc2,
     &     nqhalfgc0,nqhalfgc1,nqhalfgc2,
     &     u0,u1,u2,
     &     qtemp0,qtemp1,qtemp2,
     &     qhalf0)
c
      implicit none
include(TOP_SRCDIR/src/fortran/const.i)dnl
c
c     Input.
c
      INTEGER ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2

      INTEGER nugc0,nugc1,nugc2
      INTEGER nqhalfgc0,nqhalfgc1,nqhalfgc2

      REAL dx1,dx2,dt

      REAL u0(FACE3d0VECG(ifirst,ilast,nugc))
      REAL u1(FACE3d1VECG(ifirst,ilast,nugc))
      REAL u2(FACE3d2VECG(ifirst,ilast,nugc))

      REAL qtemp0(FACE3d0VECG(ifirst,ilast,nqhalfgc))
      REAL qtemp1(FACE3d1VECG(ifirst,ilast,nqhalfgc))
      REAL qtemp2(FACE3d2VECG(ifirst,ilast,nqhalfgc))
c
c     Input/Output.
c
      REAL qhalf0(FACE3d0VECG(ifirst,ilast,nqhalfgc))
c
c     Local variables.
c
      INTEGER ic0,ic1,ic2
      REAL Qy,Qyupwind,Qz,Qzupwind,qL_diff,qR_diff
      REAL vtan,vtanupwind,wtan,wtanupwind
      REAL vDywQz,wDzvQy
c
c     Add transverse derivitives by taking centered differences of
c     temporary predicted values.
c
c     This computation DOES include full corner transport upwinding.
c
      do ic2 = ifirst2,ilast2
         do ic1 = ifirst1,ilast1

            vtan = 0.5d0*(u1(ic1,ic2,ifirst0-1)+u1(ic1+1,ic2,ifirst0-1))
            wtan = 0.5d0*(u2(ic2,ifirst0-1,ic1)+u2(ic2+1,ifirst0-1,ic1))

!           vtan = fourth*fourth*
!    &        ( 9.d0*(u1(ic1  ,ic2,ifirst0-1)+u1(ic1+1,ic2,ifirst0-1))
!    &        - 1.d0*(u1(ic1-1,ic2,ifirst0-1)+u1(ic1+2,ic2,ifirst0-1)) )
!           wtan = fourth*fourth*
!    &        ( 9.d0*(u2(ic2  ,ifirst0-1,ic1)+u2(ic2+1,ifirst0-1,ic1))
!    &        - 1.d0*(u2(ic2-1,ifirst0-1,ic1)+u2(ic2+2,ifirst0-1,ic1)) )

            do ic0 = ifirst0-1,ilast0
               Qy = qtemp1(ic1+1,ic2,ic0)-qtemp1(ic1,ic2,ic0)
               Qz = qtemp2(ic2+1,ic0,ic1)-qtemp2(ic2,ic0,ic1)

               if ( vtan.gt.0.d0 ) then
                  wtanupwind =
     &                 0.5d0*(u2(ic2,ic0,ic1-1)+u2(ic2+1,ic0,ic1-1))
                  Qzupwind =
     &                 qtemp2(ic2+1,ic0,ic1-1)-qtemp2(ic2,ic0,ic1-1)

                  vDywQz = vtan*(wtan*Qz-wtanupwind*Qzupwind)
               else
                  wtanupwind =
     &                 0.5d0*(u2(ic2,ic0,ic1+1)+u2(ic2+1,ic0,ic1+1))
                  Qzupwind =
     &                 qtemp2(ic2+1,ic0,ic1+1)-qtemp2(ic2,ic0,ic1+1)

                  vDywQz = vtan*(wtanupwind*Qzupwind-wtan*Qz)
               endif

               if ( wtan.gt.0.d0 ) then
                  vtanupwind =
     &                 0.5d0*(u1(ic1,ic2-1,ic0)+u1(ic1+1,ic2-1,ic0))
                  Qyupwind =
     &                 qtemp1(ic1+1,ic2-1,ic0)-qtemp1(ic1,ic2-1,ic0)

                  wDzvQy = wtan*(vtan*Qy-vtanupwind*Qyupwind)
               else
                  vtanupwind =
     &                 0.5d0*(u1(ic1,ic2+1,ic0)+u1(ic1+1,ic2+1,ic0))
                  Qyupwind =
     &                 qtemp1(ic1+1,ic2+1,ic0)-qtemp1(ic1,ic2+1,ic0)

                  wDzvQy = wtan*(vtanupwind*Qyupwind-vtan*Qy)
               endif

               qL_diff =
     &              - 0.5d0*dt*vtan*Qy/dx1
     &              - 0.5d0*dt*wtan*Qz/dx2
     &              + sixth*(dt**2.0)*(wDzvQy+vDywQz)/(dx1*dx2)

               vtan = 0.5d0*(u1(ic1,ic2,ic0+1)+u1(ic1+1,ic2,ic0+1))
               wtan = 0.5d0*(u2(ic2,ic0+1,ic1)+u2(ic2+1,ic0+1,ic1))

!              vtan = fourth*fourth*
!    &              ( 9.d0*(u1(ic1  ,ic2,ic0+1)+u1(ic1+1,ic2,ic0+1))
!    &              - 1.d0*(u1(ic1-1,ic2,ic0+1)+u1(ic1+2,ic2,ic0+1)) )
!              wtan = fourth*fourth*
!    &              ( 9.d0*(u2(ic2  ,ic0+1,ic1)+u2(ic2+1,ic0+1,ic1))
!    &              - 1.d0*(u2(ic2-1,ic0+1,ic1)+u2(ic2+2,ic0+1,ic1)) )

               Qy = qtemp1(ic1+1,ic2,ic0+1)-qtemp1(ic1,ic2,ic0+1)
               Qz = qtemp2(ic2+1,ic0+1,ic1)-qtemp2(ic2,ic0+1,ic1)

               if ( vtan.gt.0.d0 ) then
                  wtanupwind =
     &                 0.5d0*(u2(ic2,ic0+1,ic1-1)+u2(ic2+1,ic0+1,ic1-1))
                  Qzupwind =
     &                 qtemp2(ic2+1,ic0+1,ic1-1)-qtemp2(ic2,ic0+1,ic1-1)

                  vDywQz = vtan*(wtan*Qz-wtanupwind*Qzupwind)
               else
                  wtanupwind =
     &                 0.5d0*(u2(ic2,ic0+1,ic1+1)+u2(ic2+1,ic0+1,ic1+1))
                  Qzupwind =
     &                 qtemp2(ic2+1,ic0+1,ic1+1)-qtemp2(ic2,ic0+1,ic1+1)

                  vDywQz = vtan*(wtanupwind*Qzupwind-wtan*Qz)
               endif

               if ( wtan.gt.0.d0 ) then
                  vtanupwind =
     &                 0.5d0*(u1(ic1,ic2-1,ic0+1)+u1(ic1+1,ic2-1,ic0+1))
                  Qyupwind =
     &                 qtemp1(ic1+1,ic2-1,ic0+1)-qtemp1(ic1,ic2-1,ic0+1)

                  wDzvQy = wtan*(vtan*Qy-vtanupwind*Qyupwind)
               else
                  vtanupwind =
     &                 0.5d0*(u1(ic1,ic2+1,ic0+1)+u1(ic1+1,ic2+1,ic0+1))
                  Qyupwind =
     &                 qtemp1(ic1+1,ic2+1,ic0+1)-qtemp1(ic1,ic2+1,ic0+1)

                  wDzvQy = wtan*(vtanupwind*Qyupwind-vtan*Qy)
               endif

               qR_diff =
     &              - 0.5d0*dt*vtan*Qy/dx1
     &              - 0.5d0*dt*wtan*Qz/dx2
     &              + sixth*(dt**2.0)*(wDzvQy+vDywQz)/(dx1*dx2)

               qhalf0(ic0+1,ic1,ic2) = qtemp0(ic0+1,ic1,ic2) +
     &              0.5d0*(qL_diff+qR_diff)+
     &              sign(1.d0,u0(ic0+1,ic1,ic2))*0.5d0*(qL_diff-qR_diff)
            enddo

         enddo
      enddo
c
      return
      end
c
      subroutine godunov_grad_fix3d(
     &     gradtype,
     &     ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &     ngradgc0,ngradgc1,ngradgc2,
     &     nqhalfgc0,nqhalfgc1,nqhalfgc2,
     &     grad0,grad1,grad2,
     &     qhalf0)
c
      implicit none
include(TOP_SRCDIR/src/fortran/const.i)dnl
c
c     Input.
c
      INTEGER ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2

      INTEGER ngradgc0,ngradgc1,ngradgc2
      INTEGER nqhalfgc0,nqhalfgc1,nqhalfgc2

      INTEGER gradtype

      REAL grad0(FACE3d0VECG(ifirst,ilast,ngradgc))
      REAL grad1(FACE3d1VECG(ifirst,ilast,ngradgc))
      REAL grad2(FACE3d2VECG(ifirst,ilast,ngradgc))
c
c     Input/Output.
c
      REAL qhalf0(FACE3d0VECG(ifirst,ilast,nqhalfgc))
c
c     Local variables.
c
      INTEGER ic0,ic1,ic2
      REAL g
c
c     Fix predicted values to account for the inclusion of the gradient
c     of a scalar to enforce incompressibility.
c
      do ic2 = ifirst2,ilast2
         do ic1 = ifirst1,ilast1
            do ic0 = ifirst0-1,ilast0
               if ( gradtype.eq.0 ) then
                  g = grad0(ic0+1,ic1,ic2)
               elseif ( gradtype.eq.1 ) then
                  g = 0.25d0*(
     &                 grad1(ic1,ic2,ic0  ) + grad1(ic1+1,ic2,ic0  ) +
     &                 grad1(ic1,ic2,ic0+1) + grad1(ic1+1,ic2,ic0+1) )
               else
                  g = 0.25d0*(
     &                 grad2(ic2,ic0  ,ic1) + grad2(ic2+1,ic0  ,ic1) +
     &                 grad2(ic2,ic0+1,ic1) + grad2(ic2+1,ic0+1,ic1) )
               endif

               qhalf0(ic0+1,ic1,ic2) = qhalf0(ic0+1,ic1,ic2)-g
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
