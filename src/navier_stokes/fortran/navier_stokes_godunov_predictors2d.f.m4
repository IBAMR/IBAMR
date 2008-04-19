dnl Process this file with m4 to produce FORTRAN source code
define(NDIM,2)dnl
define(REAL,`double precision')dnl
define(INTEGER,`integer')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim2d.i)dnl
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the sign of the input, returning zero if the absolute
c     value of x is less than a tolerance epsilon.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      REAL function sign_eps(x)
c
      implicit none
c
c     Constants.
c
      REAL EPSILON
      PARAMETER(EPSILON=1.0d-8)
c
c     Input.
c
      REAL x
c
c     Compute the sign of the input, returning zero if the absolute
c     value of x is less than a tolerance epsilon.
c
      if (dabs(x) .le. EPSILON) then
         sign_eps =  0.d0
      elseif (x  .ge. EPSILON) then
         sign_eps = +1.d0
      else
         sign_eps = -1.d0
      endif
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
c     i.e. Q satisfies an advection equation not in conservation form.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine navier_stokes_godunov_predict2d(
     &     dx,dt,
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
c
c     Input.
c
      INTEGER ifirst0,ilast0,ifirst1,ilast1

      INTEGER nQgc0,nQgc1
      INTEGER nFgc0,nFgc1

      INTEGER nugc0,nugc1
      INTEGER nqhalfgc0,nqhalfgc1

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
      call navier_stokes_godunov_predict_normal2d( ! predict values on the x-faces
     &     dx(0),dt,
     &     ifirst0,ilast0,ifirst1,ilast1,
     &     nQgc0,nQgc1,
     &     nFgc0,nFgc1,
     &     Q,
     &     F,
     &     nugc0,nugc1,
     &     nqhalfgc0,nqhalfgc1,
     &     u0,
     &     qtemp0)

      call navier_stokes_godunov_predict_normal2d( ! predict values on the y-faces
     &     dx(1),dt,
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
      call navier_stokes_godunov_transverse_fix2d( ! update values on the x-faces
     &     dx(1),dt,
     &     ifirst0,ilast0,ifirst1,ilast1,
     &     nugc0,nugc1,
     &     nqhalfgc0,nqhalfgc1,
     &     u0,u1,
     &     qtemp0,qtemp1,
     &     qhalf0)

      call navier_stokes_godunov_transverse_fix2d( ! update values on the y-faces
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
c
      subroutine navier_stokes_godunov_predict_normal2d(
     &     dx0,dt,
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
c
c     Functions.
c
      REAL sign_eps
c
c     Input.
c
      INTEGER ifirst0,ilast0,ifirst1,ilast1

      INTEGER nQgc0,nQgc1
      INTEGER nFgc0,nFgc1

      INTEGER nugc0,nugc1
      INTEGER nqhalfgc0,nqhalfgc1

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
c     Centered differences are used to approximate normal derivatives.
c     Transverse derivatives are NOT included.
c
      do ic1 = ifirst1-1,ilast1+1

!        Qx = 0.5d0*(Q(ifirst0-1+1,ic1)-Q(ifirst0-1-1,ic1))
         Qx = twothird*(Q(ifirst0-1+1,ic1)-Q(ifirst0-1-1,ic1))
     &        - sixth*half*(Q(ifirst0-1+2,ic1)-Q(ifirst0-1-2,ic1))

         unorm = 0.5d0*(u0(ifirst0-1,ic1)+u0(ifirst0-1+1,ic1))
!        unorm = fourth*fourth*
!    &        ( 9.d0*(u0(ifirst0-1  ,ic1)+u0(ifirst0-1+1,ic1))
!    &        - 1.d0*(u0(ifirst0-1-1,ic1)+u0(ifirst0-1+2,ic1)) )

         do ic0 = ifirst0-1,ilast0
            qL = Q(ic0  ,ic1)
     &           + 0.5d0*(1.d0-unorm*dt/dx0)*Qx
     &           + 0.5d0*dt*F(ic0  ,ic1)

!           Qx = 0.5d0*(Q(ic0+1+1,ic1)-Q(ic0+1-1,ic1))
            Qx = twothird*(Q(ic0+1+1,ic1)-Q(ic0+1-1,ic1))
     &           - sixth*half*(Q(ic0+1+2,ic1)-Q(ic0+1-2,ic1))

            unorm = 0.5d0*(u0(ic0+1,ic1)+u0(ic0+2,ic1))
!           unorm = fourth*fourth*
!    &           ( 9.d0*(u0(ic0+1,ic1)+u0(ic0+2,ic1))
!    &           - 1.d0*(u0(ic0  ,ic1)+u0(ic0+3,ic1)) )

            qR = Q(ic0+1,ic1)
     &           - 0.5d0*(1.d0+unorm*dt/dx0)*Qx
     &           + 0.5d0*dt*F(ic0+1,ic1)

            qhalf0(ic0+1,ic1) =
     &           0.5d0*(qL+qR)+
     &           sign_eps(u0(ic0+1,ic1))*0.5d0*(qL-qR)
         enddo

      enddo
c
      return
      end
c
      subroutine navier_stokes_godunov_transverse_fix2d(
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
c     Functions.
c
      REAL sign_eps
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

         vtan = 0.5d0*(u1(ic1,ifirst0-1)+u1(ic1+1,ifirst0-1))
!        vtan = fourth*fourth*
!    &        ( 9.d0*(u1(ic1  ,ifirst0-1)+u1(ic1+1,ifirst0-1))
!    &        - 1.d0*(u1(ic1-1,ifirst0-1)+u1(ic1+2,ifirst0-1)) )

         do ic0 = ifirst0-1,ilast0
            Qy = qtemp1(ic1+1,ic0)-qtemp1(ic1,ic0)

            qL_diff =
     &           - 0.5d0*dt*vtan*Qy/dx1

            vtan = 0.5d0*(u1(ic1,ic0+1)+u1(ic1+1,ic0+1))
!           vtan = fourth*fourth*
!    &           ( 9.d0*(u1(ic1  ,ic0+1)+u1(ic1+1,ic0+1))
!    &           - 1.d0*(u1(ic1-1,ic0+1)+u1(ic1+2,ic0+1)) )

            Qy = qtemp1(ic1+1,ic0+1)-qtemp1(ic1,ic0+1)

            qR_diff =
     &           - 0.5d0*dt*vtan*Qy/dx1

            qhalf0(ic0+1,ic1) = qtemp0(ic0+1,ic1) +
     &           0.5d0*(qL_diff+qR_diff)+
     &           sign_eps(u0(ic0+1,ic1))*0.5d0*(qL_diff-qR_diff)
         enddo

      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
