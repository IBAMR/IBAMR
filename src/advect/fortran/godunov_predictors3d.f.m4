c
c     Copyright (c) 2002-2010, Boyce Griffith
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
c        * Neither the name of New York University nor the names of its
c          contributors may be used to endorse or promote products
c          derived from this software without specific prior written
c          permission.
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
define(NDIM,3)dnl
define(REAL,`double precision')dnl
define(INTEGER,`integer')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim3d.i)dnl
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the minmod function of two values.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      REAL function minmod(a,b)
      implicit none
      REAL a,b
      minmod = (sign(0.5d0,a)+sign(0.5d0,b))*(abs(a+b)-abs(a-b))
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the median function of three values.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      REAL function median(a,b,c)
      implicit none
      REAL a,b,c
      REAL minmod2
      median = a + minmod2(b-a,c-a)
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Make the left and right values monotone.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine monotonize(Q,Q_L,Q_R,Q_star_L,Q_star_R)
      implicit none
      REAL Q(-1:1),Q_L,Q_R,Q_L_tmp,Q_R_tmp,Q_star_L,Q_star_R
      REAL median
      Q_L_tmp = median(Q(0),Q_L,Q(-1))
      Q_R_tmp = median(Q(0),Q_R,Q(+1))
      Q_star_L = median(Q(0),Q_L_tmp,3.d0*Q(0)-2.d0*Q_R_tmp)
      Q_star_R = median(Q(0),Q_R_tmp,3.d0*Q(0)-2.d0*Q_L_tmp)
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the WENO5 interpolation of several values.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      REAL function WENO5_interp(Q)
      implicit none
      REAL Q(-2:2)
      REAL f(0:2)
      REAL IS(0:2)
      REAL omega_bar(0:2)
      REAL omega(0:2),omega_sum
      REAL alpha(0:2),alpha_sum
      INTEGER i
c     Compute the candidate interpolations.
      f(0) = (11.d0*Q( 0)-7.d0*Q(-1)+2.d0*Q(-2))/6.d0
      f(1) = ( 2.d0*Q(+1)+5.d0*Q( 0)-     Q(-1))/6.d0
      f(2) = (-1.d0*Q(+2)+5.d0*Q(+1)+2.d0*Q( 0))/6.d0
c     Compute the smoothness indicators.
      IS(0) = (13.d0/12.d0)*((Q( 0)-2.d0*Q(-1)+Q(-2))**2.d0) +
     &     0.25d0*((3.d0*Q( 0)-4.d0*Q(-1)+     Q(-2))**2.d0)
      IS(1) = (13.d0/12.d0)*((Q(+1)-2.d0*Q( 0)+Q(-1))**2.d0) +
     &     0.25d0*((     Q(+1)-                Q(-1))**2.d0)
      IS(2) = (13.d0/12.d0)*((Q(+2)-2.d0*Q(+1)+Q( 0))**2.d0) +
     &     0.25d0*((     Q(+2)-4.d0*Q(+1)+3.d0*Q( 0))**2.d0)
c     Compute the weights.
      omega_bar(0) = 0.1d0
      omega_bar(1) = 0.6d0
      omega_bar(2) = 0.3d0
      do i = 0,2
         alpha(i) = omega_bar(i)/(IS(i)+1.d-40)
      enddo
      alpha_sum = 0.d0
      do i = 0,2
         alpha_sum = alpha_sum + alpha(i)
      enddo
      do i = 0,2
         omega(i) = alpha(i)/alpha_sum
      enddo
c     Improve the accuracy of the weights (following the approach of
c     Henrick, Aslam, and Powers).
      do i = 0,2
         omega(i) = omega(i)*(omega_bar(i)+omega_bar(i)**2.d0
     &        -3.d0*omega_bar(i)*omega(i)+omega(i)**2.d0)/
     &        (omega_bar(i)**2.d0+omega(i)*(1.d0-2.d0*omega_bar(i)))
      enddo
      omega_sum = 0.d0
      do i = 0,2
         omega_sum = omega_sum + omega(i)
      enddo
      do i = 0,2
         omega(i) = omega(i)/omega_sum
      enddo
c     Compute the interpolant.
      WENO5_interp = 0.d0
      do i = 0,2
         WENO5_interp = WENO5_interp + omega(i)*f(i)
      enddo
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the sign of the input, returning zero if the absolute
c     value of x is less than a specified tolerance epsilon.
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
c          dQ/dt + u * grad Q = 0
c
c     i.e. Q satisfies an advection equation that is not in conservation
c     form.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine godunov_predict3d(
     &     dx,dt,
     &     usefullctu,
     &     ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &     nQgc0,nQgc1,nQgc2,
     &     Q,Qscratch1,Qscratch2,
     &     dQ,Q_L,Q_R,
     &     nugc0,nugc1,nugc2,
     &     nqhalfgc0,nqhalfgc1,nqhalfgc2,
     &     u0,u1,u2,
     &     qtemp0,qtemp1,qtemp2,
     &     qhalf0,qhalf1,qhalf2)
c
      implicit none
c
c     Input.
c
      INTEGER ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2

      INTEGER nQgc0,nQgc1,nQgc2

      INTEGER nugc0,nugc1,nugc2
      INTEGER nqhalfgc0,nqhalfgc1,nqhalfgc2

      LOGICAL usefullctu

      REAL dx(0:NDIM-1), dt

      REAL Q(CELL3dVECG(ifirst,ilast,nQgc))
      REAL Qscratch1(ifirst1-nQgc1:ilast1+nQgc1,
     &               ifirst2-nQgc2:ilast2+nQgc2,
     &               ifirst0-nQgc0:ilast0+nQgc0)
      REAL Qscratch2(ifirst2-nQgc2:ilast2+nQgc2,
     &               ifirst0-nQgc0:ilast0+nQgc0,
     &               ifirst1-nQgc1:ilast1+nQgc1)
      REAL dQ(CELL3dVECG(ifirst,ilast,nQgc))
      REAL Q_L(CELL3dVECG(ifirst,ilast,nQgc))
      REAL Q_R(CELL3dVECG(ifirst,ilast,nQgc))

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
      call godunov_xsPPM7_predict_normal3d( ! predict values on the x-faces
     &     dx(0),dt,
     &     ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &     nQgc0,nQgc1,nQgc2,
     &     Q,dQ,Q_L,Q_R,
     &     nugc0,nugc1,nugc2,
     &     nqhalfgc0,nqhalfgc1,nqhalfgc2,
     &     u0,
     &     qtemp0)

      call godunov_xsPPM7_predict_normal3d( ! predict values on the y-faces
     &     dx(1),dt,
     &     ifirst1,ilast1,ifirst2,ilast2,ifirst0,ilast0,
     &     nQgc1,nQgc2,nQgc0,
     &     Qscratch1,dQ,Q_L,Q_R,
     &     nugc1,nugc2,nugc0,
     &     nqhalfgc1,nqhalfgc2,nqhalfgc0,
     &     u1,
     &     qtemp1)

      call godunov_xsPPM7_predict_normal3d( ! predict values on the z-faces
     &     dx(2),dt,
     &     ifirst2,ilast2,ifirst0,ilast0,ifirst1,ilast1,
     &     nQgc2,nQgc0,nQgc1,
     &     Qscratch2,dQ,Q_L,Q_R,
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
     &     usefullctu,
     &     ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &     nQgc0,nQgc1,nQgc2,
     &     nFgc0,nFgc1,nFgc2,
     &     Q,Qscratch1,Qscratch2,
     &     dQ,Q_L,Q_R,
     &     F,Fscratch1,Fscratch2,
     &     nugc0,nugc1,nugc2,
     &     nqhalfgc0,nqhalfgc1,nqhalfgc2,
     &     u0,u1,u2,
     &     qtemp0,qtemp1,qtemp2,
     &     qhalf0,qhalf1,qhalf2)
c
      implicit none
c
c     Input.
c
      INTEGER ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2

      INTEGER nQgc0,nQgc1,nQgc2
      INTEGER nFgc0,nFgc1,nFgc2

      INTEGER nugc0,nugc1,nugc2
      INTEGER nqhalfgc0,nqhalfgc1,nqhalfgc2

      LOGICAL usefullctu

      REAL dx(0:NDIM-1), dt

      REAL Q(CELL3dVECG(ifirst,ilast,nQgc))
      REAL Qscratch1(ifirst1-nQgc1:ilast1+nQgc1,
     &               ifirst2-nQgc2:ilast2+nQgc2,
     &               ifirst0-nQgc0:ilast0+nQgc0)
      REAL Qscratch2(ifirst2-nQgc2:ilast2+nQgc2,
     &               ifirst0-nQgc0:ilast0+nQgc0,
     &               ifirst1-nQgc1:ilast1+nQgc1)
      REAL dQ(CELL3dVECG(ifirst,ilast,nQgc))
      REAL Q_L(CELL3dVECG(ifirst,ilast,nQgc))
      REAL Q_R(CELL3dVECG(ifirst,ilast,nQgc))

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
      call godunov_xsPPM7_predict_normal_with_source3d( ! predict values on the x-faces
     &     dx(0),dt,
     &     ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &     nQgc0,nQgc1,nQgc2,
     &     nFgc0,nFgc1,nFgc2,
     &     Q,dQ,Q_L,Q_R,
     &     F,
     &     nugc0,nugc1,nugc2,
     &     nqhalfgc0,nqhalfgc1,nqhalfgc2,
     &     u0,
     &     qtemp0)

      call godunov_xsPPM7_predict_normal_with_source3d( ! predict values on the y-faces
     &     dx(1),dt,
     &     ifirst1,ilast1,ifirst2,ilast2,ifirst0,ilast0,
     &     nQgc1,nQgc2,nQgc0,
     &     nFgc1,nFgc2,nFgc0,
     &     Qscratch1,dQ,Q_L,Q_R,
     &     Fscratch1,
     &     nugc1,nugc2,nugc0,
     &     nqhalfgc1,nqhalfgc2,nqhalfgc0,
     &     u1,
     &     qtemp1)

      call godunov_xsPPM7_predict_normal_with_source3d( ! predict values on the z-faces
     &     dx(2),dt,
     &     ifirst2,ilast2,ifirst0,ilast0,ifirst1,ilast1,
     &     nQgc2,nQgc0,nQgc1,
     &     nFgc2,nFgc0,nFgc1,
     &     Qscratch2,dQ,Q_L,Q_R,
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
      if     ( gradtype.eq.0 ) then
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
      subroutine godunov_xsPPM7_predict_normal3d(
     &     dx0,dt,
     &     ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &     nQgc0,nQgc1,nQgc2,
     &     Q,dQ,Q_L,Q_R,
     &     nugc0,nugc1,nugc2,
     &     nqhalfgc0,nqhalfgc1,nqhalfgc2,
     &     u0,
     &     qhalf0)
c
      implicit none
c
c     Functions.
c
      REAL median,sign_eps,WENO5_interp
c
c     Input.
c
      INTEGER ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2

      INTEGER nQgc0,nQgc1,nQgc2

      INTEGER nugc0,nugc1,nugc2
      INTEGER nqhalfgc0,nqhalfgc1,nqhalfgc2

      REAL dx0,dt

      REAL Q(CELL3dVECG(ifirst,ilast,nQgc))
      REAL dQ(CELL3dVECG(ifirst,ilast,nQgc))
      REAL Q_L(CELL3dVECG(ifirst,ilast,nQgc))
      REAL Q_R(CELL3dVECG(ifirst,ilast,nQgc))

      REAL u0(FACE3d0VECG(ifirst,ilast,nugc))
c
c     Input/Output.
c
      REAL qhalf0(FACE3d0VECG(ifirst,ilast,nqhalfgc))
c
c     Local variables.
c
      INTEGER ic0,ic1,ic2
      REAL QQ,QQ_L,QQ_R
      REAL QQ_star_L,QQ_star_R
      REAL QQ_WENO(-2:2)
      REAL QQ_WENO_L,QQ_WENO_R
      REAL QQ_4th_L,QQ_4th_R
      REAL dQQ_C,dQQ_L,dQQ_R,dQQ
      REAL unorm,nu,P0,P1,P2
      INTEGER i
c
c     Predict face centered values using the xsPPM7 scheme of Rider,
c     Greenough, and Kamm.
c
      do ic2 = ifirst2-1,ilast2+1
         do ic1 = ifirst1-1,ilast1+1
            do ic0 = ifirst0-2,ilast0+2
               dQQ_C = 0.5d0*(Q(ic0+1,ic1,ic2)-Q(ic0-1,ic1,ic2))
               dQQ_L =       (Q(ic0  ,ic1,ic2)-Q(ic0-1,ic1,ic2))
               dQQ_R =       (Q(ic0+1,ic1,ic2)-Q(ic0  ,ic1,ic2))
               if (dQQ_R*dQQ_L .gt. 1.d-12) then
                  dQQ = min(abs(dQQ_C),2.d0*abs(dQQ_L),2.d0*abs(dQQ_R))*
     c                 sign(1.d0,dQQ_C)
               else
                  dQQ = 0.d0
               endif
               dQ(ic0,ic1,ic2) = dQQ
            enddo

            do ic0 = ifirst0-1,ilast0+1
c
c     Compute a 7th order interpolation.
c
               QQ   = Q(ic0,ic1,ic2)
               QQ_L = (1.d0/420.d0)*(
     &              -   3.d0*Q(ic0+3,ic1,ic2)
     &              +  25.d0*Q(ic0+2,ic1,ic2)
     &              - 101.d0*Q(ic0+1,ic1,ic2)
     &              + 319.d0*Q(ic0  ,ic1,ic2)
     &              + 214.d0*Q(ic0-1,ic1,ic2)
     &              -  38.d0*Q(ic0-2,ic1,ic2)
     &              +   4.d0*Q(ic0-3,ic1,ic2))
               QQ_R = (1.d0/420.d0)*(
     &              -   3.d0*Q(ic0-3,ic1,ic2)
     &              +  25.d0*Q(ic0-2,ic1,ic2)
     &              - 101.d0*Q(ic0-1,ic1,ic2)
     &              + 319.d0*Q(ic0  ,ic1,ic2)
     &              + 214.d0*Q(ic0+1,ic1,ic2)
     &              -  38.d0*Q(ic0+2,ic1,ic2)
     &              +   4.d0*Q(ic0+3,ic1,ic2))
               Q_L(ic0,ic1,ic2) = QQ_L
               Q_R(ic0,ic1,ic2) = QQ_R
c
c     Check for extrema or violations of monotonicity.
c
               call monotonize(
     &              Q(ic0-1,ic1,ic2),
     &              QQ_L,QQ_R,QQ_star_L,QQ_star_R)
               if ( ((QQ_star_L-QQ_L)**2.d0 .ge. 1.d-12) .or.
     &              ((QQ_star_R-QQ_R)**2.d0 .ge. 1.d-12) ) then
                  do i = -2,2
                     QQ_WENO(i) = Q(ic0-i,ic1,ic2)
                  enddo
                  QQ_WENO_L = WENO5_interp(QQ_WENO)
                  do i = -2,2
                     QQ_WENO(i) = Q(ic0+i,ic1,ic2)
                  enddo
                  QQ_WENO_R = WENO5_interp(QQ_WENO)
                  if ( ((QQ_star_L-QQ)**2.d0 .le. 1.d-12) .or.
     &                 ((QQ_star_R-QQ)**2.d0 .le. 1.d-12) ) then
                     QQ_WENO_L = median(QQ,QQ_WENO_L,QQ_L)
                     QQ_WENO_R = median(QQ,QQ_WENO_R,QQ_R)
                     call monotonize(
     &                    Q(ic0-1,ic1,ic2),
     &                    QQ_WENO_L,QQ_WENO_R,QQ_star_L,QQ_star_R)
                  else
                     QQ_4th_L = 0.5d0*(
     &                    Q(ic0-1,ic1,ic2)+Q(ic0,ic1,ic2)) -
     &                    (1.d0/6.d0)*(
     &                    dQ(ic0,ic1,ic2)-dQ(ic0-1,ic1,ic2))
                     QQ_4th_R = 0.5d0*(
     &                    Q(ic0,ic1,ic2)+Q(ic0+1,ic1,ic2)) -
     &                    (1.d0/6.d0)*(
     &                    dQ(ic0+1,ic1,ic2)-dQ(ic0,ic1,ic2))
                     QQ_4th_L = median(QQ_4th_L,QQ_WENO_L,QQ_L)
                     QQ_4th_R = median(QQ_4th_R,QQ_WENO_R,QQ_R)
                     call monotonize(
     &                    Q(ic0-1,ic1,ic2),
     &                    QQ_4th_L,QQ_4th_R,QQ_star_L,QQ_star_R)
                  endif
                  Q_L(ic0,ic1,ic2) = median(QQ_WENO_L,QQ_star_L,QQ_L)
                  Q_R(ic0,ic1,ic2) = median(QQ_WENO_R,QQ_star_R,QQ_R)
               endif
            enddo

            do ic0 = ifirst0-1,ilast0
               unorm = 0.5d0*(u0(ic0  ,ic1,ic2)+u0(ic0+1,ic1,ic2))
               nu    = unorm*dt/dx0
               QQ        = Q  (ic0  ,ic1,ic2)
               QQ_star_L = Q_L(ic0  ,ic1,ic2)
               QQ_star_R = Q_R(ic0  ,ic1,ic2)
               P0 = 1.5d0*QQ-0.25d0*(QQ_star_L+QQ_star_R)
               P1 = QQ_star_R-QQ_star_L
               P2 = 3.d0*(QQ_star_L+QQ_star_R)-6.d0*QQ
               QQ_L = P0 + 0.5d0*P1 + 0.25d0*P2
     &              - 0.5d0*nu*P1
     &              + (-(nu/2.d0)+(nu*nu/3.d0))*P2

               unorm = 0.5d0*(u0(ic0+1,ic1,ic2)+u0(ic0+2,ic1,ic2))
               nu    = unorm*dt/dx0
               QQ        = Q  (ic0+1,ic1,ic2)
               QQ_star_L = Q_L(ic0+1,ic1,ic2)
               QQ_star_R = Q_R(ic0+1,ic1,ic2)
               P0 = 1.5d0*QQ-0.25d0*(QQ_star_L+QQ_star_R)
               P1 = QQ_star_R-QQ_star_L
               P2 = 3.d0*(QQ_star_L+QQ_star_R)-6.d0*QQ
               QQ_R = P0 - 0.5d0*P1 + 0.25d0*P2
     &              - 0.5d0*nu*P1
     &              + (+(nu/2.d0)+(nu*nu/3.d0))*P2

               qhalf0(ic0+1,ic1,ic2) =
     &              0.5d0*(QQ_L+QQ_R)+
     &              sign_eps(u0(ic0+1,ic1,ic2))*0.5d0*(QQ_L-QQ_R)
            enddo
         enddo
      enddo
c
      return
      end
c
      subroutine godunov_xsPPM7_predict_normal_with_source3d(
     &     dx0,dt,
     &     ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &     nQgc0,nQgc1,nQgc2,
     &     nFgc0,nFgc1,nFgc2,
     &     Q,dQ,Q_L,Q_R,
     &     F,
     &     nugc0,nugc1,nugc2,
     &     nqhalfgc0,nqhalfgc1,nqhalfgc2,
     &     u0,
     &     qhalf0)
c
      implicit none
c
c     Functions.
c
      REAL median,sign_eps,WENO5_interp
c
c     Input.
c
      INTEGER ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2

      INTEGER nQgc0,nQgc1,nQgc2
      INTEGER nFgc0,nFgc1,nFgc2

      INTEGER nugc0,nugc1,nugc2
      INTEGER nqhalfgc0,nqhalfgc1,nqhalfgc2

      REAL dx0,dt

      REAL Q(CELL3dVECG(ifirst,ilast,nQgc))
      REAL dQ(CELL3dVECG(ifirst,ilast,nQgc))
      REAL Q_L(CELL3dVECG(ifirst,ilast,nQgc))
      REAL Q_R(CELL3dVECG(ifirst,ilast,nQgc))
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
      REAL QQ,QQ_L,QQ_R
      REAL QQ_star_L,QQ_star_R
      REAL QQ_WENO(-2:2)
      REAL QQ_WENO_L,QQ_WENO_R
      REAL QQ_4th_L,QQ_4th_R
      REAL dQQ_C,dQQ_L,dQQ_R,dQQ
      REAL unorm,nu,P0,P1,P2
      INTEGER i
c
c     Predict face centered values using the xsPPM7 scheme of Rider,
c     Greenough, and Kamm.
c
      do ic2 = ifirst2-1,ilast2+1
         do ic1 = ifirst1-1,ilast1+1
            do ic0 = ifirst0-2,ilast0+2
               dQQ_C = 0.5d0*(Q(ic0+1,ic1,ic2)-Q(ic0-1,ic1,ic2))
               dQQ_L =       (Q(ic0  ,ic1,ic2)-Q(ic0-1,ic1,ic2))
               dQQ_R =       (Q(ic0+1,ic1,ic2)-Q(ic0  ,ic1,ic2))
               if (dQQ_R*dQQ_L .gt. 1.d-12) then
                  dQQ = min(abs(dQQ_C),2.d0*abs(dQQ_L),2.d0*abs(dQQ_R))*
     c                 sign(1.d0,dQQ_C)
               else
                  dQQ = 0.d0
               endif
               dQ(ic0,ic1,ic2) = dQQ
            enddo

            do ic0 = ifirst0-1,ilast0+1
c
c     Compute a 7th order interpolation.
c
               QQ   = Q(ic0,ic1,ic2)
               QQ_L = (1.d0/420.d0)*(
     &              -   3.d0*Q(ic0+3,ic1,ic2)
     &              +  25.d0*Q(ic0+2,ic1,ic2)
     &              - 101.d0*Q(ic0+1,ic1,ic2)
     &              + 319.d0*Q(ic0  ,ic1,ic2)
     &              + 214.d0*Q(ic0-1,ic1,ic2)
     &              -  38.d0*Q(ic0-2,ic1,ic2)
     &              +   4.d0*Q(ic0-3,ic1,ic2))
               QQ_R = (1.d0/420.d0)*(
     &              -   3.d0*Q(ic0-3,ic1,ic2)
     &              +  25.d0*Q(ic0-2,ic1,ic2)
     &              - 101.d0*Q(ic0-1,ic1,ic2)
     &              + 319.d0*Q(ic0  ,ic1,ic2)
     &              + 214.d0*Q(ic0+1,ic1,ic2)
     &              -  38.d0*Q(ic0+2,ic1,ic2)
     &              +   4.d0*Q(ic0+3,ic1,ic2))
               Q_L(ic0,ic1,ic2) = QQ_L
               Q_R(ic0,ic1,ic2) = QQ_R
c
c     Check for extrema or violations of monotonicity.
c
               call monotonize(
     &              Q(ic0-1,ic1,ic2),
     &              QQ_L,QQ_R,QQ_star_L,QQ_star_R)
               if ( ((QQ_star_L-QQ_L)**2.d0 .ge. 1.d-12) .or.
     &              ((QQ_star_R-QQ_R)**2.d0 .ge. 1.d-12) ) then
                  do i = -2,2
                     QQ_WENO(i) = Q(ic0-i,ic1,ic2)
                  enddo
                  QQ_WENO_L = WENO5_interp(QQ_WENO)
                  do i = -2,2
                     QQ_WENO(i) = Q(ic0+i,ic1,ic2)
                  enddo
                  QQ_WENO_R = WENO5_interp(QQ_WENO)
                  if ( ((QQ_star_L-QQ)**2.d0 .le. 1.d-12) .or.
     &                 ((QQ_star_R-QQ)**2.d0 .le. 1.d-12) ) then
                     QQ_WENO_L = median(QQ,QQ_WENO_L,QQ_L)
                     QQ_WENO_R = median(QQ,QQ_WENO_R,QQ_R)
                     call monotonize(
     &                    Q(ic0-1,ic1,ic2),
     &                    QQ_WENO_L,QQ_WENO_R,QQ_star_L,QQ_star_R)
                  else
                     QQ_4th_L = 0.5d0*(
     &                    Q(ic0-1,ic1,ic2)+Q(ic0,ic1,ic2)) -
     &                    (1.d0/6.d0)*(
     &                    dQ(ic0,ic1,ic2)-dQ(ic0-1,ic1,ic2))
                     QQ_4th_R = 0.5d0*(
     &                    Q(ic0,ic1,ic2)+Q(ic0+1,ic1,ic2)) -
     &                    (1.d0/6.d0)*(
     &                    dQ(ic0+1,ic1,ic2)-dQ(ic0,ic1,ic2))
                     QQ_4th_L = median(QQ_4th_L,QQ_WENO_L,QQ_L)
                     QQ_4th_R = median(QQ_4th_R,QQ_WENO_R,QQ_R)
                     call monotonize(
     &                    Q(ic0-1,ic1,ic2),
     &                    QQ_4th_L,QQ_4th_R,QQ_star_L,QQ_star_R)
                  endif
                  Q_L(ic0,ic1,ic2) = median(QQ_WENO_L,QQ_star_L,QQ_L)
                  Q_R(ic0,ic1,ic2) = median(QQ_WENO_R,QQ_star_R,QQ_R)
               endif
            enddo

            do ic0 = ifirst0-1,ilast0
               unorm = 0.5d0*(u0(ic0  ,ic1,ic2)+u0(ic0+1,ic1,ic2))
               nu    = unorm*dt/dx0
               QQ        = Q  (ic0  ,ic1,ic2)
               QQ_star_L = Q_L(ic0  ,ic1,ic2)
               QQ_star_R = Q_R(ic0  ,ic1,ic2)
               P0 = 1.5d0*QQ-0.25d0*(QQ_star_L+QQ_star_R)
               P1 = QQ_star_R-QQ_star_L
               P2 = 3.d0*(QQ_star_L+QQ_star_R)-6.d0*QQ
               QQ_L = P0 + 0.5d0*P1 + 0.25d0*P2
     &              - 0.5d0*nu*P1
     &              + (-(nu/2.d0)+(nu*nu/3.d0))*P2
     &              + 0.5d0*dt*F(ic0  ,ic1,ic2)

               unorm = 0.5d0*(u0(ic0+1,ic1,ic2)+u0(ic0+2,ic1,ic2))
               nu    = unorm*dt/dx0
               QQ        = Q  (ic0+1,ic1,ic2)
               QQ_star_L = Q_L(ic0+1,ic1,ic2)
               QQ_star_R = Q_R(ic0+1,ic1,ic2)
               P0 = 1.5d0*QQ-0.25d0*(QQ_star_L+QQ_star_R)
               P1 = QQ_star_R-QQ_star_L
               P2 = 3.d0*(QQ_star_L+QQ_star_R)-6.d0*QQ
               QQ_R = P0 - 0.5d0*P1 + 0.25d0*P2
     &              - 0.5d0*nu*P1
     &              + (+(nu/2.d0)+(nu*nu/3.d0))*P2
     &              + 0.5d0*dt*F(ic0+1,ic1,ic2)

               qhalf0(ic0+1,ic1,ic2) =
     &              0.5d0*(QQ_L+QQ_R)+
     &              sign_eps(u0(ic0+1,ic1,ic2))*0.5d0*(QQ_L-QQ_R)
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
c
c     Functions.
c
      REAL sign_eps
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
            do ic0 = ifirst0-1,ilast0
               Qy = qtemp1(ic1+1,ic2,ic0)-qtemp1(ic1,ic2,ic0)
               Qz = qtemp2(ic2+1,ic0,ic1)-qtemp2(ic2,ic0,ic1)
               qL_diff =
     &              - 0.5d0*dt*vtan*Qy/dx1
     &              - 0.5d0*dt*wtan*Qz/dx2
               vtan = 0.5d0*(u1(ic1,ic2,ic0+1)+u1(ic1+1,ic2,ic0+1))
               wtan = 0.5d0*(u2(ic2,ic0+1,ic1)+u2(ic2+1,ic0+1,ic1))
               Qy = qtemp1(ic1+1,ic2,ic0+1)-qtemp1(ic1,ic2,ic0+1)
               Qz = qtemp2(ic2+1,ic0+1,ic1)-qtemp2(ic2,ic0+1,ic1)
               qR_diff =
     &              - 0.5d0*dt*vtan*Qy/dx1
     &              - 0.5d0*dt*wtan*Qz/dx2
               qhalf0(ic0+1,ic1,ic2) = qtemp0(ic0+1,ic1,ic2) +
     &              0.5d0*(qL_diff+qR_diff)+
     &              sign_eps(u0(ic0+1,ic1,ic2))*0.5d0*(qL_diff-qR_diff)
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
      REAL sixth
      parameter (sixth=0.16666666666667d0)
c
c     Functions.
c
      REAL sign_eps
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
     &              + sixth*(dt**2.d0)*(wDzvQy+vDywQz)/(dx1*dx2)

               vtan = 0.5d0*(u1(ic1,ic2,ic0+1)+u1(ic1+1,ic2,ic0+1))
               wtan = 0.5d0*(u2(ic2,ic0+1,ic1)+u2(ic2+1,ic0+1,ic1))
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
     &              + sixth*(dt**2.d0)*(wDzvQy+vDywQz)/(dx1*dx2)

               qhalf0(ic0+1,ic1,ic2) = qtemp0(ic0+1,ic1,ic2) +
     &              0.5d0*(qL_diff+qR_diff)+
     &              sign_eps(u0(ic0+1,ic1,ic2))*0.5d0*(qL_diff-qR_diff)
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
