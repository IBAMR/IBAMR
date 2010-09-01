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
     &     Q,dQ,Q_L,Q_R,Q_6,Qscratch1,
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
      REAL dQ(CELL2dVECG(ifirst,ilast,nQgc))
      REAL Q_L(CELL2dVECG(ifirst,ilast,nQgc))
      REAL Q_R(CELL2dVECG(ifirst,ilast,nQgc))
      REAL Q_6(CELL2dVECG(ifirst,ilast,nQgc))
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
c     In this computation, transverse derivatives are not included.
c
      call godunov_xsPPM7_predict_normal_with_source2d( ! predict values on the x-faces
     &     dx(0),dt,
     &     ifirst0,ilast0,ifirst1,ilast1,
     &     nQgc0,nQgc1,
     &     nFgc0,nFgc1,
     &     Q,dQ,Q_L,Q_R,Q_6,
     &     F,
     &     nugc0,nugc1,
     &     nqhalfgc0,nqhalfgc1,
     &     u0,
     &     qtemp0)

      call godunov_xsPPM7_predict_normal_with_source2d( ! predict values on the y-faces
     &     dx(1),dt,
     &     ifirst1,ilast1,ifirst0,ilast0,
     &     nQgc1,nQgc0,
     &     nFgc1,nFgc0,
     &     Qscratch1,dQ,Q_L,Q_R,Q_6,
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
      call godunov_transverse_fix2d( ! update values on the x-faces
     &     dx(1),dt,
     &     ifirst0,ilast0,ifirst1,ilast1,
     &     nugc0,nugc1,
     &     nqhalfgc0,nqhalfgc1,
     &     u0,u1,
     &     qtemp0,qtemp1,
     &     qhalf0)

      call godunov_transverse_fix2d( ! update values on the y-faces
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
