c
c     Copyright (c) 2002-2010 Boyce Griffith
c
c     Permission is hereby granted, free of charge, to any person
c     obtaining a copy of this software and associated documentation
c     files (the "Software"), to deal in the Software without
c     restriction, including without limitation the rights to use, copy,
c     modify, merge, publish, distribute, sublicense, and/or sell copies
c     of the Software, and to permit persons to whom the Software is
c     furnished to do so, subject to the following conditions:
c
c     The above copyright notice and this permission notice shall be
c     included in all copies or substantial portions of the Software.
c
c     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
c     EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
c     MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
c     NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
c     HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
c     WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
c     OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
c     DEALINGS IN THE SOFTWARE.
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
