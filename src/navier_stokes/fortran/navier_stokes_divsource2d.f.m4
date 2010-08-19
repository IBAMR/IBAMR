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
define(NDIM,2)dnl
define(REAL,`double precision')dnl
define(INTEGER,`integer')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim2d.i)dnl
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the source term F = -U max(q,0) which must be added to the
c     momentum equation to account for momentum loss due to internal
c     sinks.
c
c     NOTE: This is the source term which corresponds to the advective
c     (i.e., nonconservative) form of the equation.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine navier_stokes_advective_divsource2d(
     &     ifirst0,ilast0,ifirst1,ilast1,
     &     nugc0,nugc1,
     &     nQgc0,nQgc1,
     &     nFgc0,nFgc1,
     &     u0,u1,
     &     Q,
     &     F)
c
      implicit none
c
c     Input.
c
      INTEGER ifirst0,ilast0,ifirst1,ilast1

      INTEGER nugc0,nugc1
      INTEGER nQgc0,nQgc1
      INTEGER nFgc0,nFgc1

      REAL u0(FACE2d0VECG(ifirst,ilast,nugc))
      REAL u1(FACE2d1VECG(ifirst,ilast,nugc))

      REAL Q(CELL2dVECG(ifirst,ilast,nQgc))
c
c     Input/Output.
c
      REAL F(CELL2dVECG(ifirst,ilast,nFgc),0:NDIM-1)
c
c     Local variables.
c
      INTEGER ic0,ic1
c
c     Compute the source term F = -U max(Q,0).
c
      do ic1 = ifirst1,ilast1
         do ic0 = ifirst0,ilast0
            F(ic0,ic1,0) =
     &           -(u0(ic0+1,ic1)+u0(ic0,ic1))*dmax1(Q(ic0,ic1),0.d0)
         enddo
      enddo

      do ic1 = ifirst1,ilast1
         do ic0 = ifirst0,ilast0
            F(ic0,ic1,1) =
     &           -(u1(ic1+1,ic0)+u1(ic1,ic0))*dmax1(Q(ic0,ic1),0.d0)
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the source term F = +U min(q,0) which must be added to the
c     momentum equation to account for momentum loss due to internal
c     sinks.
c
c     NOTE: This is the source term which corresponds to the
c     conservative form of the equation.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine navier_stokes_conservative_divsource2d(
     &     ifirst0,ilast0,ifirst1,ilast1,
     &     nugc0,nugc1,
     &     nQgc0,nQgc1,
     &     nFgc0,nFgc1,
     &     u0,u1,
     &     Q,
     &     F)
c
      implicit none
c
c     Input.
c
      INTEGER ifirst0,ilast0,ifirst1,ilast1

      INTEGER nugc0,nugc1
      INTEGER nQgc0,nQgc1
      INTEGER nFgc0,nFgc1

      REAL u0(FACE2d0VECG(ifirst,ilast,nugc))
      REAL u1(FACE2d1VECG(ifirst,ilast,nugc))

      REAL Q(CELL2dVECG(ifirst,ilast,nQgc))
c
c     Input/Output.
c
      REAL F(CELL2dVECG(ifirst,ilast,nFgc),0:NDIM-1)
c
c     Local variables.
c
      INTEGER ic0,ic1
c
c     Compute the source term F = +U min(Q,0).
c
      do ic1 = ifirst1,ilast1
         do ic0 = ifirst0,ilast0
            F(ic0,ic1,0) =
     &           +(u0(ic0+1,ic1)+u0(ic0,ic1))*dmin1(Q(ic0,ic1),0.d0)
         enddo
      enddo

      do ic1 = ifirst1,ilast1
         do ic0 = ifirst0,ilast0
            F(ic0,ic1,1) =
     &           +(u1(ic1+1,ic0)+u1(ic1,ic0))*dmin1(Q(ic0,ic1),0.d0)
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
