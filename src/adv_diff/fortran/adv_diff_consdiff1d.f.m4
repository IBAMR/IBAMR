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
define(NDIM,1)dnl
define(REAL,`double precision')dnl
define(INTEGER,`integer')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim1d.i)dnl
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Update a quantity using flux differencing.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine adv_diff_consdiff1d(
     &     dx,dt,
     &     ifirst0,ilast0,
     &     nfluxgc0,
     &     nqvalgc0,
     &     flux0,
     &     qval)
c
      implicit none
include(TOP_SRCDIR/src/fortran/const.i)dnl
c
c     Input.
c
      INTEGER ifirst0,ilast0

      INTEGER nfluxgc0
      INTEGER nqvalgc0

      REAL dx(0:NDIM-1),dt

      REAL flux0(FACE1dVECG(ifirst,ilast,nfluxgc))
c
c     Input/Output.
c
      REAL qval(CELL1dVECG(ifirst,ilast,nqvalgc))
c
c     Local variables.
c
      INTEGER ic0,d
      REAL dtdx(0:NDIM-1)
c
c     Update a quantity using flux differencing.
c     
      do d = 0,NDIM-1
         dtdx(d) = dt*dx(d)
      enddo
      
      do ic0 = ifirst0,ilast0
         qval(ic0) = 
     &        -(flux0(ic0+1)-flux0(ic0))/dtdx(0)
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Update a quantity using flux differencing but include the proper
c     source term to account for a non-discretely divergence free
c     advection velocity.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine adv_diff_consdiffwithdivsource1d(
     &     dx,dt,
     &     ifirst0,ilast0,
     &     nfluxgc0,
     &     nqfluxgc0,
     &     nufluxgc0,
     &     nqvalgc0,
     &     flux0,
     &     qflux0,
     &     uflux0,
     &     qval)
c
      implicit none
include(TOP_SRCDIR/src/fortran/const.i)dnl
c
c     Input.
c
      INTEGER ifirst0,ilast0

      INTEGER nfluxgc0
      INTEGER nqfluxgc0
      INTEGER nufluxgc0
      INTEGER nqvalgc0

      REAL dx(0:NDIM-1),dt

      REAL flux0(FACE1dVECG(ifirst,ilast,nfluxgc))

      REAL qflux0(FACE1dVECG(ifirst,ilast,nqfluxgc))

      REAL uflux0(FACE1dVECG(ifirst,ilast,nufluxgc))
c
c     Input/Output.
c
      REAL qval(CELL1dVECG(ifirst,ilast,nqvalgc))
c
c     Local variables.
c
      INTEGER ic0,d
      REAL dtdx(0:NDIM-1),divsource
c
c     Update a quantity using flux differencing.
c
      do d = 0,NDIM-1
         dtdx(d) = dt*dx(d)
      enddo
      
      do ic0 = ifirst0,ilast0
         divsource = (half/(dt**2.d0))*
     &        ( qflux0(ic0+1) + qflux0(ic0) )*
     &        ( uflux0(ic0+1) - uflux0(ic0) )/dx(0)
         
         qval(ic0) = divsource
     &        -(flux0(ic0+1)-flux0(ic0))/dtdx(0)
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
