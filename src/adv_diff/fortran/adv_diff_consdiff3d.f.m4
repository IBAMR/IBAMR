c
c     Copyright (c) 2002-2014, Boyce Griffith
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
c        * Neither the name of The University of North Carolina nor the
c          names of its contributors may be used to endorse or promote
c          products derived from this software without specific prior
c          written permission.
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
define(NDIM,3)dnl
define(REAL,`double precision')dnl
define(INTEGER,`integer')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim3d.i)dnl
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Update a quantity using flux differencing.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine adv_diff_consdiff3d(
     &     dx,dt,
     &     ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &     nfluxgc0,nfluxgc1,nfluxgc2,
     &     nqvalgc0,nqvalgc1,nqvalgc2,
     &     flux0,flux1,flux2,
     &     qval)
c
      implicit none
c
c     Input.
c
      INTEGER ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2

      INTEGER nfluxgc0,nfluxgc1,nfluxgc2
      INTEGER nqvalgc0,nqvalgc1,nqvalgc2

      REAL dx(0:NDIM-1),dt

      REAL flux0(FACE3d0VECG(ifirst,ilast,nfluxgc))
      REAL flux1(FACE3d1VECG(ifirst,ilast,nfluxgc))
      REAL flux2(FACE3d2VECG(ifirst,ilast,nfluxgc))
c
c     Input/Output.
c
      REAL qval(CELL3dVECG(ifirst,ilast,nqvalgc))
c
c     Local variables.
c
      INTEGER ic0,ic1,ic2,d
      REAL dtdx(0:NDIM-1)
c
c     Update a quantity using flux differencing.
c
      do d = 0,NDIM-1
         dtdx(d) = dt*dx(d)
      enddo

      do ic2 = ifirst2,ilast2
         do ic1 = ifirst1,ilast1
            do ic0 = ifirst0,ilast0
               qval(ic0,ic1,ic2) =
     &              -(flux0(ic0+1,ic1,ic2)-flux0(ic0,ic1,ic2))/dtdx(0)
     &              -(flux1(ic1+1,ic2,ic0)-flux1(ic1,ic2,ic0))/dtdx(1)
     &              -(flux2(ic2+1,ic0,ic1)-flux2(ic2,ic0,ic1))/dtdx(2)
            enddo
         enddo
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
      subroutine adv_diff_consdiffwithdivsource3d(
     &     dx,dt,
     &     ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &     nfluxgc0,nfluxgc1,nfluxgc2,
     &     nqfluxgc0,nqfluxgc1,nqfluxgc2,
     &     nufluxgc0,nufluxgc1,nufluxgc2,
     &     nqvalgc0,nqvalgc1,nqvalgc2,
     &     flux0,flux1,flux2,
     &     qflux0,qflux1,qflux2,
     &     uflux0,uflux1,uflux2,
     &     qval)
c
      implicit none
      REAL sixth
      parameter (sixth=0.16666666666667d0)
c
c     Input.
c
      INTEGER ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2

      INTEGER nfluxgc0,nfluxgc1,nfluxgc2
      INTEGER nqfluxgc0,nqfluxgc1,nqfluxgc2
      INTEGER nufluxgc0,nufluxgc1,nufluxgc2
      INTEGER nqvalgc0,nqvalgc1,nqvalgc2

      REAL dx(0:NDIM-1),dt

      REAL flux0(FACE3d0VECG(ifirst,ilast,nfluxgc))
      REAL flux1(FACE3d1VECG(ifirst,ilast,nfluxgc))
      REAL flux2(FACE3d2VECG(ifirst,ilast,nfluxgc))

      REAL qflux0(FACE3d0VECG(ifirst,ilast,nqfluxgc))
      REAL qflux1(FACE3d1VECG(ifirst,ilast,nqfluxgc))
      REAL qflux2(FACE3d2VECG(ifirst,ilast,nqfluxgc))

      REAL uflux0(FACE3d0VECG(ifirst,ilast,nufluxgc))
      REAL uflux1(FACE3d1VECG(ifirst,ilast,nufluxgc))
      REAL uflux2(FACE3d2VECG(ifirst,ilast,nufluxgc))
c
c     Input/Output.
c
      REAL qval(CELL3dVECG(ifirst,ilast,nqvalgc))
c
c     Local variables.
c
      INTEGER ic0,ic1,ic2,d
      REAL dtdx(0:NDIM-1),divsource
c
c     Update a quantity using flux differencing.
c
      do d = 0,NDIM-1
         dtdx(d) = dt*dx(d)
      enddo

      do ic2 = ifirst2,ilast2
         do ic1 = ifirst1,ilast1
            do ic0 = ifirst0,ilast0
               divsource = (sixth/(dt**2.d0))*
     &              ( qflux0(ic0+1,ic1,ic2) + qflux0(ic0,ic1,ic2)
     &              + qflux1(ic1+1,ic2,ic0) + qflux1(ic1,ic2,ic0)
     &              + qflux2(ic2+1,ic0,ic1) + qflux2(ic2,ic0,ic1) )*
     &              ( (uflux0(ic0+1,ic1,ic2)-uflux0(ic0,ic1,ic2))/dx(0)
     &              + (uflux1(ic1+1,ic2,ic0)-uflux1(ic1,ic2,ic0))/dx(1)
     &              + (uflux2(ic2+1,ic0,ic1)-uflux2(ic2,ic0,ic1))/dx(2))

               qval(ic0,ic1,ic2) = divsource
     &              -(flux0(ic0+1,ic1,ic2)-flux0(ic0,ic1,ic2))/dtdx(0)
     &              -(flux1(ic1+1,ic2,ic0)-flux1(ic1,ic2,ic0))/dtdx(1)
     &              -(flux2(ic2+1,ic0,ic1)-flux2(ic2,ic0,ic1))/dtdx(2)
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
