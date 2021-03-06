c ---------------------------------------------------------------------
c
c Copyright (c) 2006 - 2019 by the IBAMR developers
c All rights reserved.
c
c This file is part of IBAMR.
c
c IBAMR is free software and is distributed under the 3-clause BSD
c license. The full text of the license can be found in the file
c COPYRIGHT at the top level directory of IBAMR.
c
c ---------------------------------------------------------------------

define(NDIM,2)dnl
define(REAL,`double precision')dnl
define(INTEGER,`integer')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim2d.i)dnl
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Update a quantity using flux differencing.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine adv_diff_consdiff2d(
     &     dx,dt,
     &     ifirst0,ilast0,ifirst1,ilast1,
     &     nfluxgc0,nfluxgc1,
     &     nqvalgc0,nqvalgc1,
     &     flux0,flux1,
     &     qval)
c
      implicit none
c
c     Input.
c
      INTEGER ifirst0,ilast0,ifirst1,ilast1

      INTEGER nfluxgc0,nfluxgc1
      INTEGER nqvalgc0,nqvalgc1

      REAL dx(0:NDIM-1),dt

      REAL flux0(FACE2d0VECG(ifirst,ilast,nfluxgc))
      REAL flux1(FACE2d1VECG(ifirst,ilast,nfluxgc))
c
c     Input/Output.
c
      REAL qval(CELL2dVECG(ifirst,ilast,nqvalgc))
c
c     Local variables.
c
      INTEGER ic0,ic1,d
      REAL dtdx(0:NDIM-1)
c
c     Update a quantity using flux differencing.
c
      do d = 0,NDIM-1
         dtdx(d) = dt*dx(d)
      enddo

      do ic1 = ifirst1,ilast1
         do ic0 = ifirst0,ilast0
            qval(ic0,ic1) =
     &           -(flux0(ic0+1,ic1)-flux0(ic0,ic1))/dtdx(0)
     &           -(flux1(ic1+1,ic0)-flux1(ic1,ic0))/dtdx(1)
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
      subroutine adv_diff_consdiffwithdivsource2d(
     &     dx,dt,
     &     ifirst0,ilast0,ifirst1,ilast1,
     &     nfluxgc0,nfluxgc1,
     &     nqfluxgc0,nqfluxgc1,
     &     nufluxgc0,nufluxgc1,
     &     nqvalgc0,nqvalgc1,
     &     flux0,flux1,
     &     qflux0,qflux1,
     &     uflux0,uflux1,
     &     qval)
c
      implicit none
      REAL fourth
      parameter (fourth=0.25d0)
c
c     Input.
c
      INTEGER ifirst0,ilast0,ifirst1,ilast1

      INTEGER nfluxgc0,nfluxgc1
      INTEGER nqfluxgc0,nqfluxgc1
      INTEGER nufluxgc0,nufluxgc1
      INTEGER nqvalgc0,nqvalgc1

      REAL dx(0:NDIM-1),dt

      REAL flux0(FACE2d0VECG(ifirst,ilast,nfluxgc))
      REAL flux1(FACE2d1VECG(ifirst,ilast,nfluxgc))

      REAL qflux0(FACE2d0VECG(ifirst,ilast,nqfluxgc))
      REAL qflux1(FACE2d1VECG(ifirst,ilast,nqfluxgc))

      REAL uflux0(FACE2d0VECG(ifirst,ilast,nufluxgc))
      REAL uflux1(FACE2d1VECG(ifirst,ilast,nufluxgc))
c
c     Input/Output.
c
      REAL qval(CELL2dVECG(ifirst,ilast,nqvalgc))
c
c     Local variables.
c
      INTEGER ic0,ic1,d
      REAL dtdx(0:NDIM-1),divsource
c
c     Update a quantity using flux differencing.
c
      do d = 0,NDIM-1
         dtdx(d) = dt*dx(d)
      enddo

      do ic1 = ifirst1,ilast1
         do ic0 = ifirst0,ilast0
            divsource = (fourth/(dt**2.d0))*
     &           ( qflux0(ic0+1,ic1) + qflux0(ic0,ic1)
     &           + qflux1(ic1+1,ic0) + qflux1(ic1,ic0) )*
     &           ( (uflux0(ic0+1,ic1)-uflux0(ic0,ic1))/dx(0)
     &           + (uflux1(ic1+1,ic0)-uflux1(ic1,ic0))/dx(1) )

            qval(ic0,ic1) = divsource
     &           -(flux0(ic0+1,ic1)-flux0(ic0,ic1))/dtdx(0)
     &           -(flux1(ic1+1,ic0)-flux1(ic1,ic0))/dtdx(1)
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
