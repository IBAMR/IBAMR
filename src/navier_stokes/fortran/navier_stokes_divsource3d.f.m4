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

define(NDIM,3)dnl
define(REAL,`double precision')dnl
define(INTEGER,`integer')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim3d.i)dnl
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the source term F = -U max(q,0) that must be added to the
c     momentum equation to account for momentum loss due to internal
c     fluid sources.
c
c     NOTE: This is the source term that corresponds to the advective
c     (i.e., nonconservative) form of the equation.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine navier_stokes_adv_source3d(
     &     ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &     nugc0,nugc1,nugc2,
     &     nQgc0,nQgc1,nQgc2,
     &     nFgc0,nFgc1,nFgc2,
     &     u0,u1,u2,
     &     Q,
     &     F)
c
      implicit none
c
c     Input.
c
      INTEGER ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2

      INTEGER nugc0,nugc1,nugc2
      INTEGER nQgc0,nQgc1,nQgc2
      INTEGER nFgc0,nFgc1,nFgc2

      REAL u0(FACE3d0VECG(ifirst,ilast,nugc))
      REAL u1(FACE3d1VECG(ifirst,ilast,nugc))
      REAL u2(FACE3d2VECG(ifirst,ilast,nugc))

      REAL Q(CELL3dVECG(ifirst,ilast,nQgc))
c
c     Input/Output.
c
      REAL F(CELL3dVECG(ifirst,ilast,nFgc),0:NDIM-1)
c
c     Local variables.
c
      INTEGER ic0,ic1,ic2
c
c     Compute the source term F = -U max(Q,0).
c
      do ic2 = ifirst2,ilast2
         do ic1 = ifirst1,ilast1
            do ic0 = ifirst0,ilast0
               F(ic0,ic1,ic2,0) =
     &              -0.5d0*(u0(ic0+1,ic1,ic2)+u0(ic0,ic1,ic2))*
     &              dmax1(Q(ic0,ic1,ic2),0.d0)
            enddo
         enddo
      enddo

      do ic2 = ifirst2,ilast2
         do ic1 = ifirst1,ilast1
            do ic0 = ifirst0,ilast0
               F(ic0,ic1,ic2,1) =
     &              -0.5d0*(u1(ic1+1,ic2,ic0)+u1(ic1,ic2,ic0))*
     &              dmax1(Q(ic0,ic1,ic2),0.d0)
            enddo
         enddo
      enddo

      do ic2 = ifirst2,ilast2
         do ic1 = ifirst1,ilast1
            do ic0 = ifirst0,ilast0
               F(ic0,ic1,ic2,2) =
     &              -0.5d0*(u2(ic2+1,ic0,ic1)+u2(ic2,ic0,ic1))*
     &              dmax1(Q(ic0,ic1,ic2),0.d0)
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the source term F = +U min(q,0) that must be added to the
c     momentum equation to account for momentum loss due to internal
c     fluid sinks.
c
c     NOTE: This is the source term that corresponds to the conservative
c     form of the equation.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine navier_stokes_cons_source3d(
     &     ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &     nugc0,nugc1,nugc2,
     &     nQgc0,nQgc1,nQgc2,
     &     nFgc0,nFgc1,nFgc2,
     &     u0,u1,u2,
     &     Q,
     &     F)
c
      implicit none
c
c     Input.
c
      INTEGER ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2

      INTEGER nugc0,nugc1,nugc2
      INTEGER nQgc0,nQgc1,nQgc2
      INTEGER nFgc0,nFgc1,nFgc2

      REAL u0(FACE3d0VECG(ifirst,ilast,nugc))
      REAL u1(FACE3d1VECG(ifirst,ilast,nugc))
      REAL u2(FACE3d2VECG(ifirst,ilast,nugc))

      REAL Q(CELL3dVECG(ifirst,ilast,nQgc))
c
c     Input/Output.
c
      REAL F(CELL3dVECG(ifirst,ilast,nFgc),0:NDIM-1)
c
c     Local variables.
c
      INTEGER ic0,ic1,ic2
c
c     Compute the source term F = +U min(Q,0).
c
      do ic2 = ifirst2,ilast2
         do ic1 = ifirst1,ilast1
            do ic0 = ifirst0,ilast0
               F(ic0,ic1,ic2,0) =
     &              +0.5d0*(u0(ic0+1,ic1,ic2)+u0(ic0,ic1,ic2))*
     &              dmin1(Q(ic0,ic1,ic2),0.d0)
            enddo
         enddo
      enddo

      do ic2 = ifirst2,ilast2
         do ic1 = ifirst1,ilast1
            do ic0 = ifirst0,ilast0
               F(ic0,ic1,ic2,1) =
     &              +0.5d0*(u1(ic1+1,ic2,ic0)+u1(ic1,ic2,ic0))*
     &              dmin1(Q(ic0,ic1,ic2),0.d0)
            enddo
         enddo
      enddo

      do ic2 = ifirst2,ilast2
         do ic1 = ifirst1,ilast1
            do ic0 = ifirst0,ilast0
               F(ic0,ic1,ic2,2) =
     &              +0.5d0*(u2(ic2+1,ic0,ic1)+u2(ic2,ic0,ic1))*
     &              dmin1(Q(ic0,ic1,ic2),0.d0)
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the source term F = -0.5 U abs(q) that must be added to
c     the momentum equation to account for momentum loss due to internal
c     fluid sinks.
c
c     NOTE: This is the source term that corresponds to the
c     skew-symmetric form of the equation.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine navier_stokes_skew_sym_source3d(
     &     ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &     nugc0,nugc1,nugc2,
     &     nQgc0,nQgc1,nQgc2,
     &     nFgc0,nFgc1,nFgc2,
     &     u0,u1,u2,
     &     Q,
     &     F)
c
      implicit none
c
c     Input.
c
      INTEGER ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2

      INTEGER nugc0,nugc1,nugc2
      INTEGER nQgc0,nQgc1,nQgc2
      INTEGER nFgc0,nFgc1,nFgc2

      REAL u0(FACE3d0VECG(ifirst,ilast,nugc))
      REAL u1(FACE3d1VECG(ifirst,ilast,nugc))
      REAL u2(FACE3d2VECG(ifirst,ilast,nugc))

      REAL Q(CELL3dVECG(ifirst,ilast,nQgc))
c
c     Input/Output.
c
      REAL F(CELL3dVECG(ifirst,ilast,nFgc),0:NDIM-1)
c
c     Local variables.
c
      INTEGER ic0,ic1,ic2
c
c     Compute the source term F = +U min(Q,0).
c
      do ic2 = ifirst2,ilast2
         do ic1 = ifirst1,ilast1
            do ic0 = ifirst0,ilast0
               F(ic0,ic1,ic2,0) =
     &              -0.25d0*(u0(ic0+1,ic1,ic2)+u0(ic0,ic1,ic2))*
     &              abs(Q(ic0,ic1,ic2))
            enddo
         enddo
      enddo

      do ic2 = ifirst2,ilast2
         do ic1 = ifirst1,ilast1
            do ic0 = ifirst0,ilast0
               F(ic0,ic1,ic2,1) =
     &              -0.25d0*(u1(ic1+1,ic2,ic0)+u1(ic1,ic2,ic0))*
     &              abs(Q(ic0,ic1,ic2))
            enddo
         enddo
      enddo

      do ic2 = ifirst2,ilast2
         do ic1 = ifirst1,ilast1
            do ic0 = ifirst0,ilast0
               F(ic0,ic1,ic2,2) =
     &              -0.25d0*(u2(ic2+1,ic0,ic1)+u2(ic2,ic0,ic1))*
     &              abs(Q(ic0,ic1,ic2))
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the source term F = -U max(q,0) that must be added to the
c     momentum equation to account for momentum loss due to internal
c     fluid sources.
c
c     NOTE: This is the source term that corresponds to the advective
c     (i.e., nonconservative) form of the equation.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine navier_stokes_staggered_adv_source3d(
     &     ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &     nugc0,nugc1,nugc2,
     &     nQgc0,nQgc1,nQgc2,
     &     nfgc0,nfgc1,nfgc2,
     &     u0,u1,u2,
     &     Q,
     &     f0,f1,f2)
c
      implicit none
c
c     Input.
c
      INTEGER ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2

      INTEGER nugc0,nugc1,nugc2
      INTEGER nQgc0,nQgc1,nQgc2
      INTEGER nfgc0,nfgc1,nfgc2

      REAL u0(SIDE3d0VECG(ifirst,ilast,nugc))
      REAL u1(SIDE3d1VECG(ifirst,ilast,nugc))
      REAL u2(SIDE3d2VECG(ifirst,ilast,nugc))

      REAL Q(CELL3dVECG(ifirst,ilast,nQgc))
c
c     Input/Output.
c
      REAL f0(SIDE3d0VECG(ifirst,ilast,nfgc))
      REAL f1(SIDE3d1VECG(ifirst,ilast,nfgc))
      REAL f2(SIDE3d2VECG(ifirst,ilast,nfgc))
c
c     Local variables.
c
      INTEGER ic0,ic1,ic2
c
c     Compute the source term F = -U max(Q,0).
c
      do ic2 = ifirst2,ilast2
         do ic1 = ifirst1,ilast1
            do ic0 = ifirst0,ilast0+1
               f0(ic0,ic1,ic2) =
     &              -u0(ic0,ic1,ic2)*
     &              dmax1(0.5d0*(Q(ic0-1,ic1,ic2)+Q(ic0,ic1,ic2)),0.d0)
            enddo
         enddo
      enddo

      do ic2 = ifirst2,ilast2
         do ic1 = ifirst1,ilast1+1
            do ic0 = ifirst0,ilast0
               f1(ic0,ic1,ic2) =
     &              -u1(ic0,ic1,ic2)*
     &              dmax1(0.5d0*(Q(ic0,ic1-1,ic2)+Q(ic0,ic1,ic2)),0.d0)
            enddo
         enddo
      enddo

      do ic2 = ifirst2,ilast2+1
         do ic1 = ifirst1,ilast1
            do ic0 = ifirst0,ilast0
               f2(ic0,ic1,ic2) =
     &              -u2(ic0,ic1,ic2)*
     &              dmax1(0.5d0*(Q(ic0,ic1,ic2-1)+Q(ic0,ic1,ic2)),0.d0)
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the source term F = +U min(q,0) that must be added to the
c     momentum equation to account for momentum loss due to internal
c     fluid sinks.
c
c     NOTE: This is the source term that corresponds to the conservative
c     form of the equation.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine navier_stokes_staggered_cons_source3d(
     &     ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &     nugc0,nugc1,nugc2,
     &     nQgc0,nQgc1,nQgc2,
     &     nfgc0,nfgc1,nfgc2,
     &     u0,u1,u2,
     &     Q,
     &     f0,f1,f2)
c
      implicit none
c
c     Input.
c
      INTEGER ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2

      INTEGER nugc0,nugc1,nugc2
      INTEGER nQgc0,nQgc1,nQgc2
      INTEGER nfgc0,nfgc1,nfgc2

      REAL u0(SIDE3d0VECG(ifirst,ilast,nugc))
      REAL u1(SIDE3d1VECG(ifirst,ilast,nugc))
      REAL u2(SIDE3d2VECG(ifirst,ilast,nugc))

      REAL Q(CELL3dVECG(ifirst,ilast,nQgc))
c
c     Input/Output.
c
      REAL f0(SIDE3d0VECG(ifirst,ilast,nfgc))
      REAL f1(SIDE3d1VECG(ifirst,ilast,nfgc))
      REAL f2(SIDE3d2VECG(ifirst,ilast,nfgc))
c
c     Local variables.
c
      INTEGER ic0,ic1,ic2
c
c     Compute the source term F = +U min(Q,0).
c
      do ic2 = ifirst2,ilast2
         do ic1 = ifirst1,ilast1
            do ic0 = ifirst0,ilast0+1
               f0(ic0,ic1,ic2) =
     &              +u0(ic0,ic1,ic2)*
     &              dmin1(0.5d0*(Q(ic0-1,ic1,ic2)+Q(ic0,ic1,ic2)),0.d0)
            enddo
         enddo
      enddo

      do ic2 = ifirst2,ilast2
         do ic1 = ifirst1,ilast1+1
            do ic0 = ifirst0,ilast0
               f1(ic0,ic1,ic2) =
     &              +u1(ic0,ic1,ic2)*
     &              dmin1(0.5d0*(Q(ic0,ic1-1,ic2)+Q(ic0,ic1,ic2)),0.d0)
            enddo
         enddo
      enddo

      do ic2 = ifirst2,ilast2+1
         do ic1 = ifirst1,ilast1
            do ic0 = ifirst0,ilast0
               f2(ic0,ic1,ic2) =
     &              +u2(ic0,ic1,ic2)*
     &              dmin1(0.5d0*(Q(ic0,ic1,ic2-1)+Q(ic0,ic1,ic2)),0.d0)
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the source term F = -0.5 U abs(q) that must be added to
c     the momentum equation to account for momentum loss due to internal
c     fluid sinks.
c
c     NOTE: This is the source term that corresponds to the
c     skew-symmetric form of the equation.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine navier_stokes_staggered_skew_sym_source3d(
     &     ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &     nugc0,nugc1,nugc2,
     &     nQgc0,nQgc1,nQgc2,
     &     nfgc0,nfgc1,nfgc2,
     &     u0,u1,u2,
     &     Q,
     &     f0,f1,f2)
c
      implicit none
c
c     Input.
c
      INTEGER ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2

      INTEGER nugc0,nugc1,nugc2
      INTEGER nQgc0,nQgc1,nQgc2
      INTEGER nfgc0,nfgc1,nfgc2

      REAL u0(SIDE3d0VECG(ifirst,ilast,nugc))
      REAL u1(SIDE3d1VECG(ifirst,ilast,nugc))
      REAL u2(SIDE3d2VECG(ifirst,ilast,nugc))

      REAL Q(CELL3dVECG(ifirst,ilast,nQgc))
c
c     Input/Output.
c
      REAL f0(SIDE3d0VECG(ifirst,ilast,nfgc))
      REAL f1(SIDE3d1VECG(ifirst,ilast,nfgc))
      REAL f2(SIDE3d2VECG(ifirst,ilast,nfgc))
c
c     Local variables.
c
      INTEGER ic0,ic1,ic2
c
c     Compute the source term F = +U min(Q,0).
c
      do ic2 = ifirst2,ilast2
         do ic1 = ifirst1,ilast1
            do ic0 = ifirst0,ilast0+1
               f0(ic0,ic1,ic2) =
     &              -0.5d0*u0(ic0,ic1,ic2)*
     &              abs(0.5d0*(Q(ic0-1,ic1,ic2)+Q(ic0,ic1,ic2)))
            enddo
         enddo
      enddo

      do ic2 = ifirst2,ilast2
         do ic1 = ifirst1,ilast1+1
            do ic0 = ifirst0,ilast0
               f1(ic0,ic1,ic2) =
     &              -0.5d0*u1(ic0,ic1,ic2)*
     &              abs(0.5d0*(Q(ic0,ic1-1,ic2)+Q(ic0,ic1,ic2)))
            enddo
         enddo
      enddo

      do ic2 = ifirst2,ilast2+1
         do ic1 = ifirst1,ilast1
            do ic0 = ifirst0,ilast0
               f2(ic0,ic1,ic2) =
     &              -0.5d0*u2(ic0,ic1,ic2)*
     &              abs(0.5d0*(Q(ic0,ic1,ic2-1)+Q(ic0,ic1,ic2)))
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
