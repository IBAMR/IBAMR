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
c     Determine the timestep size corresponding to a unit CFL number.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine navier_stokes_cc_stabledt3d(
     &     dx,
     &     ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &     ngc0,ngc1,ngc2,
     &     U0,U1,U2,
     &     stabdt)
c
      implicit none
c
c     Input.
c
      INTEGER ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2

      INTEGER ngc0,ngc1,ngc2

      REAL dx(0:NDIM-1)

      REAL U0(CELL3dVECG(ifirst,ilast,ngc))
      REAL U1(CELL3dVECG(ifirst,ilast,ngc))
      REAL U2(CELL3dVECG(ifirst,ilast,ngc))
c
c     Input/Output.
c
      REAL stabdt
c
c     Local variables.
c
      INTEGER i0,i1,i2,d
      REAL maxspeed(0:NDIM-1)
c
c     Determine the unit CFL number on the patch.
c
      do d = 0,NDIM-1
         maxspeed(d) = 1.0d-12  ! avoid division by zero
      enddo

      do i2 = ifirst2,ilast2
         do i1 = ifirst1,ilast1
            do i0 = ifirst0,ilast0
               maxspeed(0) = dmax1(maxspeed(0), dabs(U0(i0,i1,i2)))
            enddo
         enddo
      enddo

      do i2 = ifirst2,ilast2
         do i1 = ifirst1,ilast1
            do i0 = ifirst0,ilast0
               maxspeed(1) = dmax1(maxspeed(1), dabs(U1(i0,i1,i2)))
            enddo
         enddo
      enddo

      do i2 = ifirst2,ilast2
         do i1 = ifirst1,ilast1
            do i0 = ifirst0,ilast0
               maxspeed(2) = dmax1(maxspeed(2), dabs(U2(i0,i1,i2)))
            enddo
         enddo
      enddo

      stabdt = dmin1((dx(0)/maxspeed(0)),(dx(1)/maxspeed(1)),
     &               (dx(2)/maxspeed(2)))
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Determine the timestep size corresponding to a unit CFL number.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine navier_stokes_fc_stabledt3d(
     &     dx,
     &     ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &     ngc0,ngc1,ngc2,
     &     U0,U1,U2,
     &     stabdt)
c
      implicit none
c
c     Input.
c
      INTEGER ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2

      INTEGER ngc0,ngc1,ngc2

      REAL dx(0:NDIM-1)

      REAL U0(FACE3d0VECG(ifirst,ilast,ngc))
      REAL U1(FACE3d1VECG(ifirst,ilast,ngc))
      REAL U2(FACE3d2VECG(ifirst,ilast,ngc))
c
c     Input/Output.
c
      REAL stabdt
c
c     Local variables.
c
      INTEGER i0,i1,i2,d
      REAL maxspeed(0:NDIM-1)
c
c     Determine the unit CFL number on the patch.
c
      do d = 0,NDIM-1
         maxspeed(d) = 1.0d-12  ! avoid division by zero
      enddo

      do i2 = ifirst2,ilast2
         do i1 = ifirst1,ilast1
            do i0 = ifirst0,ilast0+1
               maxspeed(0) = dmax1(maxspeed(0), dabs(U0(i0,i1,i2)))
            enddo
         enddo
      enddo

      do i0 = ifirst0,ilast0
         do i2 = ifirst2,ilast2
            do i1 = ifirst1,ilast1+1
               maxspeed(1) = dmax1(maxspeed(1), dabs(U1(i1,i2,i0)))
            enddo
         enddo
      enddo

      do i1 = ifirst1,ilast1
         do i0 = ifirst0,ilast0
            do i2 = ifirst2,ilast2+1
               maxspeed(2) = dmax1(maxspeed(2), dabs(U2(i2,i0,i1)))
            enddo
         enddo
      enddo

      stabdt = dmin1((dx(0)/maxspeed(0)),(dx(1)/maxspeed(1)),
     &               (dx(2)/maxspeed(2)))
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Determine the timestep size corresponding to a unit CFL number.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine navier_stokes_sc_stabledt3d(
     &     dx,
     &     ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &     ngc0,ngc1,ngc2,
     &     U0,U1,U2,
     &     stabdt)
c
      implicit none
c
c     Input.
c
      INTEGER ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2

      INTEGER ngc0,ngc1,ngc2

      REAL dx(0:NDIM-1)

      REAL U0(SIDE3d0VECG(ifirst,ilast,ngc))
      REAL U1(SIDE3d1VECG(ifirst,ilast,ngc))
      REAL U2(SIDE3d2VECG(ifirst,ilast,ngc))
c
c     Input/Output.
c
      REAL stabdt
c
c     Local variables.
c
      INTEGER i0,i1,i2,d
      REAL maxspeed(0:NDIM-1)
c
c     Determine the unit CFL number on the patch.
c
      do d = 0,NDIM-1
         maxspeed(d) = 1.0d-12  ! avoid division by zero
      enddo

      do i2 = ifirst2,ilast2
         do i1 = ifirst1,ilast1
            do i0 = ifirst0,ilast0+1
               maxspeed(0) = dmax1(maxspeed(0), dabs(U0(i0,i1,i2)))
            enddo
         enddo
      enddo

      do i2 = ifirst2,ilast2
         do i1 = ifirst1,ilast1+1
            do i0 = ifirst0,ilast0
               maxspeed(1) = dmax1(maxspeed(1), dabs(U1(i0,i1,i2)))
            enddo
         enddo
      enddo

      do i2 = ifirst2,ilast2+1
         do i1 = ifirst1,ilast1
            do i0 = ifirst0,ilast0
               maxspeed(2) = dmax1(maxspeed(2), dabs(U2(i0,i1,i2)))
            enddo
         enddo
      enddo

      stabdt = dmin1((dx(0)/maxspeed(0)),(dx(1)/maxspeed(1)),
     &               (dx(2)/maxspeed(2)))
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
