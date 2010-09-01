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
