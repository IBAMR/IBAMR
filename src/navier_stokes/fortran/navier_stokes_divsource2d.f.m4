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
define(NDIM,2)dnl
define(REAL,`double precision')dnl
define(INTEGER,`integer')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim2d.i)dnl
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
      subroutine navier_stokes_adv_source2d(
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
            F(ic0,ic1,0) = -0.5d0*(u0(ic0+1,ic1)+u0(ic0,ic1))*
     &           dmax1(Q(ic0,ic1),0.d0)
         enddo
      enddo

      do ic1 = ifirst1,ilast1
         do ic0 = ifirst0,ilast0
            F(ic0,ic1,1) = -0.5d0*(u1(ic1+1,ic0)+u1(ic1,ic0))*
     &           dmax1(Q(ic0,ic1),0.d0)
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
      subroutine navier_stokes_cons_source2d(
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
            F(ic0,ic1,0) = +0.5d0*(u0(ic0+1,ic1)+u0(ic0,ic1))*
     &           dmin1(Q(ic0,ic1),0.d0)
         enddo
      enddo

      do ic1 = ifirst1,ilast1
         do ic0 = ifirst0,ilast0
            F(ic0,ic1,1) = +0.5d0*(u1(ic1+1,ic0)+u1(ic1,ic0))*
     &           dmin1(Q(ic0,ic1),0.d0)
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
      subroutine navier_stokes_skew_sym_source2d(
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
c     Compute the source term F = -0.5 U abs(Q).
c
      do ic1 = ifirst1,ilast1
         do ic0 = ifirst0,ilast0
            F(ic0,ic1,0) = -0.25d0*(u0(ic0+1,ic1)+u0(ic0,ic1))*
     &           abs(Q(ic0,ic1))
         enddo
      enddo

      do ic1 = ifirst1,ilast1
         do ic0 = ifirst0,ilast0
            F(ic0,ic1,1) = -0.25d0*(u1(ic1+1,ic0)+u1(ic1,ic0))*
     &           abs(Q(ic0,ic1))
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
      subroutine navier_stokes_staggered_adv_source2d(
     &     ifirst0,ilast0,ifirst1,ilast1,
     &     nugc0,nugc1,
     &     nQgc0,nQgc1,
     &     nfgc0,nfgc1,
     &     u0,u1,
     &     Q,
     &     f0,f1)
c
      implicit none
c
c     Input.
c
      INTEGER ifirst0,ilast0,ifirst1,ilast1

      INTEGER nugc0,nugc1
      INTEGER nQgc0,nQgc1
      INTEGER nfgc0,nfgc1

      REAL u0(SIDE2d0VECG(ifirst,ilast,nugc))
      REAL u1(SIDE2d1VECG(ifirst,ilast,nugc))

      REAL Q(CELL2dVECG(ifirst,ilast,nQgc))
c
c     Input/Output.
c
      REAL f0(SIDE2d0VECG(ifirst,ilast,nfgc))
      REAL f1(SIDE2d1VECG(ifirst,ilast,nfgc))
c
c     Local variables.
c
      INTEGER ic0,ic1
c
c     Compute the source term F = -U max(Q,0).
c
      do ic1 = ifirst1,ilast1
         do ic0 = ifirst0,ilast0+1
            f0(ic0,ic1) = -u0(ic0,ic1)*
     &           dmax1(0.5d0*(Q(ic0-1,ic1)+Q(ic0,ic1)),0.d0)
         enddo
      enddo

      do ic1 = ifirst1,ilast1+1
         do ic0 = ifirst0,ilast0
            f1(ic0,ic1) = -u1(ic0,ic1)*
     &           dmax1(0.5d0*(Q(ic0,ic1-1)+Q(ic0,ic1)),0.d0)
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
      subroutine navier_stokes_staggered_cons_source2d(
     &     ifirst0,ilast0,ifirst1,ilast1,
     &     nugc0,nugc1,
     &     nQgc0,nQgc1,
     &     nfgc0,nfgc1,
     &     u0,u1,
     &     Q,
     &     f0,f1)
c
      implicit none
c
c     Input.
c
      INTEGER ifirst0,ilast0,ifirst1,ilast1

      INTEGER nugc0,nugc1
      INTEGER nQgc0,nQgc1
      INTEGER nfgc0,nfgc1

      REAL u0(SIDE2d0VECG(ifirst,ilast,nugc))
      REAL u1(SIDE2d1VECG(ifirst,ilast,nugc))

      REAL Q(CELL2dVECG(ifirst,ilast,nQgc))
c
c     Input/Output.
c
      REAL f0(SIDE2d0VECG(ifirst,ilast,nfgc))
      REAL f1(SIDE2d1VECG(ifirst,ilast,nfgc))
c
c     Local variables.
c
      INTEGER ic0,ic1
c
c     Compute the source term F = +U min(Q,0).
c
      do ic1 = ifirst1,ilast1
         do ic0 = ifirst0,ilast0+1
            f0(ic0,ic1) = +u0(ic0,ic1)*
     &           dmin1(0.5d0*(Q(ic0-1,ic1)+Q(ic0,ic1)),0.d0)
         enddo
      enddo

      do ic1 = ifirst1,ilast1+1
         do ic0 = ifirst0,ilast0
            f1(ic0,ic1) = +u1(ic0,ic1)*
     &           dmin1(0.5d0*(Q(ic0,ic1-1)+Q(ic0,ic1)),0.d0)
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
      subroutine navier_stokes_staggered_skew_sym_source2d(
     &     ifirst0,ilast0,ifirst1,ilast1,
     &     nugc0,nugc1,
     &     nQgc0,nQgc1,
     &     nfgc0,nfgc1,
     &     u0,u1,
     &     Q,
     &     f0,f1)
c
      implicit none
c
c     Input.
c
      INTEGER ifirst0,ilast0,ifirst1,ilast1

      INTEGER nugc0,nugc1
      INTEGER nQgc0,nQgc1
      INTEGER nfgc0,nfgc1

      REAL u0(SIDE2d0VECG(ifirst,ilast,nugc))
      REAL u1(SIDE2d1VECG(ifirst,ilast,nugc))

      REAL Q(CELL2dVECG(ifirst,ilast,nQgc))
c
c     Input/Output.
c
      REAL f0(SIDE2d0VECG(ifirst,ilast,nfgc))
      REAL f1(SIDE2d1VECG(ifirst,ilast,nfgc))
c
c     Local variables.
c
      INTEGER ic0,ic1
c
c     Compute the source term F = - 0.5 U abs(Q).
c
      do ic1 = ifirst1,ilast1
         do ic0 = ifirst0,ilast0+1
            f0(ic0,ic1) = -0.5d0*u0(ic0,ic1)*
     &           abs(0.5d0*(Q(ic0-1,ic1)+Q(ic0,ic1)))
         enddo
      enddo

      do ic1 = ifirst1,ilast1+1
         do ic0 = ifirst0,ilast0
            f1(ic0,ic1) = -0.5d0*u1(ic0,ic1)*
     &           abs(0.5d0*(Q(ic0,ic1-1)+Q(ic0,ic1)))
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
