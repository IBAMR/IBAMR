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
define(NDIM,1)dnl
define(REAL,`double precision')dnl
define(INTEGER,`integer')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim1d.i)dnl
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes the advective flux corresponding to a face centered value
c     and a face centered advective velocity.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine advect_flux1d(
     &     dt,
     &     ifirst0,ilast0,
     &     nugc0,
     &     nqhalfgc0,
     &     nfluxgc0,
     &     u0,
     &     qhalf0,
     &     flux0)
c
      implicit none
include(TOP_SRCDIR/src/fortran/const.i)dnl
c
c     Input.
c
      INTEGER ifirst0,ilast0

      INTEGER nugc0
      INTEGER nqhalfgc0
      INTEGER nfluxgc0

      REAL dt

      REAL u0(FACE1dVECG(ifirst,ilast,nugc))
      REAL qhalf0(FACE1dVECG(ifirst,ilast,nqhalfgc))
c
c     Input/Output.
c
      REAL flux0(FACE1dVECG(ifirst,ilast,nfluxgc))
c
c     Local variables.
c
      INTEGER ic0
c
c     Compute the time integral of the advective flux.
c
      do ic0 = ifirst0-1,ilast0
         flux0(ic0+1) = dt*u0(ic0+1)*qhalf0(ic0+1)
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Update a quantity using flux differencing.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine advect_consdiff1d(
     &     dx,
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

      REAL dx(0:NDIM-1)

      REAL flux0(FACE1dVECG(ifirst,ilast,nfluxgc))
c
c     Input/Output.
c
      REAL qval(CELL1dVECG(ifirst,ilast,nqvalgc))
c
c     Local variables.
c
      INTEGER ic0
c
c     Update a quantity using flux differencing.
c
      do ic0 = ifirst0,ilast0
         qval(ic0) = qval(ic0)
     &        -(flux0(ic0+1)-flux0(ic0))/dx(0)
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
      subroutine advect_consdiffwithdivsource1d(
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
      INTEGER ic0
      REAL divsource
c
c     Update a quantity using flux differencing.
c
      do ic0 = ifirst0,ilast0
         divsource = (half/dt)*
     &        ( qflux0(ic0+1) + qflux0(ic0) )*
     &        ( uflux0(ic0+1) - uflux0(ic0) )/dx(0)

         qval(ic0) = qval(ic0) + divsource
     &        -(flux0(ic0+1)-flux0(ic0))/dx(0)
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes the advective derivative N = [u_ADV*grad(q)].
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine advect_derivative1d(
     &     dx,
     &     ifirst0,ilast0,
     &     nuadvgc0,
     &     nqgc0,
     &     uadv0,
     &     q0,
     &     nNgc0,
     &     N)
c
      implicit none
include(TOP_SRCDIR/src/fortran/const.i)dnl
c
c     Input.
c
      INTEGER ifirst0,ilast0

      INTEGER nuadvgc0
      INTEGER nqgc0
      INTEGER nNgc0

      REAL dx(0:NDIM-1)

      REAL uadv0(FACE1dVECG(ifirst,ilast,nuadvgc))

      REAL q0(FACE1dVECG(ifirst,ilast,nqgc))
c
c     Input/Output.
c
      REAL N(CELL1dVECG(ifirst,ilast,nNgc))
c
c     Local variables.
c
      INTEGER ic0
      REAL U
      REAL Qx0
c
c     Compute (U)*grad(q).
c
      do ic0 = ifirst0,ilast0
         U = 0.5d0*(uadv0(ic0+1)+uadv0(ic0))
         Qx0 = (q0(ic0+1)-q0(ic0))/dx(0)
         N(ic0) = U*Qx0
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
