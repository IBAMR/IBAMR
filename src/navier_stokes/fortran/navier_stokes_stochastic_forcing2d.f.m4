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
c     Compute staggered-grid divergence of the stochastic stress tensor.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine navier_stokes_stochastic_stress_div2d(
     &     dx,
     &     ifirst0,ilast0,
     &     ifirst1,ilast1,
     &     scale,
     &     n_W_cc_gc0,n_W_cc_gc1,
     &     W_cc,
     &     n_W_nc_gc0,n_W_nc_gc1,
     &     W_nc,
     &     n_divW_sc_gc0,n_divW_sc_gc1,
     &     divW_sc0,divW_sc1)
c
      implicit none
c
c     Input.
c
      INTEGER ifirst0,ilast0
      INTEGER ifirst1,ilast1

      INTEGER n_W_cc_gc0,n_W_cc_gc1
      INTEGER n_W_nc_gc0,n_W_nc_gc1
      INTEGER n_divW_sc_gc0,n_divW_sc_gc1

      REAL dx(0:NDIM-1)

      REAL scale

      REAL W_cc(
     &     CELL2dVECG(ifirst,ilast,n_W_cc_gc),
     &     0:NDIM-1
     &     )

      REAL W_nc(
     &     NODE2dVECG(ifirst,ilast,n_W_nc_gc),
     &     0:1
     &     )
c
c     Output.
c
      REAL divW_sc0(
     &     SIDE2d0VECG(ifirst,ilast,n_divW_sc_gc)
     &     )
      REAL divW_sc1(
     &     SIDE2d1VECG(ifirst,ilast,n_divW_sc_gc)
     &     )
c
c     Local variables.
c
      INTEGER i0,i1
c
c     Compute divW_sc0, the x component of div W.
c
      do i1=ifirst1,ilast1
         do i0=ifirst0,ilast0+1
            divW_sc0(i0,i1) = scale*(
     &           (W_cc(i0,i1  ,0)-W_cc(i0-1,i1,0))/dx(0) +
     &           (W_nc(i0,i1+1,0)-W_nc(i0  ,i1,0))/dx(1))
         enddo
      enddo
c
c     Compute divW_sc1, the y component of div W.
c
      do i1=ifirst1,ilast1+1
         do i0=ifirst0,ilast0
            divW_sc1(i0,i1) = scale*(
     &           (W_nc(i0+1,i1,1)-W_nc(i0,i1  ,1))/dx(0) +
     &           (W_cc(i0  ,i1,1)-W_cc(i0,i1-1,1))/dx(1))
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
