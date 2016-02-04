c
c     Routines to coarsen values via the adjoint of RT0 interpolation.
c
c     Created on 16 May 2015 by Boyce Griffith
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
c     Coarsen side-centered data via the adjoint of RT0 interpolation.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine scrt0coarsen2d(
     &     U_crse0,U_crse1,U_crse_gcw,
     &     U_fine0,U_fine1,U_fine_gcw,
     &     ilowerc0,iupperc0,
     &     ilowerc1,iupperc1,
     &     ilowerf0,iupperf0,
     &     ilowerf1,iupperf1,
     &     ratio_to_coarser,
     &     coarse_box_lower,coarse_box_upper)
c
      implicit none
c
c     Input.
c
      INTEGER ilowerc0,iupperc0
      INTEGER ilowerc1,iupperc1
      INTEGER ilowerf0,iupperf0
      INTEGER ilowerf1,iupperf1
      INTEGER U_crse_gcw,U_fine_gcw

      INTEGER ratio_to_coarser(0:NDIM-1)

      INTEGER coarse_box_lower(0:NDIM-1), coarse_box_upper(0:NDIM-1)

      REAL U_fine0(
     &     SIDE2d0(ilowerf,iupperf,U_fine_gcw)
     &     )
      REAL U_fine1(
     &     SIDE2d1(ilowerf,iupperf,U_fine_gcw)
     &     )
c
c     Input/Output.
c
      REAL U_crse0(
     &     SIDE2d0(ilowerc,iupperc,U_crse_gcw)
     &     )
      REAL U_crse1(
     &     SIDE2d1(ilowerc,iupperc,U_crse_gcw)
     &     )
c
c     Local variables.
c
      INTEGER i,i_c,i_f,j,j_c,j_f
      REAL w0,w1,w_fac
c
c     Coarsen data.
c
      w_fac = 0.d0              ! lazy lazy lazy!
      do j = 0,ratio_to_coarser(1)-1
         do i = -ratio_to_coarser(0)+1,ratio_to_coarser(0)-1
            w0 = 1.d0 - dabs(dble(i))/dble(ratio_to_coarser(0))
            w1 = 1.d0
            w_fac = w_fac + w0*w1
         enddo
      enddo
      w_fac = 1.d0/w_fac

      do j_c = coarse_box_lower(1),coarse_box_upper(1)
         do i_c = coarse_box_lower(0),coarse_box_upper(0)+1
            i_f = i_c*ratio_to_coarser(0)
            j_f = j_c*ratio_to_coarser(1)
            U_crse0(i_c,j_c) = 0.d0
            do j = 0,ratio_to_coarser(1)-1
               do i = -ratio_to_coarser(0)+1,ratio_to_coarser(0)-1
                  w0 = 1.d0 - dabs(dble(i))/dble(ratio_to_coarser(0))
                  w1 = 1.d0
                  U_crse0(i_c,j_c) = U_crse0(i_c,j_c) +
     &                 w0*w1*w_fac*U_fine0(i_f+i,j_f+j)
               enddo
            enddo
         enddo
      enddo

      w_fac = 0.d0              ! lazy lazy lazy lazy!
      do j = -ratio_to_coarser(1)+1,ratio_to_coarser(1)-1
         do i = 0,ratio_to_coarser(0)-1
            w0 = 1.d0
            w1 = 1.d0 - dabs(dble(j))/dble(ratio_to_coarser(1))
            w_fac = w_fac + w0*w1
         enddo
      enddo
      w_fac = 1.d0/w_fac

      do j_c = coarse_box_lower(1),coarse_box_upper(1)+1
         do i_c = coarse_box_lower(0),coarse_box_upper(0)
            i_f = i_c*ratio_to_coarser(0)
            j_f = j_c*ratio_to_coarser(1)
            U_crse1(i_c,j_c) = 0.d0 ! coarsen data
            do j = -ratio_to_coarser(1)+1,ratio_to_coarser(1)-1
               do i = 0,ratio_to_coarser(0)-1
                  w0 = 1.d0
                  w1 = 1.d0 - dabs(dble(j))/dble(ratio_to_coarser(1))
                  U_crse1(i_c,j_c) = U_crse1(i_c,j_c) +
     &                 w0*w1*w_fac*U_fine1(i_f+i,j_f+j)
               enddo
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
