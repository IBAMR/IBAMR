c
c     Routines to refine side-centered values.
c
c     Created on 08 Feb 2011 by Boyce Griffith
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
define(coarsen_index,`dnl
         if ($1.lt.0) then
            $2=($1+1)/$3-1
         else
            $2=$1/$3
         endif
')dnl'
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Perform specialized refine operation that employs linear
c     interpolation in the normal direction and constant interpolation
c     in the tangential direction.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine cart_side_specialized_linear_refine2d(
     &     u_fine0,u_fine1,u_fine_gcw,
     &     flower0,fupper0,
     &     flower1,fupper1,
     &     u_coarse0,u_coarse1,u_coarse_gcw,
     &     clower0,cupper0,
     &     clower1,cupper1,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ratio)
c
      implicit none
c
c     Input.
c
      INTEGER u_fine_gcw
      INTEGER flower0,fupper0
      INTEGER flower1,fupper1
      INTEGER u_coarse_gcw
      INTEGER clower0,cupper0
      INTEGER clower1,cupper1
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ratio(0:NDIM-1)

      REAL u_coarse0(
     &     SIDE2d0(clower,cupper,u_coarse_gcw))
      REAL u_coarse1(
     &     SIDE2d1(clower,cupper,u_coarse_gcw))
c
c     Input/Output.
c
      REAL u_fine0(
     &     SIDE2d0(flower,fupper,u_fine_gcw))
      REAL u_fine1(
     &     SIDE2d1(flower,fupper,u_fine_gcw))
c
c     Local variables.
c
      INTEGER i0,i1
      INTEGER i_c0,i_c1
      INTEGER i_f0,i_f1
      REAL w0,w1
c
c     Refine data.
c
      do i1=ilower1,iupper1
         do i0=ilower0,iupper0+1
            coarsen_index(i0,i_c0,ratio(0))
            coarsen_index(i1,i_c1,ratio(1))
            i_f0 = i_c0*ratio(0)
            w0 = 1.d0 - dble(i0 - i_f0)/dble(ratio(0))
            u_fine0(i0,i1) = w0*u_coarse0(i_c0,i_c1)
            if (i0 .ne. i_f0) then
               w1 = 1.d0 - w0
               u_fine0(i0,i1) = u_fine0(i0,i1) +
     &              w1*u_coarse0(i_c0+1,i_c1)
            endif
         enddo
      enddo

      do i1=ilower1,iupper1+1
         do i0=ilower0,iupper0
            coarsen_index(i0,i_c0,ratio(0))
            coarsen_index(i1,i_c1,ratio(1))
            i_f1 = i_c1*ratio(1)
            w0 = 1.d0 - dble(i1 - i_f1)/dble(ratio(1))
            u_fine1(i0,i1) = w0*u_coarse1(i_c0,i_c1)
            if (i1 .ne. i_f1) then
               w1 = 1.d0 - w0
               u_fine1(i0,i1) = u_fine1(i0,i1) +
     &              w1*u_coarse1(i_c0,i_c1+1)
            endif
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
