c
c     Routines to set CF boundary values via quadratic interpolation.
c
c     Created on 12 Dec 2007 by Boyce Griffith
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
c     Compute a cell-centered quadratic interpolation in the tangential
c     direction at coarse-fine interfaces.
c
c     NOTE: This routine imposes a maximum refinement ratio of 16.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ccquadtangentialinterpolation2d(
     &     U_fine,U_fine_gcw,
     &     U_crse,U_crse_gcw,
     &     ilowerf0,iupperf0,
     &     ilowerf1,iupperf1,
     &     ilowerc0,iupperc0,
     &     ilowerc1,iupperc1,
     &     location_index,ratio_to_coarser,
     &     blowerf,bupperf)
c
      implicit none
c
c     Constants.
c
      INTEGER MAX_RATIO
      parameter (MAX_RATIO=16)  ! max refinement ratio = 16
      REAL fourth,eighth
      parameter (fourth=0.25d0)
      parameter (eighth=0.125d0)
c
c     Input.
c
      INTEGER ilowerf0,iupperf0
      INTEGER ilowerf1,iupperf1
      INTEGER ilowerc0,iupperc0
      INTEGER ilowerc1,iupperc1
      INTEGER U_fine_gcw,U_crse_gcw

      INTEGER location_index,ratio_to_coarser(0:NDIM-1)

      INTEGER blowerf(0:NDIM-1), bupperf(0:NDIM-1)

      REAL U_crse(CELL2d(ilowerc,iupperc,U_crse_gcw))
c
c     Input/Output.
c
      REAL U_fine(CELL2d(ilowerf,iupperf,U_fine_gcw))
c
c     Local variables.
c
      INTEGER i_wgt
      INTEGER i,ibeg,iend,i_c,i_p
      INTEGER j,jbeg,jend,j_c,j_p
      REAL R0,R1,t0,t1,U_intrp(-1:1)
      REAL wgt(-1:1,0:MAX_RATIO-1)
c
c     Set the values along the appropriate side.
c
      ibeg = max(blowerf(0),ilowerf0)
      iend = min(bupperf(0),iupperf0)

      jbeg = max(blowerf(1),ilowerf1)
      jend = min(bupperf(1),iupperf1)

      if ( (location_index .eq. 0) .or.
     &     (location_index .eq. 1) ) then
c
c     Set the values along the upper/lower x side of the patch.
c
         R1 = dble(ratio_to_coarser(1))
         do j_p = 0,ratio_to_coarser(1) - 1
            t1 = dble(j_p) + 0.5d0
            wgt(-1,j_p) =
     &           eighth*(4.d0*t1*t1-8.d0*t1*R1+3.d0*R1*R1)/(R1*R1)
            wgt( 0,j_p) =
     &           fourth*(-4.d0*t1*t1+4.d0*t1*R1+3.d0*R1*R1)/(R1*R1)
            wgt(+1,j_p) =
     &           -eighth*(-4.d0*t1*t1+R1*R1)/(R1*R1)
         enddo

         do j = jbeg,jend,ratio_to_coarser(1)
            if (j .lt. 0) then
               j_c = (j+1)/ratio_to_coarser(1)-1
            else
               j_c = j/ratio_to_coarser(1)
            endif

            do i = blowerf(0),bupperf(0)
               if (i .lt. 0) then
                  i_c = (i+1)/ratio_to_coarser(0)-1
               else
                  i_c = i/ratio_to_coarser(0)
               endif

               do i_wgt = -1,1
                  U_intrp(i_wgt) = U_crse(i_c,j_c+i_wgt)
               enddo

               do j_p = 0,ratio_to_coarser(1) - 1
                  U_fine(i,j+j_p) = 0.d0
                  do i_wgt = -1,1
                     U_fine(i,j+j_p) = U_fine(i,j+j_p)
     &                    + wgt(i_wgt,j_p)*U_intrp(i_wgt)
                  enddo
               enddo

            enddo

         enddo

      elseif ( (location_index .eq. 2) .or.
     &         (location_index .eq. 3) ) then
c
c     Set the values along the upper/lower y side of the patch.
c
         R0 = dble(ratio_to_coarser(0))
         do i_p = 0,ratio_to_coarser(0) - 1
            t0 = dble(i_p) + 0.5d0
            wgt(-1,i_p) =
     &           eighth*(4.d0*t0*t0-8.d0*t0*R0+3.d0*R0*R0)/(R0*R0)
            wgt( 0,i_p) =
     &           fourth*(-4.d0*t0*t0+4.d0*t0*R0+3.d0*R0*R0)/(R0*R0)
            wgt(+1,i_p) =
     &           -eighth*(-4.d0*t0*t0+R0*R0)/(R0*R0)
         enddo

         do j = blowerf(1),bupperf(1)
            if (j .lt. 0) then
               j_c = (j+1)/ratio_to_coarser(1)-1
            else
               j_c = j/ratio_to_coarser(1)
            endif

            do i = ibeg,iend,ratio_to_coarser(0)
               if (i .lt. 0) then
                  i_c = (i+1)/ratio_to_coarser(0)-1
               else
                  i_c = i/ratio_to_coarser(0)
               endif

               do i_wgt = -1,1
                  U_intrp(i_wgt) = U_crse(i_c+i_wgt,j_c)
               enddo

               do i_p = 0,ratio_to_coarser(0) - 1
                  U_fine(i+i_p,j) = 0.d0
                  do i_wgt = -1,1
                     U_fine(i+i_p,j) = U_fine(i+i_p,j)
     &                    + wgt(i_wgt,i_p)*U_intrp(i_wgt)
                  enddo
               enddo

            enddo

         enddo

      endif
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute a cell-centered quadratic interpolation in the normal
c     direction at coarse-fine interfaces.
c
c     NOTE: This routine imposes a maximum ghost cell width of 8.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ccquadnormalinterpolation2d(
     &     U,U_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     location_index,ratio_to_coarser,
     &     blower,bupper)
c
      implicit none
c
c     Constants.
c
      INTEGER MAX_GCW
      parameter (MAX_GCW=8)     ! max ghost cell width = 8
c
c     Input.
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER U_gcw

      INTEGER location_index,ratio_to_coarser(0:NDIM-1)

      INTEGER blower(0:NDIM-1), bupper(0:NDIM-1)
c
c     Input/Output.
c
      REAL U(CELL2d(ilower,iupper,U_gcw))
c
c     Local variables.
c
      INTEGER i,ibeg,iend,igcw,i_p
      INTEGER j,jbeg,jend,jgcw,j_p
      REAL R0,R1,t0,t1
      REAL wgt(-1:1,0:MAX_GCW-1)
c
c     Set the values along the appropriate side.
c
      ibeg = max(blower(0),ilower0)
      iend = min(bupper(0),iupper0)
      igcw = abs(bupper(0)-blower(0))+1

      jbeg = max(blower(1),ilower1)
      jend = min(bupper(1),iupper1)
      jgcw = abs(bupper(1)-blower(1))+1

      if ( (location_index .eq. 0) .or.
     &     (location_index .eq. 1) ) then
c
c     Set the values along the upper/lower x side of the patch.
c
         R0 = dble(ratio_to_coarser(0))
         do i_p = 0,igcw-1
            t0 = dble(i_p) + 0.5d0
            wgt(-1,i_p) =
     &           -0.5d0*(2.d0*R0*t0+R0-4.d0*t0*t0-2.d0*t0)/(R0+3.d0)
            wgt( 0,i_p) =
     &           0.5d0*(2.d0*R0*t0+3.d0*R0-4.d0*t0*t0-6.d0*t0)/(R0+1.d0)
            wgt(+1,i_p) =
     &           (4.d0*t0*t0+8.d0*t0+3.d0)/(R0+1.d0)/(R0+3.d0)
         enddo

         if (location_index .eq. 0) then
            do j = jbeg,jend
               do i = 0,igcw-1
                  U(bupper(0)-i,j) =
     &                 wgt(+1,i)*U(bupper(0)-i,j)+
     &                 wgt( 0,i)*U(bupper(0)+1,j)+
     &                 wgt(-1,i)*U(bupper(0)+2,j)
               enddo
            enddo
         else
            do j = jbeg,jend
               do i = 0,igcw-1
                  U(blower(0)+i,j) =
     &                 wgt(+1,i)*U(blower(0)+i,j)+
     &                 wgt( 0,i)*U(blower(0)-1,j)+
     &                 wgt(-1,i)*U(blower(0)-2,j)
               enddo
            enddo
         endif

      elseif ( (location_index .eq. 2) .or.
     &         (location_index .eq. 3) ) then
c
c     Set the values along the upper/lower y side of the patch.
c
         R1 = dble(ratio_to_coarser(1))
         do j_p = 0,jgcw-1
            t1 = dble(j_p) + 0.5d0
            wgt(-1,j_p) =
     &           -0.5d0*(2.d0*R1*t1+R1-4.d0*t1*t1-2.d0*t1)/(R1+3.d0)
            wgt( 0,j_p) =
     &           0.5d0*(2.d0*R1*t1+3.d0*R1-4.d0*t1*t1-6.d0*t1)/(R1+1.d0)
            wgt(+1,j_p) =
     &           (4.d0*t1*t1+8.d0*t1+3.d0)/(R1+1.d0)/(R1+3.d0)
         enddo

         if (location_index .eq. 2) then
            do j = 0,jgcw-1
               do i = ibeg,iend
                  U(i,bupper(1)-j) =
     &                 wgt(+1,j)*U(i,bupper(1)-j)+
     &                 wgt( 0,j)*U(i,bupper(1)+1)+
     &                 wgt(-1,j)*U(i,bupper(1)+2)
               enddo
            enddo
         else
            do j = 0,jgcw-1
               do i = ibeg,iend
                  U(i,blower(1)+j) =
     &                 wgt(+1,j)*U(i,blower(1)+j)+
     &                 wgt( 0,j)*U(i,blower(1)-1)+
     &                 wgt(-1,j)*U(i,blower(1)-2)
               enddo
            enddo
         endif

      endif
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute a side-centered quadratic interpolation in the tangential
c     direction at coarse-fine interfaces.
c
c     NOTE: This routine imposes a maximum refinement ratio of 16.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine scquadtangentialinterpolation2d(
     &     U_fine0,U_fine1,U_fine_gcw,
     &     U_crse0,U_crse1,U_crse_gcw,
     &     indicator0,indicator1,indicator_gcw,
     &     ilowerf0,iupperf0,
     &     ilowerf1,iupperf1,
     &     ilowerc0,iupperc0,
     &     ilowerc1,iupperc1,
     &     location_index,ratio_to_coarser,
     &     blowerf,bupperf)
c
      implicit none
c
c     Constants.
c
      INTEGER MAX_RATIO
      parameter (MAX_RATIO=16)  ! max refinement ratio = 16
      REAL half,fourth,sixth,eighth
      parameter (half=0.5d0)
      parameter (fourth=0.25d0)
      parameter (sixth=0.16666666666667d0)
      parameter (eighth=0.125d0)
c
c     Input.
c
      INTEGER ilowerf0,iupperf0
      INTEGER ilowerf1,iupperf1
      INTEGER ilowerc0,iupperc0
      INTEGER ilowerc1,iupperc1
      INTEGER U_fine_gcw,U_crse_gcw,indicator_gcw

      INTEGER location_index,ratio_to_coarser(0:NDIM-1)

      INTEGER blowerf(0:NDIM-1), bupperf(0:NDIM-1)

      REAL U_crse0(
     &     SIDE2d0(ilowerc,iupperc,U_crse_gcw)
     &     )
      REAL U_crse1(
     &     SIDE2d1(ilowerc,iupperc,U_crse_gcw)
     &     )

      INTEGER indicator0(
     &     SIDE2d0(ilowerf,iupperf,indicator_gcw)
     &     )
      INTEGER indicator1(
     &     SIDE2d1(ilowerf,iupperf,indicator_gcw)
     &     )
c
c     Input/Output.
c
      REAL U_fine0(
     &     SIDE2d0(ilowerf,iupperf,U_fine_gcw)
     &     )
      REAL U_fine1(
     &     SIDE2d1(ilowerf,iupperf,U_fine_gcw)
     &     )
c
c     Local variables.
c
      INTEGER i_wgt,j_wgt
      INTEGER i,ibeg,iend,i_c,i_f,i_p
      INTEGER j,jbeg,jend,j_c,j_f,j_p
      REAL R0,R1,t0,t1,U_intrp(-1:2)
      REAL quad_wgt(-1:1,0:MAX_RATIO)
      REAL cubic_wgt(-1:2,0:MAX_RATIO)
c
c     Set the values along the appropriate side.
c
      ibeg = max(blowerf(0),ilowerf0)
      iend = min(bupperf(0),iupperf0)

      jbeg = max(blowerf(1),ilowerf1)
      jend = min(bupperf(1),iupperf1)

      if ( (location_index .eq. 0) .or.
     &     (location_index .eq. 1) ) then
c
c     Set the values along the upper/lower x side of the patch.
c
         R1 = dble(ratio_to_coarser(1))
         do j_p = 0,ratio_to_coarser(1) - 1
            t1 = dble(j_p) + 0.5d0
            quad_wgt(-1,j_p) =
     &           eighth*(4.d0*t1*t1-8.d0*t1*R1+3.d0*R1*R1)/(R1*R1)
            quad_wgt( 0,j_p) =
     &           fourth*(-4.d0*t1*t1+4.d0*t1*R1+3.d0*R1*R1)/(R1*R1)
            quad_wgt(+1,j_p) =
     &           -eighth*(-4.d0*t1*t1+R1*R1)/(R1*R1)
c$$$            quartic_wgt(-2,j_p) = 1.d0/384.d0*(16.d0*t1*t1*t1*t1-64.d0*
c$$$     &           t1*t1*t1*R1+56.d0*R1*R1*t1*t1+16.d0*R1*R1*R1*t1-15.d0*
c$$$     &           R1*R1*R1*R1)/(R1*R1*R1*R1)
c$$$            quartic_wgt(-1,j_p) = -1.d0/96.d0*(16.d0*t1*t1*t1*t1-48.d0*
c$$$     &           t1*t1*t1*R1-16.d0*R1*R1*t1*t1+108.d0*R1*R1*R1*t1-45.d0*
c$$$     &           R1*R1*R1*R1)/(R1*R1*R1*R1)
c$$$            quartic_wgt( 0,j_p) = 1.d0/64.d0*(16.d0*t1*t1*t1*t1-32.d0*
c$$$     &           t1*t1*t1*R1-56.d0*R1*R1*t1*t1+72.d0*R1*R1*R1*t1+45.d0*
c$$$     &           R1*R1*R1*R1)/(R1*R1*R1*R1)
c$$$            quartic_wgt( 1,j_p) = -1.d0/96.d0*(-64.d0*R1*R1*t1*t1+4.d0*
c$$$     &           R1*R1*R1*t1-16.d0*t1*t1*t1*R1+16.d0*t1*t1*t1*t1+15.d0*
c$$$     &           R1*R1*R1*R1)/(R1*R1*R1*R1)
c$$$            quartic_wgt( 2,j_p) = 1.d0/384.d0*(16.d0*t1*t1*t1*t1-40.d0*
c$$$     &           R1*R1*t1*t1+9.d0*R1*R1*R1*R1)/(R1*R1*R1*R1)
         enddo
         do j_p = 0,ratio_to_coarser(1)
            t1 = dble(j_p)
            cubic_wgt(-1,j_p) =
     &           -sixth*t1*(t1*t1-3.d0*t1*R1+2.d0*R1*R1)/(R1*R1*R1)
            cubic_wgt( 0,j_p) =
     &           half*(t1*t1*t1-2.d0*t1*t1*R1-t1*R1*R1+2.d0*R1*R1*R1)/
     &           (R1*R1*R1)
            cubic_wgt(+1,j_p) =
     &           half*t1*(-t1*t1+t1*R1+2.d0*R1*R1)/(R1*R1*R1)
            cubic_wgt(+2,j_p) =
     &           -sixth*t1*(-t1*t1+R1*R1)/(R1*R1*R1)
         enddo

         do j = jbeg,jend,ratio_to_coarser(1)
            if (j .lt. 0) then
               j_c = (j+1)/ratio_to_coarser(1)-1
            else
               j_c = j/ratio_to_coarser(1)
            endif

         do i = blowerf(0),bupperf(0)
            if (i .lt. 0) then
               i_c = (i+1)/ratio_to_coarser(0)-1
            else
               i_c = i/ratio_to_coarser(0)
            endif

c     Shift i_c and i_f for interpolating the normal component along the
c     coarse-fine interface.

            if (location_index .eq. 0) then
               i_c = i_c
               i_f = i
            else
               i_c = i_c+1
               i_f = i  +1
            endif

c     Perform quadratic interpolation in the y-direction.

            do j_wgt = -1,1
               U_intrp(j_wgt) = U_crse0(i_c,j_c+j_wgt)
            enddo

            do j_p = 0,ratio_to_coarser(1) - 1
               if (.not.(indicator0(i_f,j+j_p).eq.1)) then
                  U_fine0(i_f,j+j_p) = 0.d0
                  do j_wgt = -1,1
                     U_fine0(i_f,j+j_p) = U_fine0(i_f,j+j_p) +
     &                    quad_wgt(j_wgt,j_p)*U_intrp(j_wgt)
                  enddo
               endif
            enddo

c     Shift i_c and i_f for interpolating the tangential component along
c     the coarse-fine interface.

            if (location_index .eq. 0) then
               i_c = i_c
               i_f = i_f
            else
               i_c = i_c-1
               i_f = i_f-1
            endif

c     Perform cubic interpolation in the y-direction.

            do j_wgt = -1,2
               U_intrp(j_wgt) = U_crse1(i_c,j_c+j_wgt)
            enddo

            do j_p = 0,ratio_to_coarser(1)
               if (.not.(indicator1(i_f,j+j_p).eq.1)) then
                  U_fine1(i_f,j+j_p) = 0.d0
                  do j_wgt = -1,2
                     U_fine1(i_f,j+j_p) = U_fine1(i_f,j+j_p) +
     &                    cubic_wgt(j_wgt,j_p)*U_intrp(j_wgt)
                  enddo
               endif
            enddo

         enddo
         enddo

      elseif ( (location_index .eq. 2) .or.
     &         (location_index .eq. 3) ) then
c
c     Set the values along the upper/lower y side of the patch.
c
         R0 = dble(ratio_to_coarser(0))
         do i_p = 0,ratio_to_coarser(0) - 1
            t0 = dble(i_p) + 0.5d0
            quad_wgt(-1,i_p) =
     &           eighth*(4.d0*t0*t0-8.d0*t0*R0+3.d0*R0*R0)/(R0*R0)
            quad_wgt( 0,i_p) =
     &           fourth*(-4.d0*t0*t0+4.d0*t0*R0+3.d0*R0*R0)/(R0*R0)
            quad_wgt(+1,i_p) =
     &           -eighth*(-4.d0*t0*t0+R0*R0)/(R0*R0)
c$$$            quartic_wgt(-2,i_p) = 1.d0/384.d0*(16.d0*t0*t0*t0*t0-64.d0*
c$$$     &           t0*t0*t0*R0+56.d0*R0*R0*t0*t0+16.d0*R0*R0*R0*t0-15.d0*
c$$$     &           R0*R0*R0*R0)/(R0*R0*R0*R0)
c$$$            quartic_wgt(-1,i_p) = -1.d0/96.d0*(16.d0*t0*t0*t0*t0-48.d0*
c$$$     &           t0*t0*t0*R0-16.d0*R0*R0*t0*t0+108.d0*R0*R0*R0*t0-45.d0*
c$$$     &           R0*R0*R0*R0)/(R0*R0*R0*R0)
c$$$            quartic_wgt( 0,i_p) = 1.d0/64.d0*(16.d0*t0*t0*t0*t0-32.d0*
c$$$     &           t0*t0*t0*R0-56.d0*R0*R0*t0*t0+72.d0*R0*R0*R0*t0+45.d0*
c$$$     &           R0*R0*R0*R0)/(R0*R0*R0*R0)
c$$$            quartic_wgt( 1,i_p) = -1.d0/96.d0*(-64.d0*R0*R0*t0*t0+4.d0*
c$$$     &           R0*R0*R0*t0-16.d0*t0*t0*t0*R0+16.d0*t0*t0*t0*t0+15.d0*
c$$$     &           R0*R0*R0*R0)/(R0*R0*R0*R0)
c$$$            quartic_wgt( 2,i_p) = 1.d0/384.d0*(16.d0*t0*t0*t0*t0-40.d0*
c$$$     &           R0*R0*t0*t0+9.d0*R0*R0*R0*R0)/(R0*R0*R0*R0)
         enddo
         do i_p = 0,ratio_to_coarser(0)
            t0 = dble(i_p)
            cubic_wgt(-1,i_p) =
     &           -sixth*t0*(t0*t0-3.d0*t0*R0+2.d0*R0*R0)/(R0*R0*R0)
            cubic_wgt( 0,i_p) =
     &           half*(t0*t0*t0-2.d0*t0*t0*R0-t0*R0*R0+2.d0*R0*R0*R0)/
     &           (R0*R0*R0)
            cubic_wgt(+1,i_p) =
     &           half*t0*(-t0*t0+t0*R0+2.d0*R0*R0)/(R0*R0*R0)
            cubic_wgt(+2,i_p) =
     &           -sixth*t0*(-t0*t0+R0*R0)/(R0*R0*R0)
         enddo

         do j = blowerf(1),bupperf(1)
            if (j .lt. 0) then
               j_c = (j+1)/ratio_to_coarser(1)-1
            else
               j_c = j/ratio_to_coarser(1)
            endif

         do i = ibeg,iend,ratio_to_coarser(0)
            if (i .lt. 0) then
               i_c = (i+1)/ratio_to_coarser(0)-1
            else
               i_c = i/ratio_to_coarser(0)
            endif

c     Shift j_c and j_f for interpolating the normal component along the
c     coarse-fine interface.

            if (location_index .eq. 2) then
               j_c = j_c
               j_f = j
            else
               j_c = j_c+1
               j_f = j  +1
            endif

c     Perform quadratic interpolation in the x-direction.

            do i_wgt = -1,1
               U_intrp(i_wgt) = U_crse1(i_c+i_wgt,j_c)
            enddo

            do i_p = 0,ratio_to_coarser(0) - 1
               if (.not.(indicator1(i+i_p,j_f).eq.1)) then
                  U_fine1(i+i_p,j_f) = 0.d0
                  do i_wgt = -1,1
                     U_fine1(i+i_p,j_f) = U_fine1(i+i_p,j_f) +
     &                    quad_wgt(i_wgt,i_p)*U_intrp(i_wgt)
                  enddo
               endif
            enddo

c     Shift j_c and j_f for interpolating the tangential component along
c     the coarse-fine interface.

            if (location_index .eq. 2) then
               j_c = j_c
               j_f = j_f
            else
               j_c = j_c-1
               j_f = j_f-1
            endif

c     Perform cubic interpolation in the x-direction.

            do i_wgt = -1,2
               U_intrp(i_wgt) = U_crse0(i_c+i_wgt,j_c)
            enddo

            do i_p = 0,ratio_to_coarser(0)
               if (.not.(indicator0(i+i_p,j_f).eq.1)) then
                  U_fine0(i+i_p,j_f) = 0.d0
                  do i_wgt = -1,2
                     U_fine0(i+i_p,j_f) = U_fine0(i+i_p,j_f) +
     &                    cubic_wgt(i_wgt,i_p)*U_intrp(i_wgt)
                  enddo
               endif
            enddo

         enddo
         enddo

      endif
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute a side-centered quadratic interpolation in the normal
c     direction at coarse-fine interfaces.
c
c     NOTE: This routine imposes a maximum ghost cell width of 8.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine scquadnormalinterpolation2d(
     &     U0,U1,U_gcw,
     &     W0,W1,W_gcw,
     &     indicator0,indicator1,indicator_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     location_index,ratio_to_coarser,
     &     blower,bupper)
c
      implicit none
c
c     Constants.
c
      INTEGER MAX_GCW
      parameter (MAX_GCW=8)     ! max ghost cell width = 8
c
c     Input.
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER U_gcw,W_gcw,indicator_gcw

      INTEGER location_index,ratio_to_coarser(0:NDIM-1)

      INTEGER blower(0:NDIM-1), bupper(0:NDIM-1)

      REAL W0(
     &     SIDE2d0(ilower,iupper,W_gcw)
     &     )
      REAL W1(
     &     SIDE2d1(ilower,iupper,W_gcw)
     &     )

      INTEGER indicator0(
     &     SIDE2d0(ilower,iupper,indicator_gcw)
     &     )
      INTEGER indicator1(
     &     SIDE2d1(ilower,iupper,indicator_gcw)
     &     )
c
c     Input/Output.
c
      REAL U0(
     &     SIDE2d0(ilower,iupper,U_gcw)
     &     )
      REAL U1(
     &     SIDE2d1(ilower,iupper,U_gcw)
     &     )
c
c     Local variables.
c
      INTEGER i,ibeg,iend,igcw,i_p
      INTEGER j,jbeg,jend,jgcw,j_p
      REAL R0,R1,t0,t1
      REAL quad_wgt(-1:1,0:MAX_GCW-1)
c
c     Set the values along the appropriate side.
c
      ibeg = max(blower(0),ilower0)
      iend = min(bupper(0),iupper0)
      igcw = abs(bupper(0)-blower(0))+1

      jbeg = max(blower(1),ilower1)
      jend = min(bupper(1),iupper1)
      jgcw = abs(bupper(1)-blower(1))+1

      if ( (location_index .eq. 0) .or.
     &     (location_index .eq. 1) ) then
c
c     Set the values along the upper/lower x side of the patch.
c
         R0 = dble(ratio_to_coarser(0))
         do i_p = 0,igcw-1
            t0 = dble(i_p) + 1.d0
            quad_wgt(-1,i_p) = -t0*(-t0+R0)/(1.d0+R0)
            quad_wgt( 0,i_p) = (R0*t0+R0-t0*t0-t0)/R0
            quad_wgt(+1,i_p) = t0*(t0+1.d0)/R0/(1.d0+R0)
c$$$            cubic_wgt(-2,i_p) =
c$$$     &           1.d0/2.d0*(t0*R0+R0-t0*t0-t0)*t0/(R0+2.d0)
c$$$            cubic_wgt(-1,i_p) =
c$$$     &           -(t0*R0+2.d0*R0-t0*t0-2.d0*t0)*t0/(R0+1.d0)
c$$$            cubic_wgt( 0,i_p) =
c$$$     &           1.d0/2.d0*(3.d0*t0*R0+2.d0*R0+t0*t0*R0-3*t0*t0-t0*t0*t0
c$$$     &           -2.d0*t0)/R0
c$$$            cubic_wgt(+1,i_p) =
c$$$     &           t0*(3.d0*t0+2.d0+t0*t0)/(R0+2.d0)/(R0+1.d0)/R0
         enddo

         if (location_index .eq. 0) then
            do j = jbeg,jend
               do i = 0,igcw-1
                  if (.not.(indicator0(bupper(0)-i  ,j).eq.1)) then
                     U0(bupper(0)-i  ,j) =
     &                    quad_wgt(+1,i)*W0(bupper(0)-i  ,j)+
     &                    quad_wgt( 0,i)*W0(bupper(0)+1  ,j)+
     &                    quad_wgt(-1,i)*W0(bupper(0)+2  ,j)
                  endif
               enddo
            enddo
         else
            do j = jbeg,jend
               do i = 0,igcw-1
                  if (.not.(indicator0(blower(0)+i+1,j).eq.1)) then
                     U0(blower(0)+i+1,j) =
     &                    quad_wgt(+1,i)*W0(blower(0)+i+1,j)+
     &                    quad_wgt( 0,i)*W0(blower(0)-1+1,j)+
     &                    quad_wgt(-1,i)*W0(blower(0)-2+1,j)
                  endif
               enddo
            enddo
         endif

         R0 = dble(ratio_to_coarser(0))
         do i_p = 0,igcw-1
            t0 = dble(i_p) + 0.5d0
            quad_wgt(-1,i_p) =
     &           0.5d0*(-2.d0*R0*t0-R0+4.d0*t0*t0+2.d0*t0)/(R0+3.d0)
            quad_wgt( 0,i_p) =
     &           -0.5d0*(-2.d0*R0*t0-3.d0*R0+4.d0*t0*t0+6.d0*t0)/
     &           (R0+1.d0)
            quad_wgt(+1,i_p) =
     &           (4.d0*t0*t0+3.d0+8.d0*t0)/(R0+3.d0)/(R0+1.d0)
c$$$            cubic_wgt(-2,i_p) =
c$$$     &           1.d0/8.d0*(8.d0*t0*R0+4.d0*t0*t0*R0+3.d0*R0-16.d0*t0*t0
c$$$     &           -8.d0*t0*t0*t0-6.d0*t0)/(R0+5.d0)
c$$$            cubic_wgt(-1,i_p) =
c$$$     &           -1.d0/4.d0*(4.d0*t0*t0*R0+5.d0*R0+12.d0*t0*R0-24.d0*t0*
c$$$     &           t0-8.d0*t0*t0*t0-10.d0*t0)/(R0+3.d0)
c$$$            cubic_wgt( 0,i_p) =
c$$$     &           1.d0/8.d0*(16.d0*t0*R0+15.d0*R0+4.d0*t0*t0*R0-8.d0*t0*
c$$$     &           t0*t0-30.d0*t0-32.d0*t0*t0)/(R0+1.d0)
c$$$            cubic_wgt(+1,i_p) =
c$$$     &           (8.d0*t0*t0*t0+36.d0*t0*t0+46.d0*t0+15.d0)/(R0+5.d0)/
c$$$     &           (R0+3.d0)/(R0+1.d0)
         enddo

         if (location_index .eq. 0) then
            do j = jbeg,jend+1
               do i = 0,igcw-1
                  if (.not.(indicator1(bupper(0)-i,j).eq.1)) then
                     U1(bupper(0)-i,j) =
     &                    quad_wgt(+1,i)*W1(bupper(0)-i,j)+
     &                    quad_wgt( 0,i)*W1(bupper(0)+1,j)+
     &                    quad_wgt(-1,i)*W1(bupper(0)+2,j)
                  endif
               enddo
            enddo
         else
            do j = jbeg,jend+1
               do i = 0,igcw-1
                  if (.not.(indicator1(blower(0)+i,j).eq.1)) then
                     U1(blower(0)+i,j) =
     &                    quad_wgt(+1,i)*W1(blower(0)+i,j)+
     &                    quad_wgt( 0,i)*W1(blower(0)-1,j)+
     &                    quad_wgt(-1,i)*W1(blower(0)-2,j)
                  endif
               enddo
            enddo
         endif

      elseif ( (location_index .eq. 2) .or.
     &         (location_index .eq. 3) ) then
c
c     Set the values along the upper/lower y side of the patch.
c
         R1 = dble(ratio_to_coarser(1))
         do j_p = 0,jgcw-1
            t1 = dble(j_p) + 1.d0
            quad_wgt(-1,j_p) = -t1*(-t1+R1)/(1.d0+R1)
            quad_wgt( 0,j_p) = (R1*t1+R1-t1*t1-t1)/R1
            quad_wgt(+1,j_p) = t1*(t1+1.d0)/R1/(1.d0+R1)
c$$$            cubic_wgt(-2,j_p) =
c$$$     &           1.d0/2.d0*(t1*R1+R1-t1*t1-t1)*t1/(R1+2.d0)
c$$$            cubic_wgt(-1,j_p) =
c$$$     &           -(t1*R1+2.d0*R1-t1*t1-2.d0*t1)*t1/(R1+1.d0)
c$$$            cubic_wgt( 0,j_p) =
c$$$     &           1.d0/2.d0*(3.d0*t1*R1+2.d0*R1+t1*t1*R1-3*t1*t1-t1*t1*t1
c$$$     &           -2.d0*t1)/R1
c$$$            cubic_wgt(+1,j_p) =
c$$$     &           t1*(3.d0*t1+2.d0+t1*t1)/(R1+2.d0)/(R1+1.d0)/R1
         enddo

         if (location_index .eq. 2) then
            do j = 0,jgcw-1
               do i = ibeg,iend
                  if (.not.(indicator1(i,bupper(1)-j  ).eq.1)) then
                     U1(i,bupper(1)-j  ) =
     &                    quad_wgt(+1,j)*W1(i,bupper(1)-j  )+
     &                    quad_wgt( 0,j)*W1(i,bupper(1)+1  )+
     &                    quad_wgt(-1,j)*W1(i,bupper(1)+2  )
                  endif
               enddo
            enddo
         else
            do j = 0,jgcw-1
               do i = ibeg,iend
                  if (.not.(indicator1(i,blower(1)+j+1).eq.1)) then
                     U1(i,blower(1)+j+1) =
     &                    quad_wgt(+1,j)*W1(i,blower(1)+j+1)+
     &                    quad_wgt( 0,j)*W1(i,blower(1)-1+1)+
     &                    quad_wgt(-1,j)*W1(i,blower(1)-2+1)
                  endif
               enddo
            enddo
         endif

         R1 = dble(ratio_to_coarser(1))
         do j_p = 0,jgcw-1
            t1 = dble(j_p) + 0.5d0
            quad_wgt(-1,j_p) =
     &           0.5d0*(-2.d0*R1*t1-R1+4.d0*t1*t1+2.d0*t1)/(R1+3.d0)
            quad_wgt( 0,j_p) =
     &           -0.5d0*(-2.d0*R1*t1-3.d0*R1+4.d0*t1*t1+6.d0*t1)/
     &           (R1+1.d0)
            quad_wgt(+1,j_p) =
     &           (4.d0*t1*t1+3.d0+8.d0*t1)/(R1+3.d0)/(R1+1.d0)
c$$$            cubic_wgt(-2,j_p) =
c$$$     &           1.d0/8.d0*(8.d0*t1*R1+4.d0*t1*t1*R1+3.d0*R1-16.d0*t1*t1
c$$$     &           -8.d0*t1*t1*t1-6.d0*t1)/(R1+5.d0)
c$$$            cubic_wgt(-1,j_p) =
c$$$     &           -1.d0/4.d0*(4.d0*t1*t1*R1+5.d0*R1+12.d0*t1*R1-24.d0*t1*
c$$$     &           t1-8.d0*t1*t1*t1-10.d0*t1)/(R1+3.d0)
c$$$            cubic_wgt( 0,j_p) =
c$$$     &           1.d0/8.d0*(16.d0*t1*R1+15.d0*R1+4.d0*t1*t1*R1-8.d0*t1*
c$$$     &           t1*t1-30.d0*t1-32.d0*t1*t1)/(R1+1.d0)
c$$$            cubic_wgt(+1,j_p) =
c$$$     &           (8.d0*t1*t1*t1+36.d0*t1*t1+46.d0*t1+15.d0)/(R1+5.d0)/
c$$$     &           (R1+3.d0)/(R1+1.d0)
         enddo

         if (location_index .eq. 2) then
            do j = 0,jgcw-1
               do i = ibeg,iend+1
                  if (.not.(indicator0(i,bupper(1)-j).eq.1)) then
                     U0(i,bupper(1)-j) =
     &                    quad_wgt(+1,j)*W0(i,bupper(1)-j)+
     &                    quad_wgt( 0,j)*W0(i,bupper(1)+1)+
     &                    quad_wgt(-1,j)*W0(i,bupper(1)+2)
                  endif
               enddo
            enddo
         else
            do j = 0,jgcw-1
               do i = ibeg,iend+1
                  if (.not.(indicator0(i,blower(1)+j).eq.1)) then
                     U0(i,blower(1)+j) =
     &                    quad_wgt(+1,j)*W0(i,blower(1)+j)+
     &                    quad_wgt( 0,j)*W0(i,blower(1)-1)+
     &                    quad_wgt(-1,j)*W0(i,blower(1)-2)
                  endif
               enddo
            enddo
         endif

      endif
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
