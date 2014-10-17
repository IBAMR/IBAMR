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
define(NDIM,3)dnl
define(REAL,`double precision')dnl
define(INTEGER,`integer')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim3d.i)dnl
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
      subroutine ccquadtangentialinterpolation3d(
     &     U_fine,U_fine_gcw,
     &     U_crse,U_crse_gcw,
     &     ilowerf0,iupperf0,
     &     ilowerf1,iupperf1,
     &     ilowerf2,iupperf2,
     &     ilowerc0,iupperc0,
     &     ilowerc1,iupperc1,
     &     ilowerc2,iupperc2,
     &     location_index,ratio_to_coarser,
     &     blowerf,bupperf)
c
      implicit none
c
c     Constants.
c
      INTEGER MAX_RATIO
      parameter (MAX_RATIO=16)  ! max refinement ratio = 16
c
c     Input.
c
      INTEGER ilowerf0,iupperf0
      INTEGER ilowerf1,iupperf1
      INTEGER ilowerf2,iupperf2
      INTEGER ilowerc0,iupperc0
      INTEGER ilowerc1,iupperc1
      INTEGER ilowerc2,iupperc2
      INTEGER U_fine_gcw,U_crse_gcw

      INTEGER location_index,ratio_to_coarser(0:NDIM-1)

      INTEGER blowerf(0:NDIM-1), bupperf(0:NDIM-1)

      REAL U_crse(CELL3d(ilowerc,iupperc,U_crse_gcw))
c
c     Input/Output.
c
      REAL U_fine(CELL3d(ilowerf,iupperf,U_fine_gcw))
c
c     Local variables.
c
      INTEGER i_wgt0,i_wgt1
      INTEGER i,ibeg,iend,i_c,i_p
      INTEGER j,jbeg,jend,j_c,j_p
      INTEGER k,kbeg,kend,k_c,k_p
      REAL R0,R1,R2,t0,t1,t2,U_intrp(-1:1,-1:1)
      REAL wgt0(-1:1,0:MAX_RATIO-1)
      REAL wgt1(-1:1,0:MAX_RATIO-1)
c
c     Set the values along the appropriate side.
c
      ibeg = max(blowerf(0),ilowerf0)
      iend = min(bupperf(0),iupperf0)

      jbeg = max(blowerf(1),ilowerf1)
      jend = min(bupperf(1),iupperf1)

      kbeg = max(blowerf(2),ilowerf2)
      kend = min(bupperf(2),iupperf2)

      if ( (location_index .eq. 0) .or.
     &     (location_index .eq. 1) ) then
c
c     Pre-compute the quadratic interpolation weights.
c
         R1 = dble(ratio_to_coarser(1))
         do j_p = 0,ratio_to_coarser(1) - 1
            t1 = dble(j_p) + 0.5d0 - 0.5d0*R1
            wgt0(-1,j_p) = 0.5d0*t1*t1/(R1*R1)-0.5d0*t1/R1
            wgt0( 0,j_p) = -t1*t1/(R1*R1)+1.d0
            wgt0(+1,j_p) = 0.5d0*t1*t1/(R1*R1)+0.5d0*t1/R1
         enddo

         R2 = dble(ratio_to_coarser(2))
         do k_p = 0,ratio_to_coarser(2) - 1
            t2 = dble(k_p) + 0.5d0 - 0.5d0*R2
            wgt1(-1,k_p) = 0.5d0*t2*t2/(R2*R2)-0.5d0*t2/R2
            wgt1( 0,k_p) = -t2*t2/(R2*R2)+1.d0
            wgt1(+1,k_p) = 0.5d0*t2*t2/(R2*R2)+0.5d0*t2/R2
         enddo
c
c     Set the values along the upper/lower x side of the patch.
c
         do k = kbeg,kend,ratio_to_coarser(2)
            if (k .lt. 0) then
               k_c = (k+1)/ratio_to_coarser(2)-1
            else
               k_c = k/ratio_to_coarser(2)
            endif

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

                  do i_wgt1 = -1,1
                     do i_wgt0 = -1,1
                        U_intrp(i_wgt0,i_wgt1) =
     &                       U_crse(i_c,j_c+i_wgt0,k_c+i_wgt1)
                     enddo
                  enddo

                  do k_p = 0,ratio_to_coarser(2) - 1
                     do j_p = 0,ratio_to_coarser(1) - 1
                        U_fine(i,j+j_p,k+k_p) = 0.d0
                        do i_wgt1 = -1,1
                           do i_wgt0 = -1,1
                              U_fine(i,j+j_p,k+k_p) =
     &                             U_fine(i,j+j_p,k+k_p) +
     &                             wgt0(i_wgt0,j_p)*wgt1(i_wgt1,k_p)*
     &                             U_intrp(i_wgt0,i_wgt1)
                           enddo
                        enddo
                     enddo
                  enddo

               enddo

            enddo

         enddo

      elseif ( (location_index .eq. 2) .or.
     &         (location_index .eq. 3) ) then
c
c     Pre-compute the quadratic interpolation weights.
c
         R0 = dble(ratio_to_coarser(0))
         do i_p = 0,ratio_to_coarser(0) - 1
            t0 = dble(i_p) + 0.5d0 - 0.5d0*R0
            wgt0(-1,i_p) = 0.5d0*t0*t0/(R0*R0)-0.5d0*t0/R0
            wgt0( 0,i_p) = -t0*t0/(R0*R0)+1.d0
            wgt0(+1,i_p) = 0.5d0*t0*t0/(R0*R0)+0.5d0*t0/R0
         enddo

         R2 = dble(ratio_to_coarser(2))
         do k_p = 0,ratio_to_coarser(2) - 1
            t2 = dble(k_p) + 0.5d0 - 0.5d0*R2
            wgt1(-1,k_p) = 0.5d0*t2*t2/(R2*R2)-0.5d0*t2/R2
            wgt1( 0,k_p) = -t2*t2/(R2*R2)+1.d0
            wgt1(+1,k_p) = 0.5d0*t2*t2/(R2*R2)+0.5d0*t2/R2
         enddo
c
c     Set the values along the upper/lower y side of the patch.
c
         do k = kbeg,kend,ratio_to_coarser(2)
            if (k .lt. 0) then
               k_c = (k+1)/ratio_to_coarser(2)-1
            else
               k_c = k/ratio_to_coarser(2)
            endif

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

                  do i_wgt1 = -1,1
                     do i_wgt0 = -1,1
                        U_intrp(i_wgt0,i_wgt1) =
     &                       U_crse(i_c+i_wgt0,j_c,k_c+i_wgt1)
                     enddo
                  enddo

                  do k_p = 0,ratio_to_coarser(2) - 1
                     do i_p = 0,ratio_to_coarser(0) - 1
                        U_fine(i+i_p,j,k+k_p) = 0.d0
                        do i_wgt1 = -1,1
                           do i_wgt0 = -1,1
                              U_fine(i+i_p,j,k+k_p) =
     &                             U_fine(i+i_p,j,k+k_p) +
     &                             wgt0(i_wgt0,i_p)*wgt1(i_wgt1,k_p)*
     &                             U_intrp(i_wgt0,i_wgt1)
                           enddo
                        enddo
                     enddo
                  enddo

               enddo

            enddo

         enddo

      elseif ( (location_index .eq. 4) .or.
     &         (location_index .eq. 5) ) then
c
c     Pre-compute the quadratic interpolation weights.
c
         R0 = dble(ratio_to_coarser(0))
         do i_p = 0,ratio_to_coarser(0) - 1
            t0 = dble(i_p) + 0.5d0 - 0.5d0*R0
            wgt0(-1,i_p) = 0.5d0*t0*t0/(R0*R0)-0.5d0*t0/R0
            wgt0( 0,i_p) = -t0*t0/(R0*R0)+1.d0
            wgt0(+1,i_p) = 0.5d0*t0*t0/(R0*R0)+0.5d0*t0/R0
         enddo

         R1 = dble(ratio_to_coarser(1))
         do j_p = 0,ratio_to_coarser(1) - 1
            t1 = dble(j_p) + 0.5d0 - 0.5d0*R1
            wgt1(-1,j_p) = 0.5d0*t1*t1/(R1*R1)-0.5d0*t1/R1
            wgt1( 0,j_p) = -t1*t1/(R1*R1)+1.d0
            wgt1(+1,j_p) = 0.5d0*t1*t1/(R1*R1)+0.5d0*t1/R1
         enddo
c
c     Set the values along the upper/lower z side of the patch.
c
         do k = blowerf(2),bupperf(2)
            if (k .lt. 0) then
               k_c = (k+1)/ratio_to_coarser(2)-1
            else
               k_c = k/ratio_to_coarser(2)
            endif

            do j = jbeg,jend,ratio_to_coarser(1)
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

                  do i_wgt1 = -1,1
                     do i_wgt0 = -1,1
                        U_intrp(i_wgt0,i_wgt1) =
     &                       U_crse(i_c+i_wgt0,j_c+i_wgt1,k_c)
                     enddo
                  enddo

                  do j_p = 0,ratio_to_coarser(1) - 1
                     do i_p = 0,ratio_to_coarser(0) - 1
                        U_fine(i+i_p,j+j_p,k) = 0.d0
                        do i_wgt1 = -1,1
                           do i_wgt0 = -1,1
                              U_fine(i+i_p,j+j_p,k) =
     &                             U_fine(i+i_p,j+j_p,k) +
     &                             wgt0(i_wgt0,i_p)*wgt1(i_wgt1,j_p)*
     &                             U_intrp(i_wgt0,i_wgt1)
                           enddo
                        enddo
                     enddo
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
      subroutine ccquadnormalinterpolation3d(
     &     U,U_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
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
      INTEGER ilower2,iupper2
      INTEGER U_gcw

      INTEGER location_index,ratio_to_coarser(0:NDIM-1)

      INTEGER blower(0:NDIM-1), bupper(0:NDIM-1)
c
c     Input/Output.
c
      REAL U(CELL3d(ilower,iupper,U_gcw))
c
c     Local variables.
c
      INTEGER i,ibeg,iend,igcw,i_p
      INTEGER j,jbeg,jend,jgcw,j_p
      INTEGER k,kbeg,kend,kgcw,k_p
      REAL R0,R1,R2,t0,t1,t2
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

      kbeg = max(blower(2),ilower2)
      kend = min(bupper(2),iupper2)
      kgcw = abs(bupper(2)-blower(2))+1

      if ( (location_index .eq. 0) .or.
     &     (location_index .eq. 1) ) then
c
c     Pre-compute the quadratic interpolation weights.
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
c
c     Set the values along the upper/lower x side of the patch.
c
         if (location_index .eq. 0) then
            do k = kbeg,kend
               do j = jbeg,jend
                  do i = 0,igcw-1
                     U(bupper(0)-i,j,k) =
     &                    wgt(+1,i)*U(bupper(0)-i,j,k)+
     &                    wgt( 0,i)*U(bupper(0)+1,j,k)+
     &                    wgt(-1,i)*U(bupper(0)+2,j,k)
                  enddo
               enddo
            enddo
         else
            do k = kbeg,kend
               do j = jbeg,jend
                  do i = 0,igcw-1
                     U(blower(0)+i,j,k) =
     &                    wgt(+1,i)*U(blower(0)+i,j,k)+
     &                    wgt( 0,i)*U(blower(0)-1,j,k)+
     &                    wgt(-1,i)*U(blower(0)-2,j,k)
                  enddo
               enddo
            enddo
         endif

      elseif ( (location_index .eq. 2) .or.
     &         (location_index .eq. 3) ) then
c
c     Pre-compute the quadratic interpolation weights.
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
c
c     Set the values along the upper/lower y side of the patch.
c
         if (location_index .eq. 2) then
            do k = kbeg,kend
               do j = 0,jgcw-1
                  do i = ibeg,iend
                     U(i,bupper(1)-j,k) =
     &                    wgt(+1,j)*U(i,bupper(1)-j,k)+
     &                    wgt( 0,j)*U(i,bupper(1)+1,k)+
     &                    wgt(-1,j)*U(i,bupper(1)+2,k)
                  enddo
               enddo
            enddo
         else
            do k = kbeg,kend
               do j = 0,jgcw-1
                  do i = ibeg,iend
                     U(i,blower(1)+j,k) =
     &                    wgt(+1,j)*U(i,blower(1)+j,k)+
     &                    wgt( 0,j)*U(i,blower(1)-1,k)+
     &                    wgt(-1,j)*U(i,blower(1)-2,k)
                  enddo
               enddo
            enddo
         endif

      elseif ( (location_index .eq. 4) .or.
     &         (location_index .eq. 5) ) then
c
c     Pre-compute the quadratic interpolation weights.
c
         R2 = dble(ratio_to_coarser(2))
         do k_p = 0,kgcw-1
            t2 = dble(k_p) + 0.5d0
            wgt(-1,k_p) =
     &           -0.5d0*(2.d0*R2*t2+R2-4.d0*t2*t2-2.d0*t2)/(R2+3.d0)
            wgt( 0,k_p) =
     &           0.5d0*(2.d0*R2*t2+3.d0*R2-4.d0*t2*t2-6.d0*t2)/(R2+1.d0)
            wgt(+1,k_p) =
     &           (4.d0*t2*t2+8.d0*t2+3.d0)/(R2+1.d0)/(R2+3.d0)
         enddo
c
c     Set the values along the upper/lower z side of the patch.
c
         if (location_index .eq. 4) then
            do k = 0,kgcw-1
               do j = jbeg,jend
                  do i = ibeg,iend
                     U(i,j,bupper(2)-k) =
     &                    wgt(+1,k)*U(i,j,bupper(2)-k)+
     &                    wgt( 0,k)*U(i,j,bupper(2)+1)+
     &                    wgt(-1,k)*U(i,j,bupper(2)+2)
                  enddo
               enddo
            enddo
         else
            do k = 0,kgcw-1
               do j = jbeg,jend
                  do i = ibeg,iend
                     U(i,j,blower(2)+k) =
     &                    wgt(+1,k)*U(i,j,blower(2)+k)+
     &                    wgt( 0,k)*U(i,j,blower(2)-1)+
     &                    wgt(-1,k)*U(i,j,blower(2)-2)
                  enddo
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
      subroutine scquadtangentialinterpolation3d(
     &     U_fine0,U_fine1,U_fine2,U_fine_gcw,
     &     U_crse0,U_crse1,U_crse2,U_crse_gcw,
     &     indicator0,indicator1,indicator2,indicator_gcw,
     &     ilowerf0,iupperf0,
     &     ilowerf1,iupperf1,
     &     ilowerf2,iupperf2,
     &     ilowerc0,iupperc0,
     &     ilowerc1,iupperc1,
     &     ilowerc2,iupperc2,
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
      INTEGER ilowerf2,iupperf2
      INTEGER ilowerc0,iupperc0
      INTEGER ilowerc1,iupperc1
      INTEGER ilowerc2,iupperc2
      INTEGER U_fine_gcw,U_crse_gcw,indicator_gcw

      INTEGER location_index,ratio_to_coarser(0:NDIM-1)

      INTEGER blowerf(0:NDIM-1), bupperf(0:NDIM-1)

      REAL U_crse0(
     &     SIDE3d0(ilowerc,iupperc,U_crse_gcw)
     &     )
      REAL U_crse1(
     &     SIDE3d1(ilowerc,iupperc,U_crse_gcw)
     &     )
      REAL U_crse2(
     &     SIDE3d2(ilowerc,iupperc,U_crse_gcw)
     &     )

      INTEGER indicator0(
     &     SIDE3d0(ilowerf,iupperf,indicator_gcw)
     &     )
      INTEGER indicator1(
     &     SIDE3d1(ilowerf,iupperf,indicator_gcw)
     &     )
      INTEGER indicator2(
     &     SIDE3d2(ilowerf,iupperf,indicator_gcw)
     &     )
c
c     Input/Output.
c
      REAL U_fine0(
     &     SIDE3d0(ilowerf,iupperf,U_fine_gcw)
     &     )
      REAL U_fine1(
     &     SIDE3d1(ilowerf,iupperf,U_fine_gcw)
     &     )
      REAL U_fine2(
     &     SIDE3d2(ilowerf,iupperf,U_fine_gcw)
     &     )
c
c     Local variables.
c
      INTEGER i_wgt,j_wgt,k_wgt
      INTEGER i,ibeg,iend,i_c,i_f,i_p
      INTEGER j,jbeg,jend,j_c,j_f,j_p
      INTEGER k,kbeg,kend,k_c,k_f,k_p
      REAL R0,R1,R2,t0,t1,t2,U_intrp(-1:2,-1:2)
      REAL quad_wgt0(-1:1,0:MAX_RATIO-1)
      REAL quad_wgt1(-1:1,0:MAX_RATIO-1)
      REAL quad_wgt2(-1:1,0:MAX_RATIO-1)
      REAL cubic_wgt0(-1:2,0:MAX_RATIO)
      REAL cubic_wgt1(-1:2,0:MAX_RATIO)
      REAL cubic_wgt2(-1:2,0:MAX_RATIO)
c
c     Set the values along the appropriate side.
c
      ibeg = max(blowerf(0),ilowerf0)
      iend = min(bupperf(0),iupperf0)

      jbeg = max(blowerf(1),ilowerf1)
      jend = min(bupperf(1),iupperf1)

      kbeg = max(blowerf(2),ilowerf2)
      kend = min(bupperf(2),iupperf2)

      if ( (location_index .eq. 0) .or.
     &     (location_index .eq. 1) ) then
c
c     Pre-compute the quadratic and cubic interpolation weights.
c
         R1 = dble(ratio_to_coarser(1))
         do j_p = 0,ratio_to_coarser(1) - 1
            t1 = dble(j_p) + 0.5d0
            quad_wgt1(-1,j_p) =
     &           eighth*(4.d0*t1*t1-8.d0*t1*R1+3.d0*R1*R1)/(R1*R1)
            quad_wgt1( 0,j_p) =
     &           fourth*(-4.d0*t1*t1+4.d0*t1*R1+3.d0*R1*R1)/(R1*R1)
            quad_wgt1(+1,j_p) =
     &           -eighth*(-4.d0*t1*t1+R1*R1)/(R1*R1)
         enddo
         do j_p = 0,ratio_to_coarser(1)
            t1 = dble(j_p)
            cubic_wgt1(-1,j_p) =
     &           -sixth*t1*(t1*t1-3.d0*t1*R1+2.d0*R1*R1)/(R1*R1*R1)
            cubic_wgt1( 0,j_p) =
     &           half*(t1*t1*t1-2.d0*t1*t1*R1-t1*R1*R1+2.d0*R1*R1*R1)/
     &           (R1*R1*R1)
            cubic_wgt1(+1,j_p) =
     &           half*t1*(-t1*t1+t1*R1+2.d0*R1*R1)/(R1*R1*R1)
            cubic_wgt1(+2,j_p) =
     &           -sixth*t1*(-t1*t1+R1*R1)/(R1*R1*R1)
         enddo

         R2 = dble(ratio_to_coarser(2))
         do k_p = 0,ratio_to_coarser(2) - 1
            t2 = dble(k_p) + 0.5d0
            quad_wgt2(-1,k_p) =
     &           eighth*(4.d0*t2*t2-8.d0*t2*R2+3.d0*R2*R2)/(R2*R2)
            quad_wgt2( 0,k_p) =
     &           fourth*(-4.d0*t2*t2+4.d0*t2*R2+3.d0*R2*R2)/(R2*R2)
            quad_wgt2(+1,k_p) =
     &           -eighth*(-4.d0*t2*t2+R2*R2)/(R2*R2)
         enddo
         do k_p = 0,ratio_to_coarser(2)
            t2 = dble(k_p)
            cubic_wgt2(-1,k_p) =
     &           -sixth*t2*(t2*t2-3.d0*t2*R2+2.d0*R2*R2)/(R2*R2*R2)
            cubic_wgt2( 0,k_p) =
     &           half*(t2*t2*t2-2.d0*t2*t2*R2-t2*R2*R2+2.d0*R2*R2*R2)/
     &           (R2*R2*R2)
            cubic_wgt2(+1,k_p) =
     &           half*t2*(-t2*t2+t2*R2+2.d0*R2*R2)/(R2*R2*R2)
            cubic_wgt2(+2,k_p) =
     &           -sixth*t2*(-t2*t2+R2*R2)/(R2*R2*R2)
         enddo
c
c     Set the values along the upper/lower x side of the patch.
c
         do k = kbeg,kend,ratio_to_coarser(2)
            if (k .lt. 0) then
               k_c = (k+1)/ratio_to_coarser(2)-1
            else
               k_c = k/ratio_to_coarser(2)
            endif

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

c     Perform quadratic interpolation in the y- and z-directions.

            do k_wgt = -1,1
               do j_wgt = -1,1
                  U_intrp(j_wgt,k_wgt)=U_crse0(i_c,j_c+j_wgt,k_c+k_wgt)
               enddo
            enddo

            do k_p = 0,ratio_to_coarser(2) - 1
               do j_p = 0,ratio_to_coarser(1) - 1
                  if (.not.(indicator0(i_f,j+j_p,k+k_p).eq.1)) then
                     U_fine0(i_f,j+j_p,k+k_p) = 0.d0
                     do k_wgt = -1,1
                        do j_wgt = -1,1
                           U_fine0(i_f,j+j_p,k+k_p) =
     &                          U_fine0(i_f,j+j_p,k+k_p) +
     &                          quad_wgt1(j_wgt,j_p)*
     &                          quad_wgt2(k_wgt,k_p)*
     &                          U_intrp(j_wgt,k_wgt)
                        enddo
                     enddo
                  endif
               enddo
            enddo

c     Shift i_c and i_f for interpolating the tangential components
c     along the coarse-fine interface.

            if (location_index .eq. 0) then
               i_c = i_c
               i_f = i_f
            else
               i_c = i_c-1
               i_f = i_f-1
            endif

c     Perform cubic interpolation in the y-direction and quadratic
c     interpolation in the z-direction.

            do k_wgt = -1,1
               do j_wgt = -1,2
                  U_intrp(j_wgt,k_wgt)=U_crse1(i_c,j_c+j_wgt,k_c+k_wgt)
               enddo
            enddo

            do k_p = 0,ratio_to_coarser(2) - 1
               do j_p = 0,ratio_to_coarser(1)
                  if (.not.(indicator1(i_f,j+j_p,k+k_p).eq.1)) then
                     U_fine1(i_f,j+j_p,k+k_p) = 0.d0
                     do k_wgt = -1,1
                        do j_wgt = -1,2
                           U_fine1(i_f,j+j_p,k+k_p) =
     &                          U_fine1(i_f,j+j_p,k+k_p) +
     &                          cubic_wgt1(j_wgt,j_p)*
     &                          quad_wgt2( k_wgt,k_p)*
     &                          U_intrp(j_wgt,k_wgt)
                        enddo
                     enddo
                  endif
               enddo
            enddo

c     Perform quadratic interpolation in the y-direction and cubic
c     interpolation in the z-direction.

            do k_wgt = -1,2
               do j_wgt = -1,1
                  U_intrp(j_wgt,k_wgt)=U_crse2(i_c,j_c+j_wgt,k_c+k_wgt)
               enddo
            enddo

            do k_p = 0,ratio_to_coarser(2)
               do j_p = 0,ratio_to_coarser(1) - 1
                  if (.not.(indicator2(i_f,j+j_p,k+k_p).eq.1)) then
                     U_fine2(i_f,j+j_p,k+k_p) = 0.d0
                     do k_wgt = -1,2
                        do j_wgt = -1,1
                           U_fine2(i_f,j+j_p,k+k_p) =
     &                          U_fine2(i_f,j+j_p,k+k_p) +
     &                          quad_wgt1( j_wgt,j_p)*
     &                          cubic_wgt2(k_wgt,k_p)*
     &                          U_intrp(j_wgt,k_wgt)
                        enddo
                     enddo
                  endif
               enddo
            enddo

         enddo
         enddo
         enddo

      elseif ( (location_index .eq. 2) .or.
     &         (location_index .eq. 3) ) then
c
c     Pre-compute the quadratic and cubic interpolation weights.
c
         R0 = dble(ratio_to_coarser(0))
         do i_p = 0,ratio_to_coarser(0) - 1
            t0 = dble(i_p) + 0.5d0
            quad_wgt0(-1,i_p) =
     &           eighth*(4.d0*t0*t0-8.d0*t0*R0+3.d0*R0*R0)/(R0*R0)
            quad_wgt0( 0,i_p) =
     &           fourth*(-4.d0*t0*t0+4.d0*t0*R0+3.d0*R0*R0)/(R0*R0)
            quad_wgt0(+1,i_p) =
     &           -eighth*(-4.d0*t0*t0+R0*R0)/(R0*R0)
         enddo
         do i_p = 0,ratio_to_coarser(0)
            t0 = dble(i_p)
            cubic_wgt0(-1,i_p) =
     &           -sixth*t0*(t0*t0-3.d0*t0*R0+2.d0*R0*R0)/(R0*R0*R0)
            cubic_wgt0( 0,i_p) =
     &           half*(t0*t0*t0-2.d0*t0*t0*R0-t0*R0*R0+2.d0*R0*R0*R0)/
     &           (R0*R0*R0)
            cubic_wgt0(+1,i_p) =
     &           half*t0*(-t0*t0+t0*R0+2.d0*R0*R0)/(R0*R0*R0)
            cubic_wgt0(+2,i_p) =
     &           -sixth*t0*(-t0*t0+R0*R0)/(R0*R0*R0)
         enddo

         R2 = dble(ratio_to_coarser(2))
         do k_p = 0,ratio_to_coarser(2) - 1
            t2 = dble(k_p) + 0.5d0
            quad_wgt2(-1,k_p) =
     &           eighth*(4.d0*t2*t2-8.d0*t2*R2+3.d0*R2*R2)/(R2*R2)
            quad_wgt2( 0,k_p) =
     &           fourth*(-4.d0*t2*t2+4.d0*t2*R2+3.d0*R2*R2)/(R2*R2)
            quad_wgt2(+1,k_p) =
     &           -eighth*(-4.d0*t2*t2+R2*R2)/(R2*R2)
         enddo
         do k_p = 0,ratio_to_coarser(2)
            t2 = dble(k_p)
            cubic_wgt2(-1,k_p) =
     &           -sixth*t2*(t2*t2-3.d0*t2*R2+2.d0*R2*R2)/(R2*R2*R2)
            cubic_wgt2( 0,k_p) =
     &           half*(t2*t2*t2-2.d0*t2*t2*R2-t2*R2*R2+2.d0*R2*R2*R2)/
     &           (R2*R2*R2)
            cubic_wgt2(+1,k_p) =
     &           half*t2*(-t2*t2+t2*R2+2.d0*R2*R2)/(R2*R2*R2)
            cubic_wgt2(+2,k_p) =
     &           -sixth*t2*(-t2*t2+R2*R2)/(R2*R2*R2)
         enddo
c
c     Set the values along the upper/lower y side of the patch.
c
         do k = kbeg,kend,ratio_to_coarser(2)
            if (k .lt. 0) then
               k_c = (k+1)/ratio_to_coarser(2)-1
            else
               k_c = k/ratio_to_coarser(2)
            endif

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

c     Perform quadratic interpolation in the x- and z-directions.

            do k_wgt = -1,1
               do i_wgt = -1,1
                  U_intrp(i_wgt,k_wgt)=U_crse1(i_c+i_wgt,j_c,k_c+k_wgt)
               enddo
            enddo

            do k_p = 0,ratio_to_coarser(2) - 1
               do i_p = 0,ratio_to_coarser(0) - 1
                  if (.not.(indicator1(i+i_p,j_f,k+k_p).eq.1)) then
                     U_fine1(i+i_p,j_f,k+k_p) = 0.d0
                     do k_wgt = -1,1
                        do i_wgt = -1,1
                           U_fine1(i+i_p,j_f,k+k_p) =
     &                          U_fine1(i+i_p,j_f,k+k_p) +
     &                          quad_wgt0(i_wgt,i_p)*
     &                          quad_wgt2(k_wgt,k_p)*
     &                          U_intrp(i_wgt,k_wgt)
                        enddo
                     enddo
                  endif
               enddo
            enddo

c     Shift j_c and j_f for interpolating the tangential components
c     along the coarse-fine interface.

            if (location_index .eq. 2) then
               j_c = j_c
               j_f = j_f
            else
               j_c = j_c-1
               j_f = j_f-1
            endif

c     Perform cubic interpolation in the x-direction and quadratic
c     interpolation in the z-direction.

            do k_wgt = -1,1
               do i_wgt = -1,2
                  U_intrp(i_wgt,k_wgt)=U_crse0(i_c+i_wgt,j_c,k_c+k_wgt)
               enddo
            enddo

            do k_p = 0,ratio_to_coarser(2) - 1
               do i_p = 0,ratio_to_coarser(0)
                  if (.not.(indicator0(i+i_p,j_f,k+k_p).eq.1)) then
                     U_fine0(i+i_p,j_f,k+k_p) = 0.d0
                     do k_wgt = -1,1
                        do i_wgt = -1,2
                           U_fine0(i+i_p,j_f,k+k_p) =
     &                          U_fine0(i+i_p,j_f,k+k_p) +
     &                          cubic_wgt0(i_wgt,i_p)*
     &                          quad_wgt2( k_wgt,k_p)*
     &                          U_intrp(i_wgt,k_wgt)
                        enddo
                     enddo
                  endif
               enddo
            enddo

c     Perform quadratic interpolation in the x-direction and cubic
c     interpolation in the z-direction.

            do k_wgt = -1,2
               do i_wgt = -1,1
                  U_intrp(i_wgt,k_wgt)=U_crse2(i_c+i_wgt,j_c,k_c+k_wgt)
               enddo
            enddo

            do k_p = 0,ratio_to_coarser(2)
               do i_p = 0,ratio_to_coarser(0) - 1
                  if (.not.(indicator2(i+i_p,j_f,k+k_p).eq.1)) then
                     U_fine2(i+i_p,j_f,k+k_p) = 0.d0
                     do k_wgt = -1,2
                        do i_wgt = -1,1
                           U_fine2(i+i_p,j_f,k+k_p) =
     &                          U_fine2(i+i_p,j_f,k+k_p) +
     &                          quad_wgt0( i_wgt,i_p)*
     &                          cubic_wgt2(k_wgt,k_p)*
     &                          U_intrp(i_wgt,k_wgt)
                        enddo
                     enddo
                  endif
               enddo
            enddo

         enddo
         enddo
         enddo

      elseif ( (location_index .eq. 4) .or.
     &         (location_index .eq. 5) ) then
c
c     Pre-compute the quadratic and cubic interpolation weights.
c
         R0 = dble(ratio_to_coarser(0))
         do i_p = 0,ratio_to_coarser(0) - 1
            t0 = dble(i_p) + 0.5d0
            quad_wgt0(-1,i_p) =
     &           eighth*(4.d0*t0*t0-8.d0*t0*R0+3.d0*R0*R0)/(R0*R0)
            quad_wgt0( 0,i_p) =
     &           fourth*(-4.d0*t0*t0+4.d0*t0*R0+3.d0*R0*R0)/(R0*R0)
            quad_wgt0(+1,i_p) =
     &           -eighth*(-4.d0*t0*t0+R0*R0)/(R0*R0)
         enddo
         do i_p = 0,ratio_to_coarser(0)
            t0 = dble(i_p)
            cubic_wgt0(-1,i_p) =
     &           -sixth*t0*(t0*t0-3.d0*t0*R0+2.d0*R0*R0)/(R0*R0*R0)
            cubic_wgt0( 0,i_p) =
     &           half*(t0*t0*t0-2.d0*t0*t0*R0-t0*R0*R0+2.d0*R0*R0*R0)/
     &           (R0*R0*R0)
            cubic_wgt0(+1,i_p) =
     &           half*t0*(-t0*t0+t0*R0+2.d0*R0*R0)/(R0*R0*R0)
            cubic_wgt0(+2,i_p) =
     &           -sixth*t0*(-t0*t0+R0*R0)/(R0*R0*R0)
         enddo

         R1 = dble(ratio_to_coarser(1))
         do j_p = 0,ratio_to_coarser(1) - 1
            t1 = dble(j_p) + 0.5d0
            quad_wgt1(-1,j_p) =
     &           eighth*(4.d0*t1*t1-8.d0*t1*R1+3.d0*R1*R1)/(R1*R1)
            quad_wgt1( 0,j_p) =
     &           fourth*(-4.d0*t1*t1+4.d0*t1*R1+3.d0*R1*R1)/(R1*R1)
            quad_wgt1(+1,j_p) =
     &           -eighth*(-4.d0*t1*t1+R1*R1)/(R1*R1)
         enddo
         do j_p = 0,ratio_to_coarser(1)
            t1 = dble(j_p)
            cubic_wgt1(-1,j_p) =
     &           -sixth*t1*(t1*t1-3.d0*t1*R1+2.d0*R1*R1)/(R1*R1*R1)
            cubic_wgt1( 0,j_p) =
     &           half*(t1*t1*t1-2.d0*t1*t1*R1-t1*R1*R1+2.d0*R1*R1*R1)/
     &           (R1*R1*R1)
            cubic_wgt1(+1,j_p) =
     &           half*t1*(-t1*t1+t1*R1+2.d0*R1*R1)/(R1*R1*R1)
            cubic_wgt1(+2,j_p) =
     &           -sixth*t1*(-t1*t1+R1*R1)/(R1*R1*R1)
         enddo
c
c     Set the values along the upper/lower z side of the patch.
c
         do k = blowerf(2),bupperf(2)
            if (k .lt. 0) then
               k_c = (k+1)/ratio_to_coarser(2)-1
            else
               k_c = k/ratio_to_coarser(2)
            endif

         do j = jbeg,jend,ratio_to_coarser(1)
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

c     Shift k_c and k_f for interpolating the normal component along the
c     coarse-fine interface.

            if (location_index .eq. 4) then
               k_c = k_c
               k_f = k
            else
               k_c = k_c+1
               k_f = k  +1
            endif

c     Perform quadratic interpolation in the x- and y-directions.

            do j_wgt = -1,1
               do i_wgt = -1,1
                  U_intrp(i_wgt,j_wgt)=U_crse2(i_c+i_wgt,j_c+j_wgt,k_c)
               enddo
            enddo

            do j_p = 0,ratio_to_coarser(1) - 1
               do i_p = 0,ratio_to_coarser(0) - 1
                  if (.not.(indicator2(i+i_p,j+j_p,k_f).eq.1)) then
                     U_fine2(i+i_p,j+j_p,k_f) = 0.d0
                     do j_wgt = -1,1
                        do i_wgt = -1,1
                           U_fine2(i+i_p,j+j_p,k_f) =
     &                          U_fine2(i+i_p,j+j_p,k_f) +
     &                          quad_wgt0(i_wgt,i_p)*
     &                          quad_wgt1(j_wgt,j_p)*
     &                          U_intrp(i_wgt,j_wgt)
                        enddo
                     enddo
                  endif
               enddo
            enddo

c     Shift k_c and k_f for interpolating the tangential components
c     along the coarse-fine interface.

            if (location_index .eq. 4) then
               k_c = k_c
               k_f = k_f
            else
               k_c = k_c-1
               k_f = k_f-1
            endif

c     Perform cubic interpolation in the x-direction and quadratic
c     interpolation in the y-direction.

            do j_wgt = -1,1
               do i_wgt = -1,2
                  U_intrp(i_wgt,j_wgt)=U_crse0(i_c+i_wgt,j_c+j_wgt,k_c)
               enddo
            enddo

            do j_p = 0,ratio_to_coarser(1) - 1
               do i_p = 0,ratio_to_coarser(0)
                  if (.not.(indicator0(i+i_p,j+j_p,k_f).eq.1)) then
                     U_fine0(i+i_p,j+j_p,k_f) = 0.d0
                     do j_wgt = -1,1
                        do i_wgt = -1,2
                           U_fine0(i+i_p,j+j_p,k_f) =
     &                          U_fine0(i+i_p,j+j_p,k_f) +
     &                          cubic_wgt0(i_wgt,i_p)*
     &                          quad_wgt1( j_wgt,j_p)*
     &                          U_intrp(i_wgt,j_wgt)
                        enddo
                     enddo
                  endif
               enddo
            enddo

c     Perform quadratic interpolation in the x-direction and cubic
c     interpolation in the y-direction.

            do j_wgt = -1,2
               do i_wgt = -1,1
                  U_intrp(i_wgt,j_wgt)=U_crse1(i_c+i_wgt,j_c+j_wgt,k_c)
               enddo
            enddo

            do j_p = 0,ratio_to_coarser(1)
               do i_p = 0,ratio_to_coarser(0) - 1
                  if (.not.(indicator1(i+i_p,j+j_p,k_f).eq.1)) then
                     U_fine1(i+i_p,j+j_p,k_f) = 0.d0
                     do j_wgt = -1,2
                        do i_wgt = -1,1
                           U_fine1(i+i_p,j+j_p,k_f) =
     &                          U_fine1(i+i_p,j+j_p,k_f) +
     &                          quad_wgt0( i_wgt,i_p)*
     &                          cubic_wgt1(j_wgt,j_p)*
     &                          U_intrp(i_wgt,j_wgt)
                        enddo
                     enddo
                  endif
               enddo
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
c     Compute a side-centered quadratic interpolation in the normal
c     direction at coarse-fine interfaces.
c
c     NOTE: This routine imposes a maximum ghost cell width of 8.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine scquadnormalinterpolation3d(
     &     U0,U1,U2,U_gcw,
     &     W0,W1,W2,W_gcw,
     &     indicator0,indicator1,indicator2,indicator_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
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
      INTEGER ilower2,iupper2
      INTEGER U_gcw,W_gcw,indicator_gcw

      INTEGER location_index,ratio_to_coarser(0:NDIM-1)

      INTEGER blower(0:NDIM-1), bupper(0:NDIM-1)

      REAL W0(
     &     SIDE3d0(ilower,iupper,W_gcw)
     &     )
      REAL W1(
     &     SIDE3d1(ilower,iupper,W_gcw)
     &     )
      REAL W2(
     &     SIDE3d2(ilower,iupper,W_gcw)
     &     )

      INTEGER indicator0(
     &     SIDE3d0(ilower,iupper,indicator_gcw)
     &     )
      INTEGER indicator1(
     &     SIDE3d1(ilower,iupper,indicator_gcw)
     &     )
      INTEGER indicator2(
     &     SIDE3d2(ilower,iupper,indicator_gcw)
     &     )
c
c     Input/Output.
c
      REAL U0(
     &     SIDE3d0(ilower,iupper,U_gcw)
     &     )
      REAL U1(
     &     SIDE3d1(ilower,iupper,U_gcw)
     &     )
      REAL U2(
     &     SIDE3d2(ilower,iupper,U_gcw)
     &     )
c
c     Local variables.
c
      INTEGER i,ibeg,iend,igcw,i_p
      INTEGER j,jbeg,jend,jgcw,j_p
      INTEGER k,kbeg,kend,kgcw,k_p
      REAL R0,R1,R2,t0,t1,t2
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

      kbeg = max(blower(2),ilower2)
      kend = min(bupper(2),iupper2)
      kgcw = abs(bupper(2)-blower(2))+1

      if ( (location_index .eq. 0) .or.
     &     (location_index .eq. 1) ) then
c
c     Set the values along the upper/lower x side of the patch.
c
         R0 = dble(ratio_to_coarser(0))
         do i_p = 0,igcw-1
            t0 = dble(i_p)+1.d0
            quad_wgt(-1,i_p) = -t0*(-t0+R0)/(1.d0+R0)
            quad_wgt( 0,i_p) = (R0*t0+R0-t0*t0-t0)/R0
            quad_wgt(+1,i_p) = t0*(t0+1.d0)/R0/(1.d0+R0)
         enddo

         if (location_index .eq. 0) then
            do k = kbeg,kend
               do j = jbeg,jend
                  do i = 0,igcw-1
                     if (.not.(indicator0(bupper(0)-i  ,j,k).eq.1)) then
                        U0(bupper(0)-i  ,j,k) =
     &                       quad_wgt(+1,i)*W0(bupper(0)-i  ,j,k)+
     &                       quad_wgt( 0,i)*W0(bupper(0)+1  ,j,k)+
     &                       quad_wgt(-1,i)*W0(bupper(0)+2  ,j,k)
                     endif
                  enddo
               enddo
            enddo
         else
            do k = kbeg,kend
               do j = jbeg,jend
                  do i = 0,igcw-1
                     if (.not.(indicator0(blower(0)+i+1,j,k).eq.1)) then
                        U0(blower(0)+i+1,j,k) =
     &                       quad_wgt(+1,i)*W0(blower(0)+i+1,j,k)+
     &                       quad_wgt( 0,i)*W0(blower(0)-1+1,j,k)+
     &                       quad_wgt(-1,i)*W0(blower(0)-2+1,j,k)
                     endif
                  enddo
               enddo
            enddo
         endif

         R0 = dble(ratio_to_coarser(0))
         do i_p = 0,igcw-1
            t0 = dble(i_p)+0.5d0
            quad_wgt(-1,i_p) =
     &           0.5d0*(-2.d0*R0*t0-R0+4.d0*t0*t0+2.d0*t0)/(R0+3.d0)
            quad_wgt( 0,i_p) =
     &           -0.5d0*(-2.d0*R0*t0-3.d0*R0+4.d0*t0*t0+6.d0*t0)/
     &           (R0+1.d0)
            quad_wgt(+1,i_p) =
     &           (4.d0*t0*t0+3.d0+8.d0*t0)/(R0+3.d0)/(R0+1.d0)
         enddo

         if (location_index .eq. 0) then
            do k = kbeg,kend
               do j = jbeg,jend+1
                  do i = 0,igcw-1
                     if (.not.(indicator1(bupper(0)-i,j,k).eq.1)) then
                        U1(bupper(0)-i,j,k) =
     &                       quad_wgt(+1,i)*W1(bupper(0)-i,j,k)+
     &                       quad_wgt( 0,i)*W1(bupper(0)+1,j,k)+
     &                       quad_wgt(-1,i)*W1(bupper(0)+2,j,k)
                     endif
                  enddo
               enddo
            enddo
            do k = kbeg,kend+1
               do j = jbeg,jend
                  do i = 0,igcw-1
                     if (.not.(indicator2(bupper(0)-i,j,k).eq.1)) then
                        U2(bupper(0)-i,j,k) =
     &                       quad_wgt(+1,i)*W2(bupper(0)-i,j,k)+
     &                       quad_wgt( 0,i)*W2(bupper(0)+1,j,k)+
     &                       quad_wgt(-1,i)*W2(bupper(0)+2,j,k)
                     endif
                  enddo
               enddo
            enddo
         else
            do k = kbeg,kend
               do j = jbeg,jend+1
                  do i = 0,igcw-1
                     if (.not.(indicator1(blower(0)+i,j,k).eq.1)) then
                        U1(blower(0)+i,j,k) =
     &                       quad_wgt(+1,i)*W1(blower(0)+i,j,k)+
     &                       quad_wgt( 0,i)*W1(blower(0)-1,j,k)+
     &                       quad_wgt(-1,i)*W1(blower(0)-2,j,k)
                     endif
                  enddo
               enddo
            enddo
            do k = kbeg,kend+1
               do j = jbeg,jend
                  do i = 0,igcw-1
                     if (.not.(indicator2(blower(0)+i,j,k).eq.1)) then
                        U2(blower(0)+i,j,k) =
     &                       quad_wgt(+1,i)*W2(blower(0)+i,j,k)+
     &                       quad_wgt( 0,i)*W2(blower(0)-1,j,k)+
     &                       quad_wgt(-1,i)*W2(blower(0)-2,j,k)
                     endif
                  enddo
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
            t1 = dble(j_p)+1.d0
            quad_wgt(-1,j_p) = -t1*(-t1+R1)/(1.d0+R1)
            quad_wgt( 0,j_p) = (R1*t1+R1-t1*t1-t1)/R1
            quad_wgt(+1,j_p) = t1*(t1+1.d0)/R1/(1.d0+R1)
         enddo

         if (location_index .eq. 2) then
            do k = kbeg,kend
               do j = 0,jgcw-1
                  do i = ibeg,iend
                     if (.not.(indicator1(i,bupper(1)-j  ,k).eq.1)) then
                        U1(i,bupper(1)-j  ,k) =
     &                       quad_wgt(+1,j)*W1(i,bupper(1)-j  ,k)+
     &                       quad_wgt( 0,j)*W1(i,bupper(1)+1  ,k)+
     &                       quad_wgt(-1,j)*W1(i,bupper(1)+2  ,k)
                     endif
                  enddo
               enddo
            enddo
         else
            do k = kbeg,kend
               do j = 0,jgcw-1
                  do i = ibeg,iend
                     if (.not.(indicator1(i,bupper(1)+j+1,k).eq.1)) then
                        U1(i,bupper(1)+j+1,k) =
     &                       quad_wgt(+1,j)*W1(i,bupper(1)+j+1,k)+
     &                       quad_wgt( 0,j)*W1(i,bupper(1)-1+1,k)+
     &                       quad_wgt(-1,j)*W1(i,bupper(1)-2+1,k)
                     endif
                  enddo
               enddo
            enddo
         endif

         R1 = dble(ratio_to_coarser(1))
         do j_p = 0,jgcw-1
            t1 = dble(j_p)+0.5d0
            quad_wgt(-1,j_p) =
     &           0.5d0*(-2.d0*R1*t1-R1+4.d0*t1*t1+2.d0*t1)/(R1+3.d0)
            quad_wgt( 0,j_p) =
     &           -0.5d0*(-2.d0*R1*t1-3.d0*R1+4.d0*t1*t1+6.d0*t1)/
     &           (R1+1.d0)
            quad_wgt(+1,j_p) =
     &           (4.d0*t1*t1+3.d0+8.d0*t1)/(R1+3.d0)/(R1+1.d0)
         enddo

         if (location_index .eq. 2) then
            do k = kbeg,kend
               do j = 0,jgcw-1
                  do i = ibeg,iend+1
                     if (.not.(indicator0(i,bupper(1)-j,k).eq.1)) then
                        U0(i,bupper(1)-j,k) =
     &                       quad_wgt(+1,j)*W0(i,bupper(1)-j,k)+
     &                       quad_wgt( 0,j)*W0(i,bupper(1)+1,k)+
     &                       quad_wgt(-1,j)*W0(i,bupper(1)+2,k)
                     endif
                  enddo
               enddo
            enddo
            do k = kbeg,kend+1
               do j = 0,jgcw-1
                  do i = ibeg,iend
                     if (.not.(indicator2(i,bupper(1)-j,k).eq.1)) then
                        U2(i,bupper(1)-j,k) =
     &                       quad_wgt(+1,j)*W2(i,bupper(1)-j,k)+
     &                       quad_wgt( 0,j)*W2(i,bupper(1)+1,k)+
     &                       quad_wgt(-1,j)*W2(i,bupper(1)+2,k)
                     endif
                  enddo
               enddo
            enddo
         else
            do k = kbeg,kend
               do j = 0,jgcw-1
                  do i = ibeg,iend+1
                     if (.not.(indicator0(i,blower(1)+j,k).eq.1)) then
                        U0(i,blower(1)+j,k) =
     &                       quad_wgt(+1,j)*W0(i,blower(1)+j,k)+
     &                       quad_wgt( 0,j)*W0(i,blower(1)-1,k)+
     &                       quad_wgt(-1,j)*W0(i,blower(1)-2,k)
                     endif
                  enddo
               enddo
            enddo
            do k = kbeg,kend+1
               do j = 0,jgcw-1
                  do i = ibeg,iend
                     if (.not.(indicator2(i,blower(1)+j,k).eq.1)) then
                        U2(i,blower(1)+j,k) =
     &                       quad_wgt(+1,j)*W2(i,blower(1)+j,k)+
     &                       quad_wgt( 0,j)*W2(i,blower(1)-1,k)+
     &                       quad_wgt(-1,j)*W2(i,blower(1)-2,k)
                     endif
                  enddo
               enddo
            enddo
         endif

      elseif ( (location_index .eq. 4) .or.
     &         (location_index .eq. 5) ) then
c
c     Set the values along the upper/lower z side of the patch.
c
         R2 = dble(ratio_to_coarser(2))
         do k_p = 0,kgcw-1
            t2 = dble(k_p)+1.d0
            quad_wgt(-1,k_p) = -t2*(-t2+R2)/(1.d0+R2)
            quad_wgt( 0,k_p) = (R2*t2+R2-t2*t2-t2)/R2
            quad_wgt(+1,k_p) = t2*(t2+1.d0)/R2/(1.d0+R2)
         enddo

         if (location_index .eq. 4) then
            do k = 0,kgcw-1
               do j = jbeg,jend
                  do i = ibeg,iend
                     if (.not.(indicator2(i,j,bupper(2)-k  ).eq.1)) then
                        U2(i,j,bupper(2)-k  ) =
     &                       quad_wgt(+1,k)*W2(i,j,bupper(2)-k  )+
     &                       quad_wgt( 0,k)*W2(i,j,bupper(2)+1  )+
     &                       quad_wgt(-1,k)*W2(i,j,bupper(2)+2  )
                     endif
                  enddo
               enddo
            enddo
         else
            do k = 0,kgcw-1
               do j = jbeg,jend
                  do i = ibeg,iend
                     if (.not.(indicator2(i,j,bupper(2)+k+1).eq.1)) then
                        U2(i,j,bupper(2)+k+1) =
     &                       quad_wgt(+1,k)*W2(i,j,bupper(2)+k+1)+
     &                       quad_wgt( 0,k)*W2(i,j,bupper(2)-1+1)+
     &                       quad_wgt(-1,k)*W2(i,j,bupper(2)-2+1)
                     endif
                  enddo
               enddo
            enddo
         endif

         R2 = dble(ratio_to_coarser(2))
         do k_p = 0,kgcw-1
            t2 = dble(k_p)+0.5d0
            quad_wgt(-1,k_p) =
     &           0.5d0*(-2.d0*R2*t2-R2+4.d0*t2*t2+2.d0*t2)/(R2+3.d0)
            quad_wgt( 0,k_p) =
     &           -0.5d0*(-2.d0*R2*t2-3.d0*R2+4.d0*t2*t2+6.d0*t2)/
     &           (R2+1.d0)
            quad_wgt(+1,k_p) =
     &           (4.d0*t2*t2+3.d0+8.d0*t2)/(R2+3.d0)/(R2+1.d0)
         enddo

         if (location_index .eq. 4) then
            do k = 0,kgcw-1
               do j = jbeg,jend
                  do i = ibeg,iend+1
                     if (.not.(indicator0(i,j,bupper(2)-k).eq.1)) then
                        U0(i,j,bupper(2)-k) =
     &                       quad_wgt(+1,k)*W0(i,j,bupper(2)-k)+
     &                       quad_wgt( 0,k)*W0(i,j,bupper(2)+1)+
     &                       quad_wgt(-1,k)*W0(i,j,bupper(2)+2)
                     endif
                  enddo
               enddo
            enddo
            do k = 0,kgcw-1
               do j = jbeg,jend+1
                  do i = ibeg,iend
                     if (.not.(indicator1(i,j,bupper(2)-k).eq.1)) then
                        U1(i,j,bupper(2)-k) =
     &                       quad_wgt(+1,k)*W1(i,j,bupper(2)-k)+
     &                       quad_wgt( 0,k)*W1(i,j,bupper(2)+1)+
     &                       quad_wgt(-1,k)*W1(i,j,bupper(2)+2)
                     endif
                  enddo
               enddo
            enddo
         else
            do k = 0,kgcw-1
               do j = jbeg,jend
                  do i = ibeg,iend+1
                     if (.not.(indicator0(i,j,blower(2)+k).eq.1)) then
                        U0(i,j,blower(2)+k) =
     &                       quad_wgt(+1,k)*W0(i,j,blower(2)+k)+
     &                       quad_wgt( 0,k)*W0(i,j,blower(2)-1)+
     &                       quad_wgt(-1,k)*W0(i,j,blower(2)-2)
                     endif
                  enddo
               enddo
            enddo
            do k = 0,kgcw-1
               do j = jbeg,jend+1
                  do i = ibeg,iend
                     if (.not.(indicator1(i,j,blower(2)+k).eq.1)) then
                        U1(i,j,blower(2)+k) =
     &                       quad_wgt(+1,k)*W1(i,j,blower(2)+k)+
     &                       quad_wgt( 0,k)*W1(i,j,blower(2)-1)+
     &                       quad_wgt(-1,k)*W1(i,j,blower(2)-2)
                     endif
                  enddo
               enddo
            enddo
         endif

      endif
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
