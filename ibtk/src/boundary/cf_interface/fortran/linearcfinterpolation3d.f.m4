c ---------------------------------------------------------------------
c
c Copyright (c) 2014 - 2019 by the IBAMR developers
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
c     Compute a cell-centered linear interpolation at coarse-fine
c     interfaces.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine cclinearnormalinterpolation3d(
     &     U,U_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
     &     location_index,ratio,
     &     blower,bupper)
c
      implicit none
c
c     Input.
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2
      INTEGER U_gcw

      INTEGER location_index,ratio

      INTEGER blower(0:NDIM-1), bupper(0:NDIM-1)
c
c     Input/Output.
c
      REAL U(CELL3d(ilower,iupper,U_gcw))
c
c     Local variables.
c
      INTEGER i,ibeg,iend,i_p,i_q
      INTEGER j,jbeg,jend,j_p,j_q
      INTEGER k,kbeg,kend,k_p,k_q
      REAL R
c
c     Set the values along the appropriate side.
c
      ibeg = max(blower(0),ilower0)
      iend = min(bupper(0),iupper0)

      jbeg = max(blower(1),ilower1)
      jend = min(bupper(1),iupper1)

      kbeg = max(blower(2),ilower2)
      kend = min(bupper(2),iupper2)

      R = dble(ratio)

      if ( (location_index .eq. 0) .or.
     &     (location_index .eq. 1) ) then
c
c     Set the values along the upper/lower x side of the patch.
c
         if (location_index .eq. 0) then
            do j = jbeg,jend,ratio
               do j_p = j,j+ratio-1
                  j_q = j+(ratio-1+j-j_p)
                  do k = kbeg,kend,ratio
                     do k_p = k,k+ratio-1
                        k_q = k+(ratio-1+k-k_p)
                        U(ilower0-1,j_p,k_p) = (
     &                       2.d0*U(ilower0-1,j_p,k_p)
     &                       +  R*U(ilower0  ,j_p,k_p)
     &                       -    U(ilower0  ,j_q,k_q))/(R+1.d0)
                     enddo
                  enddo
               enddo
            enddo
         else
            do j = jbeg,jend,ratio
               do j_p = j,j+ratio-1
                  j_q = j+(ratio-1+j-j_p)
                  do k = kbeg,kend,ratio
                     do k_p = k,k+ratio-1
                        k_q = k+(ratio-1+k-k_p)
                        U(iupper0+1,j_p,k_p) = (
     &                       2.d0*U(iupper0+1,j_p,k_p)
     &                       +  R*U(iupper0  ,j_p,k_p)
     &                       -    U(iupper0  ,j_q,k_q))/(R+1.d0)
                     enddo
                  enddo
               enddo
            enddo
         endif

      elseif ( (location_index .eq. 2) .or.
     &         (location_index .eq. 3) ) then
c
c     Set the values along the upper/lower y side of the patch.
c
         if (location_index .eq. 2) then
            do i = ibeg,iend,ratio
               do i_p = i,i+ratio-1
                  i_q = i+(ratio-1+i-i_p)
                  do k = kbeg,kend,ratio
                     do k_p = k,k+ratio-1
                        k_q = k+(ratio-1+k-k_p)
                        U(i_p,ilower1-1,k_p) = (
     &                       2.d0*U(i_p,ilower1-1,k_p)
     &                       +  R*U(i_p,ilower1  ,k_p)
     &                       -    U(i_q,ilower1  ,k_q))/(R+1.d0)
                     enddo
                  enddo
               enddo
            enddo
         else
            do i = ibeg,iend,ratio
               do i_p = i,i+ratio-1
                  i_q = i+(ratio-1+i-i_p)
                  do k = kbeg,kend,ratio
                     do k_p = k,k+ratio-1
                        k_q = k+(ratio-1+k-k_p)
                        U(i_p,iupper1+1,k_p) = (
     &                       2.d0*U(i_p,iupper1+1,k_p)
     &                       +  R*U(i_p,iupper1  ,k_p)
     &                       -    U(i_q,iupper1  ,k_q))/(R+1.d0)
                     enddo
                  enddo
               enddo
            enddo
         endif

      elseif ( (location_index .eq. 4) .or.
     &         (location_index .eq. 5) ) then
c
c     Set the values along the upper/lower z side of the patch.
c
         if (location_index .eq. 4) then
            do i = ibeg,iend,ratio
               do i_p = i,i+ratio-1
                  i_q = i+(ratio-1+i-i_p)
                  do j = jbeg,jend,ratio
                     do j_p = j,j+ratio-1
                        j_q = j+(ratio-1+j-j_p)
                        U(i_p,j_p,ilower2-1) = (
     &                       2.d0*U(i_p,j_p,ilower2-1)
     &                       +  R*U(i_p,j_p,ilower2  )
     &                       -    U(i_q,j_q,ilower2  ))/(R+1.d0)
                     enddo
                  enddo
               enddo
            enddo
         else
            do i = ibeg,iend,ratio
               do i_p = i,i+ratio-1
                  i_q = i+(ratio-1+i-i_p)
                  do j = jbeg,jend,ratio
                     do j_p = j,j+ratio-1
                        j_q = j+(ratio-1+j-j_p)
                        U(i_p,j_p,iupper2+1) = (
     &                       2.d0*U(i_p,j_p,iupper2+1)
     &                       +  R*U(i_p,j_p,iupper2  )
     &                       -    U(i_q,j_q,iupper2  ))/(R+1.d0)
                     enddo
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
