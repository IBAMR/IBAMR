c ---------------------------------------------------------------------
c
c Copyright (c) 2026 by the IBAMR developers
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
define(coarsen_index,`dnl
if ($1.lt.0) then
            $2=($1+1)/$4-1
         else
            $2=$1/$4
         endif
         $3=$2*$4
')dnl'
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Coarsen cell-centered data via linear interpolation.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine cart_cell_linear_coarsen3d(
     &     u_c,u_c_gcw,
     &     u_f,u_f_gcw,
     &     clower0,cupper0,
     &     clower1,cupper1,
     &     clower2,cupper2,
     &     flower0,fupper0,
     &     flower1,fupper1,
     &     flower2,fupper2,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
     &     ratio)
c
      implicit none
c
c     Input.
c
      INTEGER u_c_gcw,u_f_gcw
      INTEGER clower0,cupper0,clower1,cupper1,clower2,cupper2
      INTEGER flower0,fupper0,flower1,fupper1,flower2,fupper2
      INTEGER ilower0,iupper0,ilower1,iupper1,ilower2,iupper2
      INTEGER ratio(0:NDIM-1)

      REAL u_f(CELL3d(flower,fupper,u_f_gcw))
c
c     Input/Output.
c
      REAL u_c(CELL3d(clower,cupper,u_c_gcw))

c
c     Local variables.
c
      INTEGER i_c0,i_c1,i_c2,i0,i1,i2
      INTEGER j_c0a,j_c0b,j_c1a,j_c1b,j_c2a,j_c2b
      INTEGER i_f0,i_f1,i_f2
      REAL xi0,xi1,xi2
      REAL w0a,w0b,w1a,w1b,w2a,w2b,w0,w1,w2,scale
c
c     Coarsen data.
c
      scale = 1.d0/dble(ratio(0)*ratio(1)*ratio(2))
      do i_c2=ilower2,iupper2
         do i_c1=ilower1,iupper1
            do i_c0=ilower0,iupper0
               u_c(i_c0,i_c1,i_c2) = 0.d0
               do i2=i_c2*ratio(2)-1,i_c2*ratio(2)+ratio(2)
                  coarsen_index(i2,j_c2b,i_f2,ratio(2))
                  xi2 = (dble(i2-i_f2)+0.5d0)/dble(ratio(2)) - 0.5d0
                  if (clower2 .eq. cupper2) then
                     j_c2a = j_c2b
                     j_c2b = j_c2a
                     w2a = 1.d0
                     w2b = 0.d0
                  else if (xi2 .le. 0.d0) then
                     j_c2a = j_c2b-1
                     w2a = -xi2
                     w2b = 1.d0+xi2
                  else
                     j_c2a = j_c2b
                     j_c2b = j_c2b+1
                     w2a = 1.d0-xi2
                     w2b = xi2
                  endif
                  w2 = 0.d0
                  if (i_c2 .eq. j_c2a) w2 = w2 + w2a
                  if (i_c2 .eq. j_c2b) w2 = w2 + w2b
                  if (w2 .eq. 0.d0) cycle

                  do i1=i_c1*ratio(1)-1,i_c1*ratio(1)+ratio(1)
                     coarsen_index(i1,j_c1b,i_f1,ratio(1))
                     xi1 = (dble(i1-i_f1)+0.5d0)/dble(ratio(1)) - 0.5d0
                     if (clower1 .eq. cupper1) then
                        j_c1a = j_c1b
                        j_c1b = j_c1a
                        w1a = 1.d0
                        w1b = 0.d0
                     else if (xi1 .le. 0.d0) then
                        j_c1a = j_c1b-1
                        w1a = -xi1
                        w1b = 1.d0+xi1
                     else
                        j_c1a = j_c1b
                        j_c1b = j_c1b+1
                        w1a = 1.d0-xi1
                        w1b = xi1
                     endif
                     w1 = 0.d0
                     if (i_c1 .eq. j_c1a) w1 = w1 + w1a
                     if (i_c1 .eq. j_c1b) w1 = w1 + w1b
                     if (w1 .eq. 0.d0) cycle

                     do i0=i_c0*ratio(0)-1,i_c0*ratio(0)+ratio(0)
                        coarsen_index(i0,j_c0b,i_f0,ratio(0))
                        xi0 = (dble(i0-i_f0)+0.5d0)/dble(ratio(0))
     &                       - 0.5d0
                        if (clower0 .eq. cupper0) then
                           j_c0a = j_c0b
                           j_c0b = j_c0a
                           w0a = 1.d0
                           w0b = 0.d0
                        else if (xi0 .le. 0.d0) then
                           j_c0a = j_c0b-1
                           w0a = -xi0
                           w0b = 1.d0+xi0
                        else
                           j_c0a = j_c0b
                           j_c0b = j_c0b+1
                           w0a = 1.d0-xi0
                           w0b = xi0
                        endif
                        w0 = 0.d0
                        if (i_c0 .eq. j_c0a) w0 = w0 + w0a
                        if (i_c0 .eq. j_c0b) w0 = w0 + w0b
                        if (w0 .ne. 0.d0) then
                           u_c(i_c0,i_c1,i_c2) = u_c(i_c0,i_c1,i_c2) +
     &                          scale*w0*w1*w2*u_f(i0,i1,i2)
                        endif
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Coarsen cell-centered boundary data via linear interpolation.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c

      subroutine cart_cell_linear_coarsen_bdry3d(
     &     u_c,u_c_gcw,
     &     u_f,u_f_gcw,
     &     clower0,cupper0,
     &     clower1,cupper1,
     &     clower2,cupper2,
     &     flower0,fupper0,
     &     flower1,fupper1,
     &     flower2,fupper2,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
     &     ratio,
     &     touches_regular_bdry,
     &     patch_ilower0,patch_iupper0,
     &     patch_ilower1,patch_iupper1,
     &     patch_ilower2,patch_iupper2,
     &     bdry_normal_axis)
c
      implicit none
c
c     Input.
c
      INTEGER u_c_gcw,u_f_gcw
      INTEGER clower0,cupper0,clower1,cupper1,clower2,cupper2
      INTEGER flower0,fupper0,flower1,fupper1,flower2,fupper2
      INTEGER ilower0,iupper0,ilower1,iupper1,ilower2,iupper2
      INTEGER ratio(0:NDIM-1)
      INTEGER touches_regular_bdry(0:2*NDIM-1)
      INTEGER patch_ilower0,patch_iupper0,patch_ilower1,patch_iupper1
      INTEGER patch_ilower2,patch_iupper2
      INTEGER bdry_normal_axis

      REAL u_f(CELL3d(flower,fupper,u_f_gcw))
c
c     Input/Output.
c
      REAL u_c(CELL3d(clower,cupper,u_c_gcw))

c
c     Local variables.
c
      INTEGER i_c0,i_c1,i_c2,i0,i1,i2
      INTEGER j_c0a,j_c0b,j_c1a,j_c1b,j_c2a,j_c2b
      INTEGER i_f0,i_f1,i_f2
      REAL xi0,xi1,xi2
      REAL w0a,w0b,w1a,w1b,w2a,w2b,w0,w1,w2,scale
c
c     Coarsen data.
c
      scale = 1.d0/dble(ratio(0)*ratio(1)*ratio(2))
      do i_c2=ilower2,iupper2
         do i_c1=ilower1,iupper1
            do i_c0=ilower0,iupper0
               u_c(i_c0,i_c1,i_c2) = 0.d0
               do i2=i_c2*ratio(2)-1,i_c2*ratio(2)+ratio(2)
                  coarsen_index(i2,j_c2b,i_f2,ratio(2))
                  xi2 = (dble(i2-i_f2)+0.5d0)/dble(ratio(2)) - 0.5d0
                  if (clower2 .eq. cupper2) then
                     j_c2a = j_c2b
                     j_c2b = j_c2a
                     w2a = 1.d0
                     w2b = 0.d0
                  else if (touches_regular_bdry(4) .ne. 0 .and.
     &                     j_c2b .eq. patch_ilower2 .and.
     &                     xi2 .le. 0.d0) then
                     j_c2a = j_c2b
                     j_c2b = j_c2b+1
                     w2a = 1.d0-xi2
                     w2b = xi2
                  else if (touches_regular_bdry(5) .ne. 0 .and.
     &                     j_c2b .eq. patch_iupper2 .and.
     &                     xi2 .ge. 0.d0) then
                     j_c2a = j_c2b-1
                     w2a = -xi2
                     w2b = 1.d0+xi2
                  else if (xi2 .le. 0.d0) then
                     j_c2a = j_c2b-1
                     w2a = -xi2
                     w2b = 1.d0+xi2
                  else
                     j_c2a = j_c2b
                     j_c2b = j_c2b+1
                     w2a = 1.d0-xi2
                     w2b = xi2
                  endif
                  w2 = 0.d0
                  if (i_c2 .eq. j_c2a) w2 = w2 + w2a
                  if (i_c2 .eq. j_c2b) w2 = w2 + w2b
                  if (w2 .eq. 0.d0) cycle

                  do i1=i_c1*ratio(1)-1,i_c1*ratio(1)+ratio(1)
                     coarsen_index(i1,j_c1b,i_f1,ratio(1))
                     xi1 = (dble(i1-i_f1)+0.5d0)/dble(ratio(1)) - 0.5d0
                     if (clower1 .eq. cupper1) then
                        j_c1a = j_c1b
                        j_c1b = j_c1a
                        w1a = 1.d0
                        w1b = 0.d0
                     else if (touches_regular_bdry(2) .ne. 0 .and.
     &                        j_c1b .eq. patch_ilower1 .and.
     &                        xi1 .le. 0.d0) then
                        j_c1a = j_c1b
                        j_c1b = j_c1b+1
                        w1a = 1.d0-xi1
                        w1b = xi1
                     else if (touches_regular_bdry(3) .ne. 0 .and.
     &                        j_c1b .eq. patch_iupper1 .and.
     &                        xi1 .ge. 0.d0) then
                        j_c1a = j_c1b-1
                        w1a = -xi1
                        w1b = 1.d0+xi1
                     else if (xi1 .le. 0.d0) then
                        j_c1a = j_c1b-1
                        w1a = -xi1
                        w1b = 1.d0+xi1
                     else
                        j_c1a = j_c1b
                        j_c1b = j_c1b+1
                        w1a = 1.d0-xi1
                        w1b = xi1
                     endif
                     w1 = 0.d0
                     if (i_c1 .eq. j_c1a) w1 = w1 + w1a
                     if (i_c1 .eq. j_c1b) w1 = w1 + w1b
                     if (w1 .eq. 0.d0) cycle

                     do i0=i_c0*ratio(0)-1,i_c0*ratio(0)+ratio(0)
                        coarsen_index(i0,j_c0b,i_f0,ratio(0))
                        xi0 = (dble(i0-i_f0)+0.5d0)/dble(ratio(0))
     &                       - 0.5d0
                        if (clower0 .eq. cupper0) then
                           j_c0a = j_c0b
                           j_c0b = j_c0a
                           w0a = 1.d0
                           w0b = 0.d0
                        else if (touches_regular_bdry(0) .ne. 0 .and.
     &                           j_c0b .eq. patch_ilower0 .and.
     &                           xi0 .le. 0.d0) then
                           j_c0a = j_c0b
                           j_c0b = j_c0b+1
                           w0a = 1.d0-xi0
                           w0b = xi0
                        else if (touches_regular_bdry(1) .ne. 0 .and.
     &                           j_c0b .eq. patch_iupper0 .and.
     &                           xi0 .ge. 0.d0) then
                           j_c0a = j_c0b-1
                           w0a = -xi0
                           w0b = 1.d0+xi0
                        else if (xi0 .le. 0.d0) then
                           j_c0a = j_c0b-1
                           w0a = -xi0
                           w0b = 1.d0+xi0
                        else
                           j_c0a = j_c0b
                           j_c0b = j_c0b+1
                           w0a = 1.d0-xi0
                           w0b = xi0
                        endif
                        w0 = 0.d0
                        if (i_c0 .eq. j_c0a) w0 = w0 + w0a
                        if (i_c0 .eq. j_c0b) w0 = w0 + w0b
                        if (w0 .ne. 0.d0) then
                           u_c(i_c0,i_c1,i_c2) = u_c(i_c0,i_c1,i_c2)
     &                          + scale*w0*w1*w2*u_f(i0,i1,i2)
                        endif
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
c
      return
      end
