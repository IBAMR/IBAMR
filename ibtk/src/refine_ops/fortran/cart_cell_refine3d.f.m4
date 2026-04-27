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
c     Refine cell-centered data via linear interpolation.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine cart_cell_linear_refine3d(
     &     u_f,u_f_gcw,
     &     flower0,fupper0,
     &     flower1,fupper1,
     &     flower2,fupper2,
     &     u_c,u_c_gcw,
     &     clower0,cupper0,
     &     clower1,cupper1,
     &     clower2,cupper2,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
     &     ratio,
     &     touches_regular_bdry)
c
      implicit none
c
c     Input.
c
      INTEGER u_f_gcw,u_c_gcw
      INTEGER flower0,fupper0,flower1,fupper1,flower2,fupper2
      INTEGER clower0,cupper0,clower1,cupper1,clower2,cupper2
      INTEGER ilower0,iupper0,ilower1,iupper1,ilower2,iupper2
      INTEGER ratio(0:NDIM-1)
      INTEGER touches_regular_bdry(0:2*NDIM-1)

      REAL u_c(CELL3d(clower,cupper,u_c_gcw))
c
c     Output.
c
      REAL u_f(CELL3d(flower,fupper,u_f_gcw))

c
c     Local variables.
c
      INTEGER i0,i1,i2
      INTEGER i_c0,i_c1,i_c2
      INTEGER i_f0,i_f1,i_f2
      INTEGER idx0(0:1),idx1(0:1),idx2(0:1)
      REAL xi0,xi1,xi2
      REAL w0(0:1),w1(0:1),w2(0:1)
c
c     Refine data.
c
      do i2=ilower2,iupper2
         coarsen_index(i2,i_c2,i_f2,ratio(2))
         xi2 = (dble(i2-i_f2)+0.5d0)/dble(ratio(2)) - 0.5d0
         if (clower2 .eq. cupper2) then
            idx2(0) = i_c2
            idx2(1) = i_c2
            w2(0) = 1.d0
            w2(1) = 0.d0
         else if (touches_regular_bdry(4) .ne. 0 .and.
     &            i_c2 .eq. clower2 .and. xi2 .le. 0.d0) then
            idx2(0) = i_c2
            idx2(1) = i_c2+1
            w2(0) = 1.d0-xi2
            w2(1) = xi2
         else if (touches_regular_bdry(5) .ne. 0 .and.
     &            i_c2 .eq. cupper2 .and. xi2 .ge. 0.d0) then
            idx2(0) = i_c2-1
            idx2(1) = i_c2
            w2(0) = -xi2
            w2(1) = 1.d0+xi2
         else if (xi2 .le. 0.d0) then
            idx2(0) = i_c2-1
            idx2(1) = i_c2
            w2(0) = -xi2
            w2(1) = 1.d0+xi2
         else
            idx2(0) = i_c2
            idx2(1) = i_c2+1
            w2(0) = 1.d0-xi2
            w2(1) = xi2
         endif

         do i1=ilower1,iupper1
            coarsen_index(i1,i_c1,i_f1,ratio(1))
            xi1 = (dble(i1-i_f1)+0.5d0)/dble(ratio(1)) - 0.5d0
            if (clower1 .eq. cupper1) then
               idx1(0) = i_c1
               idx1(1) = i_c1
               w1(0) = 1.d0
               w1(1) = 0.d0
            else if (touches_regular_bdry(2) .ne. 0 .and.
     &               i_c1 .eq. clower1 .and. xi1 .le. 0.d0) then
               idx1(0) = i_c1
               idx1(1) = i_c1+1
               w1(0) = 1.d0-xi1
               w1(1) = xi1
            else if (touches_regular_bdry(3) .ne. 0 .and.
     &               i_c1 .eq. cupper1 .and. xi1 .ge. 0.d0) then
               idx1(0) = i_c1-1
               idx1(1) = i_c1
               w1(0) = -xi1
               w1(1) = 1.d0+xi1
            else if (xi1 .le. 0.d0) then
               idx1(0) = i_c1-1
               idx1(1) = i_c1
               w1(0) = -xi1
               w1(1) = 1.d0+xi1
            else
               idx1(0) = i_c1
               idx1(1) = i_c1+1
               w1(0) = 1.d0-xi1
               w1(1) = xi1
            endif

            do i0=ilower0,iupper0
               coarsen_index(i0,i_c0,i_f0,ratio(0))
               xi0 = (dble(i0-i_f0)+0.5d0)/dble(ratio(0)) - 0.5d0
               if (clower0 .eq. cupper0) then
                  idx0(0) = i_c0
                  idx0(1) = i_c0
                  w0(0) = 1.d0
                  w0(1) = 0.d0
               else if (touches_regular_bdry(0) .ne. 0 .and.
     &                  i_c0 .eq. clower0 .and. xi0 .le. 0.d0) then
                  idx0(0) = i_c0
                  idx0(1) = i_c0+1
                  w0(0) = 1.d0-xi0
                  w0(1) = xi0
               else if (touches_regular_bdry(1) .ne. 0 .and.
     &                  i_c0 .eq. cupper0 .and. xi0 .ge. 0.d0) then
                  idx0(0) = i_c0-1
                  idx0(1) = i_c0
                  w0(0) = -xi0
                  w0(1) = 1.d0+xi0
               else if (xi0 .le. 0.d0) then
                  idx0(0) = i_c0-1
                  idx0(1) = i_c0
                  w0(0) = -xi0
                  w0(1) = 1.d0+xi0
               else
                  idx0(0) = i_c0
                  idx0(1) = i_c0+1
                  w0(0) = 1.d0-xi0
                  w0(1) = xi0
               endif

               u_f(i0,i1,i2) =
     &              w0(0)*w1(0)*w2(0)*u_c(idx0(0),idx1(0),idx2(0)) +
     &              w0(1)*w1(0)*w2(0)*u_c(idx0(1),idx1(0),idx2(0)) +
     &              w0(0)*w1(1)*w2(0)*u_c(idx0(0),idx1(1),idx2(0)) +
     &              w0(1)*w1(1)*w2(0)*u_c(idx0(1),idx1(1),idx2(0)) +
     &              w0(0)*w1(0)*w2(1)*u_c(idx0(0),idx1(0),idx2(1)) +
     &              w0(1)*w1(0)*w2(1)*u_c(idx0(1),idx1(0),idx2(1)) +
     &              w0(0)*w1(1)*w2(1)*u_c(idx0(0),idx1(1),idx2(1)) +
     &              w0(1)*w1(1)*w2(1)*u_c(idx0(1),idx1(1),idx2(1))
            enddo
         enddo
      enddo
c
      return
      end
