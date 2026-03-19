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

define(NDIM,2)dnl
define(REAL,`double precision')dnl
define(INTEGER,`integer')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim2d.i)dnl
define(coarsen_index,`dnl
if ($1.lt.0) then
            $2=($1+1)/$4-1
         else
            $2=$1/$4
         endif
         $3=$2*$4
')dnl'
c
      subroutine cart_cell_linear_refine2d(
     &     u_f,u_f_gcw,
     &     flower0,fupper0,
     &     flower1,fupper1,
     &     u_c,u_c_gcw,
     &     clower0,cupper0,
     &     clower1,cupper1,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ratio,
     &     touches_regular_bdry)
c
      implicit none
c
      INTEGER u_f_gcw,u_c_gcw
      INTEGER flower0,fupper0,flower1,fupper1
      INTEGER clower0,cupper0,clower1,cupper1
      INTEGER ilower0,iupper0,ilower1,iupper1
      INTEGER ratio(0:NDIM-1)
      INTEGER touches_regular_bdry(0:2*NDIM-1)

      REAL u_f(CELL2d(flower,fupper,u_f_gcw))
      REAL u_c(CELL2d(clower,cupper,u_c_gcw))

      INTEGER i0,i1
      INTEGER i_c0,i_c1
      INTEGER i_f0,i_f1
      INTEGER idx0(0:1),idx1(0:1)
      REAL xi0,xi1
      REAL w0(0:1),w1(0:1)
c
      do i1=ilower1,iupper1
         coarsen_index(i1,i_c1,i_f1,ratio(1))
         xi1 = (dble(i1-i_f1)+0.5d0)/dble(ratio(1)) - 0.5d0
         if (clower1 .eq. cupper1) then
            idx1(0) = i_c1
            idx1(1) = i_c1
            w1(0) = 1.d0
            w1(1) = 0.d0
         else if (touches_regular_bdry(2) .ne. 0 .and.
     &            i_c1 .eq. clower1 .and. xi1 .le. 0.d0) then
            idx1(0) = i_c1
            idx1(1) = i_c1+1
            w1(0) = 1.d0-xi1
            w1(1) = xi1
         else if (touches_regular_bdry(3) .ne. 0 .and.
     &            i_c1 .eq. cupper1 .and. xi1 .ge. 0.d0) then
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
     &               i_c0 .eq. clower0 .and. xi0 .le. 0.d0) then
               idx0(0) = i_c0
               idx0(1) = i_c0+1
               w0(0) = 1.d0-xi0
               w0(1) = xi0
            else if (touches_regular_bdry(1) .ne. 0 .and.
     &               i_c0 .eq. cupper0 .and. xi0 .ge. 0.d0) then
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

            u_f(i0,i1) =
     &           w0(0)*w1(0)*u_c(idx0(0),idx1(0)) +
     &           w0(1)*w1(0)*u_c(idx0(1),idx1(0)) +
     &           w0(0)*w1(1)*u_c(idx0(0),idx1(1)) +
     &           w0(1)*w1(1)*u_c(idx0(1),idx1(1))
         enddo
      enddo
c
      return
      end
