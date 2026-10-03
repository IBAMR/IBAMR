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
include(TOP_SRCDIR/ibtk/src/refine_ops/fortran/cart_cell_linear_helpers.f.m4)dnl
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Restrict cell-centered data with the cell-volume adjoint of linear
c     refinement.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine cart_cell_linear_coarsen2d(
     &     u_c,u_c_gcw,
     &     u_f,u_f_gcw,
     &     clower0,cupper0,
     &     clower1,cupper1,
     &     flower0,fupper0,
     &     flower1,fupper1,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ratio,
     &     touches_regular_bdry,
     &     patch_ilower0,patch_iupper0,
     &     patch_ilower1,patch_iupper1)
c
      implicit none
c
c     Input.
c
      INTEGER u_c_gcw,u_f_gcw
      INTEGER clower0,cupper0,clower1,cupper1
      INTEGER flower0,fupper0,flower1,fupper1
      INTEGER ilower0,iupper0,ilower1,iupper1
      INTEGER ratio(0:NDIM-1)
      INTEGER touches_regular_bdry(0:2*NDIM-1)
      INTEGER patch_ilower0,patch_iupper0,patch_ilower1,patch_iupper1

      REAL u_f(CELL2d(flower,fupper,u_f_gcw))
c
c     Input/Output.
c
      REAL u_c(CELL2d(clower,cupper,u_c_gcw))

c
c     Local variables.
c
      INTEGER i_c0,i_c1,i0,i1
      INTEGER j_c0,j_c1,idx0(0:1),idx1(0:1)
      INTEGER i_f0,i_f1
      INTEGER lower(0:NDIM-1),upper(0:NDIM-1)
      REAL xi0,xi1
      REAL w0(0:1),w1(0:1),wt0,wt1,scale
c
c     Coarsen data.
c
      scale = 1.d0/dble(ratio(0)*ratio(1))
      do i_c1=ilower1,iupper1
         do i_c0=ilower0,iupper0
            u_c(i_c0,i_c1) = 0.d0
            lower(0) = i_c0*ratio(0)-ratio(0)/2
            upper(0) = (i_c0+1)*ratio(0)-1+ratio(0)/2
            if (touches_regular_bdry(0) .ne. 0) then
               lower(0) = max(lower(0),
     &                          patch_ilower0*ratio(0))
            endif
            if (touches_regular_bdry(1) .ne. 0) then
               upper(0) = min(upper(0),
     &                          (patch_iupper0+1)*ratio(0)-1)
            endif
            lower(1) = i_c1*ratio(1)-ratio(1)/2
            upper(1) = (i_c1+1)*ratio(1)-1+ratio(1)/2
            if (touches_regular_bdry(2) .ne. 0) then
               lower(1) = max(lower(1),
     &                          patch_ilower1*ratio(1))
            endif
            if (touches_regular_bdry(3) .ne. 0) then
               upper(1) = min(upper(1),
     &                          (patch_iupper1+1)*ratio(1)-1)
            endif
            do i1=lower(1),upper(1)
               coarsen_index(i1,j_c1,i_f1,ratio(1))
               xi1 = (dble(i1-i_f1)+0.5d0)/dble(ratio(1)) - 0.5d0
               linear_stencil(j_c1,xi1,patch_ilower1,patch_iupper1,clower1,cupper1,touches_regular_bdry(2),touches_regular_bdry(3),1)dnl
               wt1 = 0.d0
               if (i_c1 .eq. idx1(0)) wt1 = wt1 + w1(0)
               if (i_c1 .eq. idx1(1)) wt1 = wt1 + w1(1)
               if (wt1 .eq. 0.d0) cycle

               do i0=lower(0),upper(0)
                  coarsen_index(i0,j_c0,i_f0,ratio(0))
                  xi0 = (dble(i0-i_f0)+0.5d0)/dble(ratio(0)) - 0.5d0
                  linear_stencil(j_c0,xi0,patch_ilower0,patch_iupper0,clower0,cupper0,touches_regular_bdry(0),touches_regular_bdry(1),0)dnl
                  wt0 = 0.d0
                  if (i_c0 .eq. idx0(0)) wt0 = wt0 + w0(0)
                  if (i_c0 .eq. idx0(1)) wt0 = wt0 + w0(1)
                  if (wt0 .ne. 0.d0) then
                     u_c(i_c0,i_c1) = u_c(i_c0,i_c1) +
     &                    scale*wt0*wt1*u_f(i0,i1)
                  endif
               enddo
            enddo
         enddo
      enddo
c
      return
      end
