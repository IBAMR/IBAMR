c ---------------------------------------------------------------------
c
c Copyright (c) 2011 - 2019 by the IBAMR developers
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
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Apply the divergence- and gradient-preserving correction to values
c     refined from the next coarser level of the patch hierarchy.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine div_preserving_correction2d(
     &     u0,u1,u_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     correction_box_ilower0,correction_box_iupper0,
     &     correction_box_ilower1,correction_box_iupper1,
     &     ratio,dx_fine)
c
      implicit none
c
c     Input.
c
      INTEGER u_gcw
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER correction_box_ilower0,correction_box_iupper0
      INTEGER correction_box_ilower1,correction_box_iupper1
      INTEGER ratio(0:NDIM-1)
      REAL    dx_fine(0:NDIM-1)
c
c     Input/Output.
c
      REAL u0(SIDE2d0(ilower,iupper,u_gcw))
      REAL u1(SIDE2d1(ilower,iupper,u_gcw))
c
c     Local variables.
c
      INTEGER d,i0,i1,i,j
      REAL u(-1:1,-1:1),u_xx
      REAL v(-1:1,-1:1),v_yy
      REAL dx,dy
c
c     Apply the divergence- and curl-preserving corrections.
c
      do d = 0,NDIM-1
         if ( .not.(ratio(d).eq.2) ) then
            print *,'error: invalid refinement ratio'
            call abort
         endif
      enddo

      dx = dx_fine(0)           ! NOTE: These values are not used in the
      dy = dx_fine(1)           ! 2D code, but are in the 3D routine.

      do i1=correction_box_ilower1,correction_box_iupper1,ratio(1)
         do i0=correction_box_ilower0,correction_box_iupper0,ratio(0)
            u(-1,-1) = u0(i0  ,i1  )
            u( 1,-1) = u0(i0+2,i1  )
            u(-1, 1) = u0(i0  ,i1+1)
            u( 1, 1) = u0(i0+2,i1+1)

            v(-1,-1) = u1(i0  ,i1  )
            v( 1,-1) = u1(i0+1,i1  )
            v(-1, 1) = u1(i0  ,i1+2)
            v( 1, 1) = u1(i0+1,i1+2)

            u_xx = 0.25d0*(v(-1,-1)-v(-1, 1)-v( 1,-1)+v(1,1))
            v_yy = 0.25d0*(u(-1,-1)-u( 1,-1)-u(-1, 1)+u(1,1))

            do j = -1,1,2
               u(0,j) = 0.5d0*(u(1,j)+u(-1,j))+u_xx
            enddo

            u0(i0+1,i1  ) = u(0,-1)
            u0(i0+1,i1+1) = u(0, 1)

            do i = -1,1,2
               v(i,0) = 0.5d0*(v(i,1)+v(i,-1))+v_yy
            enddo

            u1(i0  ,i1+1) = v(-1,0)
            u1(i0+1,i1+1) = v( 1,0)
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
