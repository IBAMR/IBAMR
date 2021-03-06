c ---------------------------------------------------------------------
c
c Copyright (c) 2019 - 2019 by the IBAMR developers
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes cell centered Oldroyd-B type Convective Operator
c
c     where u is vector valued face centered velocity
c     and tau is the symmetric stress tensor cell centered
c     c_data is convective term u.grad(tau)
c     returns s_data
c     computes grad(u) using centered differencing
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine upper_convective_op2d
     &        (dx, u_data_0, u_data_1,
     &        u_gcw, s_data, s_gcw, rhs_data, rhs_gcw,
     &        c_data, c_gcw, r_data, r_gcw,
     &        ilower0, iupper0, ilower1, iupper1)
      implicit none
      INTEGER ilower0, iupper0
      INTEGER ilower1, iupper1
      REAL dx(0:1)

c
c     Velocity Data
c
      INTEGER u_gcw
      REAL u_data_0(SIDE2d0(ilower,iupper,u_gcw))
      REAL u_data_1(SIDE2d1(ilower,iupper,u_gcw))
c
c     Tensor Data
c
      INTEGER s_gcw
      REAL s_data(CELL2d(ilower,iupper,s_gcw),0:2)
c
c     RHS Data
c
      INTEGER rhs_gcw
      REAL rhs_data(CELL2d(ilower,iupper,rhs_gcw),0:2)
c
c     Return Data
c
      INTEGER r_gcw
      REAL r_data(CELL2d(ilower,iupper,r_gcw),0:2)
c
c     Convective Data
c
      INTEGER c_gcw
      REAL c_data(CELL2d(ilower,iupper,c_gcw),0:2)

      INTEGER i0, i1
      REAL du_dx, dv_dx
      REAL du_dy, dv_dy
      REAL scale_ux, scale_uy
      REAL scale_vx, scale_vy
      REAL qxx, qyy, qxy

      scale_vx = 1.d0/(4.d0*dx(0))
      scale_vy = 1.d0/dx(1)
      scale_ux = 1.d0/dx(0)
      scale_uy = 1.d0/(4.d0*dx(1))

      do i1 = ilower1, iupper1
         do i0 = ilower0, iupper0
!           2nd order approximations to derivatives
           du_dx = scale_ux*(u_data_0(i0+1,i1)-u_data_0(i0,i1))
           du_dy = scale_uy*(u_data_0(i0+1,i1+1)+u_data_0(i0,i1+1)
     &              -u_data_0(i0+1,i1-1)-u_data_0(i0,i1-1))
           dv_dy = scale_vy*(u_data_1(i0,i1+1)-u_data_1(i0,i1))
           dv_dx = scale_vx*(u_data_1(i0+1,i1+1)+u_data_1(i0+1,i1)
     &              -u_data_1(i0-1,i1) - u_data_1(i0-1,i1+1))
            qxx = s_data(i0,i1,0)
            qyy = s_data(i0,i1,1)
            qxy = s_data(i0,i1,2)

            r_data(i0,i1,0) = c_data(i0,i1,0)
     &        - 2.d0*du_dx*qxx - 2.d0*du_dy*qxy
     &        - rhs_data(i0,i1,0)
            r_data(i0,i1,1) = c_data(i0,i1,1)
     &        - 2.d0*dv_dx*qxy - 2.d0*dv_dy*qyy
     &        - rhs_data(i0,i1,1)
            r_data(i0,i1,2) = c_data(i0,i1,2)
     &        - qyy*du_dy - qxx*dv_dx
     &        - rhs_data(i0,i1,2)
         enddo
      enddo
      end subroutine
