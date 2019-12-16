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

define(NDIM,3)dnl
define(REAL,`double precision')dnl
define(INTEGER,`integer')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim3d.i)dnl
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes cell centered Oldroyd-B type Convective Operator
c     for the square root evolution equation
c
c     where u is vector valued face centered velocity
c     and tau is symmetric sqrt tensor valued cell centered
c     c_data is u.grad(tau)
c     computes grad(u) using centered differences
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine sqrt_upper_convective_op3d
     &        (dx, u_data_0, u_data_1, u_data_2,
     &        u_gcw, s_data, s_gcw, rhs_data, rhs_gcw,
     &        c_data, c_gcw, r_data, r_gcw,
     &        ilower0, iupper0, ilower1,
     &        iupper1, ilower2, iupper2)

      implicit none
      INTEGER ilower0, iupper0
      INTEGER ilower1, iupper1
      INTEGER ilower2, iupper2
      REAL dx(0:2)
c
c    Velocity Data
c
      INTEGER u_gcw
      REAL u_data_0(SIDE3d0(ilower,iupper,u_gcw))
      REAL u_data_1(SIDE3d1(ilower,iupper,u_gcw))
      REAL u_data_2(SIDE3d2(ilower,iupper,u_gcw))
c
c    Tensor Data
c
      INTEGER s_gcw
      REAL s_data(CELL3d(ilower,iupper,s_gcw),0:5)
c
c    RHS Data
c
      INTEGER rhs_gcw
      REAL rhs_data(CELL3d(ilower,iupper,rhs_gcw),0:5)
c
c    Convec Data
c
      INTEGER c_gcw
      REAL c_data(CELL3d(ilower,iupper,c_gcw),0:5)
c
c    Return Data
c
      INTEGER r_gcw
      REAL r_data(CELL3d(ilower,iupper,r_gcw),0:5)

      INTEGER i0, i1, i2
      REAL du_dx, dv_dx, dw_dx
      REAL du_dy, dv_dy, dw_dy
      REAL du_dz, dv_dz, dw_dz
      REAL scale_ux, scale_uy, scale_uz
      REAL scale_vx, scale_vy, scale_vz
      REAL scale_wx, scale_wy, scale_wz
      REAL qxx_ij, qyy_ij, qzz_ij
      REAL qyz_ij, qxz_ij, qxy_ij
      REAL det
      REAL g12, g23, g13

      scale_ux = 1.d0/dx(0)
      scale_uy = 1.d0/(4.d0*dx(1))
      scale_uz = 1.d0/(4.d0*dx(2))
      scale_vx = 1.d0/(4.d0*dx(0))
      scale_vy = 1.d0/dx(1)
      scale_vz = 1.d0/(4.d0*dx(2))
      scale_wx = 1.d0/(4.d0*dx(0))
      scale_wy = 1.d0/(4.d0*dx(1))
      scale_wz = 1.d0/dx(2)

      do i2 = ilower2, iupper2
        do i1 = ilower1, iupper1
          do i0 = ilower0, iupper0
            du_dx = 
     &        scale_ux*(u_data_0(i0+1,i1,i2)-u_data_0(i0,i1,i2))
            du_dy = 
     &        scale_uy*(u_data_0(i0+1,i1+1,i2)+u_data_0(i0,i1+1,i2)
     &              -u_data_0(i0+1,i1-1,i2)-u_data_0(i0,i1-1,i2))
            du_dz = 
     &        scale_uz*(u_data_0(i0+1,i1,i2+1)+u_data_0(i0,i1,i2+1)
     &              -u_data_0(i0+1,i1,i2-1)-u_data_0(i0,i1,i2-1))
            dv_dy = 
     &        scale_vy*(u_data_1(i0,i1+1,i2)-u_data_1(i0,i1,i2))
            dv_dx = 
     &        scale_vx*(u_data_1(i0+1,i1+1,i2)+u_data_1(i0+1,i1,i2)
     &              -u_data_1(i0-1,i1+1,i2) - u_data_1(i0-1,i1,i2))
            dv_dz = 
     &        scale_vz*(u_data_1(i0,i1+1,i2+1)+u_data_1(i0,i1,i2+1)
     &              -u_data_1(i0,i1+1,i2-1) - u_data_1(i0,i1,i2-1))
            dw_dx = 
     &        scale_wx*(u_data_2(i0+1,i1,i2+1)+u_data_2(i0+1,i1,i2)
     &              -u_data_2(i0-1,i1,i2+1)-u_data_2(i0-1,i1,i2))
            dw_dy = 
     &        scale_wy*(u_data_2(i0,i1+1,i2+1)+u_data_2(i0,i1+1,i2)
     &              -u_data_2(i0,i1-1,i2+1)-u_data_2(i0,i1-1,i2))
            dw_dz = 
     &        scale_wz*(u_data_2(i0,i1,i2+1)-u_data_2(i0,i1,i2))

            qxx_ij = s_data(i0,i1,i2,0)
            qyy_ij = s_data(i0,i1,i2,1)
            qzz_ij = s_data(i0,i1,i2,2)
            qyz_ij = s_data(i0,i1,i2,3)
            qxz_ij = s_data(i0,i1,i2,4)
            qxy_ij = s_data(i0,i1,i2,5)

            det = qxx_ij*qyy_ij*qzz_ij+
     &       2.d0*qxy_ij*qxz_ij*qyz_ij-
     &       qxz_ij**2*qyy_ij-
     &       qyz_ij**2*qxx_ij-
     &       qxy_ij**2*qzz_ij

            call find_G_vals(qxx_ij, qyy_ij, qzz_ij,
     &                       qyz_ij, qxz_ij, qxy_ij,
     &                       du_dx, du_dy, du_dz,
     &                       dv_dx, dv_dy, dv_dz,
     &                       dw_dx, dw_dy, dw_dz,
     &                       g12, g13, g23)

            r_data(i0,i1,i2,0) = 
     &        c_data(i0,i1,i2,0)
     &        - qxx_ij*du_dx - qxy_ij*du_dy - qxz_ij*du_dz
     &        - g12*qxy_ij - g13*qxz_ij
     &        -0.5d0/det*(rhs_data(i0,i1,i2,0)*(
     &        qyy_ij*qzz_ij-qyz_ij**2)+rhs_data(i0,i1,i2,4)*(
     &        qxy_ij*qyz_ij-qxz_ij*qyy_ij)+
     &        rhs_data(i0,i1,i2,5)*(qxz_ij*qyz_ij-qxy_ij*qzz_ij))
            r_data(i0,i1,i2,1) = 
     &        c_data(i0,i1,i2,1)
     &        - qxy_ij*dv_dx - qyy_ij*dv_dy - qyz_ij*dv_dz
     &        + g12*qxy_ij - g23*qyz_ij
     &        -0.5d0/det*(rhs_data(i0,i1,i2,1)*(
     &        qxx_ij*qzz_ij-qxz_ij**2)+rhs_data(i0,i1,i2,3)*(
     &        qxy_ij*qxz_ij-qxx_ij*qyz_ij)+
     &        rhs_data(i0,i1,i2,5)*(qxz_ij*qyz_ij-qxy_ij*qzz_ij))
            r_data(i0,i1,i2,2) = 
     &        c_data(i0,i1,i2,2)
     &        - qxz_ij*dw_dx - qyz_ij*dw_dy - qzz_ij*dw_dz
     &        + g13*qxz_ij + g23*qyz_ij
     &        -0.5d0/det*(rhs_data(i0,i1,i2,2)*(
     &        qxx_ij*qyy_ij-qxy_ij**2)+rhs_data(i0,i1,i2,3)*(
     &        qxy_ij*qxz_ij-qxx_ij*qyz_ij)+
     &        rhs_data(i0,i1,i2,4)*(qxy_ij*qyz_ij-qxz_ij*qyy_ij))
            r_data(i0,i1,i2,3) = 
     &        c_data(i0,i1,i2,3)
     &        - qxy_ij*dw_dx - qyy_ij*dw_dy - qyz_ij*dw_dz
     &        + g12*qxz_ij - g23*qzz_ij
     &        -0.5d0/det*(rhs_data(i0,i1,i2,2)*(
     &        qxy_ij*qxz_ij-qxx_ij*qyz_ij)+rhs_data(i0,i1,i2,3)
     &        *(qxx_ij*qzz_ij-qxz_ij**2)+rhs_data(i0,i1,i2,4)*(
     &        qxz_ij*qyz_ij-qxy_ij*qzz_ij))
            r_data(i0,i1,i2,4) = 
     &        c_data(i0,i1,i2,4)
     &        - qxx_ij*dw_dx - qxy_ij*dw_dy - qxz_ij*dw_dz
     &        - g12*qyz_ij - g13*qzz_ij
     &        -0.5d0/det*(rhs_data(i0,i1,i2,2)*(
     &        qxy_ij*qyz_ij-qxz_ij*qyy_ij)+rhs_data(i0,i1,i2,3)
     &        *(qxz_ij*qyz_ij-qxy_ij*qzz_ij)+rhs_data(i0,i1,i2,4)
     &        *(qyy_ij*qzz_ij-qyz_ij**2))
            r_data(i0,i1,i2,5) = 
     &        c_data(i0,i1,i2,5)
     &        - qxx_ij*dv_dx - qxy_ij*dv_dy - qxz_ij*dv_dz
     &        - g12*qyy_ij - g13*qyz_ij
     &        -0.5d0/det*(rhs_data(i0,i1,i2,1)*(
     &        qxz_ij*qyz_ij-qxy_ij*qzz_ij)+rhs_data(i0,i1,i2,3)
     &        *(qxy_ij*qyz_ij-qxz_ij*qyy_ij)+rhs_data(i0,i1,i2,5)
     &        *(qyy_ij*qzz_ij-qyz_ij**2))
          enddo
        enddo
      enddo
      end subroutine

      subroutine find_G_vals(q11,q22,q33,
     &           q23, q13, q12,
     &           dudx, dudy, dudz,
     &           dvdx, dvdy, dvdz,
     &           dwdx, dwdy, dwdz,
     &           g12, g13, g23)
      implicit none
c     INPUTS
      REAL q11, q22, q33
      REAL q23, q13, q12
      REAL dudx, dudy, dudz
      REAL dvdx, dvdy, dvdz
      REAL dwdx, dwdy, dwdz
c     TO RETURN
      REAL g12, g13, g23
c     SCRATCH
      REAL A(3,3), b(3)
      INTEGER ipiv(3), info

      A(1,1) = q11 + q22; A(1,2) = q23; A(1,3) = -q13
      A(2,1) = q23; A(2,2) = q11 + q33; A(2,3) = q12
      A(3,1) = -q13; A(3,2) = q12; A(3,3) = q22+q33

      b(1) = q12*dudx-q11*dvdx+q22*dudy-q12*dvdy+q23*dudz-q13*dvdz
      b(2) = q13*dudx-q11*dwdx+q33*dudz-q13*dwdz+q23*dudy-q12*dwdy
      b(3) = q13*dvdx-q12*dwdx+q23*dvdy-q22*dwdy+q33*dvdz-q23*dwdz

      call dgesv(3, 1, A, 3, ipiv, b, 3, info)

      if (info /= 0) then
        print *, 'ERROR IN DGESV'
      endif

      g12 = b(1)
      g13 = b(2)
      g23 = b(3)

      endsubroutine
