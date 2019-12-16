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
c     for the log-evolution equation
c
c     where u is vector valued face centered velocity
c     and tau is symmetric log tensor valued cell centered
c     c_data is convective terms u.grad(tau)
c     computes grad(u) using centered differences
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine log_upper_convective_op3d
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
      REAL qxx, qyy, qzz
      REAL qyz, qxz, qxy

      REAL L(3,3), vecs(3,3)
      REAL vals(3)
      REAL convec_vals(3,3)
      REAL sigma(0:5)
      REAL temp(15)

      INTEGER info

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

            qxx = s_data(i0,i1,i2,0)
            qyy = s_data(i0,i1,i2,1)
            qzz = s_data(i0,i1,i2,2)
            qyz = s_data(i0,i1,i2,3)
            qxz = s_data(i0,i1,i2,4)
            qxy = s_data(i0,i1,i2,5)

            L(1,1) = du_dx; L(1,2) = du_dy; L(1,3) = du_dz
            L(2,1) = dv_dx; L(2,2) = dv_dy; L(2,3) = dv_dz
            L(3,1) = dw_dx; L(3,2) = dw_dy; L(3,3) = dw_dz

            vecs(1,1) = qxx; vecs(1,2) = qxy; vecs(1,3) = qxz
            vecs(2,1) = qxy; vecs(2,2) = qyy; vecs(2,3) = qyz
            vecs(3,1) = qxz; vecs(3,2) = qyz; vecs(3,3) = qzz

            call dsyev('V','U',3,vecs,3,vals,temp,15,info)
            if (info /= 0) then
              print *, "ERROR IN DSYEV!!"
            endif
            call log_sum(L,vecs,vals,3,convec_vals)

            sigma(0) = vecs(1,1)**2*exp(-vals(1))
     &          +vecs(1,2)**2*exp(-vals(2))
     &          +vecs(1,3)**2*exp(-vals(3))
            sigma(1) = vecs(2,1)**2*exp(-vals(1))
     &          +vecs(2,2)**2*exp(-vals(2))
     &          +vecs(2,3)**2*exp(-vals(3))
            sigma(2) = vecs(3,1)**2*exp(-vals(1))
     &          +vecs(3,2)**2*exp(-vals(2))
     &          +vecs(3,3)**2*exp(-vals(3))
            sigma(3) = vecs(2,1)*vecs(3,1)*exp(-vals(1))
     &          +vecs(2,2)*vecs(3,2)*exp(-vals(2))
     &          +vecs(2,3)*vecs(3,3)*exp(-vals(3))
            sigma(4) = vecs(1,1)*vecs(3,1)*exp(-vals(1))
     &          +vecs(1,2)*vecs(3,2)*exp(-vals(2))
     &          +vecs(1,3)*vecs(3,3)*exp(-vals(3))
            sigma(5) = vecs(1,1)*vecs(2,1)*exp(-vals(1))
     &          +vecs(1,2)*vecs(2,2)*exp(-vals(2))
     &          +vecs(1,3)*vecs(2,3)*exp(-vals(3))

            r_data(i0,i1,i2,0) = 
     &        c_data(i0,i1,i2,0)
     &        - convec_vals(1,1) - (
     &        sigma(0)*rhs_data(i0,i1,i2,0) +
     &        sigma(4)*rhs_data(i0,i1,i2,4) +
     &        sigma(5)*rhs_data(i0,i1,i2,5))

            r_data(i0,i1,i2,1) =
     &        c_data(i0,i1,i2,1)
     &        - convec_vals(2,2) - (
     &        sigma(1)*rhs_data(i0,i1,i2,1) +
     &        sigma(3)*rhs_data(i0,i1,i2,3) +
     &        sigma(5)*rhs_data(i0,i1,i2,5))

            r_data(i0,i1,i2,2) =
     &        c_data(i0,i1,i2,2)
     &        - convec_vals(3,3) - (
     &        sigma(2)*rhs_data(i0,i1,i2,2) +
     &        sigma(3)*rhs_data(i0,i1,i2,3) +
     &        sigma(4)*rhs_data(i0,i1,i2,4))

            r_data(i0,i1,i2,3) =
     &        c_data(i0,i1,i2,3)
     &        - convec_vals(2,3) - (
     &        sigma(1)*rhs_data(i0,i1,i2,3) +
     &        sigma(3)*rhs_data(i0,i1,i2,2) +
     &        sigma(5)*rhs_data(i0,i1,i2,4))

            r_data(i0,i1,i2,4) =
     &        c_data(i0,i1,i2,4)
     &        - convec_vals(1,3) - (
     &        sigma(0)*rhs_data(i0,i1,i2,4) +
     &        sigma(4)*rhs_data(i0,i1,i2,2) +
     &        sigma(5)*rhs_data(i0,i1,i2,3))

            r_data(i0,i1,i2,5) =
     &        c_data(i0,i1,i2,5)
     &        - convec_vals(1,2) - (
     &        sigma(0)*rhs_data(i0,i1,i2,5) +
     &        sigma(4)*rhs_data(i0,i1,i2,3) +
     &        sigma(5)*rhs_data(i0,i1,i2,1))
          enddo
        enddo
      enddo
      end subroutine

      subroutine log_sum(L,vecs,vals,d,to_ret)

      implicit none
      INTEGER d
      REAL L(d,d)
      REAL vecs(d,d)
      REAL vals(d)
      REAL to_ret(d,d)
      REAL Lij, Lji

      REAL exp_vals(d)
      INTEGER i, j, k, kk, ii, jj

      do i=1,d
        do j=1,d
          to_ret(i,j) = 0.d0
        enddo
        exp_vals(i) = exp(vals(i))
      enddo

      do k=1,d
        do kk=1,d
          do i=1,d
            Lij = 0
            do ii=1,d
              do jj=1,d
                Lij = Lij + vecs(jj,i)*L(jj,ii)*vecs(ii,i)
              enddo
            enddo
            to_ret(k,kk) = to_ret(k,kk) +
     &        2*Lij*vecs(k,i)*vecs(kk,i)
          enddo
        enddo
      enddo
      do k=1,d
        do kk=1,d
          do i=1,d; do j=1,d
            if (i /= j) then
              Lij = 0.d0; Lji = 0.d0
              do ii=1,d
                do jj=1,d
                  Lij = Lij + vecs(jj,i)*L(jj,ii)*vecs(ii,j)
                  Lji = Lji + vecs(jj,j)*L(jj,ii)*vecs(ii,i)
                enddo
              enddo
              if(abs(vals(k)-vals(kk))<1d-10) then
                  to_ret(k,kk) = to_ret(k,kk) +
     &              (Lji+Lij)*vecs(k,i)*vecs(kk,j)
              else
                  to_ret(k,kk) = to_ret(k,kk) +
     &              (vals(i)-vals(j))/(exp_vals(i)-exp_vals(j))
     &              *(Lij*exp_vals(j)+Lji*exp_vals(i))
     &              *vecs(k,i)*vecs(kk,j)
              endif
            endif
          enddo; enddo
        enddo
      enddo
      endsubroutine
