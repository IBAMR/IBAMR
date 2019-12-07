c Filename: log_upper_convective_op2d.f.m4
c Created on 23 Oct 2019 by Aaron Barrett
c
c Copyright (c) 2002-2019, Boyce Griffith
c All rights reserved.
c
c Redistribution and use in source and binary forms, with or without
c modification, are permitted provided that the following conditions are met:
c
c    * Redistributions of source code must retain the above copyright notice,
c      this list of conditions and the following disclaimer.
c
c    * Redistributions in binary form must reproduce the above copyright
c      notice, this list of conditions and the following disclaimer in the
c      documentation and/or other materials provided with the distribution.
c
c    * Neither the name of The University of North Carolina nor the names of
c      its contributors may be used to endorse or promote products derived from
c      this software without specific prior written permission.
c
c THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
c AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
c IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
c ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
c LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
c CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
c SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
c INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
c CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
c ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
c POSSIBILITY OF SUCH DAMAGE.

define(NDIM,2)dnl
define(REAL,`double precision')dnl
define(INTEGER,`integer')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim2d.i)dnl
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       Computes cell centered Oldroyd-B Convective Operator
c
c       where u is vector valued face centered velocity
c       and tau is symmetric log tensor valued cell centered
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine log_upper_convective_op2d
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
      REAL q00, q01, q11
      REAL du_dx, dv_dx
      REAL du_dy, dv_dy
      REAL scale0_q, scale1_q
      REAL scale_ux, scale_uy
      REAL scale_vx, scale_vy

      REAL L(2,2), vecs(2,2)
      REAL vals(2)
      REAL convec_vals(2,2)
      REAL sigma(0:2)
      REAL temp(5)
      INTEGER info

      scale_vx = 1.d0/(4.d0*dx(0))
      scale_vy = 1.d0/dx(1)
      scale_ux = 1.d0/dx(0)
      scale_uy = 1.d0/(4.d0*dx(1))
      scale0_q = 1.d0/(2.d0*dx(0))
      scale1_q = 1.d0/(2.d0*dx(1))
      do i1 = ilower1, iupper1
        do i0 = ilower0, iupper0
          q00 = s_data(i0,i1,0)
          q01 = s_data(i0,i1,2)
          q11 = s_data(i0,i1,1)
          du_dx = scale_ux*(u_data_0(i0+1,i1)-u_data_0(i0,i1))
          du_dy = scale_uy*(u_data_0(i0+1,i1+1)+u_data_0(i0,i1+1)
     &            -u_data_0(i0+1,i1-1)-u_data_0(i0,i1-1))
          dv_dy = scale_vy*(u_data_1(i0,i1+1)-u_data_1(i0,i1))
          dv_dx = scale_vx*(u_data_1(i0+1,i1+1)+u_data_1(i0+1,i1)
     &            -u_data_1(i0-1,i1) - u_data_1(i0-1,i1+1))

        L(1,1) = du_dx; L(1,2) = du_dy
        L(2,1) = dv_dx; L(2,2) = dv_dy

        vecs(1,1) = q00; vecs(1,2) = q01
        vecs(2,1) = q01; vecs(2,2) = q11

        call dsyev('V','U',2,vecs,2,vals,temp,5,info)
        if (info /= 0) then
          print *, "\nERROR IN DSYEV!!\n\n"
        endif

        call log_sum(L,vecs,vals,2,convec_vals)

        sigma(0) = vecs(1,1)**2*exp(-vals(1))
     &    + vecs(1,2)**2*exp(-vals(2))
        sigma(1) = vecs(2,1)**2*exp(-vals(1))
     &    + vecs(2,2)**2*exp(-vals(2))
        sigma(2) = vecs(1,1)*vecs(2,1)*exp(-vals(1))
     &    + vecs(1,2)*vecs(2,2)*exp(-vals(2))

        r_data(i0,i1,0) =
     &    c_data(i0,i1,0)
     &    - convec_vals(1,1) - (
     &    sigma(0)*rhs_data(i0,i1,0) +
     &    sigma(2)*rhs_data(i0,i1,2))
        r_data(i0,i1,1) =
     &    c_data(i0,i1,1)
     &    - convec_vals(2,2) - (
     &    sigma(1)*rhs_data(i0,i1,1) +
     &    sigma(2)*rhs_data(i0,i1,2))
        r_data(i0,i1,2) =
     &    c_data(i0,i1,2)
     &    - convec_vals(1,2) - (
     &    sigma(0)*rhs_data(i0,i1,2) +
     &    sigma(2)*rhs_data(i0,i1,1))
        enddo
      enddo
      end subroutine

      subroutine log_sum(L,vecs,vals,d,to_ret)

      implicit none
      INTEGER d
      REAL L(d,d)
      REAL Lij, Lji
      REAL vecs(d,d)
      REAL vals(d)
      REAL to_ret(d,d)

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
     &        2.d0*Lij*vecs(k,i)*vecs(kk,i)
          enddo
        enddo
      enddo

      do k=1,d
        do kk=1,d
          do i=1,d; do j=1,d
            if (i /= j) then
              Lij = 0
              Lji = 0
              do ii=1,d
                do jj=1,d
                  Lij = Lij + vecs(jj,i)*L(jj,ii)*vecs(ii,j)
                  Lji = Lji + vecs(jj,j)*L(jj,ii)*vecs(ii,i)
                enddo
              enddo
              if(abs(exp_vals(i)-exp_vals(j))<1d-12) then
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
