c
c     Copyright (c) 2002-2017, Boyce Griffith
c     All rights reserved.
c
c     Redistribution and use in source and binary forms, with or without
c     modification, are permitted provided that the following conditions
c     are met:
c
c        * Redistributions of source code must retain the above
c          copyright notice, this list of conditions and the following
c          disclaimer.
c
c        * Redistributions in binary form must reproduce the above
c          copyright notice, this list of conditions and the following
c          disclaimer in the documentation and/or other materials
c          provided with the distribution.
c
c        * Neither the name of The University of North Carolina nor the
c          names of its contributors may be used to endorse or promote
c          products derived from this software without specific prior
c          written permission.
c
c     THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
c     CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
c     INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
c     MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
c     DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS
c     BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
c     EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
c     TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
c     DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
c     ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
c     TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
c     THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
c     SUCH DAMAGE.
c
define(NDIM,2)dnl
define(REAL,`double precision')dnl
define(INTEGER,`integer')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim2d.i)dnl
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes r = u.grad(q)
c
c     where u is vector valued face centered velocity
c     q is cell centered with depth d
c     returns r_data at cell centeres
c     computes grad(q) using weno + wave propagation
c     interpolation coefficients and weights must be provided
c     currently only works for interp orders 3 (k=2) and 5 (k=3)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine adv_diff_wp_convective_op2d(
     &            q_data, q_gcw,
     &            u_data_0, u_data_1, u_gcw,
     &            r_data, r_gcw, d,
     &            ilower0, ilower1,
     &            iupper0, iupper1,
     &            dx, k)

      implicit none
      INTEGER ilower0, iupper0
      INTEGER ilower1, iupper1
      INTEGER d

      INTEGER q_gcw
      REAL q_data(CELL2d(ilower,iupper,q_gcw),0:(d-1))

      REAL s_data_0(FACE2d0(ilower,iupper,0),0:1)
      REAL s_data_1(FACE2d1(ilower,iupper,0),0:1)

      INTEGER u_gcw
      REAL u_data_0(FACE2d0(ilower,iupper,u_gcw))
      REAL u_data_1(FACE2d1(ilower,iupper,u_gcw))

      INTEGER r_gcw
      REAL r_data(CELL2d(ilower,iupper,r_gcw),0:(d-1))

      REAL dx(0:1)

      INTEGER k, j
      INTEGER i0, i1


      do j=0,(d-1)
      call reconstruct_data_on_patch_2d(q_data(:,:,j), q_gcw,
     &             s_data_0, s_data_1, 0,
     &             ilower0, ilower1, iupper0, iupper1, k)

      do i1 = ilower1, iupper1
        do i0 = ilower0, iupper0
         r_data(i0,i1,j) =
     &     1.d0/dx(0)*(max(u_data_0(i0,i1),0.d0)*
     &     (s_data_0(i0,i1,1)-s_data_0(i0,i1,0))
     &     + min(u_data_0(i0+1,i1),0.d0)*
     &     (s_data_0(i0+1,i1,1)-s_data_0(i0+1,i1,0))
     &     + 0.5d0*(u_data_0(i0+1,i1)+u_data_0(i0,i1))*
     &     (s_data_0(i0+1,i1,0)-s_data_0(i0,i1,1)))

         r_data(i0,i1,j) = r_data(i0,i1,j) +
     &     1.d0/dx(1)*(max(u_data_1(i1,i0),0.d0)*
     &     (s_data_1(i1,i0,1)-s_data_1(i1,i0,0))
     &     + min(u_data_1(i1+1,i0),0.d0)*
     &     (s_data_1(i1+1,i0,1)-s_data_1(i1+1,i0,0))
     &     + 0.5d0*(u_data_1(i1+1,i0)+u_data_1(i1,i0))*
     &     (s_data_1(i1+1,i0,0)-s_data_1(i1,i0,1)))
        enddo
      enddo
      enddo
      end subroutine

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       Reconstructs data on patches using a weno scheme
c       the convex and interpolation weights must be supplied
c
c       q_data is cell centered with depth 1
c       r_data_* are face centered with depth 2
c         and return the values reconstructed from each side
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine reconstruct_data_on_patch_2d(q_data, q_gcw,
     &            r_data_0, r_data_1, r_gcw,
     &            ilower0, ilower1, iupper0, iupper1, k)

       INTEGER ilower0, iupper0
       INTEGER ilower1, iupper1

       INTEGER q_gcw
       REAL q_data(CELL2d(ilower,iupper,q_gcw))

       INTEGER r_gcw
       REAL r_data_0(FACE2d0(ilower,iupper,r_gcw),0:1)
       REAL r_data_1(FACE2d1(ilower,iupper,r_gcw),0:1)


       INTEGER i0, i1

       REAL WENO5_interp
c
c     Prevent compiler warning about unused variables.
c
      k = k

c     X DIRECTION
      do i1 = ilower1, iupper1
        do i0=ilower0,iupper0+1
           r_data_0(i0,i1,1) = WENO5_interp(q_data(i0+2:i0-2:-1,i1))
           r_data_0(i0,i1,0) = WENO5_interp(q_data(i0-3:i0+1,i1))
         enddo
       enddo

c      Y DIRECTION
       do i1 = ilower1,iupper1+1
         do i0 = ilower0,iupper0
           r_data_1(i1,i0,1) = WENO5_interp(q_data(i0,i1+2:i1-2:-1))
           r_data_1(i1,i0,0) = WENO5_interp(q_data(i0,i1-3:i1+1))
         enddo
       enddo
       end subroutine
