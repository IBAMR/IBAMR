c
c     Copyright (c) 2002-2014, Boyce Griffith
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
define(NDIM,3)dnl
define(REAL,`double precision')dnl
define(INTEGER,`integer')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim3d.i)dnl
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the divergence form of the convection term corresponding
c     to the given velocity field.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine navier_stokes_staggered_div_derivative3d(
     &     dx,
     &     ifirst0,ilast0,
     &     ifirst1,ilast1,
     &     ifirst2,ilast2,
     &     n_U_gc0,n_U_gc1,n_U_gc2,
     &     U0,U1,U2,
     &     n_N_gc0,n_N_gc1,n_N_gc2,
     &     N0,N1,N2)
c
      implicit none
c
c     Input.
c
      INTEGER ifirst0,ilast0
      INTEGER ifirst1,ilast1
      INTEGER ifirst2,ilast2

      INTEGER n_U_gc0,n_U_gc1,n_U_gc2
      INTEGER n_N_gc0,n_N_gc1,n_N_gc2

      REAL dx(0:NDIM-1)

      REAL U0(
     &     SIDE3d0VECG(ifirst,ilast,n_U_gc)
     &     )
      REAL U1(
     &     SIDE3d1VECG(ifirst,ilast,n_U_gc)
     &     )
      REAL U2(
     &     SIDE3d2VECG(ifirst,ilast,n_U_gc)
     &     )
c
c     Output.
c
      REAL N0(
     &     SIDE3d0VECG(ifirst,ilast,n_N_gc)
     &     )
      REAL N1(
     &     SIDE3d1VECG(ifirst,ilast,n_N_gc)
     &     )
      REAL N2(
     &     SIDE3d2VECG(ifirst,ilast,n_N_gc)
     &     )
c
c     Local variables.
c
      INTEGER i0,i1,i2
      REAL D,D_x,D_y,D_z

c
c     Compute N = div(UU).
c
      do i2 = ifirst2,ilast2
         do i1 = ifirst1,ilast1
            do i0 = ifirst0,ilast0+1
               D_x = (0.25d0/dx(0))*(
     &              (u0(i0+1,i1,i2)+u0(i0,i1,i2))*
     &              (u0(i0+1,i1,i2)+u0(i0,i1,i2)) -
     &              (u0(i0,i1,i2)+u0(i0-1,i1,i2))*
     &              (u0(i0,i1,i2)+u0(i0-1,i1,i2)))

               D_y = (0.25d0/dx(1))*(
     &              (u1(i0,i1+1,i2)+u1(i0-1,i1+1,i2))*
     &              (u0(i0,i1+1,i2)+u0(i0,i1,i2)) -
     &              (u1(i0,i1,i2)+u1(i0-1,i1,i2))*
     &              (u0(i0,i1,i2)+u0(i0,i1-1,i2)))

               D_z = (0.25d0/dx(2))*(
     &              (u2(i0,i1,i2+1)+u2(i0-1,i1,i2+1))*
     &              (u0(i0,i1,i2+1)+u0(i0,i1,i2)) -
     &              (u2(i0,i1,i2)+u2(i0-1,i1,i2))*
     &              (u0(i0,i1,i2)+u0(i0,i1,i2-1)))

               D = D_x + D_y + D_z

               N0(i0,i1,i2) = D
            enddo
         enddo
      enddo

      do i2 = ifirst2,ilast2
         do i1 = ifirst1,ilast1+1
            do i0 = ifirst0,ilast0
               D_x = (0.25d0/dx(0))*(
     &              (u0(i0+1,i1,i2)+u0(i0+1,i1-1,i2))*
     &              (u1(i0+1,i1,i2)+u1(i0,i1,i2)) -
     &              (u0(i0  ,i1,i2)+u0(i0  ,i1-1,i2))*
     &              (u1(i0,i1,i2)+u1(i0-1,i1,i2)))

               D_y = (0.25d0/dx(1))*(
     &              (u1(i0,i1+1,i2)+u1(i0,i1,i2))*
     &              (u1(i0,i1+1,i2)+u1(i0,i1,i2)) -
     &              (u1(i0,i1,i2)+u1(i0,i1-1,i2))*
     &              (u1(i0,i1,i2)+u1(i0,i1-1,i2)))

               D_z = (0.25d0/dx(2))*(
     &              (u2(i0,i1,i2+1)+u2(i0,i1-1,i2+1))*
     &              (u1(i0,i1,i2+1)+u1(i0,i1,i2)) -
     &              (u2(i0,i1,i2)+u2(i0,i1-1,i2))*
     &              (u1(i0,i1,i2)+u1(i0,i1,i2-1)))

               D = D_x + D_y + D_z

               N1(i0,i1,i2) = D
            enddo
         enddo
      enddo

      do i2 = ifirst2,ilast2+1
         do i1 = ifirst1,ilast1
            do i0 = ifirst0,ilast0
               D_x = (0.25d0/dx(0))*(
     &              (u0(i0+1,i1,i2)+u0(i0+1,i1,i2-1))*
     &              (u2(i0+1,i1,i2)+u2(i0,i1,i2)) -
     &              (u0(i0  ,i1,i2)+u0(i0  ,i1,i2-1))*
     &              (u2(i0,i1,i2)+u2(i0-1,i1,i2)))

               D_y = (0.25d0/dx(1))*(
     &              (u1(i0,i1+1,i2)+u1(i0,i1+1,i2-1))*
     &              (u2(i0,i1+1,i2)+u2(i0,i1,i2)) -
     &              (u1(i0,i1,i2)+u1(i0,i1,i2-1))*
     &              (u2(i0,i1,i2)+u2(i0,i1-1,i2)))

               D_z = (0.25d0/dx(2))*(
     &              (u2(i0,i1,i2+1)+u2(i0,i1,i2))*
     &              (u2(i0,i1,i2+1)+u2(i0,i1,i2)) -
     &              (u2(i0,i1,i2)+u2(i0,i1,i2-1))*
     &              (u2(i0,i1,i2)+u2(i0,i1,i2-1)))

               D = D_x + D_y + D_z

               N2(i0,i1,i2) = D
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the advective form of the convection term corresponding to
c     the given velocity field.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine navier_stokes_staggered_adv_derivative3d(
     &     dx,
     &     ifirst0,ilast0,
     &     ifirst1,ilast1,
     &     ifirst2,ilast2,
     &     n_U_gc0,n_U_gc1,n_U_gc2,
     &     U0,U1,U2,
     &     n_N_gc0,n_N_gc1,n_N_gc2,
     &     N0,N1,N2)
c
      implicit none
c
c     Input.
c
      INTEGER ifirst0,ilast0
      INTEGER ifirst1,ilast1
      INTEGER ifirst2,ilast2

      INTEGER n_U_gc0,n_U_gc1,n_U_gc2
      INTEGER n_N_gc0,n_N_gc1,n_N_gc2

      REAL dx(0:NDIM-1)

      REAL U0(
     &     SIDE3d0VECG(ifirst,ilast,n_U_gc)
     &     )
      REAL U1(
     &     SIDE3d1VECG(ifirst,ilast,n_U_gc)
     &     )
      REAL U2(
     &     SIDE3d2VECG(ifirst,ilast,n_U_gc)
     &     )
c
c     Output.
c
      REAL N0(
     &     SIDE3d0VECG(ifirst,ilast,n_N_gc)
     &     )
      REAL N1(
     &     SIDE3d1VECG(ifirst,ilast,n_N_gc)
     &     )
      REAL N2(
     &     SIDE3d2VECG(ifirst,ilast,n_N_gc)
     &     )
c
c     Local variables.
c
      INTEGER i0,i1,i2
      REAL A,A_x,A_y,A_z
c
c     Compute N = (U*grad)U.
c
      do i2 = ifirst2,ilast2 
         do i1 = ifirst1,ilast1
            do i0 = ifirst0,ilast0+1
               A_x =  (0.25d0/dx(0))*(
     &              (u0(i0+1,i1,i2)+u0(i0,i1,i2))*
     &              (u0(i0+1,i1,i2)-u0(i0,i1,i2)) +
     &              (u0(i0,i1,i2)+u0(i0-1,i1,i2))*
     &              (u0(i0,i1,i2)-u0(i0-1,i1,i2)))

               A_y = (0.25d0/dx(1))*(
     &              (u1(i0,i1+1,i2)+u1(i0-1,i1+1,i2))*
     &              (u0(i0,i1+1,i2)-u0(i0,i1,i2)) +
     &              (u1(i0,i1,i2  )+u1(i0-1,i1,i2  ))*
     &              (u0(i0,i1,i2)-u0(i0,i1-1,i2)))

               A_z = (0.25d0/dx(2))*(
     &              (u2(i0,i1,i2+1)+u2(i0-1,i1,i2+1))*
     &              (u0(i0,i1,i2+1)-u0(i0,i1,i2)) +
     &              (u2(i0,i1,i2  )+u2(i0-1,i1,i2  ))*(
     &              u0(i0,i1,i2)-u0(i0,i1,i2-1)))

               A = A_x + A_y + A_z

               N0(i0,i1,i2) = A
            enddo
         enddo
      enddo

      do i2 = ifirst2,ilast2
         do i1 = ifirst1,ilast1+1
            do i0 = ifirst0,ilast0
               A_x = (0.25d0/dx(0))*(
     &              (u0(i0+1,i1,i2)+u0(i0+1,i1-1,i2))*
     &              (u1(i0+1,i1,i2)-u1(i0,i1,i2)) +
     &              (u0(i0  ,i1,i2)+u0(i0  ,i1-1,i2))*
     &              (u1(i0,i1,i2)-u1(i0-1,i1,i2)))

               A_y =  (0.25d0/dx(1))*(
     &              (u1(i0,i1+1,i2)+u1(i0,i1,i2))*
     &              (u1(i0,i1+1,i2)-u1(i0,i1,i2)) +
     &              (u1(i0,i1,i2)+u1(i0,i1-1,i2))*
     &              (u1(i0,i1,i2)-u1(i0,i1-1,i2)))

               A_z =  (0.25d0/dx(2))*(
     &              (u2(i0,i1,i2+1)+u2(i0,i1-1,i2+1))*
     &              (u1(i0,i1,i2+1)-u1(i0,i1,i2)) +
     &              (u2(i0,i1,i2)+u2(i0,i1-1,i2))*
     &              (u1(i0,i1,i2)-u1(i0,i1,i2-1)))

               A = A_x + A_y + A_z

               N1(i0,i1,i2) = A
            enddo
         enddo
      enddo

      do i2 = ifirst2,ilast2+1
         do i1 = ifirst1,ilast1
            do i0 = ifirst0,ilast0
               A_x = (0.25d0/dx(0))*(
     &              (u0(i0+1,i1,i2)+u0(i0+1,i1,i2-1))*
     &              (u2(i0+1,i1,i2)-u2(i0,i1,i2)) +
     &              (u0(i0  ,i1,i2)+u0(i0  ,i1,i2-1))*
     &              (u2(i0,i1,i2)-u2(i0-1,i1,i2)))

               A_y =  (0.25d0/dx(1))*(
     &              (u1(i0,i1+1,i2)+u1(i0,i1+1,i2-1))*
     &              (u2(i0,i1+1,i2)-u2(i0,i1,i2)) +
     &              (u1(i0,i1,i2)+u1(i0,i1,i2-1))*
     &              (u2(i0,i1,i2)-u2(i0,i1-1,i2)))

               A_z =  (0.25d0/dx(2))*(
     &              (u2(i0,i1,i2+1)+u2(i0,i1,i2))*
     &              (u2(i0,i1,i2+1)-u2(i0,i1,i2)) +
     &              (u2(i0,i1,i2)+u2(i0,i1,i2-1))*
     &              (u2(i0,i1,i2)-u2(i0,i1,i2-1)))

               A = A_x + A_y + A_z

               N2(i0,i1,i2) = A
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the skew-symmetric form of the convection term
c     corresponding to the given velocity field.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine navier_stokes_staggered_skew_sym_derivative3d(
     &     dx,
     &     ifirst0,ilast0,
     &     ifirst1,ilast1,
     &     ifirst2,ilast2,
     &     n_U_gc0,n_U_gc1,n_U_gc2,
     &     U0,U1,U2,
     &     n_N_gc0,n_N_gc1,n_N_gc2,
     &     N0,N1,N2)
c
      implicit none
c
c     Input.
c
      INTEGER ifirst0,ilast0
      INTEGER ifirst1,ilast1
      INTEGER ifirst2,ilast2

      INTEGER n_U_gc0,n_U_gc1,n_U_gc2
      INTEGER n_N_gc0,n_N_gc1,n_N_gc2

      REAL dx(0:NDIM-1)

      REAL U0(
     &     SIDE3d0VECG(ifirst,ilast,n_U_gc)
     &     )
      REAL U1(
     &     SIDE3d1VECG(ifirst,ilast,n_U_gc)
     &     )
      REAL U2(
     &     SIDE3d2VECG(ifirst,ilast,n_U_gc)
     &     )
c
c     Output.
c
      REAL N0(
     &     SIDE3d0VECG(ifirst,ilast,n_N_gc)
     &     )
      REAL N1(
     &     SIDE3d1VECG(ifirst,ilast,n_N_gc)
     &     )
      REAL N2(
     &     SIDE3d2VECG(ifirst,ilast,n_N_gc)
     &     )
c
c     Local variables.
c
      INTEGER i0,i1,i2
      REAL D,D_x,D_y,D_z
      REAL A,A_x,A_y,A_z

c
c     Compute N = 0.5*(div(UU) + (U*grad)U).
c
      do i2 = ifirst2,ilast2
         do i1 = ifirst1,ilast1
            do i0 = ifirst0,ilast0+1
               D_x = (0.25d0/dx(0))*(
     &              (u0(i0+1,i1,i2)+u0(i0,i1,i2))*
     &              (u0(i0+1,i1,i2)+u0(i0,i1,i2)) -
     &              (u0(i0,i1,i2)+u0(i0-1,i1,i2))*
     &              (u0(i0,i1,i2)+u0(i0-1,i1,i2)))

               D_y = (0.25d0/dx(1))*(
     &              (u1(i0,i1+1,i2)+u1(i0-1,i1+1,i2))*
     &              (u0(i0,i1+1,i2)+u0(i0,i1,i2)) -
     &              (u1(i0,i1,i2)+u1(i0-1,i1,i2))*
     &              (u0(i0,i1,i2)+u0(i0,i1-1,i2)))

               D_z = (0.25d0/dx(2))*(
     &              (u2(i0,i1,i2+1)+u2(i0-1,i1,i2+1))*
     &              (u0(i0,i1,i2+1)+u0(i0,i1,i2)) -
     &              (u2(i0,i1,i2)+u2(i0-1,i1,i2))*
     &              (u0(i0,i1,i2)+u0(i0,i1,i2-1)))

               D = D_x + D_y + D_z

               A_x =  (0.25d0/dx(0))*(
     &              (u0(i0+1,i1,i2)+u0(i0,i1,i2))*
     &              (u0(i0+1,i1,i2)-u0(i0,i1,i2)) +
     &              (u0(i0,i1,i2)+u0(i0-1,i1,i2))*
     &              (u0(i0,i1,i2)-u0(i0-1,i1,i2)))

               A_y = (0.25d0/dx(1))*(
     &              (u1(i0,i1+1,i2)+u1(i0-1,i1+1,i2))*
     &              (u0(i0,i1+1,i2)-u0(i0,i1,i2)) +
     &              (u1(i0,i1,i2  )+u1(i0-1,i1,i2  ))*
     &              (u0(i0,i1,i2)-u0(i0,i1-1,i2)))

               A_z = (0.25d0/dx(2))*(
     &              (u2(i0,i1,i2+1)+u2(i0-1,i1,i2+1))*
     &              (u0(i0,i1,i2+1)-u0(i0,i1,i2)) +
     &              (u2(i0,i1,i2  )+u2(i0-1,i1,i2  ))*(
     &              u0(i0,i1,i2)-u0(i0,i1,i2-1)))

               A = A_x + A_y + A_z

               N0(i0,i1,i2) = 0.5d0*(D+A)
            enddo
         enddo
      enddo

      do i2 = ifirst2,ilast2
         do i1 = ifirst1,ilast1+1
            do i0 = ifirst0,ilast0
               D_x = (0.25d0/dx(0))*(
     &              (u0(i0+1,i1,i2)+u0(i0+1,i1-1,i2))*
     &              (u1(i0+1,i1,i2)+u1(i0,i1,i2)) -
     &              (u0(i0  ,i1,i2)+u0(i0  ,i1-1,i2))*
     &              (u1(i0,i1,i2)+u1(i0-1,i1,i2)))

               D_y = (0.25d0/dx(1))*(
     &              (u1(i0,i1+1,i2)+u1(i0,i1,i2))*
     &              (u1(i0,i1+1,i2)+u1(i0,i1,i2)) -
     &              (u1(i0,i1,i2)+u1(i0,i1-1,i2))*
     &              (u1(i0,i1,i2)+u1(i0,i1-1,i2)))

               D_z = (0.25d0/dx(2))*(
     &              (u2(i0,i1,i2+1)+u2(i0,i1-1,i2+1))*
     &              (u1(i0,i1,i2+1)+u1(i0,i1,i2)) -
     &              (u2(i0,i1,i2)+u2(i0,i1-1,i2))*
     &              (u1(i0,i1,i2)+u1(i0,i1,i2-1)))

               D = D_x + D_y + D_z

               A_x = (0.25d0/dx(0))*(
     &              (u0(i0+1,i1,i2)+u0(i0+1,i1-1,i2))*
     &              (u1(i0+1,i1,i2)-u1(i0,i1,i2)) +
     &              (u0(i0  ,i1,i2)+u0(i0  ,i1-1,i2))*
     &              (u1(i0,i1,i2)-u1(i0-1,i1,i2)))

               A_y =  (0.25d0/dx(1))*(
     &              (u1(i0,i1+1,i2)+u1(i0,i1,i2))*
     &              (u1(i0,i1+1,i2)-u1(i0,i1,i2)) +
     &              (u1(i0,i1,i2)+u1(i0,i1-1,i2))*
     &              (u1(i0,i1,i2)-u1(i0,i1-1,i2)))

               A_z =  (0.25d0/dx(2))*(
     &              (u2(i0,i1,i2+1)+u2(i0,i1-1,i2+1))*
     &              (u1(i0,i1,i2+1)-u1(i0,i1,i2)) +
     &              (u2(i0,i1,i2)+u2(i0,i1-1,i2))*
     &              (u1(i0,i1,i2)-u1(i0,i1,i2-1)))

               A = A_x + A_y + A_z

               N1(i0,i1,i2) = 0.5d0*(D+A)
            enddo
         enddo
      enddo

      do i2 = ifirst2,ilast2+1
         do i1 = ifirst1,ilast1
            do i0 = ifirst0,ilast0
               D_x = (0.25d0/dx(0))*(
     &              (u0(i0+1,i1,i2)+u0(i0+1,i1,i2-1))*
     &              (u2(i0+1,i1,i2)+u2(i0,i1,i2)) -
     &              (u0(i0  ,i1,i2)+u0(i0  ,i1,i2-1))*
     &              (u2(i0,i1,i2)+u2(i0-1,i1,i2)))

               D_y = (0.25d0/dx(1))*(
     &              (u1(i0,i1+1,i2)+u1(i0,i1+1,i2-1))*
     &              (u2(i0,i1+1,i2)+u2(i0,i1,i2)) -
     &              (u1(i0,i1,i2)+u1(i0,i1,i2-1))*
     &              (u2(i0,i1,i2)+u2(i0,i1-1,i2)))

               D_z = (0.25d0/dx(2))*(
     &              (u2(i0,i1,i2+1)+u2(i0,i1,i2))*
     &              (u2(i0,i1,i2+1)+u2(i0,i1,i2)) -
     &              (u2(i0,i1,i2)+u2(i0,i1,i2-1))*
     &              (u2(i0,i1,i2)+u2(i0,i1,i2-1)))

               D = D_x + D_y + D_z

               A_x = (0.25d0/dx(0))*(
     &              (u0(i0+1,i1,i2)+u0(i0+1,i1,i2-1))*
     &              (u2(i0+1,i1,i2)-u2(i0,i1,i2)) +
     &              (u0(i0  ,i1,i2)+u0(i0  ,i1,i2-1))*
     &              (u2(i0,i1,i2)-u2(i0-1,i1,i2)))

               A_y =  (0.25d0/dx(1))*(
     &              (u1(i0,i1+1,i2)+u1(i0,i1+1,i2-1))*
     &              (u2(i0,i1+1,i2)-u2(i0,i1,i2)) +
     &              (u1(i0,i1,i2)+u1(i0,i1,i2-1))*
     &              (u2(i0,i1,i2)-u2(i0,i1-1,i2)))

               A_z =  (0.25d0/dx(2))*(
     &              (u2(i0,i1,i2+1)+u2(i0,i1,i2))*
     &              (u2(i0,i1,i2+1)-u2(i0,i1,i2)) +
     &              (u2(i0,i1,i2)+u2(i0,i1,i2-1))*
     &              (u2(i0,i1,i2)-u2(i0,i1,i2-1)))

               A = A_x + A_y + A_z

               N2(i0,i1,i2) = 0.5d0*(D+A)
            enddo
         enddo
      enddo
c
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
