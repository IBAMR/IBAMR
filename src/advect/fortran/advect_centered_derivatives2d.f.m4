c ---------------------------------------------------------------------
c
c Copyright (c) 2014 - 2019 by the IBAMR developers
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
c     Compute the divergence form of the convection term corresponding
c     to the given velocity field and concentration.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine advect_centered_side_div_derivative2d(
     &     dx,
     &     ifirst0,ilast0,
     &     ifirst1,ilast1,
     &     n_U_gc0,n_U_gc1,
     &     U0,U1,
     &     n_C_gc0,n_C_gc1,
     &     C,
     &     n_N_gc0,n_N_gc1,
     &     N)
c
      implicit none
c
c     Input.
c
      INTEGER ifirst0,ilast0
      INTEGER ifirst1,ilast1

      INTEGER n_U_gc0,n_U_gc1
      INTEGER n_C_gc0,n_C_gc1
      INTEGER n_N_gc0,n_N_gc1

      REAL dx(0:NDIM-1)

      REAL U0(
     &     SIDE2d0VECG(ifirst,ilast,n_U_gc)
     &     )
      REAL U1(
     &     SIDE2d1VECG(ifirst,ilast,n_U_gc)
     &     )

      REAL C(
     &     CELL2dVECG(ifirst,ilast,n_C_gc)
     &     )
c
c     Output.
c
      REAL N(
     &     CELL2dVECG(ifirst,ilast,n_N_gc)
     &     )
c
c     Local variables.
c
      INTEGER i0,i1
      REAL D,D_x,D_y
c
c     Compute N = div(UC).
c
      do i1 = ifirst1,ilast1
         do i0 = ifirst0,ilast0
            D_x = (0.5d0/dx(0))*(
     &           u0(i0+1,i1)*((C(i0+1,i1)+C(i0  ,i1))) -
     &           u0(i0  ,i1)*((C(i0  ,i1)+C(i0-1,i1))))

            D_y = (0.5d0/dx(1))*(
     &           u1(i0,i1+1)*((C(i0,i1+1)+C(i0,i1  ))) -
     &           u1(i0,i1  )*((C(i0,i1  )+C(i0,i1-1))))

            D = D_x + D_y

            N(i0,i1) = D
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the advective form of the convection term corresponding to
c     the given velocity field and concentration.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine advect_centered_side_adv_derivative2d(
     &     dx,
     &     ifirst0,ilast0,
     &     ifirst1,ilast1,
     &     n_U_gc0,n_U_gc1,
     &     U0,U1,
     &     n_C_gc0,n_C_gc1,
     &     C,
     &     n_N_gc0,n_N_gc1,
     &     N)
c
      implicit none
c
c     Input.
c
      INTEGER ifirst0,ilast0
      INTEGER ifirst1,ilast1

      INTEGER n_U_gc0,n_U_gc1
      INTEGER n_C_gc0,n_C_gc1
      INTEGER n_N_gc0,n_N_gc1

      REAL dx(0:NDIM-1)

      REAL U0(
     &     SIDE2d0VECG(ifirst,ilast,n_U_gc)
     &     )
      REAL U1(
     &     SIDE2d1VECG(ifirst,ilast,n_U_gc)
     &     )

      REAL C(
     &     CELL2dVECG(ifirst,ilast,n_C_gc)
     &     )
c
c     Output.
c
      REAL N(
     &     CELL2dVECG(ifirst,ilast,n_N_gc)
     &     )
c
c     Local variables.
c
      INTEGER i0,i1
      REAL A,A_x,A_y
c
c     Compute N = (U*grad)C.
c
      do i1 = ifirst1,ilast1
         do i0 = ifirst0,ilast0
            A_x = (0.5d0/dx(0))*(
     &           u0(i0+1,i1)*((C(i0+1,i1)-C(i0  ,i1))) +
     &           u0(i0  ,i1)*((C(i0  ,i1)-C(i0-1,i1))))

            A_y = (0.5d0/dx(1))*(
     &           u1(i0,i1+1)*((C(i0,i1+1)-C(i0,i1  ))) +
     &           u1(i0,i1  )*((C(i0,i1  )-C(i0,i1-1))))

            A = A_x + A_y

            N(i0,i1) = A
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the skew-symmetric form of the convection term
c     corresponding to the given velocity field and concentration.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine advect_centered_side_skew_sym_derivative2d(
     &     dx,
     &     ifirst0,ilast0,
     &     ifirst1,ilast1,
     &     n_U_gc0,n_U_gc1,
     &     U0,U1,
     &     n_C_gc0,n_C_gc1,
     &     C,
     &     n_N_gc0,n_N_gc1,
     &     N)
c
      implicit none
c
c     Input.
c
      INTEGER ifirst0,ilast0
      INTEGER ifirst1,ilast1

      INTEGER n_U_gc0,n_U_gc1
      INTEGER n_C_gc0,n_C_gc1
      INTEGER n_N_gc0,n_N_gc1

      REAL dx(0:NDIM-1)

      REAL U0(
     &     SIDE2d0VECG(ifirst,ilast,n_U_gc)
     &     )
      REAL U1(
     &     SIDE2d1VECG(ifirst,ilast,n_U_gc)
     &     )

      REAL C(
     &     CELL2dVECG(ifirst,ilast,n_C_gc)
     &     )
c
c     Output.
c
      REAL N(
     &     CELL2dVECG(ifirst,ilast,n_N_gc)
     &     )
c
c     Local variables.
c
      INTEGER i0,i1
      REAL A,A_x,A_y
      REAL D,D_x,D_y
c
c     Compute N = 0.5*(div(UC) + (U*grad)C).
c
      do i1 = ifirst1,ilast1
         do i0 = ifirst0,ilast0
            D_x = (0.5d0/dx(0))*(
     &           u0(i0+1,i1)*((C(i0+1,i1)+C(i0  ,i1))) -
     &           u0(i0  ,i1)*((C(i0  ,i1)+C(i0-1,i1))))

            D_y = (0.5d0/dx(1))*(
     &           u1(i0,i1+1)*((C(i0,i1+1)+C(i0,i1  ))) -
     &           u1(i0,i1  )*((C(i0,i1  )+C(i0,i1-1))))

            D = D_x + D_y

            A_x = (0.5d0/dx(0))*(
     &           u0(i0+1,i1)*((C(i0+1,i1)-C(i0  ,i1))) +
     &           u0(i0  ,i1)*((C(i0  ,i1)-C(i0-1,i1))))

            A_y = (0.5d0/dx(1))*(
     &           u1(i0,i1+1)*((C(i0,i1+1)-C(i0,i1  ))) +
     &           u1(i0,i1  )*((C(i0,i1  )-C(i0,i1-1))))

            A = A_x + A_y

            N(i0,i1) = 0.5d0*(D+A)
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the divergence form of the convection term corresponding
c     to the given velocity field and concentration.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine advect_centered_face_div_derivative2d(
     &     dx,
     &     ifirst0,ilast0,
     &     ifirst1,ilast1,
     &     n_U_gc0,n_U_gc1,
     &     U0,U1,
     &     n_C_gc0,n_C_gc1,
     &     C,
     &     n_N_gc0,n_N_gc1,
     &     N)
c
      implicit none
c
c     Input.
c
      INTEGER ifirst0,ilast0
      INTEGER ifirst1,ilast1

      INTEGER n_U_gc0,n_U_gc1
      INTEGER n_C_gc0,n_C_gc1
      INTEGER n_N_gc0,n_N_gc1

      REAL dx(0:NDIM-1)

      REAL U0(
     &     FACE2d0VECG(ifirst,ilast,n_U_gc)
     &     )
      REAL U1(
     &     FACE2d1VECG(ifirst,ilast,n_U_gc)
     &     )

      REAL C(
     &     CELL2dVECG(ifirst,ilast,n_C_gc)
     &     )
c
c     Output.
c
      REAL N(
     &     CELL2dVECG(ifirst,ilast,n_N_gc)
     &     )
c
c     Local variables.
c
      INTEGER i0,i1
      REAL D,D_x,D_y
c
c     Compute N = div(UC).
c
      do i1 = ifirst1,ilast1
         do i0 = ifirst0,ilast0
            D_x = (0.5d0/dx(0))*(
     &           u0(i0+1,i1)*((C(i0+1,i1)+C(i0  ,i1))) -
     &           u0(i0  ,i1)*((C(i0  ,i1)+C(i0-1,i1))))

            D_y = (0.5d0/dx(1))*(
     &           u1(i1+1,i0)*((C(i0,i1+1)+C(i0,i1  ))) -
     &           u1(i1  ,i0)*((C(i0,i1  )+C(i0,i1-1))))

            D = D_x + D_y

            N(i0,i1) = D
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the advective form of the convection term corresponding to
c     the given velocity field and concentration.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine advect_centered_face_adv_derivative2d(
     &     dx,
     &     ifirst0,ilast0,
     &     ifirst1,ilast1,
     &     n_U_gc0,n_U_gc1,
     &     U0,U1,
     &     n_C_gc0,n_C_gc1,
     &     C,
     &     n_N_gc0,n_N_gc1,
     &     N)
c
      implicit none
c
c     Input.
c
      INTEGER ifirst0,ilast0
      INTEGER ifirst1,ilast1

      INTEGER n_U_gc0,n_U_gc1
      INTEGER n_C_gc0,n_C_gc1
      INTEGER n_N_gc0,n_N_gc1

      REAL dx(0:NDIM-1)

      REAL U0(
     &     FACE2d0VECG(ifirst,ilast,n_U_gc)
     &     )
      REAL U1(
     &     FACE2d1VECG(ifirst,ilast,n_U_gc)
     &     )

      REAL C(
     &     CELL2dVECG(ifirst,ilast,n_C_gc)
     &     )
c
c     Output.
c
      REAL N(
     &     CELL2dVECG(ifirst,ilast,n_N_gc)
     &     )
c
c     Local variables.
c
      INTEGER i0,i1
      REAL A,A_x,A_y
c
c     Compute N = (U*grad)C.
c
      do i1 = ifirst1,ilast1
         do i0 = ifirst0,ilast0
            A_x = (0.5d0/dx(0))*(
     &           u0(i0+1,i1)*((C(i0+1,i1)-C(i0  ,i1))) +
     &           u0(i0  ,i1)*((C(i0  ,i1)-C(i0-1,i1))))

            A_y = (0.5d0/dx(1))*(
     &           u1(i1+1,i0)*((C(i0,i1+1)-C(i0,i1  ))) +
     &           u1(i1  ,i0)*((C(i0,i1  )-C(i0,i1-1))))

            A = A_x + A_y

            N(i0,i1) = A
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the skew-symmetric form of the convection term
c     corresponding to the given velocity field and concentration.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine advect_centered_face_skew_sym_derivative2d(
     &     dx,
     &     ifirst0,ilast0,
     &     ifirst1,ilast1,
     &     n_U_gc0,n_U_gc1,
     &     U0,U1,
     &     n_C_gc0,n_C_gc1,
     &     C,
     &     n_N_gc0,n_N_gc1,
     &     N)
c
      implicit none
c
c     Input.
c
      INTEGER ifirst0,ilast0
      INTEGER ifirst1,ilast1

      INTEGER n_U_gc0,n_U_gc1
      INTEGER n_C_gc0,n_C_gc1
      INTEGER n_N_gc0,n_N_gc1

      REAL dx(0:NDIM-1)

      REAL U0(
     &     FACE2d0VECG(ifirst,ilast,n_U_gc)
     &     )
      REAL U1(
     &     FACE2d1VECG(ifirst,ilast,n_U_gc)
     &     )

      REAL C(
     &     CELL2dVECG(ifirst,ilast,n_C_gc)
     &     )
c
c     Output.
c
      REAL N(
     &     CELL2dVECG(ifirst,ilast,n_N_gc)
     &     )
c
c     Local variables.
c
      INTEGER i0,i1
      REAL A,A_x,A_y
      REAL D,D_x,D_y
c
c     Compute N = 0.5*(div(UC) + (U*grad)C).
c
      do i1 = ifirst1,ilast1
         do i0 = ifirst0,ilast0
            D_x = (0.5d0/dx(0))*(
     &           u0(i0+1,i1)*((C(i0+1,i1)+C(i0  ,i1))) -
     &           u0(i0  ,i1)*((C(i0  ,i1)+C(i0-1,i1))))

            D_y = (0.5d0/dx(1))*(
     &           u1(i1+1,i0)*((C(i0,i1+1)+C(i0,i1  ))) -
     &           u1(i1  ,i0)*((C(i0,i1  )+C(i0,i1-1))))

            D = D_x + D_y

            A_x = (0.5d0/dx(0))*(
     &           u0(i0+1,i1)*((C(i0+1,i1)-C(i0  ,i1))) +
     &           u0(i0  ,i1)*((C(i0  ,i1)-C(i0-1,i1))))

            A_y = (0.5d0/dx(1))*(
     &           u1(i1+1,i0)*((C(i0,i1+1)-C(i0,i1  ))) +
     &           u1(i1  ,i0)*((C(i0,i1  )-C(i0,i1-1))))

            A = A_x + A_y

            N(i0,i1) = 0.5d0*(D+A)
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
