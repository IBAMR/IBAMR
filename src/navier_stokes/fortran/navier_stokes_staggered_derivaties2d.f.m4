define(NDIM,2)dnl
define(REAL,`double precision')dnl
define(INTEGER,`integer')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim2d.i)dnl
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the skew-symmetric advective derivative of the given
c     velocity field.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine navier_stokes_staggered_skew_symmetric_derivative2d(
     &     dx,
     &     ifirst0,ilast0,
     &     ifirst1,ilast1,
     &     n_U_gc0,n_U_gc1,
     &     U0,U1,
     &     n_N_gc0,n_N_gc1,
     &     N0,N1)
c
      implicit none
c
c     Input.
c
      INTEGER ifirst0,ilast0
      INTEGER ifirst1,ilast1

      INTEGER n_U_gc0,n_U_gc1
      INTEGER n_N_gc0,n_N_gc1

      REAL dx(0:NDIM-1)

      REAL U0(
     &     SIDE2d0VECG(ifirst,ilast,n_U_gc)
     &     )
      REAL U1(
     &     SIDE2d1VECG(ifirst,ilast,n_U_gc)
     &     )
c
c     Output.
c
      REAL N0(
     &     SIDE2d0VECG(ifirst,ilast,n_N_gc)
     &     )
      REAL N1(
     &     SIDE2d1VECG(ifirst,ilast,n_N_gc)
     &     )
c
c     Local variables.
c
      INTEGER i0,i1
      INTEGER gc0,gc1
      REAL A,A_x,A_y
      REAL D,D_x,D_y
c
c     Compute N = 0.5*(div(UU) + (U*grad)U).
c
      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0+1
            D_x = (0.25d0/dx(0))*(
     &           (u0(i0,i1)+u0(i0+1,i1))**2.d0 -
     &           (u0(i0,i1)+u0(i0-1,i1))**2.d0)

            D_y = (0.25d0/dx(1))*(
     &           (u0(i0,i1)+u0(i0,i1+1))*(u1(i0,i1+1)+u1(i0-1,i1+1)) -
     &           (u0(i0,i1)+u0(i0,i1-1))*(u1(i0,i1  )+u1(i0-1,i1  )))

            D = D_x + D_y

            A = D

            N(i0,i1) = 0.5d0*(D+A)
         enddo
      enddo

      do i1 = ilower1,iupper1+1
         do i0 = ilower0,iupper0
            D_x = (0.25d0/dx(1))*(
     &           (u1(i0,i1)+u1(i0+1,i1))*(u0(i0+1,i1)+u0(i0+1,i1-1)) -
     &           (u1(i0,i1)+u1(i0-1,i1))*(u0(i0  ,i1)+u0(i0  ,i1-1)))

            D_y = (0.25d0/dx(1))*(
     &           (u1(i0,i1)+u1(i0+1,i1))**2.d0 -
     &           (u1(i0,i1)+u1(i0-1,i1))**2.d0)

            D = D_x + D_y

            A = D

            N(i0,i1) = 0.5d0*(D+A)
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
