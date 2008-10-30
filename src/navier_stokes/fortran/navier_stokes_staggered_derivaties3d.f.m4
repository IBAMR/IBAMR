define(NDIM,3)dnl
define(REAL,`double precision')dnl
define(INTEGER,`integer')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim3d.i)dnl
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the skew-symmetric advective derivative of the given
c     velocity field.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine navier_stokes_staggered_skew_symmetric_derivative3d(
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
      INTEGER i0,i1
      INTEGER gc0,gc1
      REAL A,A_x,A_y
      REAL D,D_x,D_y
c
c     Compute N = 0.5*(div(UU) + (U*grad)U).
c

c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
