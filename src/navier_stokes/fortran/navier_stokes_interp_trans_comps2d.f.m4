define(NDIM,2)dnl
define(REAL,`double precision')dnl
define(INTEGER,`integer')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim2d.i)dnl
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Interpolate the transverse components of a MAC velocity field onto
c     all cell faces.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine navier_stokes_sc_interp_trans_comps2d(
     &     ifirst0,ilast0,ifirst1,ilast1,
     &     ngcU0,ngcU1,
     &     U0,U1,
     &     ngcV0,ngcV1,
     &     V0,V1)
c
      implicit none
c
c     Input.
c
      INTEGER ifirst0,ilast0,ifirst1,ilast1

      INTEGER ngcU0,ngcU1
      INTEGER ngcV0,ngcV1

      REAL dx(0:NDIM-1)

      REAL U0(SIDE2d0VECG(ifirst,ilast,ngcU))
      REAL U1(SIDE2d1VECG(ifirst,ilast,ngcU))
c
c     Input/Output.
c
      REAL V0(SIDE2d0VECG(ifirst,ilast,ngcV),0:NDIM-1)
      REAL V1(SIDE2d1VECG(ifirst,ilast,ngcV),0:NDIM-1)
c
c     Local variables.
c
      INTEGER i0,i1
c
c     Use simple averaging to compute the transverse components of the
c     velocity at each cell face.
c
      do i1 = ifirst1,ilast1
         do i0 = ifirst0,ilast0+1
            V0(i0,i1,0) = U0(i0,i1)
            V0(i0,i1,1) = 0.25d0*(U1(i0-1,i1  )+U1(i0,i1  )+
     &                            U1(i0-1,i1+1)+U1(i0,i1+1))
         enddo
      enddo

      do i1 = ifirst1,ilast1+1
         do i0 = ifirst0,ilast0
            V1(i0,i1,0) = 0.25d0*(U0(i0  ,i1-1)+U0(i0  ,i1)+
     &                            U0(i0+1,i1-1)+U0(i0+1,i1))
            V1(i0,i1,1) = U1(i0,i1)
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
