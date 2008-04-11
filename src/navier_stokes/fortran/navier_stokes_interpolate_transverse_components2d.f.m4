define(NDIM,2)dnl
define(REAL,`double precision')dnl
define(INTEGER,`integer')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim2d.i)dnl
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Determine the timestep size corresponding to a unit CFL number.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine navier_stokes_sc_interpolate_transverse_components2d(
     &     dx,
     &     ifirst0,ilast0,ifirst1,ilast1,
     &     ngc0,ngc1,
     &     U0,U1,
     &     V0,V1,
     &     stabdt)
c
      implicit none
c
c     Input.
c
      INTEGER ifirst0,ilast0,ifirst1,ilast1

      INTEGER ngc0,ngc1

      REAL dx(0:NDIM-1)

      REAL U0(SIDE2d0VECG(ifirst,ilast,ngc))
      REAL U1(SIDE2d1VECG(ifirst,ilast,ngc))

      REAL V0(SIDE2d0VECG(ifirst,ilast,ngc),0:NDIM-1)
      REAL V1(SIDE2d1VECG(ifirst,ilast,ngc),0:NDIM-1)
c
c     Input/Output.
c
      REAL stabdt
c
c     Local variables.
c
      INTEGER i0,i1,d
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

      stabdt = dmin1((dx(0)/maxspeed(0)),(dx(1)/maxspeed(1)))
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
