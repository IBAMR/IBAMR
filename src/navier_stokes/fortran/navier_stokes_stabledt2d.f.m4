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
      subroutine navier_stokes_stabledt2d(
     &     dx,
     &     ifirst0,ilast0,ifirst1,ilast1,
     &     ngc0,ngc1,
     &     F,
     &     stabdt)
c
      implicit none
include(CONSTDIR/const.i)dnl
c
c     Input.
c
      INTEGER ifirst0,ilast0,ifirst1,ilast1

      INTEGER ngc0,ngc1

      REAL dx(0:NDIM-1)

      REAL F(CELL2dVECG(ifirst,ilast,ngc),0:NDIM-1)
c
c     Input/Output.
c
      REAL stabdt
c
c     Local variables.
c
      INTEGER i0,i1
c
c     Determine the unit CFL number on the patch.
c
      do i1 = ifirst1,ilast1
         do i0 = ifirst0,ilast0
            stabdt = dmin1( stabdt,dsqrt(
     &           (2.d0*dsqrt(dx(0)**2.d0+dx(1)**2.d0))/
     &           (dsqrt(F(i0,i1,0)**2.d0+F(i0,i1,1)**2.d0)+1.d-8)) )
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
