define(NDIM,3)dnl
define(REAL,`double precision')dnl
define(INTEGER,`integer')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim3d.i)dnl
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Determine the timestep size corresponding to a unit CFL number.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine navier_stokes_stabledt3d(
     &     dx,
     &     ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &     ngc0,ngc1,ngc2,
     &     F,
     &     stabdt)
c
      implicit none
include(CONSTDIR/const.i)dnl
c
c     Input.
c
      INTEGER ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2

      INTEGER ngc0,ngc1,ngc2

      REAL dx(0:NDIM-1)

      REAL F(CELL3dVECG(ifirst,ilast,ngc),0:NDIM-1)
c
c     Input/Output.
c
      REAL stabdt
c
c     Local variables.
c
      INTEGER i0,i1,i2
c
c     Determine the unit CFL number on the patch.
c
      do i2 = ifirst2,ilast2
         do i1 = ifirst1,ilast1
            do i0 = ifirst0,ilast0
               stabdt = dmin1( stabdt,dsqrt(
     &              (2.d0*dsqrt(dx(0)**2.d0+dx(1)**2.d0+dx(2)**2.d0))/
     &              (dsqrt(F(i0,i1,i2,0)**2.d0+F(i0,i1,i2,1)**2.d0
     &              +F(i0,i1,i2,2)**2.d0)+1.d-8)) )
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
