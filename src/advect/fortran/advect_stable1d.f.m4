dnl Process this file with m4 to produce FORTRAN source code
define(NDIM,1)dnl
define(REAL,`double precision')dnl
define(INTEGER,`integer')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim1d.i)dnl
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     Determine the timestep size corresponding to a unit CFL number.
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
      subroutine advect_stabledt1d(
     &     dx,
     &     ifirst0,ilast0,
     &     ngc0,
     &     u0,
     &     stabdt)
c     
      implicit none
include(TOP_SRCDIR/src/fortran/const.i)dnl     
c     
c     Input.
c     
      INTEGER ifirst0,ilast0
      
      INTEGER ngc0
      
      REAL dx(0:NDIM-1)
      
      REAL u0(FACE1dVECG(ifirst,ilast,ngc))
c     
c     Input/Output.
c     
      REAL stabdt
c     
c     Local variables.
c     
      INTEGER i0,d
      REAL maxspeed(0:NDIM-1)
c     
c     Determine the unit CFL number on the patch.
c     
      do d = 0,NDIM-1
         maxspeed(d) = 1.d-12   ! avoid division by zero
      enddo
      
      do i0 = ifirst0,ilast0+1
         maxspeed(0) = dmax1(maxspeed(0), dabs(u0(i0)))
      enddo
      
      stabdt = dx(0)/maxspeed(0)
c     
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
