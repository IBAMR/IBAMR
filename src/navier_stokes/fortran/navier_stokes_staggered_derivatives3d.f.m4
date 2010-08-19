c
c     Copyright (c) 2002-2010 Boyce Griffith
c
c     Permission is hereby granted, free of charge, to any person
c     obtaining a copy of this software and associated documentation
c     files (the "Software"), to deal in the Software without
c     restriction, including without limitation the rights to use, copy,
c     modify, merge, publish, distribute, sublicense, and/or sell copies
c     of the Software, and to permit persons to whom the Software is
c     furnished to do so, subject to the following conditions:
c
c     The above copyright notice and this permission notice shall be
c     included in all copies or substantial portions of the Software.
c
c     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
c     EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
c     MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
c     NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
c     HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
c     WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
c     OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
c     DEALINGS IN THE SOFTWARE.
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
      subroutine navier_stokes_staggered_divergence_derivative3d(
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

c
c     Compute N = div(UU) + (U*grad)U.
c
      print *,'error: 3D case is not implemented'
      call abort
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the advection form of the convection term corresponding to
c     the given velocity field.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine navier_stokes_staggered_advection_derivative3d(
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

c
c     Compute N = (U*grad)U.
c
      print *,'error: 3D case is not implemented'
      call abort
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

c
c     Compute N = 0.5*(div(UU) + (U*grad)U).
c
      print *,'error: 3D case is not implemented'
      call abort
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
