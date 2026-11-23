c ---------------------------------------------------------------------
c
c Copyright (c) 2011 - 2025 by the IBAMR developers
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
dnl Define some macros that help us unroll interpolation loops. These assume
dnl that we are accumulating into a 2D array V from a 3D array u which are both
dnl indexed first by depth (and then by point number s, which is specified in an
dnl outer loop). Assume d, ic0, and ic1 are available as inner loop indicies.
dnl
dnl The arguments are the lower and upper bounds in y and lower and upper bounds
dnl in x.
define(INTERPOLATE_INNER_2D,
          ` do d = 0,depth-1
               V(d, s) = 0.d0
               do ic1 = $1,$2
                  do ic0 = $3,$4
                     V(d,s) = V(d,s)
     &                    +w(0,ic0-$3)
     &                    *w(1,ic1-$1)
     &                    *u(ic0,ic1,d)
                  enddo
               enddo
            enddo')
dnl Same arguments as before, but the fifth argument is the width of the
dnl stencil (e.g., 3 for bspline 3). The first branch is a hotter code path
dnl since when we are not at a boundary the number of inner loop iterations
dnl is known. Exposing this to the compiler helps generate code which speeds
dnl up the subroutine by about 25%.
define(INTERPOLATE_2D_SPECIALIZE_FIXED_WIDTH,
`   if ($2 - $1 == ($5 - 1) .and.
     &       $4 - $3 == ($5 - 1)) then
           INTERPOLATE_INNER_2D($1, ($1 + $5 - 1),
                                $3, ($3 + $5 - 1))
         else
           INTERPOLATE_INNER_2D($1, $2,
                                $3, $4)
         endif')dnl
define(SPREAD_INNER_2D,
          ` do d = 0,depth-1
               do ic1 = $1,$2
                  do ic0 = $3,$4
                     u(ic0,ic1,d) = u(ic0,ic1,d) + (
     &                    +w(0,ic0-$3)
     &                    *w(1,ic1-$1)
     &                    *V(d,s)/(dx(0)*dx(1)))
                  enddo
               enddo
            enddo')
dnl Same arguments as before, but the seventh argument is the width of the
dnl stencil (e.g., 3 for bspline 3). The first branch is a hotter code path
dnl since when we are not at a boundary the number of inner loop iterations
dnl is known. Exposing this to the compiler helps generate code which speeds
dnl up the subroutine by about 25%.
define(SPREAD_2D_SPECIALIZE_FIXED_WIDTH,
`   if ($2 - $1 == ($5 - 1) .and.
     &       $4 - $3 == ($5 - 1)) then
           SPREAD_INNER_2D($1, ($1 + $5 - 1),
                           $3, ($3 + $5 - 1))
         else
           SPREAD_INNER_2D($1, $2,
                           $3, $4)
         endif')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim2d.i)dnl

c     this is an m4 include, not a Fortran include
include(CURRENT_SRCDIR/lagrangian_delta.f.m4)

c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Interpolate u onto V at the points listed in indices using the
c     piecewise constant delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_piecewise_constant_interp2d(
     &     dx,depth,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u,
     &     indices,g,inear,ibrack,nindices,
     &     V)
c
      implicit none
c
c     Input.
c
      INTEGER depth
      INTEGER ilower0,iupper0,ilower1,iupper1
      INTEGER nugc0,nugc1
      INTEGER nindices

      INTEGER indices(0:nindices-1)
      INTEGER inear(0:NDIM-1,0:nindices-1)
      INTEGER ibrack(0:NDIM-1,0:nindices-1)

      REAL g(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1)
      REAL u(CELL2dVECG(ilower,iupper,nugc),0:depth-1)
c
c     Input/Output.
c
      REAL V(0:depth-1,0:*)
c
c     Local variables.
c
      INTEGER ic0,ic1
      INTEGER d,l,s
c
c     Prevent compiler warning about unused variables.
c
      ibrack(0,0) = ibrack(0,0)
      dx(0) = dx(0)
      g(0,0) = g(0,0)
c
c     Use the piecewise constant delta function to interpolate u onto V.
c
      do l = 0,nindices-1
         s = indices(l)
c
c     Determine the Cartesian cell in which point s is located.
c
         ic0 = inear(0,l)
         ic1 = inear(1,l)
c
c     Interpolate u onto V.
c
         do d = 0,depth-1
            V(d,s) = u(ic0,ic1,d)
         enddo
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Spread V onto u at the points listed in indices using the
c     piecewise constant delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_piecewise_constant_spread2d(
     &     dx,depth,
     &     indices,g,inear,ibrack,nindices,
     &     V,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u)
c
      implicit none
c
c     Input.
c
      INTEGER depth
      INTEGER nindices
      INTEGER ilower0,iupper0,ilower1,iupper1
      INTEGER nugc0,nugc1

      INTEGER indices(0:nindices-1)
      INTEGER inear(0:NDIM-1,0:nindices-1)
      INTEGER ibrack(0:NDIM-1,0:nindices-1)

      REAL g(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1)
      REAL V(0:depth-1,0:*)
c
c     Input/Output: spreading adds to the existing Cartesian grid values u.
c
      REAL u(CELL2dVECG(ilower,iupper,nugc),0:depth-1)
c
c     Local variables.
c
      INTEGER ic0,ic1
      INTEGER d,l,s
c
c     Prevent compiler warning about unused variables.
c
      ibrack(0,0) = ibrack(0,0)
      g(0,0) = g(0,0)
c
c     Use the piecewise constant delta function to spread V onto u.
c
      do l = 0,nindices-1
         s = indices(l)
c
c     Determine the Cartesian cell in which point s is located.
c
         ic0 = inear(0,l)
         ic1 = inear(1,l)
c
c     Spread V onto u, skipping a cell outside the ghost box of u.
c
         if ( ic0.ge.ilower0-nugc0 .and. ic0.le.iupper0+nugc0 .and.
     &        ic1.ge.ilower1-nugc1 .and. ic1.le.iupper1+nugc1 ) then
            do d = 0,depth-1
               u(ic0,ic1,d) = u(ic0,ic1,d) + V(d,s)/(dx(0)*dx(1))
            enddo
         endif
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Interpolate u onto V at the points listed in indices using the
c     discontinuous linear delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_discontinuous_linear_interp2d(
     &     dx,depth,axis,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u,
     &     indices,g,inear,ibrack,nindices,
     &     V)
c
      implicit none
c
c     Input.
c
      INTEGER depth,axis
      INTEGER ilower0,iupper0,ilower1,iupper1
      INTEGER nugc0,nugc1
      INTEGER nindices

      INTEGER indices(0:nindices-1)
      INTEGER inear(0:NDIM-1,0:nindices-1)
      INTEGER ibrack(0:NDIM-1,0:nindices-1)

      REAL g(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1)
      REAL u(CELL2dVECG(ilower,iupper,nugc),0:depth-1)
c
c     Input/Output.
c
      REAL V(0:depth-1,0:*)
c
c     Local variables.
c
      INTEGER ilower(0:NDIM-1),iupper(0:NDIM-1)
      INTEGER ic0,ic1
      INTEGER ic_center(0:NDIM-1),ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER ic_trimmed_lower(0:NDIM-1),ic_trimmed_upper(0:NDIM-1)
      INTEGER d,l,s,nugc(0:NDIM-1)

      REAL w(0:NDIM-1,0:1)
c
c     Prevent compiler warning about unused variables.
c
      dx(0) = dx(0)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1

      iupper(0) = iupper0
      iupper(1) = iupper1

      nugc(0) = nugc0
      nugc(1) = nugc1
c
c     Use the discontinuous linear delta function to interpolate u onto V.
c
      do l = 0,nindices-1
         s = indices(l)
c
c     Determine the interpolation stencils and weights.
c
         do d = 0,NDIM-1
            ic_center(d) = inear(d,l)
            if ( d.eq.axis ) then
               if ( ibrack(d,l).lt.inear(d,l) ) then
                  ic_lower(d) = ic_center(d)-1
                  ic_upper(d) = ic_center(d)
                  w(d,0) = dble(ic_center(d)-ibrack(d,l))-g(d,l)
                  w(d,1) = 1.d0 - w(d,0)
               else
                  ic_lower(d) = ic_center(d)
                  ic_upper(d) = ic_center(d)+1
                  w(d,0) = 1.d0 + dble(ic_center(d)-ibrack(d,l))-g(d,l)
                  w(d,1) = 1.d0 - w(d,0)
               endif
            else
               w(d,0) = 1.d0
               ic_lower(d) = ic_center(d)
               ic_upper(d) = ic_center(d)
            endif

            ic_trimmed_lower(d) = max(ic_lower(d),ilower(d)-nugc(d))
            ic_trimmed_upper(d) = min(ic_upper(d),iupper(d)+nugc(d))
         enddo
c
c     Interpolate u onto V.
c
         do d = 0,depth-1
            V(d,s) = 0.d0
            do ic1 = ic_trimmed_lower(1),ic_trimmed_upper(1)
               do ic0 = ic_trimmed_lower(0),ic_trimmed_upper(0)
                  V(d,s) = V(d,s)
     &                 +w(0,ic0-ic_lower(0))
     &                 *w(1,ic1-ic_lower(1))
     &                 *u(ic0,ic1,d)
               enddo
            enddo
         enddo
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Spread V onto u at the points listed in indices using the
c     discontinuous linear delta function using standard (double) precision
c     accumulation on the Cartesian grid.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_discontinuous_linear_spread2d(
     &     dx,depth,axis,
     &     indices,g,inear,ibrack,nindices,
     &     V,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u)
c
      implicit none
c
c     Input.
c
      INTEGER depth,axis
      INTEGER nindices
      INTEGER ilower0,iupper0,ilower1,iupper1
      INTEGER nugc0,nugc1

      INTEGER indices(0:nindices-1)
      INTEGER inear(0:NDIM-1,0:nindices-1)
      INTEGER ibrack(0:NDIM-1,0:nindices-1)

      REAL g(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1)
      REAL V(0:depth-1,0:*)
c
c     Input/Output: spreading adds to the existing Cartesian grid values u.
c
      REAL u(CELL2dVECG(ilower,iupper,nugc),0:depth-1)
c
c     Local variables.
c
      INTEGER ilower(0:NDIM-1),iupper(0:NDIM-1)
      INTEGER ic0,ic1
      INTEGER ic_center(0:NDIM-1),ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER ic_trimmed_lower(0:NDIM-1),ic_trimmed_upper(0:NDIM-1)
      INTEGER d,l,s,nugc(0:NDIM-1)

      REAL w(0:NDIM-1,0:1)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1

      iupper(0) = iupper0
      iupper(1) = iupper1

      nugc(0) = nugc0
      nugc(1) = nugc1
c
c     Use the discontinuous linear delta function to interpolate u onto
c     V.
c
      do l = 0,nindices-1
         s = indices(l)
c
c     Determine the interpolation stencils and weights.
c
         do d = 0,NDIM-1
            ic_center(d) = inear(d,l)
            if ( d.eq.axis ) then
               if ( ibrack(d,l).lt.inear(d,l) ) then
                  ic_lower(d) = ic_center(d)-1
                  ic_upper(d) = ic_center(d)
                  w(d,0) = dble(ic_center(d)-ibrack(d,l))-g(d,l)
                  w(d,1) = 1.d0 - w(d,0)
               else
                  ic_lower(d) = ic_center(d)
                  ic_upper(d) = ic_center(d)+1
                  w(d,0) = 1.d0 + dble(ic_center(d)-ibrack(d,l))-g(d,l)
                  w(d,1) = 1.d0 - w(d,0)
               endif
            else
               w(d,0) = 1.d0
               ic_lower(d) = ic_center(d)
               ic_upper(d) = ic_center(d)
            endif

            ic_trimmed_lower(d) = max(ic_lower(d),ilower(d)-nugc(d))
            ic_trimmed_upper(d) = min(ic_upper(d),iupper(d)+nugc(d))
         enddo
c
c     Spread V onto u.
c
         do d = 0,depth-1
            do ic1 = ic_trimmed_lower(1),ic_trimmed_upper(1)
               do ic0 = ic_trimmed_lower(0),ic_trimmed_upper(0)
                  u(ic0,ic1,d) = u(ic0,ic1,d)+(
     &                 w(0,ic0-ic_lower(0))*
     &                 w(1,ic1-ic_lower(1))*
     &                 V(d,s)/(dx(0)*dx(1)))
               enddo
            enddo
         enddo
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Interpolate u onto V at the points listed in indices using the
c     piecewise linear delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_piecewise_linear_interp2d(
     &     dx,depth,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u,
     &     indices,g,inear,ibrack,nindices,
     &     V)
c
      implicit none
c
c     Input.
c
      INTEGER depth
      INTEGER ilower0,iupper0,ilower1,iupper1
      INTEGER nugc0,nugc1
      INTEGER nindices

      INTEGER indices(0:nindices-1)
      INTEGER inear(0:NDIM-1,0:nindices-1)
      INTEGER ibrack(0:NDIM-1,0:nindices-1)

      REAL g(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1)
      REAL u(CELL2dVECG(ilower,iupper,nugc),0:depth-1)
c
c     Input/Output.
c
      REAL V(0:depth-1,0:*)
c
c     Local variables.
c
      INTEGER ilower(0:NDIM-1),iupper(0:NDIM-1)
      INTEGER ic0,ic1
      INTEGER ic_center(0:NDIM-1),ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER ic_trimmed_lower(0:NDIM-1),ic_trimmed_upper(0:NDIM-1)
      INTEGER d,l,s,nugc(0:NDIM-1)

      REAL w(0:NDIM-1,0:1)
c
c     Prevent compiler warning about unused variables.
c
      dx(0) = dx(0)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1

      iupper(0) = iupper0
      iupper(1) = iupper1

      nugc(0) = nugc0
      nugc(1) = nugc1
c
c     Use the piecewise linear delta function to interpolate u onto V.
c
      do l = 0,nindices-1
         s = indices(l)
c
c     Determine the interpolation stencils and weights.
c
         do d = 0,NDIM-1
            ic_center(d) = inear(d,l)

            if ( ibrack(d,l).lt.inear(d,l) ) then
               ic_lower(d) = ic_center(d)-1
               ic_upper(d) = ic_center(d)
               w(d,0) = dble(ic_center(d)-ibrack(d,l))-g(d,l)
               w(d,1) = 1.d0 - w(d,0)
            else
               ic_lower(d) = ic_center(d)
               ic_upper(d) = ic_center(d)+1
               w(d,0) = 1.d0 + dble(ic_center(d)-ibrack(d,l))-g(d,l)
               w(d,1) = 1.d0 - w(d,0)
            endif

            ic_trimmed_lower(d) = max(ic_lower(d),ilower(d)-nugc(d))
            ic_trimmed_upper(d) = min(ic_upper(d),iupper(d)+nugc(d))
         enddo
c
c     Interpolate u onto V.
c
         do d = 0,depth-1
            V(d,s) = 0.d0
            do ic1 = ic_trimmed_lower(1),ic_trimmed_upper(1)
               do ic0 = ic_trimmed_lower(0),ic_trimmed_upper(0)
                  V(d,s) = V(d,s)
     &                 +w(0,ic0-ic_lower(0))
     &                 *w(1,ic1-ic_lower(1))
     &                 *u(ic0,ic1,d)
               enddo
            enddo
         enddo
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Spread V onto u at the points listed in indices using the
c     piecewise linear delta function using standard (double) precision
c     accumulation on the Cartesian grid.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_piecewise_linear_spread2d(
     &     dx,depth,
     &     indices,g,inear,ibrack,nindices,
     &     V,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u)
c
      implicit none
c
c     Input.
c
      INTEGER depth
      INTEGER nindices
      INTEGER ilower0,iupper0,ilower1,iupper1
      INTEGER nugc0,nugc1

      INTEGER indices(0:nindices-1)
      INTEGER inear(0:NDIM-1,0:nindices-1)
      INTEGER ibrack(0:NDIM-1,0:nindices-1)

      REAL g(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1)
      REAL V(0:depth-1,0:*)
c
c     Input/Output: spreading adds to the existing Cartesian grid values u.
c
      REAL u(CELL2dVECG(ilower,iupper,nugc),0:depth-1)
c
c     Local variables.
c
      INTEGER ilower(0:NDIM-1),iupper(0:NDIM-1)
      INTEGER ic0,ic1
      INTEGER ic_center(0:NDIM-1),ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER ic_trimmed_lower(0:NDIM-1),ic_trimmed_upper(0:NDIM-1)
      INTEGER d,l,s,nugc(0:NDIM-1)

      REAL w(0:NDIM-1,0:1)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1

      iupper(0) = iupper0
      iupper(1) = iupper1

      nugc(0) = nugc0
      nugc(1) = nugc1
c
c     Use the piecewise linear delta function to interpolate u onto V.
c
      do l = 0,nindices-1
         s = indices(l)
c
c     Determine the interpolation stencils and weights.
c
         do d = 0,NDIM-1
            ic_center(d) = inear(d,l)

            if ( ibrack(d,l).lt.inear(d,l) ) then
               ic_lower(d) = ic_center(d)-1
               ic_upper(d) = ic_center(d)
               w(d,0) = dble(ic_center(d)-ibrack(d,l))-g(d,l)
               w(d,1) = 1.d0 - w(d,0)
            else
               ic_lower(d) = ic_center(d)
               ic_upper(d) = ic_center(d)+1
               w(d,0) = 1.d0 + dble(ic_center(d)-ibrack(d,l))-g(d,l)
               w(d,1) = 1.d0 - w(d,0)
            endif

            ic_trimmed_lower(d) = max(ic_lower(d),ilower(d)-nugc(d))
            ic_trimmed_upper(d) = min(ic_upper(d),iupper(d)+nugc(d))
         enddo
c
c     Spread V onto u.
c
         do d = 0,depth-1
            do ic1 = ic_trimmed_lower(1),ic_trimmed_upper(1)
               do ic0 = ic_trimmed_lower(0),ic_trimmed_upper(0)
                  u(ic0,ic1,d) = u(ic0,ic1,d)+(
     &                 w(0,ic0-ic_lower(0))*
     &                 w(1,ic1-ic_lower(1))*
     &                 V(d,s)/(dx(0)*dx(1)))
               enddo
            enddo
         enddo
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Interpolate u onto V at the points listed in indices using the
c     piecewise cubic delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_piecewise_cubic_interp2d(
     &     dx,depth,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u,
     &     indices,g,inear,ibrack,nindices,
     &     V)
c
      implicit none
c
c     Functions.
c
      REAL lagrangian_piecewise_cubic_delta
c
c     Input.
c
      INTEGER depth
      INTEGER ilower0,iupper0,ilower1,iupper1
      INTEGER nugc0,nugc1
      INTEGER nindices

      INTEGER indices(0:nindices-1)
      INTEGER inear(0:NDIM-1,0:nindices-1)
      INTEGER ibrack(0:NDIM-1,0:nindices-1)

      REAL g(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1)
      REAL u(CELL2dVECG(ilower,iupper,nugc),0:depth-1)
c
c     Input/Output.
c
      REAL V(0:depth-1,0:*)
c
c     Local variables.
c
      INTEGER ilower(0:NDIM-1),iupper(0:NDIM-1)
      INTEGER ic,ic0,ic1
      INTEGER ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER d,l,s,nugc(0:NDIM-1)

      REAL w(0:NDIM-1,0:3)
c
c     Prevent compiler warning about unused variables.
c
      inear(0,0) = inear(0,0)
      dx(0) = dx(0)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1

      iupper(0) = iupper0
      iupper(1) = iupper1

      nugc(0) = nugc0
      nugc(1) = nugc1
c
c     Use the piecewise cubic delta function to interpolate u onto V.
c
      do l = 0,nindices-1
         s = indices(l)
         do d=0,NDIM-1
c
c     Determine the interpolation stencil from the data point below s.
c
            ic_lower(d) = ibrack(d,l)-1
            ic_upper(d) = ibrack(d,l)+2

            ic_lower(d) = max(ic_lower(d),ilower(d)-nugc(d))
            ic_upper(d) = min(ic_upper(d),iupper(d)+nugc(d))
c
c     Compute the interpolation weights.
c
            do ic = ic_lower(d),ic_upper(d)
               w(d,ic-ic_lower(d)) =
     &              lagrangian_piecewise_cubic_delta(
     &              g(d,l)+dble(ibrack(d,l)-ic))
            enddo
         enddo
c
c     Interpolate u onto V.
c
         INTERPOLATE_2D_SPECIALIZE_FIXED_WIDTH(ic_lower(1), ic_upper(1),
                                               ic_lower(0), ic_upper(0),
                                               4)
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Spread V onto u at the points listed in indices using the
c     piecewise cubic delta function using standard (double) precision
c     accumulation on the Cartesian grid.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_piecewise_cubic_spread2d(
     &     dx,depth,
     &     indices,g,inear,ibrack,nindices,
     &     V,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u)
c
      implicit none
c
c     Functions.
c
      REAL lagrangian_piecewise_cubic_delta
c
c     Input.
c
      INTEGER depth
      INTEGER nindices
      INTEGER ilower0,iupper0,ilower1,iupper1
      INTEGER nugc0,nugc1

      INTEGER indices(0:nindices-1)
      INTEGER inear(0:NDIM-1,0:nindices-1)
      INTEGER ibrack(0:NDIM-1,0:nindices-1)

      REAL g(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1)
      REAL V(0:depth-1,0:*)
c
c     Input/Output: spreading adds to the existing Cartesian grid values u.
c
      REAL u(CELL2dVECG(ilower,iupper,nugc),0:depth-1)
c
c     Local variables.
c
      INTEGER ilower(0:NDIM-1),iupper(0:NDIM-1)
      INTEGER ic,ic0,ic1
      INTEGER ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER d,l,s,nugc(0:NDIM-1)

      REAL w(0:NDIM-1,0:3)
c
c     Prevent compiler warning about unused variables.
c
      inear(0,0) = inear(0,0)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1

      iupper(0) = iupper0
      iupper(1) = iupper1

      nugc(0) = nugc0
      nugc(1) = nugc1
c
c     Use the piecewise cubic delta function to spread V onto u.
c
      do l = 0,nindices-1
         s = indices(l)
         do d = 0,NDIM-1
c
c     Determine the spreading stencil from the data point below s.
c
            ic_lower(d) = ibrack(d,l)-1
            ic_upper(d) = ibrack(d,l)+2
            ic_lower(d) = max(ic_lower(d),ilower(d)-nugc(d))
            ic_upper(d) = min(ic_upper(d),iupper(d)+nugc(d))
c
c     Compute the spreading weights.
c
            do ic = ic_lower(d),ic_upper(d)
               w(d,ic-ic_lower(d)) =
     &              lagrangian_piecewise_cubic_delta(
     &              g(d,l)+dble(ibrack(d,l)-ic))
            enddo
         enddo
c
c     Spread V onto u.
c
         SPREAD_2D_SPECIALIZE_FIXED_WIDTH(ic_lower(1), ic_upper(1),
                                          ic_lower(0), ic_upper(0),
                                          4)
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Interpolate u onto V at the points listed in indices using the IB
c     3-point delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_ib_3_interp2d(
     &     dx,depth,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u,
     &     indices,g,inear,ibrack,nindices,
     &     V)
c
      implicit none
c
c     Functions.
c
      REAL lagrangian_ib_3_delta
c
c     Input.
c
      INTEGER depth
      INTEGER ilower0,iupper0,ilower1,iupper1
      INTEGER nugc0,nugc1
      INTEGER nindices

      INTEGER indices(0:nindices-1)
      INTEGER inear(0:NDIM-1,0:nindices-1)
      INTEGER ibrack(0:NDIM-1,0:nindices-1)

      REAL g(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1)
      REAL u(CELL2dVECG(ilower,iupper,nugc),0:depth-1)
c
c     Input/Output.
c
      REAL V(0:depth-1,0:*)
c
c     Local variables.
c
      INTEGER ilower(0:NDIM-1),iupper(0:NDIM-1)
      INTEGER ic,ic0,ic1
      INTEGER ic_center(0:NDIM-1),ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER d,l,s,nugc(0:NDIM-1)

      REAL w(0:NDIM-1,0:2)
c
c     Prevent compiler warning about unused variables.
c
      dx(0) = dx(0)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1

      iupper(0) = iupper0
      iupper(1) = iupper1

      nugc(0) = nugc0
      nugc(1) = nugc1
c
c     Use the IB 3-point delta function to interpolate u onto V.
c
      do l = 0,nindices-1
         s = indices(l)
         do d = 0,NDIM-1
c
c     Determine the Cartesian cell in which point s is located.
c
            ic_center(d) = inear(d,l)
c
c     Determine the interpolation stencil corresponding to the position
c     of point s within the cell.
c
            ic_lower(d) = ic_center(d)-1
            ic_upper(d) = ic_center(d)+1
            ic_lower(d) = max(ic_lower(d),ilower(d)-nugc(d))
            ic_upper(d) = min(ic_upper(d),iupper(d)+nugc(d))
c
c     Compute the interpolation weights.
c
            do ic = ic_lower(d),ic_upper(d)
               w(d,ic-ic_lower(d)) =
     &              lagrangian_ib_3_delta(
     &              g(d,l)+dble(ibrack(d,l)-ic))
            enddo
         enddo
c
c     Interpolate u onto V.
c
         INTERPOLATE_2D_SPECIALIZE_FIXED_WIDTH(ic_lower(1), ic_upper(1),
                                               ic_lower(0), ic_upper(0),
                                               3)
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Spread V onto u at the points listed in indices using the IB
c     3-point delta function using standard (double) precision
c     accumulation on the Cartesian grid.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_ib_3_spread2d(
     &     dx,depth,
     &     indices,g,inear,ibrack,nindices,
     &     V,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u)
c
      implicit none
c
c     Functions.
c
      REAL lagrangian_ib_3_delta
c
c     Input.
c
      INTEGER depth
      INTEGER nindices
      INTEGER ilower0,iupper0,ilower1,iupper1
      INTEGER nugc0,nugc1

      INTEGER indices(0:nindices-1)
      INTEGER inear(0:NDIM-1,0:nindices-1)
      INTEGER ibrack(0:NDIM-1,0:nindices-1)

      REAL g(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1)
      REAL V(0:depth-1,0:*)
c
c     Input/Output: spreading adds to the existing Cartesian grid values u.
c
      REAL u(CELL2dVECG(ilower,iupper,nugc),0:depth-1)
c
c     Local variables.
c
      INTEGER ilower(0:NDIM-1),iupper(0:NDIM-1)
      INTEGER ic,ic0,ic1
      INTEGER ic_center(0:NDIM-1),ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER d,l,s,nugc(0:NDIM-1)

      REAL w(0:NDIM-1,0:2)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1

      iupper(0) = iupper0
      iupper(1) = iupper1

      nugc(0) = nugc0
      nugc(1) = nugc1
c
c     Use the IB 3-point delta function to spread V onto u.
c
      do l = 0,nindices-1
         s = indices(l)
         do d = 0,NDIM-1
c
c     Determine the Cartesian cell in which point s is located.
c
            ic_center(d) = inear(d,l)
c
c     Determine the spreading stencil corresponding to the position of
c     point s within the cell.
c
            ic_lower(d) = ic_center(d)-1
            ic_upper(d) = ic_center(d)+1
            ic_lower(d) = max(ic_lower(d),ilower(d)-nugc(d))
            ic_upper(d) = min(ic_upper(d),iupper(d)+nugc(d))
c
c     Compute the spreading weights.
c
            do ic = ic_lower(d),ic_upper(d)
               w(d, ic-ic_lower(d)) =
     &              lagrangian_ib_3_delta(
     &              g(d,l)+dble(ibrack(d,l)-ic))
            enddo
         enddo
c
c     Spread V onto u.
c
         SPREAD_2D_SPECIALIZE_FIXED_WIDTH(ic_lower(1), ic_upper(1),
                                          ic_lower(0), ic_upper(0),
                                          3)
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Interpolate u onto V at the points listed in indices using the IB
c     4-point delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_ib_4_interp2d(
     &     dx,depth,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u,
     &     indices,g,inear,ibrack,nindices,
     &     V)
c
      implicit none
c
c     Input.
c
      INTEGER depth
      INTEGER ilower0,iupper0,ilower1,iupper1
      INTEGER nugc0,nugc1
      INTEGER nindices

      INTEGER indices(0:nindices-1)
      INTEGER inear(0:NDIM-1,0:nindices-1)
      INTEGER ibrack(0:NDIM-1,0:nindices-1)

      REAL g(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1)
      REAL u(CELL2dVECG(ilower,iupper,nugc),0:depth-1)
c
c     Input/Output.
c
      REAL V(0:depth-1,0:*)
c
c     Local variables.
c
      INTEGER i0,i1,ic0,ic1
      INTEGER ig_lower(0:NDIM-1),ig_upper(0:NDIM-1)
      INTEGER ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER istart0,istop0,istart1,istop1
      INTEGER d,l,s

      REAL q0,q1,r0,r1
      REAL w0(0:3),w1(0:3)
      REAL w(0:3,0:3),wy
c
c     Prevent compiler warning about unused variables.
c
      inear(0,0) = inear(0,0)
      dx(0) = dx(0)
c
c     Compute the extents of the ghost box.
c
      ig_lower(0) = ilower0-nugc0
      ig_lower(1) = ilower1-nugc1
      ig_upper(0) = iupper0+nugc0
      ig_upper(1) = iupper1+nugc1
c
c     Use the IB 4-point delta function to interpolate u onto V.
c
      do l = 0,nindices-1
         s = indices(l)
c
c     Determine the interpolation stencil corresponding to the position
c     of point s within the cell and compute the interpolation weights.
c
         ic_lower(0) = ibrack(0,l)-1
         ic_upper(0) = ic_lower(0) + 3
         r0 = g(0,l)
         q0 = sqrt(1.d0+4.d0*r0*(1.d0-r0))
         w0(0) = 0.125d0*(3.d0-2.d0*r0-q0)
         w0(1) = 0.125d0*(3.d0-2.d0*r0+q0)
         w0(2) = 0.125d0*(1.d0+2.d0*r0+q0)
         w0(3) = 0.125d0*(1.d0+2.d0*r0-q0)

         ic_lower(1) = ibrack(1,l)-1
         ic_upper(1) = ic_lower(1) + 3
         r1 = g(1,l)
         q1 = sqrt(1.d0+4.d0*r1*(1.d0-r1))
         w1(0) = 0.125d0*(3.d0-2.d0*r1-q1)
         w1(1) = 0.125d0*(3.d0-2.d0*r1+q1)
         w1(2) = 0.125d0*(1.d0+2.d0*r1+q1)
         w1(3) = 0.125d0*(1.d0+2.d0*r1-q1)
c
c     Compute the tensor product of the interpolation weights.
c
         do i1 = 0,3
            wy = w1(i1)
            do i0 = 0,3
               w(i0,i1) = w0(i0)*wy
            enddo
         enddo
c
c     Interpolate u onto V.
c
         istart0 =   max(ig_lower(0)-ic_lower(0),0)
         istop0  = 3-max(ic_upper(0)-ig_upper(0),0)
         istart1 =   max(ig_lower(1)-ic_lower(1),0)
         istop1  = 3-max(ic_upper(1)-ig_upper(1),0)
         do d = 0,depth-1
            V(d,s) = 0.d0
            do i1 = istart1,istop1
               ic1 = ic_lower(1)+i1
               do i0 = istart0,istop0
                  ic0 = ic_lower(0)+i0
                  V(d,s) = V(d,s) + w(i0,i1)*u(ic0,ic1,d)
               enddo
            enddo
         enddo
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Spread V onto u at the points listed in indices using the IB
c     4-point delta function using standard (double) precision
c     accumulation on the Cartesian grid.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_ib_4_spread2d(
     &     dx,depth,
     &     indices,g,inear,ibrack,nindices,
     &     V,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u)
c
      implicit none
c
c     Input.
c
      INTEGER depth
      INTEGER nindices
      INTEGER ilower0,iupper0,ilower1,iupper1
      INTEGER nugc0,nugc1

      INTEGER indices(0:nindices-1)
      INTEGER inear(0:NDIM-1,0:nindices-1)
      INTEGER ibrack(0:NDIM-1,0:nindices-1)

      REAL g(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1)
      REAL V(0:depth-1,0:*)
c
c     Input/Output: spreading adds to the existing Cartesian grid values u.
c
      REAL u(CELL2dVECG(ilower,iupper,nugc),0:depth-1)
c
c     Local variables.
c
      INTEGER i0,i1,ic0,ic1
      INTEGER ig_lower(0:NDIM-1),ig_upper(0:NDIM-1)
      INTEGER ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER istart0,istop0,istart1,istop1
      INTEGER d,l,s

      REAL q0,q1,r0,r1
      REAL w0(0:3),w1(0:3)
      REAL w(0:3,0:3),wy
c
c     Prevent compiler warning about unused variables.
c
      inear(0,0) = inear(0,0)
c
c     Compute the extents of the ghost box.
c
      ig_lower(0) = ilower0-nugc0
      ig_lower(1) = ilower1-nugc1
      ig_upper(0) = iupper0+nugc0
      ig_upper(1) = iupper1+nugc1
c
c     Use the IB 4-point delta function to spread V onto u.
c
      do l = 0,nindices-1
         s = indices(l)
c
c     Determine the interpolation stencil corresponding to the position
c     of point s within the cell and compute the interpolation weights.
c
         ic_lower(0) = ibrack(0,l)-1
         ic_upper(0) = ic_lower(0) + 3
         r0 = g(0,l)
         q0 = sqrt(1.d0+4.d0*r0*(1.d0-r0))
         w0(0) = 0.125d0*(3.d0-2.d0*r0-q0)
         w0(1) = 0.125d0*(3.d0-2.d0*r0+q0)
         w0(2) = 0.125d0*(1.d0+2.d0*r0+q0)
         w0(3) = 0.125d0*(1.d0+2.d0*r0-q0)

         ic_lower(1) = ibrack(1,l)-1
         ic_upper(1) = ic_lower(1) + 3
         r1 = g(1,l)
         q1 = sqrt(1.d0+4.d0*r1*(1.d0-r1))
         w1(0) = 0.125d0*(3.d0-2.d0*r1-q1)
         w1(1) = 0.125d0*(3.d0-2.d0*r1+q1)
         w1(2) = 0.125d0*(1.d0+2.d0*r1+q1)
         w1(3) = 0.125d0*(1.d0+2.d0*r1-q1)
c
c     Compute the tensor product of the scaled interpolation weights.
c
         do i1 = 0,3
            wy = w1(i1)/(dx(0)*dx(1))
            do i0 = 0,3
               w(i0,i1) = w0(i0)*wy
            enddo
         enddo
c
c     Spread V onto u.
c
         istart0 =   max(ig_lower(0)-ic_lower(0),0)
         istop0  = 3-max(ic_upper(0)-ig_upper(0),0)
         istart1 =   max(ig_lower(1)-ic_lower(1),0)
         istop1  = 3-max(ic_upper(1)-ig_upper(1),0)
         do d = 0,depth-1
            do i1 = istart1,istop1
               ic1 = ic_lower(1)+i1
               do i0 = istart0,istop0
                  ic0 = ic_lower(0)+i0
                  u(ic0,ic1,d) = u(ic0,ic1,d) + w(i0,i1)*V(d,s)
               enddo
            enddo
         enddo
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Interpolate u onto V at the points listed in indices using the IB
c     4-point delta function that has been broadened to have a support
c     of 8 meshwidths.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_ib_4_w8_interp2d(
     &     dx,depth,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u,
     &     indices,g,inear,ibrack,nindices,
     &     V)
c
      implicit none
c
c     Input.
c
      INTEGER depth
      INTEGER ilower0,iupper0,ilower1,iupper1
      INTEGER nugc0,nugc1
      INTEGER nindices

      INTEGER indices(0:nindices-1)
      INTEGER inear(0:NDIM-1,0:nindices-1)
      INTEGER ibrack(0:NDIM-1,0:nindices-1)

      REAL g(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1)
      REAL u(CELL2dVECG(ilower,iupper,nugc),0:depth-1)
c
c     Input/Output.
c
      REAL V(0:depth-1,0:*)
c
c     Local variables.
c
      INTEGER i0,i1,ic0,ic1
      INTEGER ig_lower(0:NDIM-1),ig_upper(0:NDIM-1)
      INTEGER ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER istart0,istop0,istart1,istop1
      INTEGER d,l,s

      REAL q0,q1,r0,r1
      REAL w0(0:7),w1(0:7)
      REAL w(0:7,0:7),wy
c
c     Prevent compiler warning about unused variables.
c
      inear(0,0) = inear(0,0)
      dx(0) = dx(0)
c
c     Compute the extents of the ghost box.
c
      ig_lower(0) = ilower0-nugc0
      ig_lower(1) = ilower1-nugc1
      ig_upper(0) = iupper0+nugc0
      ig_upper(1) = iupper1+nugc1
c
c     Use a broadened version of the IB 4-point delta function to
c     interpolate u onto V.
c
      do l = 0,nindices-1
         s = indices(l)
c
c     Determine the interpolation stencil corresponding to the position
c     of point s within the cell and compute the interpolation weights.
c
         ic_lower(0) = ibrack(0,l)-3
         ic_upper(0) = ic_lower(0) + 7
         r0 = 0.5d0*g(0,l)
         q0 = sqrt(1.d0+4.d0*r0*(1.d0-r0))
         w0(1) = 0.0625d0*(3.d0-2.d0*r0-q0)
         w0(3) = 0.0625d0*(3.d0-2.d0*r0+q0)
         w0(5) = 0.0625d0*(1.d0+2.d0*r0+q0)
         w0(7) = 0.0625d0*(1.d0+2.d0*r0-q0)
         r0 = r0+0.5d0
         q0 = sqrt(1.d0+4.d0*r0*(1.d0-r0))
         w0(0) = 0.0625d0*(3.d0-2.d0*r0-q0)
         w0(2) = 0.0625d0*(3.d0-2.d0*r0+q0)
         w0(4) = 0.0625d0*(1.d0+2.d0*r0+q0)
         w0(6) = 0.0625d0*(1.d0+2.d0*r0-q0)

         ic_lower(1) = ibrack(1,l)-3
         ic_upper(1) = ic_lower(1) + 7
         r1 = 0.5d0*g(1,l)
         q1 = sqrt(1.d0+4.d0*r1*(1.d0-r1))
         w1(1) = 0.0625d0*(3.d0-2.d0*r1-q1)
         w1(3) = 0.0625d0*(3.d0-2.d0*r1+q1)
         w1(5) = 0.0625d0*(1.d0+2.d0*r1+q1)
         w1(7) = 0.0625d0*(1.d0+2.d0*r1-q1)
         r1 = r1+0.5d0
         q1 = sqrt(1.d0+4.d0*r1*(1.d0-r1))
         w1(0) = 0.0625d0*(3.d0-2.d0*r1-q1)
         w1(2) = 0.0625d0*(3.d0-2.d0*r1+q1)
         w1(4) = 0.0625d0*(1.d0+2.d0*r1+q1)
         w1(6) = 0.0625d0*(1.d0+2.d0*r1-q1)
c
c     Compute the tensor product of the interpolation weights.
c
         do i1 = 0,7
            wy = w1(i1)
            do i0 = 0,7
               w(i0,i1) = w0(i0)*wy
            enddo
         enddo
c
c     Interpolate u onto V.
c
         istart0 =   max(ig_lower(0)-ic_lower(0),0)
         istop0  = 7-max(ic_upper(0)-ig_upper(0),0)
         istart1 =   max(ig_lower(1)-ic_lower(1),0)
         istop1  = 7-max(ic_upper(1)-ig_upper(1),0)
         do d = 0,depth-1
            V(d,s) = 0.d0
            do i1 = istart1,istop1
               ic1 = ic_lower(1)+i1
               do i0 = istart0,istop0
                  ic0 = ic_lower(0)+i0
                  V(d,s) = V(d,s) + w(i0,i1)*u(ic0,ic1,d)
               enddo
            enddo
         enddo
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Spread V onto u at the points listed in indices using the IB
c     4-point delta function that has been broadened to have a support
c     of 8 meshwidths.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_ib_4_w8_spread2d(
     &     dx,depth,
     &     indices,g,inear,ibrack,nindices,
     &     V,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u)
c
      implicit none
c
c     Input.
c
      INTEGER depth
      INTEGER nindices
      INTEGER ilower0,iupper0,ilower1,iupper1
      INTEGER nugc0,nugc1

      INTEGER indices(0:nindices-1)
      INTEGER inear(0:NDIM-1,0:nindices-1)
      INTEGER ibrack(0:NDIM-1,0:nindices-1)

      REAL g(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1)
      REAL V(0:depth-1,0:*)
c
c     Input/Output: spreading adds to the existing Cartesian grid values u.
c
      REAL u(CELL2dVECG(ilower,iupper,nugc),0:depth-1)
c
c     Local variables.
c
      INTEGER i0,i1,ic0,ic1
      INTEGER ig_lower(0:NDIM-1),ig_upper(0:NDIM-1)
      INTEGER ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER istart0,istop0,istart1,istop1
      INTEGER d,l,s

      REAL q0,q1,r0,r1
      REAL w0(0:7),w1(0:7)
      REAL w(0:7,0:7),wy
c
c     Prevent compiler warning about unused variables.
c
      inear(0,0) = inear(0,0)
c
c     Compute the extents of the ghost box.
c
      ig_lower(0) = ilower0-nugc0
      ig_lower(1) = ilower1-nugc1
      ig_upper(0) = iupper0+nugc0
      ig_upper(1) = iupper1+nugc1
c
c     Use a broadened version of the IB 4-point delta function to spread
c     V onto u.
c
      do l = 0,nindices-1
         s = indices(l)
c
c     Determine the interpolation stencil corresponding to the position
c     of point s within the cell and compute the interpolation weights.
c
         ic_lower(0) = ibrack(0,l)-3
         ic_upper(0) = ic_lower(0) + 7
         r0 = 0.5d0*g(0,l)
         q0 = sqrt(1.d0+4.d0*r0*(1.d0-r0))
         w0(1) = 0.0625d0*(3.d0-2.d0*r0-q0)
         w0(3) = 0.0625d0*(3.d0-2.d0*r0+q0)
         w0(5) = 0.0625d0*(1.d0+2.d0*r0+q0)
         w0(7) = 0.0625d0*(1.d0+2.d0*r0-q0)
         r0 = r0+0.5d0
         q0 = sqrt(1.d0+4.d0*r0*(1.d0-r0))
         w0(0) = 0.0625d0*(3.d0-2.d0*r0-q0)
         w0(2) = 0.0625d0*(3.d0-2.d0*r0+q0)
         w0(4) = 0.0625d0*(1.d0+2.d0*r0+q0)
         w0(6) = 0.0625d0*(1.d0+2.d0*r0-q0)

         ic_lower(1) = ibrack(1,l)-3
         ic_upper(1) = ic_lower(1) + 7
         r1 = 0.5d0*g(1,l)
         q1 = sqrt(1.d0+4.d0*r1*(1.d0-r1))
         w1(1) = 0.0625d0*(3.d0-2.d0*r1-q1)
         w1(3) = 0.0625d0*(3.d0-2.d0*r1+q1)
         w1(5) = 0.0625d0*(1.d0+2.d0*r1+q1)
         w1(7) = 0.0625d0*(1.d0+2.d0*r1-q1)
         r1 = r1+0.5d0
         q1 = sqrt(1.d0+4.d0*r1*(1.d0-r1))
         w1(0) = 0.0625d0*(3.d0-2.d0*r1-q1)
         w1(2) = 0.0625d0*(3.d0-2.d0*r1+q1)
         w1(4) = 0.0625d0*(1.d0+2.d0*r1+q1)
         w1(6) = 0.0625d0*(1.d0+2.d0*r1-q1)
c
c     Compute the tensor product of the scaled interpolation weights.
c
         do i1 = 0,7
            wy = w1(i1)/(dx(0)*dx(1))
            do i0 = 0,7
               w(i0,i1) = w0(i0)*wy
            enddo
         enddo
c
c     Spread V onto u.
c
         istart0 =   max(ig_lower(0)-ic_lower(0),0)
         istop0  = 7-max(ic_upper(0)-ig_upper(0),0)
         istart1 =   max(ig_lower(1)-ic_lower(1),0)
         istop1  = 7-max(ic_upper(1)-ig_upper(1),0)
         do d = 0,depth-1
            do i1 = istart1,istop1
               ic1 = ic_lower(1)+i1
               do i0 = istart0,istop0
                  ic0 = ic_lower(0)+i0
                  u(ic0,ic1,d) = u(ic0,ic1,d) + w(i0,i1)*V(d,s)
               enddo
            enddo
         enddo
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Interpolate u onto V at the points listed in indices using the
c     5-point IB delta with three continuous derivatives
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_ib_5_interp2d(
     &     dx,depth,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u,
     &     indices,g,inear,ibrack,nindices,
     &     V)
c
      implicit none
c
c     Input.
c
      INTEGER depth
      INTEGER ilower0,iupper0,ilower1,iupper1
      INTEGER nugc0,nugc1
      INTEGER nindices

      INTEGER indices(0:nindices-1)
      INTEGER inear(0:NDIM-1,0:nindices-1)
      INTEGER ibrack(0:NDIM-1,0:nindices-1)

      REAL g(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1)
      REAL u(CELL2dVECG(ilower,iupper,nugc),0:depth-1)

c
c     Input/Output.
c
      REAL V(0:depth-1,0:*)
c
c     Local variables.
c
      REAL r
      REAL phi
      REAL K
      INTEGER i0,i1,ic0,ic1
      INTEGER ig_lower(0:NDIM-1),ig_upper(0:NDIM-1)
      INTEGER ic_center(0:NDIM-1),ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER istart0,istop0,istart1,istop1
      INTEGER d,l,s

      REAL w0(0:4),w1(0:4)

      PARAMETER (K = (38.0d0 - sqrt(69.0d0))/60.0d0)
c
c     Prevent compiler warning about unused variables.
c
      dx(0) = dx(0)
c
c     Compute the extents of the ghost box.
c
      ig_lower(0) = ilower0-nugc0
      ig_lower(1) = ilower1-nugc1
      ig_upper(0) = iupper0+nugc0
      ig_upper(1) = iupper1+nugc1
c
c     Use a 5-point IB delta function to interpolate u onto V.
c
      do l = 0,nindices-1
         s = indices(l)
c
c     Determine the Cartesian cell in which point s is located.
c
         ic_center(0) = inear(0,l)
         ic_center(1) = inear(1,l)

c
c     Determine the interpolation stencil corresponding to the position
c     of point s within the cell.
c
         do d = 0,NDIM-1
            ic_lower(d) = ic_center(d)-2
            ic_upper(d) = ic_center(d)+2
         enddo
c
c     Compute the interpolation weights.
c
         ic0 = ic_center(0)
         r = g(0,l)+dble(ibrack(0,l)-ic0)
         phi = (136.0d0 - 40.0d0*K - 40.0d0*r**2 + sqrt(2.0d0)
     &          *sqrt(3123.0d0 - 6840.0d0*K + 3600.0d0*(K**2)
     &          - 12440.0d0*(r**2) + 25680.0d0*K*(r**2)
     &          - 12600.0d0*(K**2)*(r**2) + 8080.0d0*(r**4)
     &          - 8400.0d0*K*(r**4) - 1400.0d0*(r**6)))
     &          /280.0d0

         w0(0) = (1.0d0/12.0d0) * (-2.0d0 + 2.0d0*phi + 2.0d0*K +
     &                         r - 3.0d0*K*r + 2.0d0*r**2 - r**3)
         w0(1) = (1.0d0/ 6.0d0) * ( 4.0d0 - 4.0d0*phi -       K -
     &                   4.0d0*r + 3.0d0*K*r -       r**2 + r**3)
         w0(2) = phi
         w0(3) = (1.0d0/ 6.0d0) * ( 4.0d0 - 4.0d0*phi -       K +
     &                   4.0d0*r - 3.0d0*K*r -       r**2 - r**3)
         w0(4) = (1.0d0/12.0d0) * (-2.0d0 + 2.0d0*phi + 2.0d0*K -
     &                         r + 3.0d0*K*r + 2.0d0*r**2 + r**3)

         ic1 = ic_center(1)
         r = g(1,l)+dble(ibrack(1,l)-ic1)
         phi = (136.0d0 - 40.0d0*K - 40.0d0*r**2 + sqrt(2.0d0)
     &          *sqrt(3123.0d0 - 6840.0d0*K + 3600.0d0*(K**2)
     &          - 12440.0d0*(r**2) + 25680.0d0*K*(r**2)
     &          - 12600.0d0*(K**2)*(r**2) + 8080.0d0*(r**4)
     &          - 8400.0d0*K*(r**4) - 1400.0d0*(r**6)))
     &          /280.0d0

         w1(0) = (1.0d0/12.0d0) * (-2.0d0 + 2.0d0*phi + 2.0d0*K +
     &                         r - 3.0d0*K*r + 2.0d0*r**2 - r**3)
         w1(1) = (1.0d0/ 6.0d0) * ( 4.0d0 - 4.0d0*phi -       K -
     &                   4.0d0*r + 3.0d0*K*r -       r**2 + r**3)
         w1(2) = phi
         w1(3) = (1.0d0/ 6.0d0) * ( 4.0d0 - 4.0d0*phi -       K +
     &                   4.0d0*r - 3.0d0*K*r -       r**2 - r**3)
         w1(4) = (1.0d0/12.0d0) * (-2.0d0 + 2.0d0*phi + 2.0d0*K -
     &                         r + 3.0d0*K*r + 2.0d0*r**2 + r**3)

c
c     Interpolate u onto V.
c
         istart0 =   max(ig_lower(0)-ic_lower(0),0)
         istop0  = 4-max(ic_upper(0)-ig_upper(0),0)
         istart1 =   max(ig_lower(1)-ic_lower(1),0)
         istop1  = 4-max(ic_upper(1)-ig_upper(1),0)
         do d = 0,depth-1
            V(d,s) = 0.d0
            do i1 = istart1,istop1
               ic1 = ic_lower(1)+i1
               do i0 = istart0,istop0
                  ic0 = ic_lower(0)+i0
                  V(d,s) = V(d,s) + w0(i0) * w1(i1) * u(ic0,ic1,d)
               enddo
            enddo
         enddo
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Spread V onto u at the points listed in indices using the
c     5-point IB delta with three continuous derivatives
c     using standard (double) precision
c     accumulation on the Cartesian grid.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_ib_5_spread2d(
     &     dx,depth,
     &     indices,g,inear,ibrack,nindices,
     &     V,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u)
c
      implicit none
c
c     Input.
c
      INTEGER depth
      INTEGER nindices
      INTEGER ilower0,iupper0,ilower1,iupper1
      INTEGER nugc0,nugc1

      INTEGER indices(0:nindices-1)
      INTEGER inear(0:NDIM-1,0:nindices-1)
      INTEGER ibrack(0:NDIM-1,0:nindices-1)

      REAL g(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1)
      REAL V(0:depth-1,0:*)
c
c     Input/Output: spreading adds to the existing Cartesian grid values u.
c
      REAL u(CELL2dVECG(ilower,iupper,nugc),0:depth-1)
c
c     Local variables.
c
      REAL r
      REAL phi
      REAL K
      INTEGER i0,i1,ic0,ic1
      INTEGER ig_lower(0:NDIM-1),ig_upper(0:NDIM-1)
      INTEGER ic_center(0:NDIM-1),ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER istart0,istop0,istart1,istop1
      INTEGER d,l,s

      REAL w0(0:4),w1(0:4)

      PARAMETER (K = (38.0d0 - sqrt(69.0d0))/60.0d0)
c
c     Compute the extents of the ghost box.
c
      ig_lower(0) = ilower0-nugc0
      ig_lower(1) = ilower1-nugc1
      ig_upper(0) = iupper0+nugc0
      ig_upper(1) = iupper1+nugc1
c
c     Use a 5-point IB delta function to spread V onto u.
c
      do l = 0,nindices-1
         s = indices(l)
c
c     Determine the Cartesian cell in which point s is located.
c
         ic_center(0) = inear(0,l)
         ic_center(1) = inear(1,l)

c
c     Determine the spreading stencil corresponding to the position of
c     point s within the cell.
c
         do d = 0,NDIM-1
            ic_lower(d) = ic_center(d)-2
            ic_upper(d) = ic_center(d)+2
         enddo
c
c     Compute the spreading weights.
c

         ic0 = ic_center(0)
         r = g(0,l)+dble(ibrack(0,l)-ic0)
         phi = (136.0d0 - 40.0d0*K - 40.0d0*r**2 + sqrt(2.0d0)
     &          *sqrt(3123.0d0 - 6840.0d0*K + 3600.0d0*(K**2)
     &          - 12440.0d0*(r**2) + 25680.0d0*K*(r**2)
     &          - 12600.0d0*(K**2)*(r**2) + 8080.0d0*(r**4)
     &          - 8400.0d0*K*(r**4) - 1400.0d0*(r**6)))
     &          /280.0d0

         w0(0) = (1.0d0/12.0d0) * (-2.0d0 + 2.0d0*phi + 2.0d0*K +
     &                         r - 3.0d0*K*r + 2.0d0*r**2 - r**3)
         w0(1) = (1.0d0/ 6.0d0) * ( 4.0d0 - 4.0d0*phi -       K -
     &                   4.0d0*r + 3.0d0*K*r -       r**2 + r**3)
         w0(2) = phi
         w0(3) = (1.0d0/ 6.0d0) * ( 4.0d0 - 4.0d0*phi -       K +
     &                   4.0d0*r - 3.0d0*K*r -       r**2 - r**3)
         w0(4) = (1.0d0/12.0d0) * (-2.0d0 + 2.0d0*phi + 2.0d0*K -
     &                         r + 3.0d0*K*r + 2.0d0*r**2 + r**3)

         ic1 = ic_center(1)
         r = g(1,l)+dble(ibrack(1,l)-ic1)
         phi = (136.0d0 - 40.0d0*K - 40.0d0*r**2 + sqrt(2.0d0)
     &          *sqrt(3123.0d0 - 6840.0d0*K + 3600.0d0*(K**2)
     &          - 12440.0d0*(r**2) + 25680.0d0*K*(r**2)
     &          - 12600.0d0*(K**2)*(r**2) + 8080.0d0*(r**4)
     &          - 8400.0d0*K*(r**4) - 1400.0d0*(r**6)))
     &          /280.0d0

         w1(0) = (1.0d0/12.0d0) * (-2.0d0 + 2.0d0*phi + 2.0d0*K +
     &                         r - 3.0d0*K*r + 2.0d0*r**2 - r**3)
         w1(1) = (1.0d0/ 6.0d0) * ( 4.0d0 - 4.0d0*phi -       K -
     &                   4.0d0*r + 3.0d0*K*r -       r**2 + r**3)
         w1(2) = phi
         w1(3) = (1.0d0/ 6.0d0) * ( 4.0d0 - 4.0d0*phi -       K +
     &                   4.0d0*r - 3.0d0*K*r -       r**2 - r**3)
         w1(4) = (1.0d0/12.0d0) * (-2.0d0 + 2.0d0*phi + 2.0d0*K -
     &                         r + 3.0d0*K*r + 2.0d0*r**2 + r**3)

c
c     Spread V onto u.
c
         istart0 =   max(ig_lower(0)-ic_lower(0),0)
         istop0  = 4-max(ic_upper(0)-ig_upper(0),0)
         istart1 =   max(ig_lower(1)-ic_lower(1),0)
         istop1  = 4-max(ic_upper(1)-ig_upper(1),0)
         do d = 0,depth-1
            do i1 = istart1,istop1
               ic1 = ic_lower(1)+i1
               do i0 = istart0,istop0
                  ic0 = ic_lower(0)+i0
                  u(ic0,ic1,d) = u(ic0,ic1,d)+(
     &                 w0(i0) * w1(i1) * V(d,s)/(dx(0)*dx(1)))
               enddo
            enddo
         enddo
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Interpolate u onto V at the points listed in indices using the IB
c     6-point delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_ib_6_interp2d(
     &     dx,depth,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u,
     &     indices,g,inear,ibrack,nindices,
     &     V)
c
      implicit none
c
c     Input.
c
      INTEGER depth
      INTEGER ilower0,iupper0,ilower1,iupper1
      INTEGER nugc0,nugc1
      INTEGER nindices

      INTEGER indices(0:nindices-1)
      INTEGER inear(0:NDIM-1,0:nindices-1)
      INTEGER ibrack(0:NDIM-1,0:nindices-1)

      REAL g(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1)
      REAL u(CELL2dVECG(ilower,iupper,nugc),0:depth-1)
c
c     Input/Output.
c
      REAL V(0:depth-1,0:*)
c
c     Local variables.
c
      INTEGER i0,i1,ic0,ic1
      INTEGER ig_lower(0:NDIM-1),ig_upper(0:NDIM-1)
      INTEGER ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER istart0,istop0,istart1,istop1
      INTEGER d,l,s

      REAL r,alpha,beta,gamma,discr,K
      REAL pm3,pm2,pm1,p,pp1,pp2
      REAL w0(0:5),w1(0:5)
      REAL w(0:5,0:5),wy

      PARAMETER (K = (59.d0/60.d0)*(1.d0-sqrt(1.d0-(3220.d0/3481.d0))))
c
c     Prevent compiler warning about unused variables.
c
      inear(0,0) = inear(0,0)
      dx(0) = dx(0)
c
c     Compute the extents of the ghost box.
c
      ig_lower(0) = ilower0-nugc0
      ig_lower(1) = ilower1-nugc1
      ig_upper(0) = iupper0+nugc0
      ig_upper(1) = iupper1+nugc1
c
c     Use the IB 6-point delta function to interpolate u onto V.
c
      do l = 0,nindices-1
         s = indices(l)
c
c     Determine the interpolation stencil corresponding to the position
c     of point s within the cell and compute the interpolation weights.
c
         ic_lower(0) = ibrack(0,l)-2
         ic_upper(0) = ic_lower(0) + 5
         r = 1.d0 - g(0,l)

         alpha = 28.d0
         beta = (9.d0/4.d0)-(3.d0/2.d0)*(K+r**2)+((22.d0/3.d0)-7.d0*K)*r
     $        -(7.d0/3.d0)*r**3
         gamma = (1.d0/4.d0)*( ((161.d0/36.d0)-(59.d0/6.d0)*K+5.d0*K**2)
     $        *(1.d0/2.d0)*r**2 + (-(109.d0/24.d0)+5.d0*K)*(1.d0/3.d0)*r
     $        **4+ (5.d0/18.d0)*r**6 )
         discr = beta**2-4.d0*alpha*gamma

         pm3 = (-beta+sign(1.d0,(3.d0/2.d0)-K)*sqrt(discr))/(2.d0*alpha)
         pm2 =  -3.d0*pm3 - (1.d0/16.d0) + (1.d0/8.d0)*(K+r**2) + (1.d0
     $        /12.d0)*(3.d0*K-1.d0)*r + (1.d0/12.d0)*r**3
         pm1 =   2.d0*pm3 + (1.d0/4.d0) +  (1.d0/6.d0)*(4.d0-3.d0*K)*r -
     $        (1.d0/6.d0)*r**3
         p   =   2.d0*pm3 + (5.d0/8.d0)  - (1.d0/4.d0)*(K+r**2)
         pp1 =  -3.d0*pm3 + (1.d0/4.d0) -  (1.d0/6.d0)*(4.d0-3.d0*K)*r +
     $        (1.d0/6.d0)*r**3
         pp2 =        pm3 - (1.d0/16.d0) + (1.d0/8.d0)*(K+r**2) - (1.d0
     $        /12.d0)*(3.d0*K-1.d0)*r - (1.d0/12.d0)*r**3

         w0(0) = pm3
         w0(1) = pm2
         w0(2) = pm1
         w0(3) = p
         w0(4) = pp1
         w0(5) = pp2

         ic_lower(1) = ibrack(1,l)-2
         ic_upper(1) = ic_lower(1) + 5
         r = 1.d0 - g(1,l)

         alpha = 28.d0
         beta = (9.d0/4.d0)-(3.d0/2.d0)*(K+r**2)+((22.d0/3.d0)-7.d0*K)*r
     $        -(7.d0/3.d0)*r**3
         gamma = (1.d0/4.d0)*( ((161.d0/36.d0)-(59.d0/6.d0)*K+5.d0*K**2)
     $        *(1.d0/2.d0)*r**2 + (-(109.d0/24.d0)+5.d0*K)*(1.d0/3.d0)*r
     $        **4+ (5.d0/18.d0)*r**6 )
         discr = beta**2-4.d0*alpha*gamma

         pm3 = (-beta+sign(1.d0,(3.d0/2.d0)-K)*sqrt(discr))/(2.d0*alpha)
         pm2 =  -3.d0*pm3 - (1.d0/16.d0) + (1.d0/8.d0)*(K+r**2) + (1.d0
     $        /12.d0)*(3.d0*K-1.d0)*r + (1.d0/12.d0)*r**3
         pm1 =   2.d0*pm3 + (1.d0/4.d0) +  (1.d0/6.d0)*(4.d0-3.d0*K)*r -
     $        (1.d0/6.d0)*r**3
         p   =   2.d0*pm3 + (5.d0/8.d0)  - (1.d0/4.d0)*(K+r**2)
         pp1 =  -3.d0*pm3 + (1.d0/4.d0) -  (1.d0/6.d0)*(4.d0-3.d0*K)*r +
     $        (1.d0/6.d0)*r**3
         pp2 =        pm3 - (1.d0/16.d0) + (1.d0/8.d0)*(K+r**2) - (1.d0
     $        /12.d0)*(3.d0*K-1.d0)*r - (1.d0/12.d0)*r**3

         w1(0) = pm3
         w1(1) = pm2
         w1(2) = pm1
         w1(3) = p
         w1(4) = pp1
         w1(5) = pp2
c
c     Compute the tensor product of the interpolation weights.
c
         do i1 = 0,5
            wy = w1(i1)
            do i0 = 0,5
               w(i0,i1) = w0(i0)*wy
            enddo
         enddo
c
c     Interpolate u onto V.
c
         istart0 =   max(ig_lower(0)-ic_lower(0),0)
         istop0  = 5-max(ic_upper(0)-ig_upper(0),0)
         istart1 =   max(ig_lower(1)-ic_lower(1),0)
         istop1  = 5-max(ic_upper(1)-ig_upper(1),0)
         do d = 0,depth-1
            V(d,s) = 0.d0
            do i1 = istart1,istop1
               ic1 = ic_lower(1)+i1
               do i0 = istart0,istop0
                  ic0 = ic_lower(0)+i0
                  V(d,s) = V(d,s) + w(i0,i1)*u(ic0,ic1,d)
               enddo
            enddo
         enddo
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Spread V onto u at the points listed in indices using the IB
c     6-point delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_ib_6_spread2d(
     &     dx,depth,
     &     indices,g,inear,ibrack,nindices,
     &     V,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u)
c
      implicit none
c
c     Input.
c
      INTEGER depth
      INTEGER nindices
      INTEGER ilower0,iupper0,ilower1,iupper1
      INTEGER nugc0,nugc1

      INTEGER indices(0:nindices-1)
      INTEGER inear(0:NDIM-1,0:nindices-1)
      INTEGER ibrack(0:NDIM-1,0:nindices-1)

      REAL g(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1)
      REAL V(0:depth-1,0:*)
c
c     Input/Output: spreading adds to the existing Cartesian grid values u.
c
      REAL u(CELL2dVECG(ilower,iupper,nugc),0:depth-1)
c
c     Local variables.
c
      INTEGER i0,i1,ic0,ic1
      INTEGER ig_lower(0:NDIM-1),ig_upper(0:NDIM-1)
      INTEGER ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER istart0,istop0,istart1,istop1
      INTEGER d,l,s

      REAL r,alpha,beta,gamma,discr,K
      REAL pm3,pm2,pm1,p,pp1,pp2
      REAL w0(0:5),w1(0:5)
      REAL w(0:5,0:5),wy

      PARAMETER (K = (59.d0/60.d0)*(1.d0-sqrt(1.d0-(3220.d0/3481.d0))))
c
c     Prevent compiler warning about unused variables.
c
      inear(0,0) = inear(0,0)
c
c     Compute the extents of the ghost box.
c
      ig_lower(0) = ilower0-nugc0
      ig_lower(1) = ilower1-nugc1
      ig_upper(0) = iupper0+nugc0
      ig_upper(1) = iupper1+nugc1
c
c     Use the IB 6-point delta function to spread V onto u.
c
      do l = 0,nindices-1
         s = indices(l)
c
c     Determine the interpolation stencil corresponding to the position
c     of point s within the cell and compute the interpolation weights.
c
         ic_lower(0) = ibrack(0,l)-2
         ic_upper(0) = ic_lower(0) + 5
         r = 1.d0 - g(0,l)

         alpha = 28.d0
         beta = (9.d0/4.d0)-(3.d0/2.d0)*(K+r**2)+((22.d0/3.d0)-7.d0*K)*r
     $        -(7.d0/3.d0)*r**3
         gamma = (1.d0/4.d0)*( ((161.d0/36.d0)-(59.d0/6.d0)*K+5.d0*K**2)
     $        *(1.d0/2.d0)*r**2 + (-(109.d0/24.d0)+5.d0*K)*(1.d0/3.d0)*r
     $        **4+ (5.d0/18.d0)*r**6 )
         discr = beta**2-4.d0*alpha*gamma

         pm3 = (-beta+sign(1.d0,(3.d0/2.d0)-K)*sqrt(discr))/(2.d0*alpha)
         pm2 =  -3.d0*pm3 - (1.d0/16.d0) + (1.d0/8.d0)*(K+r**2) + (1.d0
     $        /12.d0)*(3.d0*K-1.d0)*r + (1.d0/12.d0)*r**3
         pm1 =   2.d0*pm3 + (1.d0/4.d0) +  (1.d0/6.d0)*(4.d0-3.d0*K)*r -
     $        (1.d0/6.d0)*r**3
         p   =   2.d0*pm3 + (5.d0/8.d0)  - (1.d0/4.d0)*(K+r**2)
         pp1 =  -3.d0*pm3 + (1.d0/4.d0) -  (1.d0/6.d0)*(4.d0-3.d0*K)*r +
     $        (1.d0/6.d0)*r**3
         pp2 =        pm3 - (1.d0/16.d0) + (1.d0/8.d0)*(K+r**2) - (1.d0
     $        /12.d0)*(3.d0*K-1.d0)*r - (1.d0/12.d0)*r**3

         w0(0) = pm3
         w0(1) = pm2
         w0(2) = pm1
         w0(3) = p
         w0(4) = pp1
         w0(5) = pp2

         ic_lower(1) = ibrack(1,l)-2
         ic_upper(1) = ic_lower(1) + 5
         r = 1.d0 - g(1,l)

         alpha = 28.d0
         beta = (9.d0/4.d0)-(3.d0/2.d0)*(K+r**2)+((22.d0/3.d0)-7.d0*K)*r
     $        -(7.d0/3.d0)*r**3
         gamma = (1.d0/4.d0)*( ((161.d0/36.d0)-(59.d0/6.d0)*K+5.d0*K**2)
     $        *(1.d0/2.d0)*r**2 + (-(109.d0/24.d0)+5.d0*K)*(1.d0/3.d0)*r
     $        **4+ (5.d0/18.d0)*r**6 )
         discr = beta**2-4.d0*alpha*gamma

         pm3 = (-beta+sign(1.d0,(3.d0/2.d0)-K)*sqrt(discr))/(2.d0*alpha)
         pm2 =  -3.d0*pm3 - (1.d0/16.d0) + (1.d0/8.d0)*(K+r**2) + (1.d0
     $        /12.d0)*(3.d0*K-1.d0)*r + (1.d0/12.d0)*r**3
         pm1 =   2.d0*pm3 + (1.d0/4.d0) +  (1.d0/6.d0)*(4.d0-3.d0*K)*r -
     $        (1.d0/6.d0)*r**3
         p   =   2.d0*pm3 + (5.d0/8.d0)  - (1.d0/4.d0)*(K+r**2)
         pp1 =  -3.d0*pm3 + (1.d0/4.d0) -  (1.d0/6.d0)*(4.d0-3.d0*K)*r +
     $        (1.d0/6.d0)*r**3
         pp2 =        pm3 - (1.d0/16.d0) + (1.d0/8.d0)*(K+r**2) - (1.d0
     $        /12.d0)*(3.d0*K-1.d0)*r - (1.d0/12.d0)*r**3

         w1(0) = pm3
         w1(1) = pm2
         w1(2) = pm1
         w1(3) = p
         w1(4) = pp1
         w1(5) = pp2
c
c     Compute the tensor product of the scaled interpolation weights.
c
         do i1 = 0,5
            wy = w1(i1)/(dx(0)*dx(1))
            do i0 = 0,5
               w(i0,i1) = w0(i0)*wy
            enddo
         enddo
c
c     Spread V onto u.
c
         istart0 =   max(ig_lower(0)-ic_lower(0),0)
         istop0  = 5-max(ic_upper(0)-ig_upper(0),0)
         istart1 =   max(ig_lower(1)-ic_lower(1),0)
         istop1  = 5-max(ic_upper(1)-ig_upper(1),0)
         do d = 0,depth-1
            do i1 = istart1,istop1
               ic1 = ic_lower(1)+i1
               do i0 = istart0,istop0
                  ic0 = ic_lower(0)+i0
                  u(ic0,ic1,d) = u(ic0,ic1,d) + w(i0,i1)*V(d,s)
               enddo
            enddo
         enddo
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Interpolate u onto V at the points listed in indices using the
c     3-point B-spline delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_bspline_3_interp2d(
     &     dx,depth,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u,
     &     indices,g,inear,ibrack,nindices,
     &     V)
c
      implicit none
c
c     Functions.
c
      REAL lagrangian_bspline_3_delta
c
c     Input.
c
      INTEGER depth
      INTEGER ilower0,iupper0,ilower1,iupper1
      INTEGER nugc0,nugc1
      INTEGER nindices

      INTEGER indices(0:nindices-1)
      INTEGER inear(0:NDIM-1,0:nindices-1)
      INTEGER ibrack(0:NDIM-1,0:nindices-1)

      REAL g(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1)
      REAL u(CELL2dVECG(ilower,iupper,nugc),0:depth-1)
c
c     Input/Output.
c
      REAL V(0:depth-1,0:*)
c
c     Local variables.
c
      INTEGER ilower(0:NDIM-1),iupper(0:NDIM-1)
      INTEGER ic,ic0,ic1
      INTEGER ic_center(0:NDIM-1),ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER d,l,s,nugc(0:NDIM-1)

      REAL w(0:NDIM-1,0:2)
c
c     Prevent compiler warning about unused variables.
c
      dx(0) = dx(0)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1

      iupper(0) = iupper0
      iupper(1) = iupper1

      nugc(0) = nugc0
      nugc(1) = nugc1
c
c     Use a 3-point B-spline function to interpolate u onto V.
c
      do l = 0,nindices-1
         s = indices(l)
         do d=0,NDIM-1
c
c     Determine the Cartesian cell in which point s is located.
c
            ic_center(d) = inear(d,l)
c
c     Determine the interpolation stencil corresponding to the position
c     of point s within the cell.
c
            ic_lower(d) = ic_center(d)-1
            ic_upper(d) = ic_center(d)+1
            ic_lower(d) = max(ic_lower(d),ilower(d)-nugc(d))
            ic_upper(d) = min(ic_upper(d),iupper(d)+nugc(d))
c
c     Compute the interpolation weights.
c
            do ic = ic_lower(d),ic_upper(d)
               w(d,ic-ic_lower(d)) =
     &              lagrangian_bspline_3_delta(
     &              g(d,l)+dble(ibrack(d,l)-ic))
            enddo
         enddo
c
c     Interpolate u onto V.
c
         INTERPOLATE_2D_SPECIALIZE_FIXED_WIDTH(ic_lower(1), ic_upper(1),
                                               ic_lower(0), ic_upper(0),
                                               3)
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Spread V onto u at the points listed in indices using the 3-point
c     B-spline delta function using standard (double) precision
c     accumulation on the Cartesian grid.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_bspline_3_spread2d(
     &     dx,depth,
     &     indices,g,inear,ibrack,nindices,
     &     V,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u)
c
      implicit none
c
c     Functions.
c
      REAL lagrangian_bspline_3_delta
c
c     Input.
c
      INTEGER depth
      INTEGER nindices
      INTEGER ilower0,iupper0,ilower1,iupper1
      INTEGER nugc0,nugc1

      INTEGER indices(0:nindices-1)
      INTEGER inear(0:NDIM-1,0:nindices-1)
      INTEGER ibrack(0:NDIM-1,0:nindices-1)

      REAL g(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1)
      REAL V(0:depth-1,0:*)
c
c     Input/Output: spreading adds to the existing Cartesian grid values u.
c
      REAL u(CELL2dVECG(ilower,iupper,nugc),0:depth-1)
c
c     Local variables.
c
      INTEGER ilower(0:NDIM-1),iupper(0:NDIM-1)
      INTEGER ic,ic0,ic1
      INTEGER ic_center(0:NDIM-1),ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER d,l,s,nugc(0:NDIM-1)

      REAL w(0:NDIM-1,0:2)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1

      iupper(0) = iupper0
      iupper(1) = iupper1

      nugc(0) = nugc0
      nugc(1) = nugc1
c
c     Use a 3-point B-spline function to spread V onto u.
c
      do l = 0,nindices-1
         s = indices(l)
         do d=0,NDIM-1
c
c     Determine the Cartesian cell in which point s is located.
c
            ic_center(d) = inear(d,l)
c
c     Determine the spreading stencil corresponding to the position of
c     point s within the cell.
c
            ic_lower(d) = ic_center(d)-1
            ic_upper(d) = ic_center(d)+1
            ic_lower(d) = max(ic_lower(d),ilower(d)-nugc(d))
            ic_upper(d) = min(ic_upper(d),iupper(d)+nugc(d))
c
c     Compute the spreading weights.
c
            do ic = ic_lower(d),ic_upper(d)
               w(d,ic-ic_lower(d)) =
     &              lagrangian_bspline_3_delta(
     &              g(d,l)+dble(ibrack(d,l)-ic))
            enddo
         enddo
c
c     Spread V onto u.
c
         SPREAD_2D_SPECIALIZE_FIXED_WIDTH(ic_lower(1), ic_upper(1),
                                          ic_lower(0), ic_upper(0),
                                          3)
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Interpolate u onto V at the points listed in indices using the
c     4-point B-spline delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_bspline_4_interp2d(
     &     dx,depth,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u,
     &     indices,g,inear,ibrack,nindices,
     &     V)
c
      implicit none
c
c     Functions.
c
      REAL lagrangian_bspline_4_delta
c
c     Input.
c
      INTEGER depth
      INTEGER ilower0,iupper0,ilower1,iupper1
      INTEGER nugc0,nugc1
      INTEGER nindices

      INTEGER indices(0:nindices-1)
      INTEGER inear(0:NDIM-1,0:nindices-1)
      INTEGER ibrack(0:NDIM-1,0:nindices-1)

      REAL g(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1)
      REAL u(CELL2dVECG(ilower,iupper,nugc),0:depth-1)
c
c     Input/Output.
c
      REAL V(0:depth-1,0:*)
c
c     Local variables.
c
      INTEGER ilower(0:NDIM-1),iupper(0:NDIM-1)
      INTEGER ic,ic0,ic1
      INTEGER ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER d,l,s,nugc(0:NDIM-1)

      REAL w(0:NDIM-1,0:3)
c
c     Prevent compiler warning about unused variables.
c
      inear(0,0) = inear(0,0)
      dx(0) = dx(0)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1

      iupper(0) = iupper0
      iupper(1) = iupper1

      nugc(0) = nugc0
      nugc(1) = nugc1
c
c     Use a 4-point B-spline function to interpolate u onto V.
c
      do l = 0,nindices-1
         s = indices(l)
         do d=0,NDIM-1
c
c     Determine the interpolation stencil from the data point below s.
c
            ic_lower(d) = ibrack(d,l)-1
            ic_upper(d) = ibrack(d,l)+2
            ic_lower(d) = max(ic_lower(d),ilower(d)-nugc(d))
            ic_upper(d) = min(ic_upper(d),iupper(d)+nugc(d))
c
c     Compute the interpolation weights.
c
            do ic = ic_lower(d),ic_upper(d)
               w(d,ic-ic_lower(d)) =
     &              lagrangian_bspline_4_delta(
     &              g(d,l)+dble(ibrack(d,l)-ic))
            enddo
         enddo
c
c     Interpolate u onto V.
c
         INTERPOLATE_2D_SPECIALIZE_FIXED_WIDTH(ic_lower(1), ic_upper(1),
                                               ic_lower(0), ic_upper(0),
                                               4)
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Spread V onto u at the points listed in indices using the 4-point
c     B-spline delta function using standard (double) precision
c     accumulation on the Cartesian grid.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_bspline_4_spread2d(
     &     dx,depth,
     &     indices,g,inear,ibrack,nindices,
     &     V,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u)
c
      implicit none
c
c     Functions.
c
      REAL lagrangian_bspline_4_delta
c
c     Input.
c
      INTEGER depth
      INTEGER nindices
      INTEGER ilower0,iupper0,ilower1,iupper1
      INTEGER nugc0,nugc1

      INTEGER indices(0:nindices-1)
      INTEGER inear(0:NDIM-1,0:nindices-1)
      INTEGER ibrack(0:NDIM-1,0:nindices-1)

      REAL g(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1)
      REAL V(0:depth-1,0:*)
c
c     Input/Output: spreading adds to the existing Cartesian grid values u.
c
      REAL u(CELL2dVECG(ilower,iupper,nugc),0:depth-1)
c
c     Local variables.
c
      INTEGER ilower(0:NDIM-1),iupper(0:NDIM-1)
      INTEGER ic,ic0,ic1
      INTEGER ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER d,l,s,nugc(0:NDIM-1)

      REAL w(0:NDIM-1,0:3)
c
c     Prevent compiler warning about unused variables.
c
      inear(0,0) = inear(0,0)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1

      iupper(0) = iupper0
      iupper(1) = iupper1

      nugc(0) = nugc0
      nugc(1) = nugc1
c
c     Use a 4-point B-spline function to spread V onto u.
c
      do l = 0,nindices-1
         s = indices(l)
         do d=0,NDIM-1
c
c     Determine the spreading stencil from the data point below s.
c
            ic_lower(d) = ibrack(d,l)-1
            ic_upper(d) = ibrack(d,l)+2
            ic_lower(d) = max(ic_lower(d),ilower(d)-nugc(d))
            ic_upper(d) = min(ic_upper(d),iupper(d)+nugc(d))
c
c     Compute the spreading weights.
c
            do ic = ic_lower(d),ic_upper(d)
               w(d,ic-ic_lower(d)) =
     &              lagrangian_bspline_4_delta(
     &              g(d,l)+dble(ibrack(d,l)-ic))
            enddo
         enddo
c
c     Spread V onto u.
c
         SPREAD_2D_SPECIALIZE_FIXED_WIDTH(ic_lower(1), ic_upper(1),
                                          ic_lower(0), ic_upper(0),
                                          4)
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Interpolate u onto V at the points listed in indices using the
c     5-point B-spline delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_bspline_5_interp2d(
     &     dx,depth,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u,
     &     indices,g,inear,ibrack,nindices,
     &     V)
c
      implicit none
c
c     Functions.
c
      REAL lagrangian_bspline_5_delta
c
c     Input.
c
      INTEGER depth
      INTEGER ilower0,iupper0,ilower1,iupper1
      INTEGER nugc0,nugc1
      INTEGER nindices

      INTEGER indices(0:nindices-1)
      INTEGER inear(0:NDIM-1,0:nindices-1)
      INTEGER ibrack(0:NDIM-1,0:nindices-1)

      REAL g(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1)
      REAL u(CELL2dVECG(ilower,iupper,nugc),0:depth-1)
c
c     Input/Output.
c
      REAL V(0:depth-1,0:*)
c
c     Local variables.
c
      INTEGER ilower(0:NDIM-1),iupper(0:NDIM-1)
      INTEGER ic,ic0,ic1
      INTEGER ic_center(0:NDIM-1),ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER d,l,s,nugc(0:NDIM-1)

      REAL w(0:NDIM-1,0:4)
c
c     Prevent compiler warning about unused variables.
c
      dx(0) = dx(0)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1

      iupper(0) = iupper0
      iupper(1) = iupper1

      nugc(0) = nugc0
      nugc(1) = nugc1
c
c     Use a 5-point B-spline function to interpolate u onto V.
c
      do l = 0,nindices-1
         s = indices(l)
         do d=0,NDIM-1
c
c     Determine the Cartesian cell in which point s is located.
c
            ic_center(d) = inear(d,l)
c
c     Determine the interpolation stencil corresponding to the position
c     of point s within the cell.
c
            ic_lower(d) = ic_center(d)-2
            ic_upper(d) = ic_center(d)+2
            ic_lower(d) = max(ic_lower(d),ilower(d)-nugc(d))
            ic_upper(d) = min(ic_upper(d),iupper(d)+nugc(d))
c
c     Compute the interpolation weights.
c
            do ic = ic_lower(d),ic_upper(d)
               w(d,ic-ic_lower(d)) =
     &              lagrangian_bspline_5_delta(
     &              g(d,l)+dble(ibrack(d,l)-ic))
            enddo
         enddo
c
c     Interpolate u onto V.
c
         INTERPOLATE_2D_SPECIALIZE_FIXED_WIDTH(ic_lower(1), ic_upper(1),
                                               ic_lower(0), ic_upper(0),
                                               5)
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Spread V onto u at the points listed in indices using the 5-point
c     B-spline delta function using standard (double) precision
c     accumulation on the Cartesian grid.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_bspline_5_spread2d(
     &     dx,depth,
     &     indices,g,inear,ibrack,nindices,
     &     V,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u)
c
      implicit none
c
c     Functions.
c
      REAL lagrangian_bspline_5_delta
c
c     Input.
c
      INTEGER depth
      INTEGER nindices
      INTEGER ilower0,iupper0,ilower1,iupper1
      INTEGER nugc0,nugc1

      INTEGER indices(0:nindices-1)
      INTEGER inear(0:NDIM-1,0:nindices-1)
      INTEGER ibrack(0:NDIM-1,0:nindices-1)

      REAL g(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1)
      REAL V(0:depth-1,0:*)
c
c     Input/Output: spreading adds to the existing Cartesian grid values u.
c
      REAL u(CELL2dVECG(ilower,iupper,nugc),0:depth-1)
c
c     Local variables.
c
      INTEGER ilower(0:NDIM-1),iupper(0:NDIM-1)
      INTEGER ic,ic0,ic1
      INTEGER ic_center(0:NDIM-1),ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER d,l,s,nugc(0:NDIM-1)

      REAL w(0:NDIM-1,0:4)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1

      iupper(0) = iupper0
      iupper(1) = iupper1

      nugc(0) = nugc0
      nugc(1) = nugc1
c
c     Use a 5-point B-spline function to spread V onto u.
c
      do l = 0,nindices-1
         s = indices(l)
         do d = 0,NDIM-1
c
c     Determine the Cartesian cell in which point s is located.
c
            ic_center(d) = inear(d,l)
c
c     Determine the spreading stencil corresponding to the position of
c     point s within the cell.
c
            ic_lower(d) = ic_center(d)-2
            ic_upper(d) = ic_center(d)+2
            ic_lower(d) = max(ic_lower(d),ilower(d)-nugc(d))
            ic_upper(d) = min(ic_upper(d),iupper(d)+nugc(d))
c
c     Compute the spreading weights.
c
            do ic = ic_lower(d),ic_upper(d)
               w(d,ic-ic_lower(d)) =
     &              lagrangian_bspline_5_delta(
     &              g(d,l)+dble(ibrack(d,l)-ic))
            enddo
         enddo
c
c     Spread V onto u.
c
         SPREAD_2D_SPECIALIZE_FIXED_WIDTH(ic_lower(1), ic_upper(1),
                                          ic_lower(0), ic_upper(0),
                                          5)
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Interpolate u onto V at the points listed in indices using the
c     6-point B-spline delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_bspline_6_interp2d(
     &     dx,depth,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u,
     &     indices,g,inear,ibrack,nindices,
     &     V)
c
      implicit none
c
c     Functions.
c
      REAL lagrangian_bspline_6_delta
c
c     Input.
c
      INTEGER depth
      INTEGER ilower0,iupper0,ilower1,iupper1
      INTEGER nugc0,nugc1
      INTEGER nindices

      INTEGER indices(0:nindices-1)
      INTEGER inear(0:NDIM-1,0:nindices-1)
      INTEGER ibrack(0:NDIM-1,0:nindices-1)

      REAL g(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1)
      REAL u(CELL2dVECG(ilower,iupper,nugc),0:depth-1)
c
c     Input/Output.
c
      REAL V(0:depth-1,0:*)
c
c     Local variables.
c
      INTEGER ilower(0:NDIM-1),iupper(0:NDIM-1)
      INTEGER ic,ic0,ic1
      INTEGER ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER d,l,s,nugc(0:NDIM-1)

      REAL w(0:NDIM-1,0:5)
c
c     Prevent compiler warning about unused variables.
c
      inear(0,0) = inear(0,0)
      dx(0) = dx(0)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1

      iupper(0) = iupper0
      iupper(1) = iupper1

      nugc(0) = nugc0
      nugc(1) = nugc1
c
c     Use a 6-point B-spline function to interpolate u onto V.
c
      do l = 0,nindices-1
         s = indices(l)
         do d=0,NDIM-1
c
c     Determine the interpolation stencil from the data point below s.
c
            ic_lower(d) = ibrack(d,l)-2
            ic_upper(d) = ibrack(d,l)+3
            ic_lower(d) = max(ic_lower(d),ilower(d)-nugc(d))
            ic_upper(d) = min(ic_upper(d),iupper(d)+nugc(d))
c
c     Compute the interpolation weights.
c
            do ic = ic_lower(d),ic_upper(d)
               w(d,ic-ic_lower(d)) =
     &              lagrangian_bspline_6_delta(
     &              g(d,l)+dble(ibrack(d,l)-ic))
            enddo
         enddo
c
c     Interpolate u onto V.
c
         INTERPOLATE_2D_SPECIALIZE_FIXED_WIDTH(ic_lower(1), ic_upper(1),
                                               ic_lower(0), ic_upper(0),
                                               6)
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Spread V onto u at the points listed in indices using the 6-point
c     B-spline delta function using standard (double) precision
c     accumulation on the Cartesian grid.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_bspline_6_spread2d(
     &     dx,depth,
     &     indices,g,inear,ibrack,nindices,
     &     V,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u)
c
      implicit none
c
c     Functions.
c
      REAL lagrangian_bspline_6_delta
c
c     Input.
c
      INTEGER depth
      INTEGER nindices
      INTEGER ilower0,iupper0,ilower1,iupper1
      INTEGER nugc0,nugc1

      INTEGER indices(0:nindices-1)
      INTEGER inear(0:NDIM-1,0:nindices-1)
      INTEGER ibrack(0:NDIM-1,0:nindices-1)

      REAL g(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1)
      REAL V(0:depth-1,0:*)
c
c     Input/Output: spreading adds to the existing Cartesian grid values u.
c
      REAL u(CELL2dVECG(ilower,iupper,nugc),0:depth-1)
c
c     Local variables.
c
      INTEGER ilower(0:NDIM-1),iupper(0:NDIM-1)
      INTEGER ic,ic0,ic1
      INTEGER ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER d,l,s,nugc(0:NDIM-1)

      REAL w(0:NDIM-1,0:5)
c
c     Prevent compiler warning about unused variables.
c
      inear(0,0) = inear(0,0)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1

      iupper(0) = iupper0
      iupper(1) = iupper1

      nugc(0) = nugc0
      nugc(1) = nugc1
c
c     Use a 6-point B-spline function to spread V onto u.
c
      do l = 0,nindices-1
         s = indices(l)
         do d=0,NDIM-1
c
c     Determine the spreading stencil from the data point below s.
c
            ic_lower(d) = ibrack(d,l)-2
            ic_upper(d) = ibrack(d,l)+3
            ic_lower(d) = max(ic_lower(d),ilower(d)-nugc(d))
            ic_upper(d) = min(ic_upper(d),iupper(d)+nugc(d))
c
c     Compute the spreading weights.
c
            do ic = ic_lower(d),ic_upper(d)
               w(d,ic-ic_lower(d)) =
     &              lagrangian_bspline_6_delta(
     &              g(d,l)+dble(ibrack(d,l)-ic))
            enddo
         enddo
c
c     Spread V onto u.
c
         SPREAD_2D_SPECIALIZE_FIXED_WIDTH(ic_lower(1), ic_upper(1),
                                          ic_lower(0), ic_upper(0),
                                          6)
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Interpolate u onto V at the points listed in indices using the
c     composite 2-point/3-point B-spline delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_composite_bspline_23_interp2d(
     &     dx,depth,axis,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u,
     &     indices,g,inear,ibrack,nindices,
     &     V)
c
      implicit none
c
c     Functions.
c
      REAL lagrangian_piecewise_linear_delta
      REAL lagrangian_bspline_3_delta
c
c     Input.
c
      INTEGER depth,axis
      INTEGER ilower0,iupper0,ilower1,iupper1
      INTEGER nugc0,nugc1
      INTEGER nindices

      INTEGER indices(0:nindices-1)
      INTEGER inear(0:NDIM-1,0:nindices-1)
      INTEGER ibrack(0:NDIM-1,0:nindices-1)

      REAL g(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1)
      REAL u(CELL2dVECG(ilower,iupper,nugc),0:depth-1)
c
c     Input/Output.
c
      REAL V(0:depth-1,0:*)
c
c     Local variables.
c
      INTEGER ilower(0:NDIM-1),iupper(0:NDIM-1)
      INTEGER ic,ic0,ic1
      INTEGER ic_center(0:NDIM-1),ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER d,l,s,nugc(0:NDIM-1)

      REAL w(0:NDIM-1,0:2)
c
c     Prevent compiler warning about unused variables.
c
      dx(0) = dx(0)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1

      iupper(0) = iupper0
      iupper(1) = iupper1

      nugc(0) = nugc0
      nugc(1) = nugc1
c
c     Use a composite B-spline function to interpolate u onto V.
c
      do l = 0,nindices-1
         s = indices(l)
         do d=0,NDIM-1
c
c     Determine the Cartesian cell in which point s is located.
c
            ic_center(d) = inear(d,l)
c
c     Determine the interpolation stencil corresponding to the position
c     of point s within the cell.
c
            ic_lower(d) = ic_center(d)-1
            ic_upper(d) = ic_center(d)+1
            ic_lower(d) = max(ic_lower(d),ilower(d)-nugc(d))
            ic_upper(d) = min(ic_upper(d),iupper(d)+nugc(d))
c
c     Compute the interpolation weights.
c
            do ic = ic_lower(d),ic_upper(d)
               if (d .eq. axis) then
                  w(d,ic-ic_lower(d)) =
     &                 lagrangian_piecewise_linear_delta(
     &                 g(d,l)+dble(ibrack(d,l)-ic))
               else
                  w(d,ic-ic_lower(d)) =
     &                 lagrangian_bspline_3_delta(
     &                 g(d,l)+dble(ibrack(d,l)-ic))
               endif
            enddo
         enddo
c
c     Interpolate u onto V.
c
         INTERPOLATE_2D_SPECIALIZE_FIXED_WIDTH(ic_lower(1), ic_upper(1),
                                               ic_lower(0), ic_upper(0),
                                               3)
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Spread V onto u at the points listed in indices using the composite
c     2-point/3-point B-spline delta function using standard (double)
c     precision accumulation on the Cartesian grid.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_composite_bspline_23_spread2d(
     &     dx,depth,axis,
     &     indices,g,inear,ibrack,nindices,
     &     V,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u)
c
      implicit none
c
c     Functions.
c
      REAL lagrangian_piecewise_linear_delta
      REAL lagrangian_bspline_3_delta
c
c     Input.
c
      INTEGER depth,axis
      INTEGER nindices
      INTEGER ilower0,iupper0,ilower1,iupper1
      INTEGER nugc0,nugc1

      INTEGER indices(0:nindices-1)
      INTEGER inear(0:NDIM-1,0:nindices-1)
      INTEGER ibrack(0:NDIM-1,0:nindices-1)

      REAL g(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1)
      REAL V(0:depth-1,0:*)
c
c     Input/Output: spreading adds to the existing Cartesian grid values u.
c
      REAL u(CELL2dVECG(ilower,iupper,nugc),0:depth-1)
c
c     Local variables.
c
      INTEGER ilower(0:NDIM-1),iupper(0:NDIM-1)
      INTEGER ic,ic0,ic1
      INTEGER ic_center(0:NDIM-1),ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER d,l,s,nugc(0:NDIM-1)

      REAL w(0:NDIM-1,0:2)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1

      iupper(0) = iupper0
      iupper(1) = iupper1

      nugc(0) = nugc0
      nugc(1) = nugc1
c
c     Use a composite B-spline function to spread V onto u.
c
      do l = 0,nindices-1
         s = indices(l)
         do d=0,NDIM-1
c
c     Determine the Cartesian cell in which point s is located.
c
            ic_center(d) = inear(d,l)
c
c     Determine the spreading stencil corresponding to the position of
c     point s within the cell.
c
            ic_lower(d) = ic_center(d)-1
            ic_upper(d) = ic_center(d)+1
            ic_lower(d) = max(ic_lower(d),ilower(d)-nugc(d))
            ic_upper(d) = min(ic_upper(d),iupper(d)+nugc(d))
c
c     Compute the spreading weights.
c
            do ic = ic_lower(d),ic_upper(d)
               if (d .eq. axis) then
                  w(d,ic-ic_lower(d)) =
     &                 lagrangian_piecewise_linear_delta(
     &                 g(d,l)+dble(ibrack(d,l)-ic))
               else
                  w(d,ic-ic_lower(d)) =
     &                 lagrangian_bspline_3_delta(
     &                 g(d,l)+dble(ibrack(d,l)-ic))
               endif
            enddo
         enddo
c
c     Spread V onto u.
c
         SPREAD_2D_SPECIALIZE_FIXED_WIDTH(ic_lower(1), ic_upper(1),
                                          ic_lower(0), ic_upper(0),
                                          3)
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Interpolate u onto V at the points listed in indices using the
c     composite 3-point/2-point B-spline delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_composite_bspline_32_interp2d(
     &     dx,depth,axis,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u,
     &     indices,g,inear,ibrack,nindices,
     &     V)
c
      implicit none
c
c     Functions.
c
      REAL lagrangian_piecewise_linear_delta
      REAL lagrangian_bspline_3_delta
c
c     Input.
c
      INTEGER depth,axis
      INTEGER ilower0,iupper0,ilower1,iupper1
      INTEGER nugc0,nugc1
      INTEGER nindices

      INTEGER indices(0:nindices-1)
      INTEGER inear(0:NDIM-1,0:nindices-1)
      INTEGER ibrack(0:NDIM-1,0:nindices-1)

      REAL g(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1)
      REAL u(CELL2dVECG(ilower,iupper,nugc),0:depth-1)
c
c     Input/Output.
c
      REAL V(0:depth-1,0:*)
c
c     Local variables.
c
      INTEGER ilower(0:NDIM-1),iupper(0:NDIM-1)
      INTEGER ic,ic0,ic1
      INTEGER ic_center(0:NDIM-1),ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER d,l,s,nugc(0:NDIM-1)

      REAL w(0:NDIM-1,0:2)
c
c     Prevent compiler warning about unused variables.
c
      dx(0) = dx(0)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1

      iupper(0) = iupper0
      iupper(1) = iupper1

      nugc(0) = nugc0
      nugc(1) = nugc1
c
c     Use a composite B-spline function to interpolate u onto V.
c
      do l = 0,nindices-1
         s = indices(l)
         do d=0,NDIM-1
c
c     Determine the Cartesian cell in which point s is located.
c
            ic_center(d) = inear(d,l)
c
c     Determine the interpolation stencil corresponding to the position
c     of point s within the cell.
c
            ic_lower(d) = ic_center(d)-1
            ic_upper(d) = ic_center(d)+1
            ic_lower(d) = max(ic_lower(d),ilower(d)-nugc(d))
            ic_upper(d) = min(ic_upper(d),iupper(d)+nugc(d))
c
c     Compute the interpolation weights.
c
            do ic = ic_lower(d),ic_upper(d)
               if (d .eq. axis) then
                  w(d,ic-ic_lower(d)) =
     &                 lagrangian_bspline_3_delta(
     &                 g(d,l)+dble(ibrack(d,l)-ic))
               else
                  w(d,ic-ic_lower(d)) =
     &                 lagrangian_piecewise_linear_delta(
     &                 g(d,l)+dble(ibrack(d,l)-ic))
               endif
            enddo
         enddo
c
c     Interpolate u onto V.
c
         INTERPOLATE_2D_SPECIALIZE_FIXED_WIDTH(ic_lower(1), ic_upper(1),
                                               ic_lower(0), ic_upper(0),
                                               3)
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Spread V onto u at the points listed in indices using the composite
c     3-point/2-point B-spline delta function using standard (double)
c     precision accumulation on the Cartesian grid.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_composite_bspline_32_spread2d(
     &     dx,depth,axis,
     &     indices,g,inear,ibrack,nindices,
     &     V,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u)
c
      implicit none
c
c     Functions.
c
      REAL lagrangian_piecewise_linear_delta
      REAL lagrangian_bspline_3_delta
c
c     Input.
c
      INTEGER depth,axis
      INTEGER nindices
      INTEGER ilower0,iupper0,ilower1,iupper1
      INTEGER nugc0,nugc1

      INTEGER indices(0:nindices-1)
      INTEGER inear(0:NDIM-1,0:nindices-1)
      INTEGER ibrack(0:NDIM-1,0:nindices-1)

      REAL g(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1)
      REAL V(0:depth-1,0:*)
c
c     Input/Output: spreading adds to the existing Cartesian grid values u.
c
      REAL u(CELL2dVECG(ilower,iupper,nugc),0:depth-1)
c
c     Local variables.
c
      INTEGER ilower(0:NDIM-1),iupper(0:NDIM-1)
      INTEGER ic,ic0,ic1
      INTEGER ic_center(0:NDIM-1),ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER d,l,s,nugc(0:NDIM-1)

      REAL w(0:NDIM-1,0:2)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1

      iupper(0) = iupper0
      iupper(1) = iupper1

      nugc(0) = nugc0
      nugc(1) = nugc1
c
c     Use a composite B-spline function to spread V onto u.
c
      do l = 0,nindices-1
         s = indices(l)
         do d=0,NDIM-1
c
c     Determine the Cartesian cell in which point s is located.
c
            ic_center(d) = inear(d,l)
c
c     Determine the spreading stencil corresponding to the position of
c     point s within the cell.
c
            ic_lower(d) = ic_center(d)-1
            ic_upper(d) = ic_center(d)+1
            ic_lower(d) = max(ic_lower(d),ilower(d)-nugc(d))
            ic_upper(d) = min(ic_upper(d),iupper(d)+nugc(d))
c
c     Compute the spreading weights.
c
            do ic = ic_lower(d),ic_upper(d)
               if (d .eq. axis) then
                  w(d,ic-ic_lower(d)) =
     &                 lagrangian_bspline_3_delta(
     &                 g(d,l)+dble(ibrack(d,l)-ic))
               else
                  w(d,ic-ic_lower(d)) =
     &                 lagrangian_piecewise_linear_delta(
     &                 g(d,l)+dble(ibrack(d,l)-ic))
               endif
            enddo
         enddo
c
c     Spread V onto u.
c
         SPREAD_2D_SPECIALIZE_FIXED_WIDTH(ic_lower(1), ic_upper(1),
                                          ic_lower(0), ic_upper(0),
                                          3)
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Interpolate u onto V at the points listed in indices using the
c     composite 3-point/4-point B-spline delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_composite_bspline_34_interp2d(
     &     dx,depth,axis,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u,
     &     indices,g,inear,ibrack,nindices,
     &     V)
c
      implicit none
c
c     Functions.
c
      REAL lagrangian_bspline_3_delta
      REAL lagrangian_bspline_4_delta
c
c     Input.
c
      INTEGER depth,axis
      INTEGER ilower0,iupper0,ilower1,iupper1
      INTEGER nugc0,nugc1
      INTEGER nindices

      INTEGER indices(0:nindices-1)
      INTEGER inear(0:NDIM-1,0:nindices-1)
      INTEGER ibrack(0:NDIM-1,0:nindices-1)

      REAL g(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1)
      REAL u(CELL2dVECG(ilower,iupper,nugc),0:depth-1)
c
c     Input/Output.
c
      REAL V(0:depth-1,0:*)
c
c     Local variables.
c
      INTEGER ilower(0:NDIM-1),iupper(0:NDIM-1)
      INTEGER ic,ic0,ic1
      INTEGER ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER d,l,s,nugc(0:NDIM-1)

      REAL w(0:NDIM-1,0:3)
c
c     Prevent compiler warning about unused variables.
c
      inear(0,0) = inear(0,0)
      dx(0) = dx(0)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1

      iupper(0) = iupper0
      iupper(1) = iupper1

      nugc(0) = nugc0
      nugc(1) = nugc1
c
c     Use a composite B-spline function to interpolate u onto V.
c
      do l = 0,nindices-1
         s = indices(l)
         do d=0,NDIM-1
c
c     Determine the interpolation stencil from the data point below s.
c
            ic_lower(d) = ibrack(d,l)-1
            ic_upper(d) = ibrack(d,l)+2
            ic_lower(d) = max(ic_lower(d),ilower(d)-nugc(d))
            ic_upper(d) = min(ic_upper(d),iupper(d)+nugc(d))
c
c     Compute the interpolation weights.
c
            do ic = ic_lower(d),ic_upper(d)
               if (d .eq. axis) then
                  w(d,ic-ic_lower(d)) =
     &                 lagrangian_bspline_3_delta(
     &                 g(d,l)+dble(ibrack(d,l)-ic))
               else
                  w(d,ic-ic_lower(d)) =
     &                 lagrangian_bspline_4_delta(
     &                 g(d,l)+dble(ibrack(d,l)-ic))
               endif
            enddo
         enddo
c
c     Interpolate u onto V.
c
         INTERPOLATE_2D_SPECIALIZE_FIXED_WIDTH(ic_lower(1), ic_upper(1),
                                               ic_lower(0), ic_upper(0),
                                               4)
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Spread V onto u at the points listed in indices using the composite
c     3-point/4-point B-spline delta function using standard (double)
c     precision accumulation on the Cartesian grid.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_composite_bspline_34_spread2d(
     &     dx,depth,axis,
     &     indices,g,inear,ibrack,nindices,
     &     V,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u)
c
      implicit none
c
c     Functions.
c
      REAL lagrangian_bspline_3_delta
      REAL lagrangian_bspline_4_delta
c
c     Input.
c
      INTEGER depth,axis
      INTEGER nindices
      INTEGER ilower0,iupper0,ilower1,iupper1
      INTEGER nugc0,nugc1

      INTEGER indices(0:nindices-1)
      INTEGER inear(0:NDIM-1,0:nindices-1)
      INTEGER ibrack(0:NDIM-1,0:nindices-1)

      REAL g(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1)
      REAL V(0:depth-1,0:*)
c
c     Input/Output: spreading adds to the existing Cartesian grid values u.
c
      REAL u(CELL2dVECG(ilower,iupper,nugc),0:depth-1)
c
c     Local variables.
c
      INTEGER ilower(0:NDIM-1),iupper(0:NDIM-1)
      INTEGER ic,ic0,ic1
      INTEGER ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER d,l,s,nugc(0:NDIM-1)

      REAL w(0:NDIM-1,0:3)
c
c     Prevent compiler warning about unused variables.
c
      inear(0,0) = inear(0,0)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1

      iupper(0) = iupper0
      iupper(1) = iupper1

      nugc(0) = nugc0
      nugc(1) = nugc1
c
c     Use a composite B-spline function to spread V onto u.
c
      do l = 0,nindices-1
         s = indices(l)
         do d=0,NDIM-1
c
c     Determine the spreading stencil from the data point below s.
c
            ic_lower(d) = ibrack(d,l)-1
            ic_upper(d) = ibrack(d,l)+2
            ic_lower(d) = max(ic_lower(d),ilower(d)-nugc(d))
            ic_upper(d) = min(ic_upper(d),iupper(d)+nugc(d))
c
c     Compute the spreading weights.
c
            do ic = ic_lower(d),ic_upper(d)
               if (d .eq. axis) then
                  w(d,ic-ic_lower(d)) =
     &                 lagrangian_bspline_3_delta(
     &                 g(d,l)+dble(ibrack(d,l)-ic))
               else
                  w(d,ic-ic_lower(d)) =
     &                 lagrangian_bspline_4_delta(
     &                 g(d,l)+dble(ibrack(d,l)-ic))
               endif
            enddo
         enddo
c
c     Spread V onto u.
c
         SPREAD_2D_SPECIALIZE_FIXED_WIDTH(ic_lower(1), ic_upper(1),
                                          ic_lower(0), ic_upper(0),
                                          4)
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Interpolate u onto V at the points listed in indices using the
c     composite 4-point/3-point B-spline delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_composite_bspline_43_interp2d(
     &     dx,depth,axis,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u,
     &     indices,g,inear,ibrack,nindices,
     &     V)
c
      implicit none
c
c     Functions.
c
      REAL lagrangian_bspline_3_delta
      REAL lagrangian_bspline_4_delta
c
c     Input.
c
      INTEGER depth,axis
      INTEGER ilower0,iupper0,ilower1,iupper1
      INTEGER nugc0,nugc1
      INTEGER nindices

      INTEGER indices(0:nindices-1)
      INTEGER inear(0:NDIM-1,0:nindices-1)
      INTEGER ibrack(0:NDIM-1,0:nindices-1)

      REAL g(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1)
      REAL u(CELL2dVECG(ilower,iupper,nugc),0:depth-1)
c
c     Input/Output.
c
      REAL V(0:depth-1,0:*)
c
c     Local variables.
c
      INTEGER ilower(0:NDIM-1),iupper(0:NDIM-1)
      INTEGER ic,ic0,ic1
      INTEGER ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER d,l,s,nugc(0:NDIM-1)

      REAL w(0:NDIM-1,0:3)
c
c     Prevent compiler warning about unused variables.
c
      inear(0,0) = inear(0,0)
      dx(0) = dx(0)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1

      iupper(0) = iupper0
      iupper(1) = iupper1

      nugc(0) = nugc0
      nugc(1) = nugc1
c
c     Use a composite B-spline function to interpolate u onto V.
c
      do l = 0,nindices-1
         s = indices(l)
         do d=0,NDIM-1
c
c     Determine the interpolation stencil from the data point below s.
c
            ic_lower(d) = ibrack(d,l)-1
            ic_upper(d) = ibrack(d,l)+2
            ic_lower(d) = max(ic_lower(d),ilower(d)-nugc(d))
            ic_upper(d) = min(ic_upper(d),iupper(d)+nugc(d))
c
c     Compute the interpolation weights.
c
            do ic = ic_lower(d),ic_upper(d)
               if (d .eq. axis) then
                  w(d,ic-ic_lower(d)) =
     &                 lagrangian_bspline_4_delta(
     &                 g(d,l)+dble(ibrack(d,l)-ic))
               else
                  w(d,ic-ic_lower(d)) =
     &                 lagrangian_bspline_3_delta(
     &                 g(d,l)+dble(ibrack(d,l)-ic))
               endif
            enddo
         enddo
c
c     Interpolate u onto V.
c
         INTERPOLATE_2D_SPECIALIZE_FIXED_WIDTH(ic_lower(1), ic_upper(1),
                                               ic_lower(0), ic_upper(0),
                                               4)
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Spread V onto u at the points listed in indices using the composite
c     4-point/3-point B-spline delta function using standard (double)
c     precision accumulation on the Cartesian grid.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_composite_bspline_43_spread2d(
     &     dx,depth,axis,
     &     indices,g,inear,ibrack,nindices,
     &     V,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u)
c
      implicit none
c
c     Functions.
c
      REAL lagrangian_bspline_3_delta
      REAL lagrangian_bspline_4_delta
c
c     Input.
c
      INTEGER depth,axis
      INTEGER nindices
      INTEGER ilower0,iupper0,ilower1,iupper1
      INTEGER nugc0,nugc1

      INTEGER indices(0:nindices-1)
      INTEGER inear(0:NDIM-1,0:nindices-1)
      INTEGER ibrack(0:NDIM-1,0:nindices-1)

      REAL g(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1)
      REAL V(0:depth-1,0:*)
c
c     Input/Output: spreading adds to the existing Cartesian grid values u.
c
      REAL u(CELL2dVECG(ilower,iupper,nugc),0:depth-1)
c
c     Local variables.
c
      INTEGER ilower(0:NDIM-1),iupper(0:NDIM-1)
      INTEGER ic,ic0,ic1
      INTEGER ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER d,l,s,nugc(0:NDIM-1)

      REAL w(0:NDIM-1,0:3)
c
c     Prevent compiler warning about unused variables.
c
      inear(0,0) = inear(0,0)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1

      iupper(0) = iupper0
      iupper(1) = iupper1

      nugc(0) = nugc0
      nugc(1) = nugc1
c
c     Use a composite B-spline function to spread V onto u.
c
      do l = 0,nindices-1
         s = indices(l)
         do d=0,NDIM-1
c
c     Determine the spreading stencil from the data point below s.
c
            ic_lower(d) = ibrack(d,l)-1
            ic_upper(d) = ibrack(d,l)+2
            ic_lower(d) = max(ic_lower(d),ilower(d)-nugc(d))
            ic_upper(d) = min(ic_upper(d),iupper(d)+nugc(d))
c
c     Compute the spreading weights.
c
            do ic = ic_lower(d),ic_upper(d)
               if (d .eq. axis) then
                  w(d,ic-ic_lower(d)) =
     &                 lagrangian_bspline_4_delta(
     &                 g(d,l)+dble(ibrack(d,l)-ic))
               else
                  w(d,ic-ic_lower(d)) =
     &                 lagrangian_bspline_3_delta(
     &                 g(d,l)+dble(ibrack(d,l)-ic))
               endif
            enddo
         enddo
c
c     Spread V onto u.
c
         SPREAD_2D_SPECIALIZE_FIXED_WIDTH(ic_lower(1), ic_upper(1),
                                          ic_lower(0), ic_upper(0),
                                          4)
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Interpolate u onto V at the points listed in indices using the
c     composite 4-point/5-point B-spline delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_composite_bspline_45_interp2d(
     &     dx,depth,axis,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u,
     &     indices,g,inear,ibrack,nindices,
     &     V)
c
      implicit none
c
c     Functions.
c
      REAL lagrangian_bspline_4_delta
      REAL lagrangian_bspline_5_delta
c
c     Input.
c
      INTEGER depth,axis
      INTEGER ilower0,iupper0,ilower1,iupper1
      INTEGER nugc0,nugc1
      INTEGER nindices

      INTEGER indices(0:nindices-1)
      INTEGER inear(0:NDIM-1,0:nindices-1)
      INTEGER ibrack(0:NDIM-1,0:nindices-1)

      REAL g(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1)
      REAL u(CELL2dVECG(ilower,iupper,nugc),0:depth-1)
c
c     Input/Output.
c
      REAL V(0:depth-1,0:*)
c
c     Local variables.
c
      INTEGER ilower(0:NDIM-1),iupper(0:NDIM-1)
      INTEGER ic,ic0,ic1
      INTEGER ic_center(0:NDIM-1),ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER d,l,s,nugc(0:NDIM-1)

      REAL w(0:NDIM-1,0:4)
c
c     Prevent compiler warning about unused variables.
c
      dx(0) = dx(0)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1

      iupper(0) = iupper0
      iupper(1) = iupper1

      nugc(0) = nugc0
      nugc(1) = nugc1
c
c     Use a composite B-spline function to interpolate u onto V.
c
      do l = 0,nindices-1
         s = indices(l)
         do d=0,NDIM-1
c
c     Determine the Cartesian cell in which point s is located.
c
            ic_center(d) = inear(d,l)
c
c     Determine the interpolation stencil corresponding to the position
c     of point s within the cell.
c
            ic_lower(d) = ic_center(d)-2
            ic_upper(d) = ic_center(d)+2
            ic_lower(d) = max(ic_lower(d),ilower(d)-nugc(d))
            ic_upper(d) = min(ic_upper(d),iupper(d)+nugc(d))
c
c     Compute the interpolation weights.
c
            do ic = ic_lower(d),ic_upper(d)
               if (d .eq. axis) then
                  w(d,ic-ic_lower(d)) =
     &                 lagrangian_bspline_4_delta(
     &                 g(d,l)+dble(ibrack(d,l)-ic))
               else
                  w(d,ic-ic_lower(d)) =
     &                 lagrangian_bspline_5_delta(
     &                 g(d,l)+dble(ibrack(d,l)-ic))
               endif
            enddo
         enddo
c
c     Interpolate u onto V.
c
         INTERPOLATE_2D_SPECIALIZE_FIXED_WIDTH(ic_lower(1), ic_upper(1),
                                               ic_lower(0), ic_upper(0),
                                               5)
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Spread V onto u at the points listed in indices using the composite
c     4-point/5-point B-spline delta function using standard (double)
c     precision accumulation on the Cartesian grid.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_composite_bspline_45_spread2d(
     &     dx,depth,axis,
     &     indices,g,inear,ibrack,nindices,
     &     V,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u)
c
      implicit none
c
c     Functions.
c
      REAL lagrangian_bspline_4_delta
      REAL lagrangian_bspline_5_delta
c
c     Input.
c
      INTEGER depth,axis
      INTEGER nindices
      INTEGER ilower0,iupper0,ilower1,iupper1
      INTEGER nugc0,nugc1

      INTEGER indices(0:nindices-1)
      INTEGER inear(0:NDIM-1,0:nindices-1)
      INTEGER ibrack(0:NDIM-1,0:nindices-1)

      REAL g(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1)
      REAL V(0:depth-1,0:*)
c
c     Input/Output: spreading adds to the existing Cartesian grid values u.
c
      REAL u(CELL2dVECG(ilower,iupper,nugc),0:depth-1)
c
c     Local variables.
c
      INTEGER ilower(0:NDIM-1),iupper(0:NDIM-1)
      INTEGER ic,ic0,ic1
      INTEGER ic_center(0:NDIM-1),ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER d,l,s,nugc(0:NDIM-1)

      REAL w(0:NDIM-1,0:4)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1

      iupper(0) = iupper0
      iupper(1) = iupper1

      nugc(0) = nugc0
      nugc(1) = nugc1
c
c     Use a composite B-spline function to spread V onto u.
c
      do l = 0,nindices-1
         s = indices(l)
         do d = 0,NDIM-1
c
c     Determine the Cartesian cell in which point s is located.
c
            ic_center(d) = inear(d,l)
c
c     Determine the spreading stencil corresponding to the position of
c     point s within the cell.
c
            ic_lower(d) = ic_center(d)-2
            ic_upper(d) = ic_center(d)+2
            ic_lower(d) = max(ic_lower(d),ilower(d)-nugc(d))
            ic_upper(d) = min(ic_upper(d),iupper(d)+nugc(d))
c
c     Compute the spreading weights.
c
            do ic = ic_lower(d),ic_upper(d)
               if (d .eq. axis) then
                  w(d,ic-ic_lower(d)) =
     &                 lagrangian_bspline_4_delta(
     &                 g(d,l)+dble(ibrack(d,l)-ic))
               else
                  w(d,ic-ic_lower(d)) =
     &                 lagrangian_bspline_5_delta(
     &                 g(d,l)+dble(ibrack(d,l)-ic))
               endif
            enddo
         enddo
c
c     Spread V onto u.
c
         SPREAD_2D_SPECIALIZE_FIXED_WIDTH(ic_lower(1), ic_upper(1),
                                          ic_lower(0), ic_upper(0),
                                          5)
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Interpolate u onto V at the points listed in indices using the
c     composite 5-point/4-point B-spline delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_composite_bspline_54_interp2d(
     &     dx,depth,axis,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u,
     &     indices,g,inear,ibrack,nindices,
     &     V)
c
      implicit none
c
c     Functions.
c
      REAL lagrangian_bspline_4_delta
      REAL lagrangian_bspline_5_delta
c
c     Input.
c
      INTEGER depth,axis
      INTEGER ilower0,iupper0,ilower1,iupper1
      INTEGER nugc0,nugc1
      INTEGER nindices

      INTEGER indices(0:nindices-1)
      INTEGER inear(0:NDIM-1,0:nindices-1)
      INTEGER ibrack(0:NDIM-1,0:nindices-1)

      REAL g(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1)
      REAL u(CELL2dVECG(ilower,iupper,nugc),0:depth-1)
c
c     Input/Output.
c
      REAL V(0:depth-1,0:*)
c
c     Local variables.
c
      INTEGER ilower(0:NDIM-1),iupper(0:NDIM-1)
      INTEGER ic,ic0,ic1
      INTEGER ic_center(0:NDIM-1),ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER d,l,s,nugc(0:NDIM-1)

      REAL w(0:NDIM-1,0:4)
c
c     Prevent compiler warning about unused variables.
c
      dx(0) = dx(0)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1

      iupper(0) = iupper0
      iupper(1) = iupper1

      nugc(0) = nugc0
      nugc(1) = nugc1
c
c     Use a composite B-spline function to interpolate u onto V.
c
      do l = 0,nindices-1
         s = indices(l)
         do d=0,NDIM-1
c
c     Determine the Cartesian cell in which point s is located.
c
            ic_center(d) = inear(d,l)
c
c     Determine the interpolation stencil corresponding to the position
c     of point s within the cell.
c
            ic_lower(d) = ic_center(d)-2
            ic_upper(d) = ic_center(d)+2
            ic_lower(d) = max(ic_lower(d),ilower(d)-nugc(d))
            ic_upper(d) = min(ic_upper(d),iupper(d)+nugc(d))
c
c     Compute the interpolation weights.
c
            do ic = ic_lower(d),ic_upper(d)
               if (d .eq. axis) then
                  w(d,ic-ic_lower(d)) =
     &                 lagrangian_bspline_5_delta(
     &                 g(d,l)+dble(ibrack(d,l)-ic))
               else
                  w(d,ic-ic_lower(d)) =
     &                 lagrangian_bspline_4_delta(
     &                 g(d,l)+dble(ibrack(d,l)-ic))
               endif
            enddo
         enddo
c
c     Interpolate u onto V.
c
         INTERPOLATE_2D_SPECIALIZE_FIXED_WIDTH(ic_lower(1), ic_upper(1),
                                               ic_lower(0), ic_upper(0),
                                               5)
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Spread V onto u at the points listed in indices using the composite
c     5-point/4-point B-spline delta function using standard (double)
c     precision accumulation on the Cartesian grid.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_composite_bspline_54_spread2d(
     &     dx,depth,axis,
     &     indices,g,inear,ibrack,nindices,
     &     V,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u)
c
      implicit none
c
c     Functions.
c
      REAL lagrangian_bspline_4_delta
      REAL lagrangian_bspline_5_delta
c
c     Input.
c
      INTEGER depth,axis
      INTEGER nindices
      INTEGER ilower0,iupper0,ilower1,iupper1
      INTEGER nugc0,nugc1

      INTEGER indices(0:nindices-1)
      INTEGER inear(0:NDIM-1,0:nindices-1)
      INTEGER ibrack(0:NDIM-1,0:nindices-1)

      REAL g(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1)
      REAL V(0:depth-1,0:*)
c
c     Input/Output: spreading adds to the existing Cartesian grid values u.
c
      REAL u(CELL2dVECG(ilower,iupper,nugc),0:depth-1)
c
c     Local variables.
c
      INTEGER ilower(0:NDIM-1),iupper(0:NDIM-1)
      INTEGER ic,ic0,ic1
      INTEGER ic_center(0:NDIM-1),ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER d,l,s,nugc(0:NDIM-1)

      REAL w(0:NDIM-1,0:4)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1

      iupper(0) = iupper0
      iupper(1) = iupper1

      nugc(0) = nugc0
      nugc(1) = nugc1
c
c     Use a composite B-spline function to spread V onto u.
c
      do l = 0,nindices-1
         s = indices(l)
         do d = 0,NDIM-1
c
c     Determine the Cartesian cell in which point s is located.
c
            ic_center(d) = inear(d,l)
c
c     Determine the spreading stencil corresponding to the position of
c     point s within the cell.
c
            ic_lower(d) = ic_center(d)-2
            ic_upper(d) = ic_center(d)+2
            ic_lower(d) = max(ic_lower(d),ilower(d)-nugc(d))
            ic_upper(d) = min(ic_upper(d),iupper(d)+nugc(d))
c
c     Compute the spreading weights.
c
            do ic = ic_lower(d),ic_upper(d)
               if (d .eq. axis) then
                  w(d,ic-ic_lower(d)) =
     &                 lagrangian_bspline_5_delta(
     &                 g(d,l)+dble(ibrack(d,l)-ic))
               else
                  w(d,ic-ic_lower(d)) =
     &                 lagrangian_bspline_4_delta(
     &                 g(d,l)+dble(ibrack(d,l)-ic))
               endif
            enddo
         enddo
c
c     Spread V onto u.
c
         SPREAD_2D_SPECIALIZE_FIXED_WIDTH(ic_lower(1), ic_upper(1),
                                          ic_lower(0), ic_upper(0),
                                          5)
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Interpolate u onto V at the points listed in indices using the
c     composite 5-point/6-point B-spline delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_composite_bspline_56_interp2d(
     &     dx,depth,axis,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u,
     &     indices,g,inear,ibrack,nindices,
     &     V)
c
      implicit none
c
c     Functions.
c
      REAL lagrangian_bspline_5_delta
      REAL lagrangian_bspline_6_delta
c
c     Input.
c
      INTEGER depth,axis
      INTEGER ilower0,iupper0,ilower1,iupper1
      INTEGER nugc0,nugc1
      INTEGER nindices

      INTEGER indices(0:nindices-1)
      INTEGER inear(0:NDIM-1,0:nindices-1)
      INTEGER ibrack(0:NDIM-1,0:nindices-1)

      REAL g(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1)
      REAL u(CELL2dVECG(ilower,iupper,nugc),0:depth-1)
c
c     Input/Output.
c
      REAL V(0:depth-1,0:*)
c
c     Local variables.
c
      INTEGER ilower(0:NDIM-1),iupper(0:NDIM-1)
      INTEGER ic,ic0,ic1
      INTEGER ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER d,l,s,nugc(0:NDIM-1)

      REAL w(0:NDIM-1,0:5)
c
c     Prevent compiler warning about unused variables.
c
      inear(0,0) = inear(0,0)
      dx(0) = dx(0)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1

      iupper(0) = iupper0
      iupper(1) = iupper1

      nugc(0) = nugc0
      nugc(1) = nugc1
c
c     Use a composite B-spline function to interpolate u onto V.
c
      do l = 0,nindices-1
         s = indices(l)
         do d=0,NDIM-1
c
c     Determine the interpolation stencil from the data point below s.
c
            ic_lower(d) = ibrack(d,l)-2
            ic_upper(d) = ibrack(d,l)+3
            ic_lower(d) = max(ic_lower(d),ilower(d)-nugc(d))
            ic_upper(d) = min(ic_upper(d),iupper(d)+nugc(d))
c
c     Compute the interpolation weights.
c
            do ic = ic_lower(d),ic_upper(d)
               if (d .eq. axis) then
                  w(d,ic-ic_lower(d)) =
     &                 lagrangian_bspline_5_delta(
     &                 g(d,l)+dble(ibrack(d,l)-ic))
               else
                  w(d,ic-ic_lower(d)) =
     &                 lagrangian_bspline_6_delta(
     &                 g(d,l)+dble(ibrack(d,l)-ic))
               endif
            enddo
         enddo
c
c     Interpolate u onto V.
c
         INTERPOLATE_2D_SPECIALIZE_FIXED_WIDTH(ic_lower(1), ic_upper(1),
                                               ic_lower(0), ic_upper(0),
                                               6)
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Spread V onto u at the points listed in indices using the composite
c     5-point/6-point B-spline delta function using standard (double)
c     precision accumulation on the Cartesian grid.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_composite_bspline_56_spread2d(
     &     dx,depth,axis,
     &     indices,g,inear,ibrack,nindices,
     &     V,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u)
c
      implicit none
c
c     Functions.
c
      REAL lagrangian_bspline_5_delta
      REAL lagrangian_bspline_6_delta
c
c     Input.
c
      INTEGER depth,axis
      INTEGER nindices
      INTEGER ilower0,iupper0,ilower1,iupper1
      INTEGER nugc0,nugc1

      INTEGER indices(0:nindices-1)
      INTEGER inear(0:NDIM-1,0:nindices-1)
      INTEGER ibrack(0:NDIM-1,0:nindices-1)

      REAL g(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1)
      REAL V(0:depth-1,0:*)
c
c     Input/Output: spreading adds to the existing Cartesian grid values u.
c
      REAL u(CELL2dVECG(ilower,iupper,nugc),0:depth-1)
c
c     Local variables.
c
      INTEGER ilower(0:NDIM-1),iupper(0:NDIM-1)
      INTEGER ic,ic0,ic1
      INTEGER ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER d,l,s,nugc(0:NDIM-1)

      REAL w(0:NDIM-1,0:5)
c
c     Prevent compiler warning about unused variables.
c
      inear(0,0) = inear(0,0)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1

      iupper(0) = iupper0
      iupper(1) = iupper1

      nugc(0) = nugc0
      nugc(1) = nugc1
c
c     Use a composite B-spline function to spread V onto u.
c
      do l = 0,nindices-1
         s = indices(l)
         do d=0,NDIM-1
c
c     Determine the spreading stencil from the data point below s.
c
            ic_lower(d) = ibrack(d,l)-2
            ic_upper(d) = ibrack(d,l)+3
            ic_lower(d) = max(ic_lower(d),ilower(d)-nugc(d))
            ic_upper(d) = min(ic_upper(d),iupper(d)+nugc(d))
c
c     Compute the spreading weights.
c
            do ic = ic_lower(d),ic_upper(d)
               if (d .eq. axis) then
                  w(d,ic-ic_lower(d)) =
     &                 lagrangian_bspline_5_delta(
     &                 g(d,l)+dble(ibrack(d,l)-ic))
               else
                  w(d,ic-ic_lower(d)) =
     &                 lagrangian_bspline_6_delta(
     &                 g(d,l)+dble(ibrack(d,l)-ic))
               endif
            enddo
         enddo
c
c     Spread V onto u.
c
         SPREAD_2D_SPECIALIZE_FIXED_WIDTH(ic_lower(1), ic_upper(1),
                                          ic_lower(0), ic_upper(0),
                                          6)
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Interpolate u onto V at the points listed in indices using the
c     composite 6-point/5-point B-spline delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_composite_bspline_65_interp2d(
     &     dx,depth,axis,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u,
     &     indices,g,inear,ibrack,nindices,
     &     V)
c
      implicit none
c
c     Functions.
c
      REAL lagrangian_bspline_5_delta
      REAL lagrangian_bspline_6_delta
c
c     Input.
c
      INTEGER depth,axis
      INTEGER ilower0,iupper0,ilower1,iupper1
      INTEGER nugc0,nugc1
      INTEGER nindices

      INTEGER indices(0:nindices-1)
      INTEGER inear(0:NDIM-1,0:nindices-1)
      INTEGER ibrack(0:NDIM-1,0:nindices-1)

      REAL g(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1)
      REAL u(CELL2dVECG(ilower,iupper,nugc),0:depth-1)
c
c     Input/Output.
c
      REAL V(0:depth-1,0:*)
c
c     Local variables.
c
      INTEGER ilower(0:NDIM-1),iupper(0:NDIM-1)
      INTEGER ic,ic0,ic1
      INTEGER ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER d,l,s,nugc(0:NDIM-1)

      REAL w(0:NDIM-1,0:5)
c
c     Prevent compiler warning about unused variables.
c
      inear(0,0) = inear(0,0)
      dx(0) = dx(0)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1

      iupper(0) = iupper0
      iupper(1) = iupper1

      nugc(0) = nugc0
      nugc(1) = nugc1
c
c     Use a composite B-spline function to interpolate u onto V.
c
      do l = 0,nindices-1
         s = indices(l)
         do d=0,NDIM-1
c
c     Determine the interpolation stencil from the data point below s.
c
            ic_lower(d) = ibrack(d,l)-2
            ic_upper(d) = ibrack(d,l)+3
            ic_lower(d) = max(ic_lower(d),ilower(d)-nugc(d))
            ic_upper(d) = min(ic_upper(d),iupper(d)+nugc(d))
c
c     Compute the interpolation weights.
c
            do ic = ic_lower(d),ic_upper(d)
               if (d .eq. axis) then
                  w(d,ic-ic_lower(d)) =
     &                 lagrangian_bspline_6_delta(
     &                 g(d,l)+dble(ibrack(d,l)-ic))
               else
                  w(d,ic-ic_lower(d)) =
     &                 lagrangian_bspline_5_delta(
     &                 g(d,l)+dble(ibrack(d,l)-ic))
               endif
            enddo
         enddo
c
c     Interpolate u onto V.
c
         INTERPOLATE_2D_SPECIALIZE_FIXED_WIDTH(ic_lower(1), ic_upper(1),
                                               ic_lower(0), ic_upper(0),
                                               6)
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Spread V onto u at the points listed in indices using the composite
c     6-point/5-point B-spline delta function using standard (double)
c     precision accumulation on the Cartesian grid.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_composite_bspline_65_spread2d(
     &     dx,depth,axis,
     &     indices,g,inear,ibrack,nindices,
     &     V,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u)
c
      implicit none
c
c     Functions.
c
      REAL lagrangian_bspline_5_delta
      REAL lagrangian_bspline_6_delta
c
c     Input.
c
      INTEGER depth,axis
      INTEGER nindices
      INTEGER ilower0,iupper0,ilower1,iupper1
      INTEGER nugc0,nugc1

      INTEGER indices(0:nindices-1)
      INTEGER inear(0:NDIM-1,0:nindices-1)
      INTEGER ibrack(0:NDIM-1,0:nindices-1)

      REAL g(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1)
      REAL V(0:depth-1,0:*)
c
c     Input/Output: spreading adds to the existing Cartesian grid values u.
c
      REAL u(CELL2dVECG(ilower,iupper,nugc),0:depth-1)
c
c     Local variables.
c
      INTEGER ilower(0:NDIM-1),iupper(0:NDIM-1)
      INTEGER ic,ic0,ic1
      INTEGER ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER d,l,s,nugc(0:NDIM-1)

      REAL w(0:NDIM-1,0:5)
c
c     Prevent compiler warning about unused variables.
c
      inear(0,0) = inear(0,0)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1

      iupper(0) = iupper0
      iupper(1) = iupper1

      nugc(0) = nugc0
      nugc(1) = nugc1
c
c     Use a composite B-spline function to spread V onto u.
c
      do l = 0,nindices-1
         s = indices(l)
         do d=0,NDIM-1
c
c     Determine the spreading stencil from the data point below s.
c
            ic_lower(d) = ibrack(d,l)-2
            ic_upper(d) = ibrack(d,l)+3
            ic_lower(d) = max(ic_lower(d),ilower(d)-nugc(d))
            ic_upper(d) = min(ic_upper(d),iupper(d)+nugc(d))
c
c     Compute the spreading weights.
c
            do ic = ic_lower(d),ic_upper(d)
               if (d .eq. axis) then
                  w(d,ic-ic_lower(d)) =
     &                 lagrangian_bspline_6_delta(
     &                 g(d,l)+dble(ibrack(d,l)-ic))
               else
                  w(d,ic-ic_lower(d)) =
     &                 lagrangian_bspline_5_delta(
     &                 g(d,l)+dble(ibrack(d,l)-ic))
               endif
            enddo
         enddo
c
c     Spread V onto u.
c
         SPREAD_2D_SPECIALIZE_FIXED_WIDTH(ic_lower(1), ic_upper(1),
                                          ic_lower(0), ic_upper(0),
                                          6)
c
c     End loop over points.
c
      enddo
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
