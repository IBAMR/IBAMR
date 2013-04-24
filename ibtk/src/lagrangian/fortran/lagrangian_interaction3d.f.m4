c
c     Copyright (c) 2002-2013, Boyce Griffith
c     All rights reserved.
c
c     Redistribution and use in source and binary forms, with or without
c     modification, are permitted provided that the following conditions
c     are met:
c
c        * Redistributions of source code must retain the above
c          copyright notice, this list of conditions and the following
c          disclaimer.
c
c        * Redistributions in binary form must reproduce the above
c          copyright notice, this list of conditions and the following
c          disclaimer in the documentation and/or other materials
c          provided with the distribution.
c
c        * Neither the name of New York University nor the names of its
c          contributors may be used to endorse or promote products
c          derived from this software without specific prior written
c          permission.
c
c     THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
c     CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
c     INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
c     MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
c     DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS
c     BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
c     EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
c     TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
c     DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
c     ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
c     TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
c     THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
c     SUCH DAMAGE.
c
define(NDIM,3)dnl
define(REAL,`double precision')dnl
define(INTEGER,`integer')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim3d.i)dnl
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Interpolate u onto V at the positions specified by X using the
c     piecewise constant delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_piecewise_constant_interp3d(
     &     dx,x_lower,x_upper,depth,
     &     ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &     nugc0,nugc1,nugc2,
     &     u,
     &     indices,Xshift,nindices,
     &     X,V)
c
      implicit none
c
c     Input.
c
      INTEGER depth
      INTEGER ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      INTEGER nugc0,nugc1,nugc2
      INTEGER nindices

      INTEGER indices(0:nindices-1)

      REAL Xshift(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1),x_lower(0:NDIM-1),x_upper(0:NDIM-1)
      REAL u(CELL3dVECG(ifirst,ilast,nugc),0:depth-1)
      REAL X(0:NDIM-1,0:*)
c
c     Input/Output.
c
      REAL V(0:depth-1,0:*)
c
c     Local variables.
c
      INTEGER ic0,ic1,ic2
      INTEGER d,l,s
c
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
c
c     Use the piecewise constant delta function to interpolate u onto V.
c
      do l = 0,nindices-1
         s = indices(l)
c
c     Determine the Cartesian cell in which X(s) is located.
c
         ic0 = NINT((X(0,s)+Xshift(0,l)-x_lower(0))/dx(0)-0.5d0)+ifirst0
         ic1 = NINT((X(1,s)+Xshift(1,l)-x_lower(1))/dx(1)-0.5d0)+ifirst1
         ic2 = NINT((X(2,s)+Xshift(2,l)-x_lower(2))/dx(2)-0.5d0)+ifirst2
c
c     Interpolate u onto V.
c
         do d = 0,depth-1
            V(d,s) = u(ic0,ic1,ic2,d)
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Spread V onto u at the positions specified by X using the
c     piecewise constant delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_piecewise_constant_spread3d(
     &     dx,x_lower,x_upper,depth,
     &     indices,Xshift,nindices,
     &     X,V,
     &     ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &     nugc0,nugc1,nugc2,
     &     u)
c
      implicit none
c
c     Input.
c
      INTEGER depth
      INTEGER nindices
      INTEGER ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      INTEGER nugc0,nugc1,nugc2

      INTEGER indices(0:nindices-1)

      REAL Xshift(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1),x_lower(0:NDIM-1),x_upper(0:NDIM-1)
      REAL u(CELL3dVECG(ifirst,ilast,nugc),0:depth-1)
      REAL X(0:NDIM-1,0:*)
c
c     Input/Output.
c
      REAL V(0:depth-1,0:*)
c
c     Local variables.
c
      INTEGER ic0,ic1,ic2
      INTEGER d,l,s
c
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
c
c     Use the piecewise constant delta function to spread V onto u.
c
      do l = 0,nindices-1
         s = indices(l)
c
c     Determine the Cartesian cell in which X(s) is located.
c
         ic0 = NINT((X(0,s)+Xshift(0,l)-x_lower(0))/dx(0)-0.5d0)+ifirst0
         ic1 = NINT((X(1,s)+Xshift(1,l)-x_lower(1))/dx(1)-0.5d0)+ifirst1
         ic2 = NINT((X(2,s)+Xshift(2,l)-x_lower(2))/dx(2)-0.5d0)+ifirst2
c
c     Spread V onto u.
c
         do d = 0,depth-1
            u(ic0,ic1,ic2,d) = u(ic0,ic1,ic2,d)
     &           + V(d,s)/(dx(0)*dx(1)*dx(2))
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Interpolate u onto V at the positions specified by X using the
c     piecewise linear delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_piecewise_linear_interp3d(
     &     dx,x_lower,x_upper,depth,
     &     ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &     nugc0,nugc1,nugc2,
     &     u,
     &     indices,Xshift,nindices,
     &     X,V)
c
      implicit none
c
c     Functions.
c
      EXTERNAL lagrangian_floor
      INTEGER lagrangian_floor
      REAL lagrangian_piecewise_linear_delta
c
c     Input.
c
      INTEGER depth
      INTEGER ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      INTEGER nugc0,nugc1,nugc2
      INTEGER nindices

      INTEGER indices(0:nindices-1)

      REAL Xshift(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1),x_lower(0:NDIM-1),x_upper(0:NDIM-1)
      REAL u(CELL3dVECG(ifirst,ilast,nugc),0:depth-1)
      REAL X(0:NDIM-1,0:*)
c
c     Input/Output.
c
      REAL V(0:depth-1,0:*)
c
c     Local variables.
c
      INTEGER ic0,ic1,ic2
      INTEGER ic_center(0:NDIM-1),ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER d,l,s

      REAL X_cell(0:NDIM-1),w0(0:1),w1(0:1),w2(0:1)
c
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
c
c     Use the piecewise linear delta function to interpolate u onto V.
c
      do l = 0,nindices-1
         s = indices(l)
c
c     Determine the Cartesian cell in which X(s) is located.
c
         ic_center(0) =
     &        lagrangian_floor((X(0,s)+Xshift(0,l)-x_lower(0))/dx(0))
     &        + ifirst0
         ic_center(1) =
     &        lagrangian_floor((X(1,s)+Xshift(1,l)-x_lower(1))/dx(1))
     &        + ifirst1
         ic_center(2) =
     &        lagrangian_floor((X(2,s)+Xshift(2,l)-x_lower(2))/dx(2))
     &        + ifirst2

         X_cell(0) = x_lower(0)+(dble(ic_center(0)-ifirst0)+0.5d0)*dx(0)
         X_cell(1) = x_lower(1)+(dble(ic_center(1)-ifirst1)+0.5d0)*dx(1)
         X_cell(2) = x_lower(2)+(dble(ic_center(2)-ifirst2)+0.5d0)*dx(2)
c
c     Determine the standard interpolation stencil corresponding to the
c     position of X(s) within the cell.
c
         do d = 0,NDIM-1
            if ( X(d,s).lt.X_cell(d) ) then
               ic_lower(d) = ic_center(d)-1
               ic_upper(d) = ic_center(d)
            else
               ic_lower(d) = ic_center(d)
               ic_upper(d) = ic_center(d)+1
            endif
         enddo

         ic_lower(0) = max(ic_lower(0),ifirst0-nugc0)
         ic_upper(0) = min(ic_upper(0),ilast0 +nugc0)

         ic_lower(1) = max(ic_lower(1),ifirst1-nugc1)
         ic_upper(1) = min(ic_upper(1),ilast1 +nugc1)

         ic_lower(2) = max(ic_lower(2),ifirst2-nugc2)
         ic_upper(2) = min(ic_upper(2),ilast2 +nugc2)
c
c     Compute the standard interpolation weights.
c
CDEC$ LOOP COUNT(2)
         do ic0 = ic_lower(0),ic_upper(0)
            X_cell(0) = x_lower(0)+(dble(ic0-ifirst0)+0.5d0)*dx(0)
            w0(ic0-ic_lower(0)) =
     &           lagrangian_piecewise_linear_delta(
     &           (X(0,s)+Xshift(0,l)-X_cell(0))/dx(0))
         enddo
CDEC$ LOOP COUNT(2)
         do ic1 = ic_lower(1),ic_upper(1)
            X_cell(1) = x_lower(1)+(dble(ic1-ifirst1)+0.5d0)*dx(1)
            w1(ic1-ic_lower(1)) =
     &           lagrangian_piecewise_linear_delta(
     &           (X(1,s)+Xshift(1,l)-X_cell(1))/dx(1))
         enddo
CDEC$ LOOP COUNT(2)
         do ic2 = ic_lower(2),ic_upper(2)
            X_cell(2) = x_lower(2)+(dble(ic2-ifirst2)+0.5d0)*dx(2)
            w2(ic2-ic_lower(2)) =
     &           lagrangian_piecewise_linear_delta(
     &           (X(2,s)+Xshift(2,l)-X_cell(2))/dx(2))
         enddo
c
c     Interpolate u onto V.
c
         do d = 0,depth-1
            V(d,s) = 0.d0
CDEC$ LOOP COUNT(2)
            do ic2 = ic_lower(2),ic_upper(2)
CDEC$ LOOP COUNT(2)
               do ic1 = ic_lower(1),ic_upper(1)
CDEC$ LOOP COUNT(2)
                  do ic0 = ic_lower(0),ic_upper(0)
                     V(d,s) = V(d,s)
     &                    +w0(ic0-ic_lower(0))
     &                    *w1(ic1-ic_lower(1))
     &                    *w2(ic2-ic_lower(2))
     &                    *u(ic0,ic1,ic2,d)
                  enddo
               enddo
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Spread V onto u at the positions specified by X using the
c     piecewise linear delta function using standard (double) precision
c     accumulation on the Cartesian grid.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_piecewise_linear_spread3d(
     &     dx,x_lower,x_upper,depth,
     &     indices,Xshift,nindices,
     &     X,V,
     &     ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &     nugc0,nugc1,nugc2,
     &     u)
c
      implicit none
c
c     Functions.
c
      EXTERNAL lagrangian_floor
      INTEGER lagrangian_floor
      REAL lagrangian_piecewise_linear_delta
c
c     Input.
c
      INTEGER depth
      INTEGER nindices
      INTEGER ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      INTEGER nugc0,nugc1,nugc2

      INTEGER indices(0:nindices-1)

      REAL Xshift(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1),x_lower(0:NDIM-1),x_upper(0:NDIM-1)
      REAL u(CELL3dVECG(ifirst,ilast,nugc),0:depth-1)
      REAL X(0:NDIM-1,0:*)
c
c     Input/Output.
c
      REAL V(0:depth-1,0:*)
c
c     Local variables.
c
      INTEGER ic0,ic1,ic2
      INTEGER ic_center(0:NDIM-1),ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER d,l,s

      REAL X_cell(0:NDIM-1),w0(0:1),w1(0:1),w2(0:1)
c
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
c
c     Use the piecewise linear delta function to spread V onto u.
c
      do l = 0,nindices-1
         s = indices(l)
c
c     Determine the Cartesian cell in which X(s) is located.
c
         ic_center(0) =
     &        lagrangian_floor((X(0,s)+Xshift(0,l)-x_lower(0))/dx(0))
     &        + ifirst0
         ic_center(1) =
     &        lagrangian_floor((X(1,s)+Xshift(1,l)-x_lower(1))/dx(1))
     &        + ifirst1
         ic_center(2) =
     &        lagrangian_floor((X(2,s)+Xshift(2,l)-x_lower(2))/dx(2))
     &        + ifirst2

         X_cell(0) = x_lower(0)+(dble(ic_center(0)-ifirst0)+0.5d0)*dx(0)
         X_cell(1) = x_lower(1)+(dble(ic_center(1)-ifirst1)+0.5d0)*dx(1)
         X_cell(2) = x_lower(2)+(dble(ic_center(2)-ifirst2)+0.5d0)*dx(2)
c
c     Determine the standard spreading stencil corresponding to the
c     position of X(s) within the cell.
c
         do d = 0,NDIM-1
            if ( X(d,s).lt.X_cell(d) ) then
               ic_lower(d) = ic_center(d)-1
               ic_upper(d) = ic_center(d)
            else
               ic_lower(d) = ic_center(d)
               ic_upper(d) = ic_center(d)+1
            endif
         enddo

         ic_lower(0) = max(ic_lower(0),ifirst0-nugc0)
         ic_upper(0) = min(ic_upper(0),ilast0 +nugc0)

         ic_lower(1) = max(ic_lower(1),ifirst1-nugc1)
         ic_upper(1) = min(ic_upper(1),ilast1 +nugc1)

         ic_lower(2) = max(ic_lower(2),ifirst2-nugc2)
         ic_upper(2) = min(ic_upper(2),ilast2 +nugc2)
c
c     Compute the standard spreading weights.
c
CDEC$ LOOP COUNT(2)
         do ic0 = ic_lower(0),ic_upper(0)
            X_cell(0) = x_lower(0)+(dble(ic0-ifirst0)+0.5d0)*dx(0)
            w0(ic0-ic_lower(0)) =
     &           lagrangian_piecewise_linear_delta(
     &           (X(0,s)+Xshift(0,l)-X_cell(0))/dx(0))
         enddo
CDEC$ LOOP COUNT(2)
         do ic1 = ic_lower(1),ic_upper(1)
            X_cell(1) = x_lower(1)+(dble(ic1-ifirst1)+0.5d0)*dx(1)
            w1(ic1-ic_lower(1)) =
     &           lagrangian_piecewise_linear_delta(
     &           (X(1,s)+Xshift(1,l)-X_cell(1))/dx(1))
         enddo
CDEC$ LOOP COUNT(2)
         do ic2 = ic_lower(2),ic_upper(2)
            X_cell(2) = x_lower(2)+(dble(ic2-ifirst2)+0.5d0)*dx(2)
            w2(ic2-ic_lower(2)) =
     &           lagrangian_piecewise_linear_delta(
     &           (X(2,s)+Xshift(2,l)-X_cell(2))/dx(2))
         enddo
c
c     Spread V onto u.
c
         do d = 0,depth-1
CDEC$ LOOP COUNT(2)
            do ic2 = ic_lower(2),ic_upper(2)
CDEC$ LOOP COUNT(2)
               do ic1 = ic_lower(1),ic_upper(1)
CDEC$ LOOP COUNT(2)
                  do ic0 = ic_lower(0),ic_upper(0)
                     u(ic0,ic1,ic2,d) = u(ic0,ic1,ic2,d)+(
     &                    w0(ic0-ic_lower(0))*
     &                    w1(ic1-ic_lower(1))*
     &                    w2(ic2-ic_lower(2))*
     &                    V(d,s)/(dx(0)*dx(1)*dx(2)))
                  enddo
               enddo
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Interpolate u onto V at the positions specified by X using the IB
c     3-point delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_ib_3_interp3d(
     &     dx,x_lower,x_upper,depth,
     &     ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &     nugc0,nugc1,nugc2,
     &     u,
     &     indices,Xshift,nindices,
     &     X,V)
c
      implicit none
c
c     Functions.
c
      EXTERNAL lagrangian_floor
      INTEGER lagrangian_floor
      REAL lagrangian_ib_3_delta
c
c     Input.
c
      INTEGER depth
      INTEGER ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      INTEGER nugc0,nugc1,nugc2
      INTEGER nindices

      INTEGER indices(0:nindices-1)

      REAL Xshift(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1),x_lower(0:NDIM-1),x_upper(0:NDIM-1)
      REAL u(CELL3dVECG(ifirst,ilast,nugc),0:depth-1)
      REAL X(0:NDIM-1,0:*)
c
c     Input/Output.
c
      REAL V(0:depth-1,0:*)
c
c     Local variables.
c
      INTEGER ic0,ic1,ic2
      INTEGER ic_center(0:NDIM-1),ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER d,l,s

      REAL X_cell(0:NDIM-1),w0(0:2),w1(0:2),w2(0:2)
c
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
c
c     Use the IB 3-point delta function to interpolate u onto V.
c
      do l = 0,nindices-1
         s = indices(l)
c
c     Determine the Cartesian cell in which X(s) is located.
c
         ic_center(0) =
     &        lagrangian_floor((X(0,s)+Xshift(0,l)-x_lower(0))/dx(0))
     &        + ifirst0
         ic_center(1) =
     &        lagrangian_floor((X(1,s)+Xshift(1,l)-x_lower(1))/dx(1))
     &        + ifirst1
         ic_center(2) =
     &        lagrangian_floor((X(2,s)+Xshift(2,l)-x_lower(2))/dx(2))
     &        + ifirst2

         X_cell(0) = x_lower(0)+(dble(ic_center(0)-ifirst0)+0.5d0)*dx(0)
         X_cell(1) = x_lower(1)+(dble(ic_center(1)-ifirst1)+0.5d0)*dx(1)
         X_cell(2) = x_lower(2)+(dble(ic_center(2)-ifirst2)+0.5d0)*dx(2)
c
c     Determine the standard interpolation stencil corresponding to the
c     position of X(s) within the cell.
c
         do d = 0,NDIM-1
            ic_lower(d) = ic_center(d)-1
            ic_upper(d) = ic_center(d)+1
         enddo

         ic_lower(0) = max(ic_lower(0),ifirst0-nugc0)
         ic_upper(0) = min(ic_upper(0),ilast0 +nugc0)

         ic_lower(1) = max(ic_lower(1),ifirst1-nugc1)
         ic_upper(1) = min(ic_upper(1),ilast1 +nugc1)

         ic_lower(2) = max(ic_lower(2),ifirst2-nugc2)
         ic_upper(2) = min(ic_upper(2),ilast2 +nugc2)
c
c     Compute the standard interpolation weights.
c
CDEC$ LOOP COUNT(3)
         do ic0 = ic_lower(0),ic_upper(0)
            X_cell(0) = x_lower(0)+(dble(ic0-ifirst0)+0.5d0)*dx(0)
            w0(ic0-ic_lower(0)) =
     &           lagrangian_ib_3_delta(
     &           (X(0,s)+Xshift(0,l)-X_cell(0))/dx(0))
         enddo
CDEC$ LOOP COUNT(3)
         do ic1 = ic_lower(1),ic_upper(1)
            X_cell(1) = x_lower(1)+(dble(ic1-ifirst1)+0.5d0)*dx(1)
            w1(ic1-ic_lower(1)) =
     &           lagrangian_ib_3_delta(
     &           (X(1,s)+Xshift(1,l)-X_cell(1))/dx(1))
         enddo
CDEC$ LOOP COUNT(3)
         do ic2 = ic_lower(2),ic_upper(2)
            X_cell(2) = x_lower(2)+(dble(ic2-ifirst2)+0.5d0)*dx(2)
            w2(ic2-ic_lower(2)) =
     &           lagrangian_ib_3_delta(
     &           (X(2,s)+Xshift(2,l)-X_cell(2))/dx(2))
         enddo
c
c     Interpolate u onto V.
c
         do d = 0,depth-1
            V(d,s) = 0.d0
CDEC$ LOOP COUNT(3)
            do ic2 = ic_lower(2),ic_upper(2)
CDEC$ LOOP COUNT(3)
               do ic1 = ic_lower(1),ic_upper(1)
CDEC$ LOOP COUNT(3)
                  do ic0 = ic_lower(0),ic_upper(0)
                     V(d,s) = V(d,s)
     &                    +w0(ic0-ic_lower(0))
     &                    *w1(ic1-ic_lower(1))
     &                    *w2(ic2-ic_lower(2))
     &                    *u(ic0,ic1,ic2,d)
                  enddo
               enddo
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Spread V onto u at the positions specified by X using the IB
c     3-point delta function using standard (double) precision
c     accumulation on the Cartesian grid.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_ib_3_spread3d(
     &     dx,x_lower,x_upper,depth,
     &     indices,Xshift,nindices,
     &     X,V,
     &     ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &     nugc0,nugc1,nugc2,
     &     u)
c
      implicit none
c
c     Functions.
c
      EXTERNAL lagrangian_floor
      INTEGER lagrangian_floor
      REAL lagrangian_ib_3_delta
c
c     Input.
c
      INTEGER depth
      INTEGER nindices
      INTEGER ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      INTEGER nugc0,nugc1,nugc2

      INTEGER indices(0:nindices-1)

      REAL Xshift(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1),x_lower(0:NDIM-1),x_upper(0:NDIM-1)
      REAL u(CELL3dVECG(ifirst,ilast,nugc),0:depth-1)
      REAL X(0:NDIM-1,0:*)
c
c     Input/Output.
c
      REAL V(0:depth-1,0:*)
c
c     Local variables.
c
      INTEGER ic0,ic1,ic2
      INTEGER ic_center(0:NDIM-1),ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER d,l,s

      REAL X_cell(0:NDIM-1),w0(0:2),w1(0:2),w2(0:2)
c
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
c
c     Use the IB 3-point delta function to spread V onto u.
c
      do l = 0,nindices-1
         s = indices(l)
c
c     Determine the Cartesian cell in which X(s) is located.
c
         ic_center(0) =
     &        lagrangian_floor((X(0,s)+Xshift(0,l)-x_lower(0))/dx(0))
     &        + ifirst0
         ic_center(1) =
     &        lagrangian_floor((X(1,s)+Xshift(1,l)-x_lower(1))/dx(1))
     &        + ifirst1
         ic_center(2) =
     &        lagrangian_floor((X(2,s)+Xshift(2,l)-x_lower(2))/dx(2))
     &        + ifirst2

         X_cell(0) = x_lower(0)+(dble(ic_center(0)-ifirst0)+0.5d0)*dx(0)
         X_cell(1) = x_lower(1)+(dble(ic_center(1)-ifirst1)+0.5d0)*dx(1)
         X_cell(2) = x_lower(2)+(dble(ic_center(2)-ifirst2)+0.5d0)*dx(2)
c
c     Determine the standard spreading stencil corresponding to the
c     position of X(s) within the cell.
c
         do d = 0,NDIM-1
            ic_lower(d) = ic_center(d)-1
            ic_upper(d) = ic_center(d)+1
         enddo

         ic_lower(0) = max(ic_lower(0),ifirst0-nugc0)
         ic_upper(0) = min(ic_upper(0),ilast0 +nugc0)

         ic_lower(1) = max(ic_lower(1),ifirst1-nugc1)
         ic_upper(1) = min(ic_upper(1),ilast1 +nugc1)

         ic_lower(2) = max(ic_lower(2),ifirst2-nugc2)
         ic_upper(2) = min(ic_upper(2),ilast2 +nugc2)
c
c     Compute the standard spreading weights.
c
CDEC$ LOOP COUNT(3)
         do ic0 = ic_lower(0),ic_upper(0)
            X_cell(0) = x_lower(0)+(dble(ic0-ifirst0)+0.5d0)*dx(0)
            w0(ic0-ic_lower(0)) =
     &           lagrangian_ib_3_delta(
     &           (X(0,s)+Xshift(0,l)-X_cell(0))/dx(0))
         enddo
CDEC$ LOOP COUNT(3)
         do ic1 = ic_lower(1),ic_upper(1)
            X_cell(1) = x_lower(1)+(dble(ic1-ifirst1)+0.5d0)*dx(1)
            w1(ic1-ic_lower(1)) =
     &           lagrangian_ib_3_delta(
     &           (X(1,s)+Xshift(1,l)-X_cell(1))/dx(1))
         enddo
CDEC$ LOOP COUNT(3)
         do ic2 = ic_lower(2),ic_upper(2)
            X_cell(2) = x_lower(2)+(dble(ic2-ifirst2)+0.5d0)*dx(2)
            w2(ic2-ic_lower(2)) =
     &           lagrangian_ib_3_delta(
     &           (X(2,s)+Xshift(2,l)-X_cell(2))/dx(2))
         enddo
c
c     Spread V onto u.
c
         do d = 0,depth-1
CDEC$ LOOP COUNT(3)
            do ic2 = ic_lower(2),ic_upper(2)
CDEC$ LOOP COUNT(3)
               do ic1 = ic_lower(1),ic_upper(1)
CDEC$ LOOP COUNT(3)
                  do ic0 = ic_lower(0),ic_upper(0)
                     u(ic0,ic1,ic2,d) = u(ic0,ic1,ic2,d)+(
     &                    w0(ic0-ic_lower(0))*
     &                    w1(ic1-ic_lower(1))*
     &                    w2(ic2-ic_lower(2))*
     &                    V(d,s)/(dx(0)*dx(1)*dx(2)))
                  enddo
               enddo
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Interpolate u onto V at the positions specified by X using the IB
c     4-point delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_ib_4_interp3d(
     &     dx,x_lower,x_upper,depth,
     &     ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &     patch_touches_lower_physical_bdry,
     &     patch_touches_upper_physical_bdry,
     &     nugc0,nugc1,nugc2,
     &     u,
     &     indices,Xshift,nindices,
     &     X,V)
c
      implicit none
c
c     Input.
c
      INTEGER depth
      INTEGER ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      INTEGER nugc0,nugc1,nugc2
      INTEGER nindices

      INTEGER patch_touches_lower_physical_bdry(0:NDIM-1)
      INTEGER patch_touches_upper_physical_bdry(0:NDIM-1)

      INTEGER indices(0:nindices-1)

      REAL Xshift(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1),x_lower(0:NDIM-1),x_upper(0:NDIM-1)
      REAL u(CELL3dVECG(ifirst,ilast,nugc),0:depth-1)
      REAL X(0:NDIM-1,0:*)
c
c     Input/Output.
c
      REAL V(0:depth-1,0:*)
c
c     Local variables.
c
      INTEGER i0,i1,i2,ic0,ic1,ic2
      INTEGER ig_lower(0:NDIM-1),ig_upper(0:NDIM-1)
      INTEGER ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER istart0,istop0,istart1,istop1,istart2,istop2
      INTEGER d,k,l,s

      REAL X_o_dx,q0,q1,q2,r0,r1,r2
      REAL w0(0:3),w1(0:3),w2(0:3),f(0:3)
      REAL w(0:3,0:3,0:3),wyz,wz

      LOGICAL account_for_phys_bdry
      LOGICAL touches_lower_bdry(0:NDIM-1)
      LOGICAL touches_upper_bdry(0:NDIM-1)
c
c     Compute the extents of the ghost box.
c
      ig_lower(0) = ifirst0-nugc0
      ig_lower(1) = ifirst1-nugc1
      ig_lower(2) = ifirst2-nugc2
      ig_upper(0) = ilast0 +nugc0
      ig_upper(1) = ilast1 +nugc1
      ig_upper(2) = ilast2 +nugc2
c
c     Determine if we need to account for physical boundaries.
c
      account_for_phys_bdry = .false.
      do d = 0,NDIM-1
         account_for_phys_bdry = account_for_phys_bdry    .or.
     &        (patch_touches_lower_physical_bdry(d).eq.1) .or.
     &        (patch_touches_upper_physical_bdry(d).eq.1)
      enddo
c
c     Use the IB 4-point delta function to interpolate u onto V, but use
c     a modified delta function near physical boundaries.
c
      do l = 0,nindices-1
         s = indices(l)
c
c     Determine the standard interpolation stencil corresponding to the
c     position of X(s) within the cell and compute the standard
c     interpolation weights.
c
         X_o_dx = (X(0,s)+Xshift(0,l)-x_lower(0))/dx(0)
         ic_lower(0) = NINT(X_o_dx)+ifirst0-2
         ic_upper(0) = ic_lower(0) + 3
         r0 = X_o_dx - ((ic_lower(0)+1-ifirst0)+0.5d0)
         q0 = sqrt(1.d0+4.d0*r0*(1.d0-r0))
         w0(0) = 0.125d0*(3.d0-2.d0*r0-q0)
         w0(1) = 0.125d0*(3.d0-2.d0*r0+q0)
         w0(2) = 0.125d0*(1.d0+2.d0*r0+q0)
         w0(3) = 0.125d0*(1.d0+2.d0*r0-q0)

         X_o_dx = (X(1,s)+Xshift(1,l)-x_lower(1))/dx(1)
         ic_lower(1) = NINT(X_o_dx)+ifirst1-2
         ic_upper(1) = ic_lower(1) + 3
         r1 = X_o_dx - ((ic_lower(1)+1-ifirst1)+0.5d0)
         q1 = sqrt(1.d0+4.d0*r1*(1.d0-r1))
         w1(0) = 0.125d0*(3.d0-2.d0*r1-q1)
         w1(1) = 0.125d0*(3.d0-2.d0*r1+q1)
         w1(2) = 0.125d0*(1.d0+2.d0*r1+q1)
         w1(3) = 0.125d0*(1.d0+2.d0*r1-q1)

         X_o_dx = (X(2,s)+Xshift(2,l)-x_lower(2))/dx(2)
         ic_lower(2) = NINT(X_o_dx)+ifirst2-2
         ic_upper(2) = ic_lower(2) + 3
         r2 = X_o_dx - ((ic_lower(2)+1-ifirst2)+0.5d0)
         q2 = sqrt(1.d0+4.d0*r2*(1.d0-r2))
         w2(0) = 0.125d0*(3.d0-2.d0*r2-q2)
         w2(1) = 0.125d0*(3.d0-2.d0*r2+q2)
         w2(2) = 0.125d0*(1.d0+2.d0*r2+q2)
         w2(3) = 0.125d0*(1.d0+2.d0*r2-q2)
c
c     When necessary, modify the interpolation stencil and weights near
c     physical boundaries.
c
         if ( account_for_phys_bdry ) then
            do d = 0,NDIM-1
               touches_lower_bdry(d) =
     &              (patch_touches_lower_physical_bdry(d).eq.1) .and.
     &              (X(d,s) - x_lower(d) .lt. 1.5d0*dx(d))
               touches_upper_bdry(d) =
     &              (patch_touches_upper_physical_bdry(d).eq.1) .and.
     &              (x_upper(d) - X(d,s) .lt. 1.5d0*dx(d))
            enddo

            if (touches_lower_bdry(0)) then
               call lagrangian_one_sided_ib_4_delta(
     &              w0,(X(0,s)-x_lower(0))/dx(0))
               ic_lower(0) = ifirst0
               ic_upper(0) = ifirst0+3
            elseif (touches_upper_bdry(0)) then
               call lagrangian_one_sided_ib_4_delta(
     &              f,(x_upper(0)-X(0,s))/dx(0))
               ic_lower(0) = ilast0-3
               ic_upper(0) = ilast0
               do k = 0,3
                  w0(3-k) = f(k)
               enddo
            endif

            if (touches_lower_bdry(1)) then
               call lagrangian_one_sided_ib_4_delta(
     &              w1,(X(1,s)-x_lower(1))/dx(1))
               ic_lower(1) = ifirst1
               ic_upper(1) = ifirst1+3
            elseif (touches_upper_bdry(1)) then
               call lagrangian_one_sided_ib_4_delta(
     &              f,(x_upper(1)-X(1,s))/dx(1))
               ic_lower(1) = ilast1-3
               ic_upper(1) = ilast1
               do k = 0,3
                  w1(3-k) = f(k)
               enddo
            endif

            if (touches_lower_bdry(2)) then
               call lagrangian_one_sided_ib_4_delta(
     &              w2,(X(2,s)-x_lower(2))/dx(2))
               ic_lower(2) = ifirst2
               ic_upper(2) = ifirst2+3
            elseif (touches_upper_bdry(2)) then
               call lagrangian_one_sided_ib_4_delta(
     &              f,(x_upper(2)-X(2,s))/dx(2))
               ic_lower(2) = ilast2-3
               ic_upper(2) = ilast2
               do k = 0,3
                  w2(3-k) = f(k)
               enddo
            endif
         endif
c
c     Compute the tensor product of the interpolation weights.
c
         do i2 = 0,3
            wz = w2(i2)
            do i1 = 0,3
               wyz = w1(i1)*wz
               do i0 = 0,3
                  w(i0,i1,i2) = w0(i0)*wyz
               enddo
            enddo
         enddo
c
c     Interpolate u onto V.
c
         if ( ic_lower(0).lt.ig_lower(0) .or.
     &        ic_lower(1).lt.ig_lower(1) .or.
     &        ic_lower(2).lt.ig_lower(2) .or.
     &        ic_upper(0).gt.ig_upper(0) .or.
     &        ic_upper(1).gt.ig_upper(1) .or.
     &        ic_upper(2).gt.ig_upper(2) ) then
            istart0 =   max(ig_lower(0)-ic_lower(0),0)
            istop0  = 3-max(ic_upper(0)-ig_upper(0),0)
            istart1 =   max(ig_lower(1)-ic_lower(1),0)
            istop1  = 3-max(ic_upper(1)-ig_upper(1),0)
            istart2 =   max(ig_lower(2)-ic_lower(2),0)
            istop2  = 3-max(ic_upper(2)-ig_upper(2),0)
            do d = 0,depth-1
               V(d,s) = 0.d0
               do i2 = istart2,istop2
                  ic2 = ic_lower(2)+i2
                  do i1 = istart1,istop1
                     ic1 = ic_lower(1)+i1
                     do i0 = istart0,istop0
                        ic0 = ic_lower(0)+i0
                        V(d,s) = V(d,s) + w(i0,i1,i2)*u(ic0,ic1,ic2,d)
                     enddo
                  enddo
               enddo
            enddo
         else
            do d = 0,depth-1
               V(d,s) = 0.d0
               do i2 = 0,3
                  ic2 = ic_lower(2)+i2
                  do i1 = 0,3
                     ic1 = ic_lower(1)+i1
                     do i0 = 0,3
                        ic0 = ic_lower(0)+i0
                        V(d,s) = V(d,s) + w(i0,i1,i2)*u(ic0,ic1,ic2,d)
                     enddo
                  enddo
               enddo
            enddo
         endif
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Spread V onto u at the positions specified by X using the IB
c     4-point delta function using standard (double) precision
c     accumulation on the Cartesian grid.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_ib_4_spread3d(
     &     dx,x_lower,x_upper,depth,
     &     indices,Xshift,nindices,
     &     X,V,
     &     ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &     patch_touches_lower_physical_bdry,
     &     patch_touches_upper_physical_bdry,
     &     nugc0,nugc1,nugc2,
     &     u)
c
      implicit none
c
c     Input.
c
      INTEGER depth
      INTEGER nindices
      INTEGER ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      INTEGER nugc0,nugc1,nugc2

      INTEGER indices(0:nindices-1)

      INTEGER patch_touches_lower_physical_bdry(0:NDIM-1)
      INTEGER patch_touches_upper_physical_bdry(0:NDIM-1)

      REAL Xshift(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1),x_lower(0:NDIM-1),x_upper(0:NDIM-1)
      REAL u(CELL3dVECG(ifirst,ilast,nugc),0:depth-1)
      REAL X(0:NDIM-1,0:*)
c
c     Input/Output.
c
      REAL V(0:depth-1,0:*)
c
c     Local variables.
c
      INTEGER i0,i1,i2,ic0,ic1,ic2
      INTEGER ig_lower(0:NDIM-1),ig_upper(0:NDIM-1)
      INTEGER ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER istart0,istop0,istart1,istop1,istart2,istop2
      INTEGER d,k,l,s

      REAL X_o_dx,q0,q1,q2,r0,r1,r2
      REAL w0(0:3),w1(0:3),w2(0:3),f(0:3)
      REAL w(0:3,0:3,0:3),wyz,wz

      LOGICAL account_for_phys_bdry
      LOGICAL touches_lower_bdry(0:NDIM-1)
      LOGICAL touches_upper_bdry(0:NDIM-1)
c
c     Compute the extents of the ghost box.
c
      ig_lower(0) = ifirst0-nugc0
      ig_lower(1) = ifirst1-nugc1
      ig_lower(2) = ifirst2-nugc2
      ig_upper(0) = ilast0 +nugc0
      ig_upper(1) = ilast1 +nugc1
      ig_upper(2) = ilast2 +nugc2
c
c     Determine if we need to account for physical boundaries.
c
      account_for_phys_bdry = .false.
      do d = 0,NDIM-1
         account_for_phys_bdry = account_for_phys_bdry    .or.
     &        (patch_touches_lower_physical_bdry(d).eq.1) .or.
     &        (patch_touches_upper_physical_bdry(d).eq.1)
      enddo
c
c     Use the IB 4-point delta function to spread V onto u, but use
c     a modified delta function near physical boundaries.
c
      do l = 0,nindices-1
         s = indices(l)
c
c     Determine the standard interpolation stencil corresponding to the
c     position of X(s) within the cell and compute the standard
c     interpolation weights.
c
         X_o_dx = (X(0,s)+Xshift(0,l)-x_lower(0))/dx(0)
         ic_lower(0) = NINT(X_o_dx)+ifirst0-2
         ic_upper(0) = ic_lower(0) + 3
         r0 = X_o_dx - ((ic_lower(0)+1-ifirst0)+0.5d0)
         q0 = sqrt(1.d0+4.d0*r0*(1.d0-r0))
         w0(0) = 0.125d0*(3.d0-2.d0*r0-q0)
         w0(1) = 0.125d0*(3.d0-2.d0*r0+q0)
         w0(2) = 0.125d0*(1.d0+2.d0*r0+q0)
         w0(3) = 0.125d0*(1.d0+2.d0*r0-q0)

         X_o_dx = (X(1,s)+Xshift(1,l)-x_lower(1))/dx(1)
         ic_lower(1) = NINT(X_o_dx)+ifirst1-2
         ic_upper(1) = ic_lower(1) + 3
         r1 = X_o_dx - ((ic_lower(1)+1-ifirst1)+0.5d0)
         q1 = sqrt(1.d0+4.d0*r1*(1.d0-r1))
         w1(0) = 0.125d0*(3.d0-2.d0*r1-q1)
         w1(1) = 0.125d0*(3.d0-2.d0*r1+q1)
         w1(2) = 0.125d0*(1.d0+2.d0*r1+q1)
         w1(3) = 0.125d0*(1.d0+2.d0*r1-q1)

         X_o_dx = (X(2,s)+Xshift(2,l)-x_lower(2))/dx(2)
         ic_lower(2) = NINT(X_o_dx)+ifirst2-2
         ic_upper(2) = ic_lower(2) + 3
         r2 = X_o_dx - ((ic_lower(2)+1-ifirst2)+0.5d0)
         q2 = sqrt(1.d0+4.d0*r2*(1.d0-r2))
         w2(0) = 0.125d0*(3.d0-2.d0*r2-q2)
         w2(1) = 0.125d0*(3.d0-2.d0*r2+q2)
         w2(2) = 0.125d0*(1.d0+2.d0*r2+q2)
         w2(3) = 0.125d0*(1.d0+2.d0*r2-q2)
c
c     When necessary, modify the interpolation stencil and weights near
c     physical boundaries.
c
         if ( account_for_phys_bdry ) then
            do d = 0,NDIM-1
               touches_lower_bdry(d) =
     &              (patch_touches_lower_physical_bdry(d).eq.1) .and.
     &              (X(d,s) - x_lower(d) .lt. 1.5d0*dx(d))
               touches_upper_bdry(d) =
     &              (patch_touches_upper_physical_bdry(d).eq.1) .and.
     &              (x_upper(d) - X(d,s) .lt. 1.5d0*dx(d))
            enddo

            if (touches_lower_bdry(0)) then
               call lagrangian_one_sided_ib_4_delta(
     &              w0,(X(0,s)-x_lower(0))/dx(0))
               ic_lower(0) = ifirst0
               ic_upper(0) = ifirst0+3
            elseif (touches_upper_bdry(0)) then
               call lagrangian_one_sided_ib_4_delta(
     &              f,(x_upper(0)-X(0,s))/dx(0))
               ic_lower(0) = ilast0-3
               ic_upper(0) = ilast0
               do k = 0,3
                  w0(3-k) = f(k)
               enddo
            endif

            if (touches_lower_bdry(1)) then
               call lagrangian_one_sided_ib_4_delta(
     &              w1,(X(1,s)-x_lower(1))/dx(1))
               ic_lower(1) = ifirst1
               ic_upper(1) = ifirst1+3
            elseif (touches_upper_bdry(1)) then
               call lagrangian_one_sided_ib_4_delta(
     &              f,(x_upper(1)-X(1,s))/dx(1))
               ic_lower(1) = ilast1-3
               ic_upper(1) = ilast1
               do k = 0,3
                  w1(3-k) = f(k)
               enddo
            endif

            if (touches_lower_bdry(2)) then
               call lagrangian_one_sided_ib_4_delta(
     &              w2,(X(2,s)-x_lower(2))/dx(2))
               ic_lower(2) = ifirst2
               ic_upper(2) = ifirst2+3
            elseif (touches_upper_bdry(2)) then
               call lagrangian_one_sided_ib_4_delta(
     &              f,(x_upper(2)-X(2,s))/dx(2))
               ic_lower(2) = ilast2-3
               ic_upper(2) = ilast2
               do k = 0,3
                  w2(3-k) = f(k)
               enddo
            endif
         endif
c
c     Compute the tensor product of the scaled interpolation weights.
c
         do i2 = 0,3
            wz = w2(i2)/(dx(0)*dx(1)*dx(2))
            do i1 = 0,3
               wyz = w1(i1)*wz
               do i0 = 0,3
                  w(i0,i1,i2) = w0(i0)*wyz
               enddo
            enddo
         enddo
c
c     Spread V onto u.
c
         if ( ic_lower(0).lt.ig_lower(0) .or.
     &        ic_lower(1).lt.ig_lower(1) .or.
     &        ic_lower(2).lt.ig_lower(2) .or.
     &        ic_upper(0).gt.ig_upper(0) .or.
     &        ic_upper(1).gt.ig_upper(1) .or.
     &        ic_upper(2).gt.ig_upper(2) ) then
            istart0 =   max(ig_lower(0)-ic_lower(0),0)
            istop0  = 3-max(ic_upper(0)-ig_upper(0),0)
            istart1 =   max(ig_lower(1)-ic_lower(1),0)
            istop1  = 3-max(ic_upper(1)-ig_upper(1),0)
            istart2 =   max(ig_lower(2)-ic_lower(2),0)
            istop2  = 3-max(ic_upper(2)-ig_upper(2),0)
            do d = 0,depth-1
               do i2 = istart2,istop2
                  ic2 = ic_lower(2)+i2
                  do i1 = istart1,istop1
                     ic1 = ic_lower(1)+i1
                     do i0 = istart0,istop0
                        ic0 = ic_lower(0)+i0
                        u(ic0,ic1,ic2,d) = u(ic0,ic1,ic2,d) +
     &                       w(i0,i1,i2)*V(d,s)
                     enddo
                  enddo
               enddo
            enddo
         else
            do d = 0,depth-1
               do i2 = 0,3
                  ic2 = ic_lower(2)+i2
                  do i1 = 0,3
                     ic1 = ic_lower(1)+i1
                     do i0 = 0,3
                        ic0 = ic_lower(0)+i0
                        u(ic0,ic1,ic2,d) = u(ic0,ic1,ic2,d) +
     &                       w(i0,i1,i2)*V(d,s)
                     enddo
                  enddo
               enddo
            enddo
         endif
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Interpolate u onto V at the positions specified by X using the
c     broadened 8-point version of the IB 4-point delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_wide8_ib_4_interp3d(
     &     dx,x_lower,x_upper,depth,
     &     ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &     nugc0,nugc1,nugc2,
     &     u,
     &     indices,Xshift,nindices,
     &     X,V)
c
      implicit none
c
c     Functions.
c
      EXTERNAL lagrangian_floor
      INTEGER lagrangian_floor
      REAL lagrangian_wide8_ib_4_delta
c
c     Input.
c
      INTEGER depth
      INTEGER ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      INTEGER nugc0,nugc1,nugc2
      INTEGER nindices

      INTEGER indices(0:nindices-1)

      REAL Xshift(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1),x_lower(0:NDIM-1),x_upper(0:NDIM-1)
      REAL u(CELL3dVECG(ifirst,ilast,nugc),0:depth-1)
      REAL X(0:NDIM-1,0:*)
c
c     Input/Output.
c
      REAL V(0:depth-1,0:*)
c
c     Local variables.
c
      INTEGER ic0,ic1,ic2
      INTEGER ic_center(0:NDIM-1),ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER d,l,s

      REAL X_cell(0:NDIM-1),w0(0:7),w1(0:7),w2(0:7)
c
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
c
c     Use the broadened (8-point) version of the IB 4-point delta
c     function to interpolate u onto V.
c
      do l = 0,nindices-1
         s = indices(l)
c
c     Determine the Cartesian cell in which X(s) is located.
c
         ic_center(0) =
     &        lagrangian_floor((X(0,s)+Xshift(0,l)-x_lower(0))/dx(0))
     &        + ifirst0
         ic_center(1) =
     &        lagrangian_floor((X(1,s)+Xshift(1,l)-x_lower(1))/dx(1))
     &        + ifirst1
         ic_center(2) =
     &        lagrangian_floor((X(2,s)+Xshift(2,l)-x_lower(2))/dx(2))
     &        + ifirst2

         X_cell(0) = x_lower(0)+(dble(ic_center(0)-ifirst0)+0.5d0)*dx(0)
         X_cell(1) = x_lower(1)+(dble(ic_center(1)-ifirst1)+0.5d0)*dx(1)
         X_cell(2) = x_lower(2)+(dble(ic_center(2)-ifirst2)+0.5d0)*dx(2)
c
c     Determine the standard interpolation stencil corresponding to the
c     position of X(s) within the cell.
c
         do d = 0,NDIM-1
            if ( X(d,s).lt.X_cell(d) ) then
               ic_lower(d) = ic_center(d)-4
               ic_upper(d) = ic_center(d)+3
            else
               ic_lower(d) = ic_center(d)-3
               ic_upper(d) = ic_center(d)+4
            endif
         enddo

         ic_lower(0) = max(ic_lower(0),ifirst0-nugc0)
         ic_upper(0) = min(ic_upper(0),ilast0 +nugc0)

         ic_lower(1) = max(ic_lower(1),ifirst1-nugc1)
         ic_upper(1) = min(ic_upper(1),ilast1 +nugc1)

         ic_lower(2) = max(ic_lower(2),ifirst2-nugc2)
         ic_upper(2) = min(ic_upper(2),ilast2 +nugc2)
c
c     Compute the standard interpolation weights.
c
C     DEC$ LOOP COUNT(8)
         do ic0 = ic_lower(0),ic_upper(0)
            X_cell(0) = x_lower(0)+(dble(ic0-ifirst0)+0.5d0)*dx(0)
            w0(ic0-ic_lower(0)) =
     &           lagrangian_wide8_ib_4_delta(
     &           (X(0,s)+Xshift(0,l)-X_cell(0))/dx(0))
         enddo
C     DEC$ LOOP COUNT(8)
         do ic1 = ic_lower(1),ic_upper(1)
            X_cell(1) = x_lower(1)+(dble(ic1-ifirst1)+0.5d0)*dx(1)
            w1(ic1-ic_lower(1)) =
     &           lagrangian_wide8_ib_4_delta(
     &           (X(1,s)+Xshift(1,l)-X_cell(1))/dx(1))
         enddo
C     DEC$ LOOP COUNT(8)
         do ic2 = ic_lower(2),ic_upper(2)
            X_cell(2) = x_lower(2)+(dble(ic2-ifirst2)+0.5d0)*dx(2)
            w2(ic2-ic_lower(2)) =
     &           lagrangian_wide8_ib_4_delta(
     &           (X(2,s)+Xshift(2,l)-X_cell(2))/dx(2))
         enddo
c
c     Interpolate u onto V.
c
         do d = 0,depth-1
            V(d,s) = 0.d0
C     DEC$ LOOP COUNT(8)
            do ic2 = ic_lower(2),ic_upper(2)
C     DEC$ LOOP COUNT(8)
               do ic1 = ic_lower(1),ic_upper(1)
C     DEC$ LOOP COUNT(8)
                  do ic0 = ic_lower(0),ic_upper(0)
                     V(d,s) = V(d,s)
     &                    +w0(ic0-ic_lower(0))
     &                    *w1(ic1-ic_lower(1))
     &                    *w2(ic2-ic_lower(2))
     &                    *u(ic0,ic1,ic2,d)
                  enddo
               enddo
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Spread V onto u at the positions specified by X using the
c     broadened 8-point version of the IB 4-point delta function using
c     standard (double) precision accumulation on the Cartesian grid.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_wide8_ib_4_spread3d(
     &     dx,x_lower,x_upper,depth,
     &     indices,Xshift,nindices,
     &     X,V,
     &     ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &     nugc0,nugc1,nugc2,
     &     u)
c
      implicit none
c
c     Functions.
c
      EXTERNAL lagrangian_floor
      INTEGER lagrangian_floor
      REAL lagrangian_wide8_ib_4_delta
c
c     Input.
c
      INTEGER depth
      INTEGER nindices
      INTEGER ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      INTEGER nugc0,nugc1,nugc2

      INTEGER indices(0:nindices-1)

      REAL Xshift(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1),x_lower(0:NDIM-1),x_upper(0:NDIM-1)
      REAL u(CELL3dVECG(ifirst,ilast,nugc),0:depth-1)
      REAL X(0:NDIM-1,0:*)
c
c     Input/Output.
c
      REAL V(0:depth-1,0:*)
c
c     Local variables.
c
      INTEGER ic0,ic1,ic2
      INTEGER ic_center(0:NDIM-1),ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER d,l,s

      REAL X_cell(0:NDIM-1),w0(0:7),w1(0:7),w2(0:7)
c
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
c
c     Use the broadened (8-point) version of the IB 4-point delta
c     function to spread V onto u.
c
      do l = 0,nindices-1
         s = indices(l)
c
c     Determine the Cartesian cell in which X(s) is located.
c
         ic_center(0) =
     &        lagrangian_floor((X(0,s)+Xshift(0,l)-x_lower(0))/dx(0))
     &        + ifirst0
         ic_center(1) =
     &        lagrangian_floor((X(1,s)+Xshift(1,l)-x_lower(1))/dx(1))
     &        + ifirst1
         ic_center(2) =
     &        lagrangian_floor((X(2,s)+Xshift(2,l)-x_lower(2))/dx(2))
     &        + ifirst2

         X_cell(0) = x_lower(0)+(dble(ic_center(0)-ifirst0)+0.5d0)*dx(0)
         X_cell(1) = x_lower(1)+(dble(ic_center(1)-ifirst1)+0.5d0)*dx(1)
         X_cell(2) = x_lower(2)+(dble(ic_center(2)-ifirst2)+0.5d0)*dx(2)
c
c     Determine the standard spreading stencil corresponding to the
c     position of X(s) within the cell.
c
         do d = 0,NDIM-1
            if ( X(d,s).lt.X_cell(d) ) then
               ic_lower(d) = ic_center(d)-4
               ic_upper(d) = ic_center(d)+3
            else
               ic_lower(d) = ic_center(d)-3
               ic_upper(d) = ic_center(d)+4
            endif
         enddo

         ic_lower(0) = max(ic_lower(0),ifirst0-nugc0)
         ic_upper(0) = min(ic_upper(0),ilast0 +nugc0)

         ic_lower(1) = max(ic_lower(1),ifirst1-nugc1)
         ic_upper(1) = min(ic_upper(1),ilast1 +nugc1)

         ic_lower(2) = max(ic_lower(2),ifirst2-nugc2)
         ic_upper(2) = min(ic_upper(2),ilast2 +nugc2)
c
c     Compute the standard spreading weights.
c
C     DEC$ LOOP COUNT(8)
         do ic0 = ic_lower(0),ic_upper(0)
            X_cell(0) = x_lower(0)+(dble(ic0-ifirst0)+0.5d0)*dx(0)
            w0(ic0-ic_lower(0)) =
     &           lagrangian_wide8_ib_4_delta(
     &           (X(0,s)+Xshift(0,l)-X_cell(0))/dx(0))
         enddo
C     DEC$ LOOP COUNT(8)
         do ic1 = ic_lower(1),ic_upper(1)
            X_cell(1) = x_lower(1)+(dble(ic1-ifirst1)+0.5d0)*dx(1)
            w1(ic1-ic_lower(1)) =
     &           lagrangian_wide8_ib_4_delta(
     &           (X(1,s)+Xshift(1,l)-X_cell(1))/dx(1))
         enddo
C     DEC$ LOOP COUNT(8)
         do ic2 = ic_lower(2),ic_upper(2)
            X_cell(2) = x_lower(2)+(dble(ic2-ifirst2)+0.5d0)*dx(2)
            w2(ic2-ic_lower(2)) =
     &           lagrangian_wide8_ib_4_delta(
     &           (X(2,s)+Xshift(2,l)-X_cell(2))/dx(2))
         enddo
c
c     Spread V onto u.
c
         do d = 0,depth-1
C     DEC$ LOOP COUNT(8)
            do ic2 = ic_lower(2),ic_upper(2)
C     DEC$ LOOP COUNT(8)
               do ic1 = ic_lower(1),ic_upper(1)
C     DEC$ LOOP COUNT(8)
                  do ic0 = ic_lower(0),ic_upper(0)
                     u(ic0,ic1,ic2,d) = u(ic0,ic1,ic2,d)+(
     &                    w0(ic0-ic_lower(0))*
     &                    w1(ic1-ic_lower(1))*
     &                    w2(ic2-ic_lower(2))*
     &                    V(d,s)/(dx(0)*dx(1)*dx(2)))
                  enddo
               enddo
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Interpolate u onto V at the positions specified by X using the
c     broadened 16-point version of the IB 4-point delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_wide16_ib_4_interp3d(
     &     dx,x_lower,x_upper,depth,
     &     ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &     nugc0,nugc1,nugc2,
     &     u,
     &     indices,Xshift,nindices,
     &     X,V)
c
      implicit none
c
c     Functions.
c
      EXTERNAL lagrangian_floor
      INTEGER lagrangian_floor
      REAL lagrangian_wide16_ib_4_delta
c
c     Input.
c
      INTEGER depth
      INTEGER ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      INTEGER nugc0,nugc1,nugc2
      INTEGER nindices

      INTEGER indices(0:nindices-1)

      REAL Xshift(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1),x_lower(0:NDIM-1),x_upper(0:NDIM-1)
      REAL u(CELL3dVECG(ifirst,ilast,nugc),0:depth-1)
      REAL X(0:NDIM-1,0:*)
c
c     Input/Output.
c
      REAL V(0:depth-1,0:*)
c
c     Local variables.
c
      INTEGER ic0,ic1,ic2
      INTEGER ic_center(0:NDIM-1),ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER d,l,s

      REAL X_cell(0:NDIM-1),w0(0:15),w1(0:15),w2(0:15)
c
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
c
c     Use the broadened (16-point) version of the IB 4-point delta
c     function to interpolate u onto V.
c
      do l = 0,nindices-1
         s = indices(l)
c
c     Determine the Cartesian cell in which X(s) is located.
c
         ic_center(0) =
     &        lagrangian_floor((X(0,s)+Xshift(0,l)-x_lower(0))/dx(0))
     &        + ifirst0
         ic_center(1) =
     &        lagrangian_floor((X(1,s)+Xshift(1,l)-x_lower(1))/dx(1))
     &        + ifirst1
         ic_center(2) =
     &        lagrangian_floor((X(2,s)+Xshift(2,l)-x_lower(2))/dx(2))
     &        + ifirst2

         X_cell(0) = x_lower(0)+(dble(ic_center(0)-ifirst0)+0.5d0)*dx(0)
         X_cell(1) = x_lower(1)+(dble(ic_center(1)-ifirst1)+0.5d0)*dx(1)
         X_cell(2) = x_lower(2)+(dble(ic_center(2)-ifirst2)+0.5d0)*dx(2)
c
c     Determine the standard interpolation stencil corresponding to the
c     position of X(s) within the cell.
c
         do d = 0,NDIM-1
            if ( X(d,s).lt.X_cell(d) ) then
               ic_lower(d) = ic_center(d)-8
               ic_upper(d) = ic_center(d)+7
            else
               ic_lower(d) = ic_center(d)-7
               ic_upper(d) = ic_center(d)+8
            endif
         enddo

         ic_lower(0) = max(ic_lower(0),ifirst0-nugc0)
         ic_upper(0) = min(ic_upper(0),ilast0 +nugc0)

         ic_lower(1) = max(ic_lower(1),ifirst1-nugc1)
         ic_upper(1) = min(ic_upper(1),ilast1 +nugc1)

         ic_lower(2) = max(ic_lower(2),ifirst2-nugc2)
         ic_upper(2) = min(ic_upper(2),ilast2 +nugc2)
c
c     Compute the standard interpolation weights.
c
C     DEC$ LOOP COUNT(16)
         do ic0 = ic_lower(0),ic_upper(0)
            X_cell(0) = x_lower(0)+(dble(ic0-ifirst0)+0.5d0)*dx(0)
            w0(ic0-ic_lower(0)) =
     &           lagrangian_wide16_ib_4_delta(
     &           (X(0,s)+Xshift(0,l)-X_cell(0))/dx(0))
         enddo
C     DEC$ LOOP COUNT(16)
         do ic1 = ic_lower(1),ic_upper(1)
            X_cell(1) = x_lower(1)+(dble(ic1-ifirst1)+0.5d0)*dx(1)
            w1(ic1-ic_lower(1)) =
     &           lagrangian_wide16_ib_4_delta(
     &           (X(1,s)+Xshift(1,l)-X_cell(1))/dx(1))
         enddo
C     DEC$ LOOP COUNT(16)
         do ic2 = ic_lower(2),ic_upper(2)
            X_cell(2) = x_lower(2)+(dble(ic2-ifirst2)+0.5d0)*dx(2)
            w2(ic2-ic_lower(2)) =
     &           lagrangian_wide16_ib_4_delta(
     &           (X(2,s)+Xshift(2,l)-X_cell(2))/dx(2))
         enddo
c
c     Interpolate u onto V.
c
         do d = 0,depth-1
            V(d,s) = 0.d0
C     DEC$ LOOP COUNT(16)
            do ic2 = ic_lower(2),ic_upper(2)
C     DEC$ LOOP COUNT(16)
               do ic1 = ic_lower(1),ic_upper(1)
C     DEC$ LOOP COUNT(16)
                  do ic0 = ic_lower(0),ic_upper(0)
                     V(d,s) = V(d,s)
     &                    +w0(ic0-ic_lower(0))
     &                    *w1(ic1-ic_lower(1))
     &                    *w2(ic2-ic_lower(2))
     &                    *u(ic0,ic1,ic2,d)
                  enddo
               enddo
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Spread V onto u at the positions specified by X using the
c     broadened 16-point version of the IB 4-point delta function using
c     standard (double) precision accumulation on the Cartesian grid.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_wide16_ib_4_spread3d(
     &     dx,x_lower,x_upper,depth,
     &     indices,Xshift,nindices,
     &     X,V,
     &     ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &     nugc0,nugc1,nugc2,
     &     u)
c
      implicit none
c
c     Functions.
c
      EXTERNAL lagrangian_floor
      INTEGER lagrangian_floor
      REAL lagrangian_wide16_ib_4_delta
c
c     Input.
c
      INTEGER depth
      INTEGER nindices
      INTEGER ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      INTEGER nugc0,nugc1,nugc2

      INTEGER indices(0:nindices-1)

      REAL Xshift(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1),x_lower(0:NDIM-1),x_upper(0:NDIM-1)
      REAL u(CELL3dVECG(ifirst,ilast,nugc),0:depth-1)
      REAL X(0:NDIM-1,0:*)
c
c     Input/Output.
c
      REAL V(0:depth-1,0:*)
c
c     Local variables.
c
      INTEGER ic0,ic1,ic2
      INTEGER ic_center(0:NDIM-1),ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER d,l,s

      REAL X_cell(0:NDIM-1),w0(0:7),w1(0:7),w2(0:7)
c
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
c
c     Use the broadened (8-point) version of the IB 4-point delta
c     function to spread V onto u.
c
      do l = 0,nindices-1
         s = indices(l)
c
c     Determine the Cartesian cell in which X(s) is located.
c
         ic_center(0) =
     &        lagrangian_floor((X(0,s)+Xshift(0,l)-x_lower(0))/dx(0))
     &        + ifirst0
         ic_center(1) =
     &        lagrangian_floor((X(1,s)+Xshift(1,l)-x_lower(1))/dx(1))
     &        + ifirst1
         ic_center(2) =
     &        lagrangian_floor((X(2,s)+Xshift(2,l)-x_lower(2))/dx(2))
     &        + ifirst2

         X_cell(0) = x_lower(0)+(dble(ic_center(0)-ifirst0)+0.5d0)*dx(0)
         X_cell(1) = x_lower(1)+(dble(ic_center(1)-ifirst1)+0.5d0)*dx(1)
         X_cell(2) = x_lower(2)+(dble(ic_center(2)-ifirst2)+0.5d0)*dx(2)
c
c     Determine the standard spreading stencil corresponding to the
c     position of X(s) within the cell.
c
         do d = 0,NDIM-1
            if ( X(d,s).lt.X_cell(d) ) then
               ic_lower(d) = ic_center(d)-4
               ic_upper(d) = ic_center(d)+3
            else
               ic_lower(d) = ic_center(d)-3
               ic_upper(d) = ic_center(d)+4
            endif
         enddo

         ic_lower(0) = max(ic_lower(0),ifirst0-nugc0)
         ic_upper(0) = min(ic_upper(0),ilast0 +nugc0)

         ic_lower(1) = max(ic_lower(1),ifirst1-nugc1)
         ic_upper(1) = min(ic_upper(1),ilast1 +nugc1)

         ic_lower(2) = max(ic_lower(2),ifirst2-nugc2)
         ic_upper(2) = min(ic_upper(2),ilast2 +nugc2)
c
c     Compute the standard spreading weights.
c
C     DEC$ LOOP COUNT(8)
         do ic0 = ic_lower(0),ic_upper(0)
            X_cell(0) = x_lower(0)+(dble(ic0-ifirst0)+0.5d0)*dx(0)
            w0(ic0-ic_lower(0)) =
     &           lagrangian_wide16_ib_4_delta(
     &           (X(0,s)+Xshift(0,l)-X_cell(0))/dx(0))
         enddo
C     DEC$ LOOP COUNT(8)
         do ic1 = ic_lower(1),ic_upper(1)
            X_cell(1) = x_lower(1)+(dble(ic1-ifirst1)+0.5d0)*dx(1)
            w1(ic1-ic_lower(1)) =
     &           lagrangian_wide16_ib_4_delta(
     &           (X(1,s)+Xshift(1,l)-X_cell(1))/dx(1))
         enddo
C     DEC$ LOOP COUNT(8)
         do ic2 = ic_lower(2),ic_upper(2)
            X_cell(2) = x_lower(2)+(dble(ic2-ifirst2)+0.5d0)*dx(2)
            w2(ic2-ic_lower(2)) =
     &           lagrangian_wide16_ib_4_delta(
     &           (X(2,s)+Xshift(2,l)-X_cell(2))/dx(2))
         enddo
c
c     Spread V onto u.
c
         do d = 0,depth-1
C     DEC$ LOOP COUNT(8)
            do ic2 = ic_lower(2),ic_upper(2)
C     DEC$ LOOP COUNT(8)
               do ic1 = ic_lower(1),ic_upper(1)
C     DEC$ LOOP COUNT(8)
                  do ic0 = ic_lower(0),ic_upper(0)
                     u(ic0,ic1,ic2,d) = u(ic0,ic1,ic2,d)+(
     &                    w0(ic0-ic_lower(0))*
     &                    w1(ic1-ic_lower(1))*
     &                    w2(ic2-ic_lower(2))*
     &                    V(d,s)/(dx(0)*dx(1)*dx(2)))
                  enddo
               enddo
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Interpolate u onto V at the positions specified by X using the IB
c     6-point delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_ib_6_interp3d(
     &     dx,x_lower,x_upper,depth,
     &     ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &     patch_touches_lower_physical_bdry,
     &     patch_touches_upper_physical_bdry,
     &     nugc0,nugc1,nugc2,
     &     u,
     &     indices,Xshift,nindices,
     &     X,V)
c
      implicit none
c
c     Functions.
c
      EXTERNAL lagrangian_floor
      INTEGER lagrangian_floor
      REAL lagrangian_ib_6_delta
c
c     Input.
c
      INTEGER depth
      INTEGER ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      INTEGER nugc0,nugc1,nugc2
      INTEGER nindices

      INTEGER patch_touches_lower_physical_bdry(0:NDIM-1)
      INTEGER patch_touches_upper_physical_bdry(0:NDIM-1)

      INTEGER indices(0:nindices-1)

      REAL Xshift(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1),x_lower(0:NDIM-1),x_upper(0:NDIM-1)
      REAL u(CELL3dVECG(ifirst,ilast,nugc),0:depth-1)
      REAL X(0:NDIM-1,0:*)
c
c     Input/Output.
c
      REAL V(0:depth-1,0:*)
c
c     Local variables.
c
      INTEGER ic0,ic1,ic2
      INTEGER ic_center(0:NDIM-1),ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER d,k,l,s

      REAL X_cell(0:NDIM-1),f(0:5),w0(0:5),w1(0:5),w2(0:5)

      LOGICAL touches_lower_bdry(0:NDIM-1)
      LOGICAL touches_upper_bdry(0:NDIM-1)
c
c     Use the IB 6-point delta function to interpolate u onto V.
c
      do l = 0,nindices-1
         s = indices(l)
c
c     Determine the Cartesian cell in which X(s) is located.
c
         ic_center(0) =
     &        lagrangian_floor((X(0,s)+Xshift(0,l)-x_lower(0))/dx(0))
     &        + ifirst0
         ic_center(1) =
     &        lagrangian_floor((X(1,s)+Xshift(1,l)-x_lower(1))/dx(1))
     &        + ifirst1
         ic_center(2) =
     &        lagrangian_floor((X(2,s)+Xshift(2,l)-x_lower(2))/dx(2))
     &        + ifirst2

         X_cell(0) = x_lower(0)+(dble(ic_center(0)-ifirst0)+0.5d0)*dx(0)
         X_cell(1) = x_lower(1)+(dble(ic_center(1)-ifirst1)+0.5d0)*dx(1)
         X_cell(2) = x_lower(2)+(dble(ic_center(2)-ifirst2)+0.5d0)*dx(2)
c
c     Determine the standard interpolation stencil corresponding to the
c     position of X(s) within the cell.
c
         do d = 0,NDIM-1
            if ( X(d,s).lt.X_cell(d) ) then
               ic_lower(d) = ic_center(d)-3
               ic_upper(d) = ic_center(d)+2
            else
               ic_lower(d) = ic_center(d)-2
               ic_upper(d) = ic_center(d)+3
            endif
         enddo

         ic_lower(0) = max(ic_lower(0),ifirst0-nugc0)
         ic_upper(0) = min(ic_upper(0),ilast0 +nugc0)

         ic_lower(1) = max(ic_lower(1),ifirst1-nugc1)
         ic_upper(1) = min(ic_upper(1),ilast1 +nugc1)

         ic_lower(2) = max(ic_lower(2),ifirst2-nugc2)
         ic_upper(2) = min(ic_upper(2),ilast2 +nugc2)
c
c     Compute the standard interpolation weights.
c
CDEC$ LOOP COUNT(6)
         do ic0 = ic_lower(0),ic_upper(0)
            X_cell(0) = x_lower(0)+(dble(ic0-ifirst0)+0.5d0)*dx(0)
            w0(ic0-ic_lower(0)) =
     &           lagrangian_ib_6_delta(
     &           (X(0,s)+Xshift(0,l)-X_cell(0))/dx(0))
         enddo
CDEC$ LOOP COUNT(6)
         do ic1 = ic_lower(1),ic_upper(1)
            X_cell(1) = x_lower(1)+(dble(ic1-ifirst1)+0.5d0)*dx(1)
            w1(ic1-ic_lower(1)) =
     &           lagrangian_ib_6_delta(
     &           (X(1,s)+Xshift(1,l)-X_cell(1))/dx(1))
         enddo
CDEC$ LOOP COUNT(6)
         do ic2 = ic_lower(2),ic_upper(2)
            X_cell(2) = x_lower(2)+(dble(ic2-ifirst2)+0.5d0)*dx(2)
            w2(ic2-ic_lower(2)) =
     &           lagrangian_ib_6_delta(
     &           (X(2,s)+Xshift(2,l)-X_cell(2))/dx(2))
         enddo
c
c     Determine whether special interpolation weights are needed to
c     handle physical boundary conditions.
c
         do d = 0,NDIM-1
            if ( (patch_touches_lower_physical_bdry(d).eq.1).and.
     &           (X(d,s) - x_lower(d) .lt. 2.5d0*dx(d)) ) then
               touches_lower_bdry(d) = .true.
            else
               touches_lower_bdry(d) = .false.
            endif
            if ( (patch_touches_upper_physical_bdry(d).eq.1).and.
     &           (x_upper(d) - X(d,s) .lt. 2.5d0*dx(d)) ) then
               touches_upper_bdry(d) = .true.
            else
               touches_upper_bdry(d) = .false.
            endif
         enddo
c
c     Modify the interpolation stencil and weights near physical
c     boundaries.
c
         if (touches_lower_bdry(0)) then
            call lagrangian_one_sided_ib_6_delta(
     &           w0,(X(0,s)-x_lower(0))/dx(0))
            ic_lower(0) = ifirst0
            ic_upper(0) = ifirst0+5
         elseif (touches_upper_bdry(0)) then
            call lagrangian_one_sided_ib_6_delta(
     &           f,(x_upper(0)-X(0,s))/dx(0))
            do k = 0,5
               w0(5-k) = f(k)
            enddo
            ic_lower(0) = ilast0-5
            ic_upper(0) = ilast0
         endif

         if (touches_lower_bdry(1)) then
            call lagrangian_one_sided_ib_6_delta(
     &           w1,(X(1,s)-x_lower(1))/dx(1))
            ic_lower(1) = ifirst1
            ic_upper(1) = ifirst1+5
         elseif (touches_upper_bdry(1)) then
            call lagrangian_one_sided_ib_6_delta(
     &           f,(x_upper(1)-X(1,s))/dx(1))
            do k = 0,5
               w1(5-k) = f(k)
            enddo
            ic_lower(1) = ilast1-5
            ic_upper(1) = ilast1
         endif

         if (touches_lower_bdry(2)) then
            call lagrangian_one_sided_ib_6_delta(
     &           w2,(X(2,s)-x_lower(2))/dx(2))
            ic_lower(2) = ifirst2
            ic_upper(2) = ifirst2+5
         elseif (touches_upper_bdry(2)) then
            call lagrangian_one_sided_ib_6_delta(
     &           f,(x_upper(2)-X(2,s))/dx(2))
            do k = 0,5
               w2(5-k) = f(k)
            enddo
            ic_lower(2) = ilast2-5
            ic_upper(2) = ilast2
         endif
c
c     Interpolate u onto V.
c
         do d = 0,depth-1
            V(d,s) = 0.d0
CDEC$ LOOP COUNT(6)
            do ic2 = ic_lower(2),ic_upper(2)
CDEC$ LOOP COUNT(6)
               do ic1 = ic_lower(1),ic_upper(1)
CDEC$ LOOP COUNT(6)
CDEC$ NOVECTOR
                  do ic0 = ic_lower(0),ic_upper(0)
                     V(d,s) = V(d,s)
     &                    +w0(ic0-ic_lower(0))
     &                    *w1(ic1-ic_lower(1))
     &                    *w2(ic2-ic_lower(2))
     &                    *u(ic0,ic1,ic2,d)
                  enddo
               enddo
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Spread V onto u at the positions specified by X using the IB
c     6-point delta function using standard (double) precision
c     accumulation on the Cartesian grid.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_ib_6_spread3d(
     &     dx,x_lower,x_upper,depth,
     &     indices,Xshift,nindices,
     &     X,V,
     &     ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &     patch_touches_lower_physical_bdry,
     &     patch_touches_upper_physical_bdry,
     &     nugc0,nugc1,nugc2,
     &     u)
c
      implicit none
c
c     Functions.
c
      EXTERNAL lagrangian_floor
      INTEGER lagrangian_floor
      REAL lagrangian_ib_6_delta
c
c     Input.
c
      INTEGER depth
      INTEGER nindices
      INTEGER ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      INTEGER nugc0,nugc1,nugc2

      INTEGER indices(0:nindices-1)

      INTEGER patch_touches_lower_physical_bdry(0:NDIM-1)
      INTEGER patch_touches_upper_physical_bdry(0:NDIM-1)

      REAL Xshift(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1),x_lower(0:NDIM-1),x_upper(0:NDIM-1)
      REAL u(CELL3dVECG(ifirst,ilast,nugc),0:depth-1)
      REAL X(0:NDIM-1,0:*)
c
c     Input/Output.
c
      REAL V(0:depth-1,0:*)
c
c     Local variables.
c
      INTEGER ic0,ic1,ic2
      INTEGER ic_center(0:NDIM-1),ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER d,k,l,s

      REAL X_cell(0:NDIM-1),f(0:5),w0(0:5),w1(0:5),w2(0:5)

      LOGICAL touches_lower_bdry(0:NDIM-1)
      LOGICAL touches_upper_bdry(0:NDIM-1)
c
c     Use the IB 6-point delta function to spread V onto u.
c
      do l = 0,nindices-1
         s = indices(l)
c
c     Determine the Cartesian cell in which X(s) is located.
c
         ic_center(0) =
     &        lagrangian_floor((X(0,s)+Xshift(0,l)-x_lower(0))/dx(0))
     &        + ifirst0
         ic_center(1) =
     &        lagrangian_floor((X(1,s)+Xshift(1,l)-x_lower(1))/dx(1))
     &        + ifirst1
         ic_center(2) =
     &        lagrangian_floor((X(2,s)+Xshift(2,l)-x_lower(2))/dx(2))
     &        + ifirst2

         X_cell(0) = x_lower(0)+(dble(ic_center(0)-ifirst0)+0.5d0)*dx(0)
         X_cell(1) = x_lower(1)+(dble(ic_center(1)-ifirst1)+0.5d0)*dx(1)
         X_cell(2) = x_lower(2)+(dble(ic_center(2)-ifirst2)+0.5d0)*dx(2)
c
c     Determine the standard spreading stencil corresponding to the
c     position of X(s) within the cell.
c
         do d = 0,NDIM-1
            if ( X(d,s).lt.X_cell(d) ) then
               ic_lower(d) = ic_center(d)-3
               ic_upper(d) = ic_center(d)+2
            else
               ic_lower(d) = ic_center(d)-2
               ic_upper(d) = ic_center(d)+3
            endif
         enddo

         ic_lower(0) = max(ic_lower(0),ifirst0-nugc0)
         ic_upper(0) = min(ic_upper(0),ilast0 +nugc0)

         ic_lower(1) = max(ic_lower(1),ifirst1-nugc1)
         ic_upper(1) = min(ic_upper(1),ilast1 +nugc1)

         ic_lower(2) = max(ic_lower(2),ifirst2-nugc2)
         ic_upper(2) = min(ic_upper(2),ilast2 +nugc2)
c
c     Compute the standard spreading weights.
c
CDEC$ LOOP COUNT(6)
         do ic0 = ic_lower(0),ic_upper(0)
            X_cell(0) = x_lower(0)+(dble(ic0-ifirst0)+0.5d0)*dx(0)
            w0(ic0-ic_lower(0)) =
     &           lagrangian_ib_6_delta(
     &           (X(0,s)+Xshift(0,l)-X_cell(0))/dx(0))
         enddo
CDEC$ LOOP COUNT(6)
         do ic1 = ic_lower(1),ic_upper(1)
            X_cell(1) = x_lower(1)+(dble(ic1-ifirst1)+0.5d0)*dx(1)
            w1(ic1-ic_lower(1)) =
     &           lagrangian_ib_6_delta(
     &           (X(1,s)+Xshift(1,l)-X_cell(1))/dx(1))
         enddo
CDEC$ LOOP COUNT(6)
         do ic2 = ic_lower(2),ic_upper(2)
            X_cell(2) = x_lower(2)+(dble(ic2-ifirst2)+0.5d0)*dx(2)
            w2(ic2-ic_lower(2)) =
     &           lagrangian_ib_6_delta(
     &           (X(2,s)+Xshift(2,l)-X_cell(2))/dx(2))
         enddo
c
c     Determine whether special spreading weights are needed to handle
c     physical boundary conditions.
c
         do d = 0,NDIM-1
            if ( (patch_touches_lower_physical_bdry(d).eq.1).and.
     &           (X(d,s) - x_lower(d) .lt. 2.5d0*dx(d)) ) then
               touches_lower_bdry(d) = .true.
            else
               touches_lower_bdry(d) = .false.
            endif
            if ( (patch_touches_upper_physical_bdry(d).eq.1).and.
     &           (x_upper(d) - X(d,s) .lt. 2.5d0*dx(d)) ) then
               touches_upper_bdry(d) = .true.
            else
               touches_upper_bdry(d) = .false.
            endif
         enddo
c
c     Modify the spreading stencil and weights near physical boundaries.
c
         if (touches_lower_bdry(0)) then
            call lagrangian_one_sided_ib_6_delta(
     &           w0,(X(0,s)-x_lower(0))/dx(0))
            ic_lower(0) = ifirst0
            ic_upper(0) = ifirst0+5
         elseif (touches_upper_bdry(0)) then
            call lagrangian_one_sided_ib_6_delta(
     &           f,(x_upper(0)-X(0,s))/dx(0))
            do k = 0,5
               w0(5-k) = f(k)
            enddo
            ic_lower(0) = ilast0-5
            ic_upper(0) = ilast0
         endif

         if (touches_lower_bdry(1)) then
            call lagrangian_one_sided_ib_6_delta(
     &           w1,(X(1,s)-x_lower(1))/dx(1))
            ic_lower(1) = ifirst1
            ic_upper(1) = ifirst1+5
         elseif (touches_upper_bdry(1)) then
            call lagrangian_one_sided_ib_6_delta(
     &           f,(x_upper(1)-X(1,s))/dx(1))
            do k = 0,5
               w1(5-k) = f(k)
            enddo
            ic_lower(1) = ilast1-5
            ic_upper(1) = ilast1
         endif

         if (touches_lower_bdry(2)) then
            call lagrangian_one_sided_ib_6_delta(
     &           w2,(X(2,s)-x_lower(2))/dx(2))
            ic_lower(2) = ifirst2
            ic_upper(2) = ifirst2+5
         elseif (touches_upper_bdry(2)) then
            call lagrangian_one_sided_ib_6_delta(
     &           f,(x_upper(2)-X(2,s))/dx(2))
            do k = 0,5
               w2(5-k) = f(k)
            enddo
            ic_lower(2) = ilast2-5
            ic_upper(2) = ilast2
         endif
c
c     Spread V onto u.
c
         do d = 0,depth-1
CDEC$ LOOP COUNT(6)
            do ic2 = ic_lower(2),ic_upper(2)
CDEC$ LOOP COUNT(6)
               do ic1 = ic_lower(1),ic_upper(1)
CDEC$ LOOP COUNT(6)
CDEC$ NOVECTOR
                  do ic0 = ic_lower(0),ic_upper(0)
                     u(ic0,ic1,ic2,d) = u(ic0,ic1,ic2,d)+(
     &                    w0(ic0-ic_lower(0))*
     &                    w1(ic1-ic_lower(1))*
     &                    w2(ic2-ic_lower(2))*
     &                    V(d,s)/(dx(0)*dx(1)*dx(2)))
                  enddo
               enddo
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
