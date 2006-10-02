define(NDIM,2)dnl
define(REAL,`double precision')dnl
define(INTEGER,`integer')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim2d.i)dnl
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Interpolate u onto V at the positions specified by X using the IB
c     4-point delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_ib4_interp2d(
     &     dx,x_lower,x_upper,depth,
     &     ifirst0,ilast0,ifirst1,ilast1,
     &     nugc0,nugc1,
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
      REAL lagrangian_ib4_delta
c
c     Input.
c
      INTEGER depth
      INTEGER ifirst0,ilast0,ifirst1,ilast1
      INTEGER nugc0,nugc1
      INTEGER nindices

      INTEGER indices(0:nindices-1)

      REAL Xshift(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1),x_lower(0:NDIM-1),x_upper(0:NDIM-1)
      REAL u(CELL2dVECG(ifirst,ilast,nugc),0:depth-1)
      REAL X(0:NDIM-1,0:*)
c
c     Input/Output.
c
      REAL V(0:depth-1,0:*)
c
c     Local variables.
c
      INTEGER ic0,ic1
      INTEGER ic_center(0:NDIM-1),ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER d,l,s

      REAL X_cell(0:NDIM-1),w0(0:3),w1(0:3)
c
      x_upper(0) = x_upper(0)   ! to prevent warnings about
      x_upper(1) = x_upper(1)   ! x_upper being unused.
c
c     Use the IB 4-point delta function to interpolate u onto V.
c
      do l = 0,nindices-1
         s = indices(l)
c
c     Determine the Cartesian cell which X(s) lies in and the
c     interpolation stencil corresponding to the position of X(s) within
c     the cell.
c
         ic_center(0) =
     &        lagrangian_floor((X(0,s)+Xshift(0,l)-x_lower(0))/dx(0))
     &        + ifirst0
         ic_center(1) =
     &        lagrangian_floor((X(1,s)+Xshift(1,l)-x_lower(1))/dx(1))
     &        + ifirst1

         X_cell(0) = x_lower(0)+(dble(ic_center(0)-ifirst0)+0.5d0)*dx(0)
         X_cell(1) = x_lower(1)+(dble(ic_center(1)-ifirst1)+0.5d0)*dx(1)

         do d = 0,NDIM-1
            if ( X(d,s).lt.X_cell(d) ) then
               ic_lower(d) = ic_center(d)-2
               ic_upper(d) = ic_center(d)+1
            else
               ic_lower(d) = ic_center(d)-1
               ic_upper(d) = ic_center(d)+2
            endif
         enddo

         ic_lower(0) = max(ic_lower(0),ifirst0-nugc0)
         ic_upper(0) = min(ic_upper(0),ilast0 +nugc0)

         ic_lower(1) = max(ic_lower(1),ifirst1-nugc1)
         ic_upper(1) = min(ic_upper(1),ilast1 +nugc1)
c
c     Compute the interpolation weights.
c
CDEC$ LOOP COUNT(4)
         do ic0 = ic_lower(0),ic_upper(0)
            X_cell(0) = x_lower(0)+(dble(ic0-ifirst0)+0.5d0)*dx(0)
            w0(ic0-ic_lower(0)) =
     &           lagrangian_ib4_delta(
     &           (X(0,s)+Xshift(0,l)-X_cell(0))/dx(0))
         enddo
CDEC$ LOOP COUNT(4)
         do ic1 = ic_lower(1),ic_upper(1)
            X_cell(1) = x_lower(1)+(dble(ic1-ifirst1)+0.5d0)*dx(1)
            w1(ic1-ic_lower(1)) =
     &           lagrangian_ib4_delta(
     &           (X(1,s)+Xshift(1,l)-X_cell(1))/dx(1))
         enddo
c
c     Interpolate u onto V.
c
         do d = 0,depth-1
            V(d,s) = 0.d0
CDEC$ LOOP COUNT(4)
            do ic1 = ic_lower(1),ic_upper(1)
CDEC$ LOOP COUNT(4)
               do ic0 = ic_lower(0),ic_upper(0)
                  V(d,s) = V(d,s)
     &                 +w0(ic0-ic_lower(0))*w1(ic1-ic_lower(1))
     &                 *u(ic0,ic1,d)
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
c     4-point delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_ib4_spread2d(
     &     dx,x_lower,x_upper,depth,
     &     indices,Xshift,nindices,
     &     X,V,
     &     ifirst0,ilast0,ifirst1,ilast1,
     &     nugc0,nugc1,
     &     u)
c
      implicit none
c
c     Functions.
c
      EXTERNAL lagrangian_floor
      INTEGER lagrangian_floor
      REAL lagrangian_ib4_delta
c
c     Input.
c
      INTEGER depth
      INTEGER nindices
      INTEGER ifirst0,ilast0,ifirst1,ilast1
      INTEGER nugc0,nugc1

      INTEGER indices(0:nindices-1)

      REAL Xshift(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1),x_lower(0:NDIM-1),x_upper(0:NDIM-1)
      REAL u(CELL2dVECG(ifirst,ilast,nugc),0:depth-1)
      REAL X(0:NDIM-1,0:*)
c
c     Input/Output.
c
      REAL V(0:depth-1,0:*)
c
c     Local variables.
c
      INTEGER ic0,ic1
      INTEGER ic_center(0:NDIM-1),ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER d,l,s

      REAL X_cell(0:NDIM-1),w0(0:3),w1(0:3)
c
      x_upper(0) = x_upper(0)   ! to prevent warnings about
      x_upper(1) = x_upper(1)   ! x_upper being unused.
c
c     Use the IB 4-point delta function to spread V onto u.
c
      do l = 0,nindices-1
         s = indices(l)
c
c     Determine the Cartesian cell which X(s) lies in and the
c     interpolation stencil corresponding to the position of X(s) within
c     the cell.
c
         ic_center(0) =
     &        lagrangian_floor((X(0,s)+Xshift(0,l)-x_lower(0))/dx(0))
     &        + ifirst0
         ic_center(1) =
     &        lagrangian_floor((X(1,s)+Xshift(1,l)-x_lower(1))/dx(1))
     &        + ifirst1

         X_cell(0) = x_lower(0)+(dble(ic_center(0)-ifirst0)+0.5d0)*dx(0)
         X_cell(1) = x_lower(1)+(dble(ic_center(1)-ifirst1)+0.5d0)*dx(1)

         do d = 0,NDIM-1
            if ( X(d,s).lt.X_cell(d) ) then
               ic_lower(d) = ic_center(d)-2
               ic_upper(d) = ic_center(d)+1
            else
               ic_lower(d) = ic_center(d)-1
               ic_upper(d) = ic_center(d)+2
            endif
         enddo

         ic_lower(0) = max(ic_lower(0),ifirst0-nugc0)
         ic_upper(0) = min(ic_upper(0),ilast0 +nugc0)

         ic_lower(1) = max(ic_lower(1),ifirst1-nugc1)
         ic_upper(1) = min(ic_upper(1),ilast1 +nugc1)
c
c     Compute the spreading weights.
c
CDEC$ LOOP COUNT(4)
         do ic0 = ic_lower(0),ic_upper(0)
            X_cell(0) = x_lower(0)+(dble(ic0-ifirst0)+0.5d0)*dx(0)
            w0(ic0-ic_lower(0)) =
     &           lagrangian_ib4_delta(
     &           (X(0,s)+Xshift(0,l)-X_cell(0))/dx(0))/dx(0)
         enddo
CDEC$ LOOP COUNT(4)
         do ic1 = ic_lower(1),ic_upper(1)
            X_cell(1) = x_lower(1)+(dble(ic1-ifirst1)+0.5d0)*dx(1)
            w1(ic1-ic_lower(1)) =
     &           lagrangian_ib4_delta(
     &           (X(1,s)+Xshift(1,l)-X_cell(1))/dx(1))/dx(1)
         enddo
c
c     Spread V onto u.
c
         do d = 0,depth-1
CDEC$ LOOP COUNT(4)
            do ic1 = ic_lower(1),ic_upper(1)
CDEC$ LOOP COUNT(4)
               do ic0 = ic_lower(0),ic_upper(0)
                  u(ic0,ic1,d) = u(ic0,ic1,d)
     &                 +w0(ic0-ic_lower(0))*w1(ic1-ic_lower(1))
     &                 *V(d,s)
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
c     4-point delta function, spread out over 8 meshwidths.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_wib4_interp2d(
     &     dx,x_lower,x_upper,depth,
     &     ifirst0,ilast0,ifirst1,ilast1,
     &     nugc0,nugc1,
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
      REAL lagrangian_wib4_delta
c
c     Input.
c
      INTEGER depth
      INTEGER ifirst0,ilast0,ifirst1,ilast1
      INTEGER nugc0,nugc1
      INTEGER nindices

      INTEGER indices(0:nindices-1)

      REAL Xshift(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1),x_lower(0:NDIM-1),x_upper(0:NDIM-1)
      REAL u(CELL2dVECG(ifirst,ilast,nugc),0:depth-1)
      REAL X(0:NDIM-1,0:*)
c
c     Input/Output.
c
      REAL V(0:depth-1,0:*)
c
c     Local variables.
c
      INTEGER ic0,ic1
      INTEGER ic_center(0:NDIM-1),ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER d,l,s

      REAL X_cell(0:NDIM-1),w0(0:7),w1(0:7)
c
      x_upper(0) = x_upper(0)   ! to prevent warnings about
      x_upper(1) = x_upper(1)   ! x_upper being unused.
c
c     Use the IB 4-point delta function to interpolate u onto V.
c
      do l = 0,nindices-1
         s = indices(l)
c
c     Determine the Cartesian cell which X(s) lies in and the
c     interpolation stencil corresponding to the position of X(s) within
c     the cell.
c
         ic_center(0) =
     &        lagrangian_floor((X(0,s)+Xshift(0,l)-x_lower(0))/dx(0))
     &        + ifirst0
         ic_center(1) =
     &        lagrangian_floor((X(1,s)+Xshift(1,l)-x_lower(1))/dx(1))
     &        + ifirst1

         X_cell(0) = x_lower(0)+(dble(ic_center(0)-ifirst0)+0.5d0)*dx(0)
         X_cell(1) = x_lower(1)+(dble(ic_center(1)-ifirst1)+0.5d0)*dx(1)

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
c
c     Compute the interpolation weights.
c
CDEC$ LOOP COUNT(8)
         do ic0 = ic_lower(0),ic_upper(0)
            X_cell(0) = x_lower(0)+(dble(ic0-ifirst0)+0.5d0)*dx(0)
            w0(ic0-ic_lower(0)) =
     &           lagrangian_wib4_delta(
     &           (X(0,s)+Xshift(0,l)-X_cell(0))/dx(0))
         enddo
CDEC$ LOOP COUNT(8)
         do ic1 = ic_lower(1),ic_upper(1)
            X_cell(1) = x_lower(1)+(dble(ic1-ifirst1)+0.5d0)*dx(1)
            w1(ic1-ic_lower(1)) =
     &           lagrangian_wib4_delta(
     &           (X(1,s)+Xshift(1,l)-X_cell(1))/dx(1))
         enddo
c
c     Interpolate u onto V.
c
         do d = 0,depth-1
            V(d,s) = 0.d0
CDEC$ LOOP COUNT(8)
            do ic1 = ic_lower(1),ic_upper(1)
CDEC$ LOOP COUNT(8)
CDEC$ NOVECTOR
               do ic0 = ic_lower(0),ic_upper(0)
                  V(d,s) = V(d,s)
     &                 +w0(ic0-ic_lower(0))*w1(ic1-ic_lower(1))
     &                 *u(ic0,ic1,d)
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
c     4-point delta function, spread out over 8 meshwidths.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_wib4_spread2d(
     &     dx,x_lower,x_upper,depth,
     &     indices,Xshift,nindices,
     &     X,V,
     &     ifirst0,ilast0,ifirst1,ilast1,
     &     nugc0,nugc1,
     &     u)
c
      implicit none
c
c     Functions.
c
      EXTERNAL lagrangian_floor
      INTEGER lagrangian_floor
      REAL lagrangian_wib4_delta
c
c     Input.
c
      INTEGER depth
      INTEGER nindices
      INTEGER ifirst0,ilast0,ifirst1,ilast1
      INTEGER nugc0,nugc1

      INTEGER indices(0:nindices-1)

      REAL Xshift(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1),x_lower(0:NDIM-1),x_upper(0:NDIM-1)
      REAL u(CELL2dVECG(ifirst,ilast,nugc),0:depth-1)
      REAL X(0:NDIM-1,0:*)
c
c     Input/Output.
c
      REAL V(0:depth-1,0:*)
c
c     Local variables.
c
      INTEGER ic0,ic1
      INTEGER ic_center(0:NDIM-1),ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER d,l,s

      REAL X_cell(0:NDIM-1),w0(0:7),w1(0:7)
c
      x_upper(0) = x_upper(0)   ! to prevent warnings about
      x_upper(1) = x_upper(1)   ! x_upper being unused.
c
c     Use the IB 4-point delta function to spread V onto u.
c
      do l = 0,nindices-1
         s = indices(l)
c
c     Determine the Cartesian cell which X(s) lies in and the
c     interpolation stencil corresponding to the position of X(s) within
c     the cell.
c
         ic_center(0) =
     &        lagrangian_floor((X(0,s)+Xshift(0,l)-x_lower(0))/dx(0))
     &        + ifirst0
         ic_center(1) =
     &        lagrangian_floor((X(1,s)+Xshift(1,l)-x_lower(1))/dx(1))
     &        + ifirst1

         X_cell(0) = x_lower(0)+(dble(ic_center(0)-ifirst0)+0.5d0)*dx(0)
         X_cell(1) = x_lower(1)+(dble(ic_center(1)-ifirst1)+0.5d0)*dx(1)

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
c
c     Compute the spreading weights.
c
CDEC$ LOOP COUNT(8)
         do ic0 = ic_lower(0),ic_upper(0)
            X_cell(0) = x_lower(0)+(dble(ic0-ifirst0)+0.5d0)*dx(0)
            w0(ic0-ic_lower(0)) =
     &           lagrangian_wib4_delta(
     &           (X(0,s)+Xshift(0,l)-X_cell(0))/dx(0))/dx(0)
         enddo
CDEC$ LOOP COUNT(8)
         do ic1 = ic_lower(1),ic_upper(1)
            X_cell(1) = x_lower(1)+(dble(ic1-ifirst1)+0.5d0)*dx(1)
            w1(ic1-ic_lower(1)) =
     &           lagrangian_wib4_delta(
     &           (X(1,s)+Xshift(1,l)-X_cell(1))/dx(1))/dx(1)
         enddo
c
c     Spread V onto u.
c
         do d = 0,depth-1
CDEC$ LOOP COUNT(8)
            do ic1 = ic_lower(1),ic_upper(1)
CDEC$ LOOP COUNT(8)
CDEC$ NOVECTOR
               do ic0 = ic_lower(0),ic_upper(0)
                  u(ic0,ic1,d) = u(ic0,ic1,d)
     &                 +w0(ic0-ic_lower(0))*w1(ic1-ic_lower(1))
     &                 *V(d,s)
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
      subroutine lagrangian_ib6_interp2d(
     &     dx,x_lower,x_upper,depth,
     &     ifirst0,ilast0,ifirst1,ilast1,
     &     nugc0,nugc1,
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
      REAL lagrangian_ib6_delta
c
c     Input.
c
      INTEGER depth
      INTEGER ifirst0,ilast0,ifirst1,ilast1
      INTEGER nugc0,nugc1
      INTEGER nindices

      INTEGER indices(0:nindices-1)

      REAL Xshift(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1),x_lower(0:NDIM-1),x_upper(0:NDIM-1)
      REAL u(CELL2dVECG(ifirst,ilast,nugc),0:depth-1)
      REAL X(0:NDIM-1,0:*)
c
c     Input/Output.
c
      REAL V(0:depth-1,0:*)
c
c     Local variables.
c
      INTEGER ic0,ic1
      INTEGER ic_center(0:NDIM-1),ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER d,l,s

      REAL X_cell(0:NDIM-1),w0(0:5),w1(0:5)
c
      x_upper(0) = x_upper(0)   ! to prevent warnings about
      x_upper(1) = x_upper(1)   ! x_upper being unused.
c
c     Use the IB 6-point delta function to interpolate u onto V.
c
      do l = 0,nindices-1
         s = indices(l)
c
c     Determine the Cartesian cell which X(s) lies in and the
c     interpolation stencil corresponding to the position of X(s) within
c     the cell.
c
         ic_center(0) =
     &        lagrangian_floor((X(0,s)+Xshift(0,l)-x_lower(0))/dx(0))
     &        + ifirst0
         ic_center(1) =
     &        lagrangian_floor((X(1,s)+Xshift(1,l)-x_lower(1))/dx(1))
     &        + ifirst1

         X_cell(0) = x_lower(0)+(dble(ic_center(0)-ifirst0)+0.5d0)*dx(0)
         X_cell(1) = x_lower(1)+(dble(ic_center(1)-ifirst1)+0.5d0)*dx(1)

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
c
c     Compute the interpolation weights.
c
CDEC$ LOOP COUNT(6)
         do ic0 = ic_lower(0),ic_upper(0)
            X_cell(0) = x_lower(0)+(dble(ic0-ifirst0)+0.5d0)*dx(0)
            w0(ic0-ic_lower(0)) =
     &           lagrangian_ib6_delta(
     &           (X(0,s)+Xshift(0,l)-X_cell(0))/dx(0))
         enddo
CDEC$ LOOP COUNT(6)
         do ic1 = ic_lower(1),ic_upper(1)
            X_cell(1) = x_lower(1)+(dble(ic1-ifirst1)+0.5d0)*dx(1)
            w1(ic1-ic_lower(1)) =
     &           lagrangian_ib6_delta(
     &           (X(1,s)+Xshift(1,l)-X_cell(1))/dx(1))
         enddo
c
c     Interpolate u onto V.
c
         do d = 0,depth-1
            V(d,s) = 0.d0
CDEC$ LOOP COUNT(6)
            do ic1 = ic_lower(1),ic_upper(1)
CDEC$ LOOP COUNT(6)
CDEC$ NOVECTOR
               do ic0 = ic_lower(0),ic_upper(0)
                  V(d,s) = V(d,s)
     &                 +w0(ic0-ic_lower(0))*w1(ic1-ic_lower(1))
     &                 *u(ic0,ic1,d)
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
c     6-point delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_ib6_spread2d(
     &     dx,x_lower,x_upper,depth,
     &     indices,Xshift,nindices,
     &     X,V,
     &     ifirst0,ilast0,ifirst1,ilast1,
     &     nugc0,nugc1,
     &     u)
c
      implicit none
c
c     Functions.
c
      EXTERNAL lagrangian_floor
      INTEGER lagrangian_floor
      REAL lagrangian_ib6_delta
c
c     Input.
c
      INTEGER depth
      INTEGER nindices
      INTEGER ifirst0,ilast0,ifirst1,ilast1
      INTEGER nugc0,nugc1

      INTEGER indices(0:nindices-1)

      REAL Xshift(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1),x_lower(0:NDIM-1),x_upper(0:NDIM-1)
      REAL u(CELL2dVECG(ifirst,ilast,nugc),0:depth-1)
      REAL X(0:NDIM-1,0:*)
c
c     Input/Output.
c
      REAL V(0:depth-1,0:*)
c
c     Local variables.
c
      INTEGER ic0,ic1
      INTEGER ic_center(0:NDIM-1),ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER d,l,s

      REAL X_cell(0:NDIM-1),w0(0:5),w1(0:5)
c
      x_upper(0) = x_upper(0)   ! to prevent warnings about
      x_upper(1) = x_upper(1)   ! x_upper being unused.
c
c     Use the IB 6-point delta function to spread V onto u.
c
      do l = 0,nindices-1
         s = indices(l)
c
c     Determine the Cartesian cell which X(s) lies in and the
c     interpolation stencil corresponding to the position of X(s) within
c     the cell.
c
         ic_center(0) =
     &        lagrangian_floor((X(0,s)+Xshift(0,l)-x_lower(0))/dx(0))
     &        + ifirst0
         ic_center(1) =
     &        lagrangian_floor((X(1,s)+Xshift(1,l)-x_lower(1))/dx(1))
     &        + ifirst1

         X_cell(0) = x_lower(0)+(dble(ic_center(0)-ifirst0)+0.5d0)*dx(0)
         X_cell(1) = x_lower(1)+(dble(ic_center(1)-ifirst1)+0.5d0)*dx(1)

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
c
c     Compute the spreading weights.
c
CDEC$ LOOP COUNT(6)
         do ic0 = ic_lower(0),ic_upper(0)
            X_cell(0) = x_lower(0)+(dble(ic0-ifirst0)+0.5d0)*dx(0)
            w0(ic0-ic_lower(0)) =
     &           lagrangian_ib6_delta(
     &           (X(0,s)+Xshift(0,l)-X_cell(0))/dx(0))/dx(0)
         enddo
CDEC$ LOOP COUNT(6)
         do ic1 = ic_lower(1),ic_upper(1)
            X_cell(1) = x_lower(1)+(dble(ic1-ifirst1)+0.5d0)*dx(1)
            w1(ic1-ic_lower(1)) =
     &           lagrangian_ib6_delta(
     &           (X(1,s)+Xshift(1,l)-X_cell(1))/dx(1))/dx(1)
         enddo
c
c     Spread V onto u.
c
         do d = 0,depth-1
CDEC$ LOOP COUNT(6)
            do ic1 = ic_lower(1),ic_upper(1)
CDEC$ LOOP COUNT(6)
CDEC$ NOVECTOR
               do ic0 = ic_lower(0),ic_upper(0)
                  u(ic0,ic1,d) = u(ic0,ic1,d)
     &                 +w0(ic0-ic_lower(0))*w1(ic1-ic_lower(1))
     &                 *V(d,s)
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
c     Interpolate u onto V at the positions specified by X using
c     piecewise cubic interpolation.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_pwcubic_interp2d(
     &     dx,x_lower,x_upper,depth,
     &     ifirst0,ilast0,ifirst1,ilast1,
     &     nugc0,nugc1,
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
      REAL lagrangian_pwcubic_delta
c
c     Input.
c
      INTEGER depth
      INTEGER ifirst0,ilast0,ifirst1,ilast1
      INTEGER nugc0,nugc1
      INTEGER nindices

      INTEGER indices(0:nindices-1)

      REAL Xshift(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1),x_lower(0:NDIM-1),x_upper(0:NDIM-1)
      REAL u(CELL2dVECG(ifirst,ilast,nugc),0:depth-1)
      REAL X(0:NDIM-1,0:*)
c
c     Input/Output.
c
      REAL V(0:depth-1,0:*)
c
c     Local variables.
c
      INTEGER ic0,ic1
      INTEGER ic_center(0:NDIM-1),ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER d,l,s

      REAL X_cell(0:NDIM-1),w0(0:3),w1(0:3)
c
      x_upper(0) = x_upper(0)   ! to prevent warnings about
      x_upper(1) = x_upper(1)   ! x_upper being unused.
c
c     Use picewise cubic interpolation to interpolate u onto V.
c
      do l = 0,nindices-1
         s = indices(l)
c
c     Determine the Cartesian cell which X(s) lies in and the
c     interpolation stencil corresponding to the position of X(s) within
c     the cell.
c
         ic_center(0) =
     &        lagrangian_floor((X(0,s)+Xshift(0,l)-x_lower(0))/dx(0))
     &        + ifirst0
         ic_center(1) =
     &        lagrangian_floor((X(1,s)+Xshift(1,l)-x_lower(1))/dx(1))
     &        + ifirst1

         X_cell(0) = x_lower(0)+(dble(ic_center(0)-ifirst0)+0.5d0)*dx(0)
         X_cell(1) = x_lower(1)+(dble(ic_center(1)-ifirst1)+0.5d0)*dx(1)

         do d = 0,NDIM-1
            if ( X(d,s).lt.X_cell(d) ) then
               ic_lower(d) = ic_center(d)-2
               ic_upper(d) = ic_center(d)+1
            else
               ic_lower(d) = ic_center(d)-1
               ic_upper(d) = ic_center(d)+2
            endif
         enddo

         ic_lower(0) = max(ic_lower(0),ifirst0-nugc0)
         ic_upper(0) = min(ic_upper(0),ilast0 +nugc0)

         ic_lower(1) = max(ic_lower(1),ifirst1-nugc1)
         ic_upper(1) = min(ic_upper(1),ilast1 +nugc1)
c
c     Compute the interpolation weights.
c
CDEC$ LOOP COUNT(4)
         do ic0 = ic_lower(0),ic_upper(0)
            X_cell(0) = x_lower(0)+(dble(ic0-ifirst0)+0.5d0)*dx(0)
            w0(ic0-ic_lower(0)) =
     &           lagrangian_pwcubic_delta(
     &           (X(0,s)+Xshift(0,l)-X_cell(0))/dx(0))
         enddo
CDEC$ LOOP COUNT(4)
         do ic1 = ic_lower(1),ic_upper(1)
            X_cell(1) = x_lower(1)+(dble(ic1-ifirst1)+0.5d0)*dx(1)
            w1(ic1-ic_lower(1)) =
     &           lagrangian_pwcubic_delta(
     &           (X(1,s)+Xshift(1,l)-X_cell(1))/dx(1))
         enddo
c
c     Interpolate u onto V.
c
         do d = 0,depth-1
            V(d,s) = 0.d0
CDEC$ LOOP COUNT(4)
            do ic1 = ic_lower(1),ic_upper(1)
CDEC$ LOOP COUNT(4)
               do ic0 = ic_lower(0),ic_upper(0)
                  V(d,s) = V(d,s)
     &                 +w0(ic0-ic_lower(0))*w1(ic1-ic_lower(1))
     &                 *u(ic0,ic1,d)
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
c     Spread V onto u at the positions specified by X using a piecewise
c     cubic delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_pwcubic_spread2d(
     &     dx,x_lower,x_upper,depth,
     &     indices,Xshift,nindices,
     &     X,V,
     &     ifirst0,ilast0,ifirst1,ilast1,
     &     nugc0,nugc1,
     &     u)
c
      implicit none
c
c     Functions.
c
      EXTERNAL lagrangian_floor
      INTEGER lagrangian_floor
      REAL lagrangian_pwcubic_delta
c
c     Input.
c
      INTEGER depth
      INTEGER nindices
      INTEGER ifirst0,ilast0,ifirst1,ilast1
      INTEGER nugc0,nugc1

      INTEGER indices(0:nindices-1)

      REAL Xshift(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1),x_lower(0:NDIM-1),x_upper(0:NDIM-1)
      REAL u(CELL2dVECG(ifirst,ilast,nugc),0:depth-1)
      REAL X(0:NDIM-1,0:*)
c
c     Input/Output.
c
      REAL V(0:depth-1,0:*)
c
c     Local variables.
c
      INTEGER ic0,ic1
      INTEGER ic_center(0:NDIM-1),ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER d,l,s

      REAL X_cell(0:NDIM-1),w0(0:3),w1(0:3)
c
      x_upper(0) = x_upper(0)   ! to prevent warnings about
      x_upper(1) = x_upper(1)   ! x_upper being unused.
c
c     Use picewise cubic interpolation to spread V onto u.
c
      do l = 0,nindices-1
         s = indices(l)
c
c     Determine the Cartesian cell which X(s) lies in and the
c     interpolation stencil corresponding to the position of X(s) within
c     the cell.
c
         ic_center(0) =
     &        lagrangian_floor((X(0,s)+Xshift(0,l)-x_lower(0))/dx(0))
     &        + ifirst0
         ic_center(1) =
     &        lagrangian_floor((X(1,s)+Xshift(1,l)-x_lower(1))/dx(1))
     &        + ifirst1

         X_cell(0) = x_lower(0)+(dble(ic_center(0)-ifirst0)+0.5d0)*dx(0)
         X_cell(1) = x_lower(1)+(dble(ic_center(1)-ifirst1)+0.5d0)*dx(1)

         do d = 0,NDIM-1
            if ( X(d,s).lt.X_cell(d) ) then
               ic_lower(d) = ic_center(d)-2
               ic_upper(d) = ic_center(d)+1
            else
               ic_lower(d) = ic_center(d)-1
               ic_upper(d) = ic_center(d)+2
            endif
         enddo

         ic_lower(0) = max(ic_lower(0),ifirst0-nugc0)
         ic_upper(0) = min(ic_upper(0),ilast0 +nugc0)

         ic_lower(1) = max(ic_lower(1),ifirst1-nugc1)
         ic_upper(1) = min(ic_upper(1),ilast1 +nugc1)
c
c     Compute the spreading weights.
c
CDEC$ LOOP COUNT(4)
         do ic0 = ic_lower(0),ic_upper(0)
            X_cell(0) = x_lower(0)+(dble(ic0-ifirst0)+0.5d0)*dx(0)
            w0(ic0-ic_lower(0)) =
     &           lagrangian_pwcubic_delta(
     &           (X(0,s)+Xshift(0,l)-X_cell(0))/dx(0))/dx(0)
         enddo
CDEC$ LOOP COUNT(4)
         do ic1 = ic_lower(1),ic_upper(1)
            X_cell(1) = x_lower(1)+(dble(ic1-ifirst1)+0.5d0)*dx(1)
            w1(ic1-ic_lower(1)) =
     &           lagrangian_pwcubic_delta(
     &           (X(1,s)+Xshift(1,l)-X_cell(1))/dx(1))/dx(1)
         enddo
c
c     Spread V onto u.
c
         do d = 0,depth-1
CDEC$ LOOP COUNT(4)
            do ic1 = ic_lower(1),ic_upper(1)
CDEC$ LOOP COUNT(4)
               do ic0 = ic_lower(0),ic_upper(0)
                  u(ic0,ic1,d) = u(ic0,ic1,d)
     &                 +w0(ic0-ic_lower(0))*w1(ic1-ic_lower(1))
     &                 *V(d,s)
               enddo
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
