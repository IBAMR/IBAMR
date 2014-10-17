c
c     Copyright (c) 2002-2014, Boyce Griffith
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
c        * Neither the name of The University of North Carolina nor the
c          names of its contributors may be used to endorse or promote
c          products derived from this software without specific prior
c          written permission.
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
define(NDIM,2)dnl
define(REAL,`double precision')dnl
define(INTEGER,`integer')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim2d.i)dnl
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Interpolate u onto V at the positions specified by X using the
c     piecewise constant delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_piecewise_constant_interp2d(
     &     dx,x_lower,x_upper,depth,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u,
     &     indices,Xshift,nindices,
     &     X,V)
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

      REAL Xshift(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1),x_lower(0:NDIM-1),x_upper(0:NDIM-1)
      REAL u(CELL2dVECG(ilower,iupper,nugc),0:depth-1)
      REAL X(0:NDIM-1,0:*)
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
      x_upper(0) = x_upper(0)
c
c     Use the piecewise constant delta function to interpolate u onto V.
c
      do l = 0,nindices-1
         s = indices(l)
c
c     Determine the Cartesian cell in which X(s) is located.
c
         ic0 = NINT((X(0,s)+Xshift(0,l)-x_lower(0))/dx(0)-0.5d0)+ilower0
         ic1 = NINT((X(1,s)+Xshift(1,l)-x_lower(1))/dx(1)-0.5d0)+ilower1
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
c     Spread V onto u at the positions specified by X using the
c     piecewise constant delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_piecewise_constant_spread2d(
     &     dx,x_lower,x_upper,depth,
     &     indices,Xshift,nindices,
     &     X,V,
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

      REAL Xshift(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1),x_lower(0:NDIM-1),x_upper(0:NDIM-1)
      REAL u(CELL2dVECG(ilower,iupper,nugc),0:depth-1)
      REAL X(0:NDIM-1,0:*)
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
      x_upper(0) = x_upper(0)
c
c     Use the piecewise constant delta function to spread V onto u.
c
      do l = 0,nindices-1
         s = indices(l)
c
c     Determine the Cartesian cell in which X(s) is located.
c
         ic0 = NINT((X(0,s)+Xshift(0,l)-x_lower(0))/dx(0)-0.5d0)+ilower0
         ic1 = NINT((X(1,s)+Xshift(1,l)-x_lower(1))/dx(1)-0.5d0)+ilower1
c
c     Spread V onto u.
c
         do d = 0,depth-1
            u(ic0,ic1,d) = u(ic0,ic1,d) + V(d,s)/(dx(0)*dx(1))
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
c     Interpolate u onto V at the positions specified by X using the
c     discontinuous linear delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_discontinuous_linear_interp2d(
     &     dx,x_lower,x_upper,depth,axis,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u,
     &     indices,Xshift,nindices,
     &     X,V)
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

      REAL Xshift(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1),x_lower(0:NDIM-1),x_upper(0:NDIM-1)
      REAL u(CELL2dVECG(ilower,iupper,nugc),0:depth-1)
      REAL X(0:NDIM-1,0:*)
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

      REAL X_cell(0:NDIM-1),X_shifted(0:NDIM-1),w(0:NDIM-1,0:1)
c
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
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
c     Shift the point by X_shift.
c
         do d = 0,NDIM-1
            X_shifted(d) = X(d,s)+Xshift(d,l)
         enddo
c
c     Determine the interpolation stencils and weights.
c
         do d = 0,NDIM-1
            ic_center(d) = ilower(d) +
     &           NINT((X_shifted(d)-x_lower(d))/dx(d)-0.5d0)
            X_cell(d) = x_lower(d) +
     &           (dble(ic_center(d)-ilower(d))+0.5d0)*dx(d)

            if ( d.eq.axis ) then
               if ( X_shifted(d).lt.X_cell(d) ) then
                  ic_lower(d) = ic_center(d)-1
                  ic_upper(d) = ic_center(d)
                  w(d,0) = (X_cell(d)-X_shifted(d))/dx(d)
                  w(d,1) = 1.d0 - w(d,0)
               else
                  ic_lower(d) = ic_center(d)
                  ic_upper(d) = ic_center(d)+1
                  w(d,0) = 1.d0 + (X_cell(d)-X_shifted(d))/dx(d)
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
c     Spread V onto u at the positions specified by X using the
c     discontinuous linear delta function using standard (double) precision
c     accumulation on the Cartesian grid.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_discontinuous_linear_spread2d(
     &     dx,x_lower,x_upper,depth,axis,
     &     indices,Xshift,nindices,
     &     X,V,
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

      REAL Xshift(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1),x_lower(0:NDIM-1),x_upper(0:NDIM-1)
      REAL u(CELL2dVECG(ilower,iupper,nugc),0:depth-1)
      REAL X(0:NDIM-1,0:*)
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

      REAL X_cell(0:NDIM-1),X_shifted(0:NDIM-1),w(0:NDIM-1,0:1)
c
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
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
c     Shift the point by X_shift.
c
         do d = 0,NDIM-1
            X_shifted(d) = X(d,s)+Xshift(d,l)
         enddo
c
c     Determine the interpolation stencils and weights.
c
         do d = 0,NDIM-1
            ic_center(d) = ilower(d) +
     &           NINT((X_shifted(d)-x_lower(d))/dx(d)-0.5d0)
            X_cell(d) = x_lower(d) +
     &           (dble(ic_center(d)-ilower(d))+0.5d0)*dx(d)

            if ( d.eq.axis ) then
               if ( X_shifted(d).lt.X_cell(d) ) then
                  ic_lower(d) = ic_center(d)-1
                  ic_upper(d) = ic_center(d)
                  w(d,0) = (X_cell(d)-X_shifted(d))/dx(d)
                  w(d,1) = 1.d0 - w(d,0)
               else
                  ic_lower(d) = ic_center(d)
                  ic_upper(d) = ic_center(d)+1
                  w(d,0) = 1.d0 + (X_cell(d)-X_shifted(d))/dx(d)
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
c     Interpolate u onto V at the positions specified by X using the
c     piecewise linear delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_piecewise_linear_interp2d(
     &     dx,x_lower,x_upper,depth,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u,
     &     indices,Xshift,nindices,
     &     X,V)
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

      REAL Xshift(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1),x_lower(0:NDIM-1),x_upper(0:NDIM-1)
      REAL u(CELL2dVECG(ilower,iupper,nugc),0:depth-1)
      REAL X(0:NDIM-1,0:*)
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

      REAL X_cell(0:NDIM-1),X_shifted(0:NDIM-1),w(0:NDIM-1,0:1)
c
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
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
c     Shift the point by X_shift.
c
         do d = 0,NDIM-1
            X_shifted(d) = X(d,s)+Xshift(d,l)
         enddo
c
c     Determine the interpolation stencils and weights.
c
         do d = 0,NDIM-1
            ic_center(d) = ilower(d) +
     &           NINT((X_shifted(d)-x_lower(d))/dx(d)-0.5d0)
            X_cell(d) = x_lower(d) +
     &           (dble(ic_center(d)-ilower(d))+0.5d0)*dx(d)

            if ( X_shifted(d).lt.X_cell(d) ) then
               ic_lower(d) = ic_center(d)-1
               ic_upper(d) = ic_center(d)
               w(d,0) = (X_cell(d)-X_shifted(d))/dx(d)
               w(d,1) = 1.d0 - w(d,0)
            else
               ic_lower(d) = ic_center(d)
               ic_upper(d) = ic_center(d)+1
               w(d,0) = 1.d0 + (X_cell(d)-X_shifted(d))/dx(d)
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
c     Spread V onto u at the positions specified by X using the
c     piecewise linear delta function using standard (double) precision
c     accumulation on the Cartesian grid.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_piecewise_linear_spread2d(
     &     dx,x_lower,x_upper,depth,
     &     indices,Xshift,nindices,
     &     X,V,
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

      REAL Xshift(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1),x_lower(0:NDIM-1),x_upper(0:NDIM-1)
      REAL u(CELL2dVECG(ilower,iupper,nugc),0:depth-1)
      REAL X(0:NDIM-1,0:*)
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

      REAL X_cell(0:NDIM-1),X_shifted(0:NDIM-1),w(0:NDIM-1,0:1)
c
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
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
c     Shift the point by X_shift.
c
         do d = 0,NDIM-1
            X_shifted(d) = X(d,s)+Xshift(d,l)
         enddo
c
c     Determine the interpolation stencils and weights.
c
         do d = 0,NDIM-1
            ic_center(d) = ilower(d) +
     &           NINT((X_shifted(d)-x_lower(d))/dx(d)-0.5d0)
            X_cell(d) = x_lower(d) +
     &           (dble(ic_center(d)-ilower(d))+0.5d0)*dx(d)

            if ( X_shifted(d).lt.X_cell(d) ) then
               ic_lower(d) = ic_center(d)-1
               ic_upper(d) = ic_center(d)
               w(d,0) = (X_cell(d)-X_shifted(d))/dx(d)
               w(d,1) = 1.d0 - w(d,0)
            else
               ic_lower(d) = ic_center(d)
               ic_upper(d) = ic_center(d)+1
               w(d,0) = 1.d0 + (X_cell(d)-X_shifted(d))/dx(d)
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
c     Interpolate u onto V at the positions specified by X using the
c     piecewise cubic delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_piecewise_cubic_interp2d(
     &     dx,x_lower,x_upper,depth,
     &     ilower0,iupper0,ilower1,iupper1,
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
      REAL lagrangian_piecewise_cubic_delta
c
c     Input.
c
      INTEGER depth
      INTEGER ilower0,iupper0,ilower1,iupper1
      INTEGER nugc0,nugc1
      INTEGER nindices

      INTEGER indices(0:nindices-1)

      REAL Xshift(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1),x_lower(0:NDIM-1),x_upper(0:NDIM-1)
      REAL u(CELL2dVECG(ilower,iupper,nugc),0:depth-1)
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
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
c
c     Use the piecewise cubic delta function to interpolate u onto V.
c
      do l = 0,nindices-1
         s = indices(l)
c
c     Determine the Cartesian cell in which X(s) is located.
c
         ic_center(0) =
     &        lagrangian_floor((X(0,s)+Xshift(0,l)-x_lower(0))/dx(0))
     &        + ilower0
         ic_center(1) =
     &        lagrangian_floor((X(1,s)+Xshift(1,l)-x_lower(1))/dx(1))
     &        + ilower1

         X_cell(0) = x_lower(0)+(dble(ic_center(0)-ilower0)+0.5d0)*dx(0)
         X_cell(1) = x_lower(1)+(dble(ic_center(1)-ilower1)+0.5d0)*dx(1)
c
c     Determine the interpolation stencil corresponding to the position
c     of X(s) within the cell.
c
         do d = 0,NDIM-1
            if ( X(d,s).lt.X_cell(d) ) then
               ic_lower(d) = ic_center(d)-2
               ic_upper(d) = ic_center(d)+1
            else
               ic_lower(d) = ic_center(d)-1
               ic_upper(d) = ic_center(d)+2
            endif
         enddo

         ic_lower(0) = max(ic_lower(0),ilower0-nugc0)
         ic_upper(0) = min(ic_upper(0),iupper0+nugc0)

         ic_lower(1) = max(ic_lower(1),ilower1-nugc1)
         ic_upper(1) = min(ic_upper(1),iupper1+nugc1)
c
c     Compute the interpolation weights.
c
         do ic0 = ic_lower(0),ic_upper(0)
            X_cell(0) = x_lower(0)+(dble(ic0-ilower0)+0.5d0)*dx(0)
            w0(ic0-ic_lower(0)) =
     &           lagrangian_piecewise_cubic_delta(
     &           (X(0,s)+Xshift(0,l)-X_cell(0))/dx(0))
         enddo

         do ic1 = ic_lower(1),ic_upper(1)
            X_cell(1) = x_lower(1)+(dble(ic1-ilower1)+0.5d0)*dx(1)
            w1(ic1-ic_lower(1)) =
     &           lagrangian_piecewise_cubic_delta(
     &           (X(1,s)+Xshift(1,l)-X_cell(1))/dx(1))
         enddo
c
c     Interpolate u onto V.
c
         do d = 0,depth-1
            V(d,s) = 0.d0
            do ic1 = ic_lower(1),ic_upper(1)
               do ic0 = ic_lower(0),ic_upper(0)
                  V(d,s) = V(d,s)
     &                 +w0(ic0-ic_lower(0))
     &                 *w1(ic1-ic_lower(1))
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
c     Spread V onto u at the positions specified by X using the
c     piecewise cubic delta function using standard (double) precision
c     accumulation on the Cartesian grid.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_piecewise_cubic_spread2d(
     &     dx,x_lower,x_upper,depth,
     &     indices,Xshift,nindices,
     &     X,V,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u)
c
      implicit none
c
c     Functions.
c
      EXTERNAL lagrangian_floor
      INTEGER lagrangian_floor
      REAL lagrangian_piecewise_cubic_delta
c
c     Input.
c
      INTEGER depth
      INTEGER nindices
      INTEGER ilower0,iupper0,ilower1,iupper1
      INTEGER nugc0,nugc1

      INTEGER indices(0:nindices-1)

      REAL Xshift(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1),x_lower(0:NDIM-1),x_upper(0:NDIM-1)
      REAL u(CELL2dVECG(ilower,iupper,nugc),0:depth-1)
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
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
c
c     Use the piecewise cubic delta function to spread V onto u.
c
      do l = 0,nindices-1
         s = indices(l)
c
c     Determine the Cartesian cell in which X(s) is located.
c
         ic_center(0) =
     &        lagrangian_floor((X(0,s)+Xshift(0,l)-x_lower(0))/dx(0))
     &        + ilower0
         ic_center(1) =
     &        lagrangian_floor((X(1,s)+Xshift(1,l)-x_lower(1))/dx(1))
     &        + ilower1

         X_cell(0) = x_lower(0)+(dble(ic_center(0)-ilower0)+0.5d0)*dx(0)
         X_cell(1) = x_lower(1)+(dble(ic_center(1)-ilower1)+0.5d0)*dx(1)
c
c     Determine the spreading stencil corresponding to the position of
c     X(s) within the cell.
c
         do d = 0,NDIM-1
            if ( X(d,s).lt.X_cell(d) ) then
               ic_lower(d) = ic_center(d)-2
               ic_upper(d) = ic_center(d)+1
            else
               ic_lower(d) = ic_center(d)-1
               ic_upper(d) = ic_center(d)+2
            endif
         enddo

         ic_lower(0) = max(ic_lower(0),ilower0-nugc0)
         ic_upper(0) = min(ic_upper(0),iupper0+nugc0)

         ic_lower(1) = max(ic_lower(1),ilower1-nugc1)
         ic_upper(1) = min(ic_upper(1),iupper1+nugc1)
c
c     Compute the spreading weights.
c
         do ic0 = ic_lower(0),ic_upper(0)
            X_cell(0) = x_lower(0)+(dble(ic0-ilower0)+0.5d0)*dx(0)
            w0(ic0-ic_lower(0)) =
     &           lagrangian_piecewise_cubic_delta(
     &           (X(0,s)+Xshift(0,l)-X_cell(0))/dx(0))
         enddo

         do ic1 = ic_lower(1),ic_upper(1)
            X_cell(1) = x_lower(1)+(dble(ic1-ilower1)+0.5d0)*dx(1)
            w1(ic1-ic_lower(1)) =
     &           lagrangian_piecewise_cubic_delta(
     &           (X(1,s)+Xshift(1,l)-X_cell(1))/dx(1))
         enddo
c
c     Spread V onto u.
c
         do d = 0,depth-1
            do ic1 = ic_lower(1),ic_upper(1)
               do ic0 = ic_lower(0),ic_upper(0)
                  u(ic0,ic1,d) = u(ic0,ic1,d)+(
     &                 w0(ic0-ic_lower(0))*
     &                 w1(ic1-ic_lower(1))*
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
c     Interpolate u onto V at the positions specified by X using the IB
c     3-point delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_ib_3_interp2d(
     &     dx,x_lower,x_upper,depth,
     &     ilower0,iupper0,ilower1,iupper1,
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
      REAL lagrangian_ib_3_delta
c
c     Input.
c
      INTEGER depth
      INTEGER ilower0,iupper0,ilower1,iupper1
      INTEGER nugc0,nugc1
      INTEGER nindices

      INTEGER indices(0:nindices-1)

      REAL Xshift(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1),x_lower(0:NDIM-1),x_upper(0:NDIM-1)
      REAL u(CELL2dVECG(ilower,iupper,nugc),0:depth-1)
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

      REAL X_cell(0:NDIM-1),w0(0:2),w1(0:2)
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
     &        + ilower0
         ic_center(1) =
     &        lagrangian_floor((X(1,s)+Xshift(1,l)-x_lower(1))/dx(1))
     &        + ilower1

         X_cell(0) = x_lower(0)+(dble(ic_center(0)-ilower0)+0.5d0)*dx(0)
         X_cell(1) = x_lower(1)+(dble(ic_center(1)-ilower1)+0.5d0)*dx(1)
c
c     Determine the interpolation stencil corresponding to the position
c     of X(s) within the cell.
c
         do d = 0,NDIM-1
            ic_lower(d) = ic_center(d)-1
            ic_upper(d) = ic_center(d)+1
         enddo

         ic_lower(0) = max(ic_lower(0),ilower0-nugc0)
         ic_upper(0) = min(ic_upper(0),iupper0+nugc0)

         ic_lower(1) = max(ic_lower(1),ilower1-nugc1)
         ic_upper(1) = min(ic_upper(1),iupper1+nugc1)
c
c     Compute the interpolation weights.
c
         do ic0 = ic_lower(0),ic_upper(0)
            X_cell(0) = x_lower(0)+(dble(ic0-ilower0)+0.5d0)*dx(0)
            w0(ic0-ic_lower(0)) =
     &           lagrangian_ib_3_delta(
     &           (X(0,s)+Xshift(0,l)-X_cell(0))/dx(0))
         enddo

         do ic1 = ic_lower(1),ic_upper(1)
            X_cell(1) = x_lower(1)+(dble(ic1-ilower1)+0.5d0)*dx(1)
            w1(ic1-ic_lower(1)) =
     &           lagrangian_ib_3_delta(
     &           (X(1,s)+Xshift(1,l)-X_cell(1))/dx(1))
         enddo
c
c     Interpolate u onto V.
c
         do d = 0,depth-1
            V(d,s) = 0.d0
            do ic1 = ic_lower(1),ic_upper(1)
               do ic0 = ic_lower(0),ic_upper(0)
                  V(d,s) = V(d,s)
     &                 +w0(ic0-ic_lower(0))
     &                 *w1(ic1-ic_lower(1))
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
c     Spread V onto u at the positions specified by X using the IB
c     3-point delta function using standard (double) precision
c     accumulation on the Cartesian grid.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_ib_3_spread2d(
     &     dx,x_lower,x_upper,depth,
     &     indices,Xshift,nindices,
     &     X,V,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
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
      INTEGER ilower0,iupper0,ilower1,iupper1
      INTEGER nugc0,nugc1

      INTEGER indices(0:nindices-1)

      REAL Xshift(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1),x_lower(0:NDIM-1),x_upper(0:NDIM-1)
      REAL u(CELL2dVECG(ilower,iupper,nugc),0:depth-1)
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

      REAL X_cell(0:NDIM-1),w0(0:2),w1(0:2)
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
     &        + ilower0
         ic_center(1) =
     &        lagrangian_floor((X(1,s)+Xshift(1,l)-x_lower(1))/dx(1))
     &        + ilower1

         X_cell(0) = x_lower(0)+(dble(ic_center(0)-ilower0)+0.5d0)*dx(0)
         X_cell(1) = x_lower(1)+(dble(ic_center(1)-ilower1)+0.5d0)*dx(1)
c
c     Determine the spreading stencil corresponding to the position of
c     X(s) within the cell.
c
         do d = 0,NDIM-1
            ic_lower(d) = ic_center(d)-1
            ic_upper(d) = ic_center(d)+1
         enddo

         ic_lower(0) = max(ic_lower(0),ilower0-nugc0)
         ic_upper(0) = min(ic_upper(0),iupper0+nugc0)

         ic_lower(1) = max(ic_lower(1),ilower1-nugc1)
         ic_upper(1) = min(ic_upper(1),iupper1+nugc1)
c
c     Compute the spreading weights.
c
         do ic0 = ic_lower(0),ic_upper(0)
            X_cell(0) = x_lower(0)+(dble(ic0-ilower0)+0.5d0)*dx(0)
            w0(ic0-ic_lower(0)) =
     &           lagrangian_ib_3_delta(
     &           (X(0,s)+Xshift(0,l)-X_cell(0))/dx(0))
         enddo

         do ic1 = ic_lower(1),ic_upper(1)
            X_cell(1) = x_lower(1)+(dble(ic1-ilower1)+0.5d0)*dx(1)
            w1(ic1-ic_lower(1)) =
     &           lagrangian_ib_3_delta(
     &           (X(1,s)+Xshift(1,l)-X_cell(1))/dx(1))
         enddo
c
c     Spread V onto u.
c
         do d = 0,depth-1
            do ic1 = ic_lower(1),ic_upper(1)
               do ic0 = ic_lower(0),ic_upper(0)
                  u(ic0,ic1,d) = u(ic0,ic1,d)+(
     &                 w0(ic0-ic_lower(0))*
     &                 w1(ic1-ic_lower(1))*
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
c     Interpolate u onto V at the positions specified by X using the IB
c     4-point delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_ib_4_interp2d(
     &     dx,x_lower,x_upper,depth,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u,
     &     indices,Xshift,nindices,
     &     X,V)
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

      REAL Xshift(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1),x_lower(0:NDIM-1),x_upper(0:NDIM-1)
      REAL u(CELL2dVECG(ilower,iupper,nugc),0:depth-1)
      REAL X(0:NDIM-1,0:*)
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
      INTEGER d,k,l,s

      REAL X_o_dx,q0,q1,r0,r1
      REAL w0(0:3),w1(0:3),f(0:3)
      REAL w(0:3,0:3),wy
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
c     of X(s) within the cell and compute the interpolation weights.
c
         X_o_dx = (X(0,s)+Xshift(0,l)-x_lower(0))/dx(0)
         ic_lower(0) = NINT(X_o_dx)+ilower0-2
         ic_upper(0) = ic_lower(0) + 3
         r0 = X_o_dx - ((ic_lower(0)+1-ilower0)+0.5d0)
         q0 = sqrt(1.d0+4.d0*r0*(1.d0-r0))
         w0(0) = 0.125d0*(3.d0-2.d0*r0-q0)
         w0(1) = 0.125d0*(3.d0-2.d0*r0+q0)
         w0(2) = 0.125d0*(1.d0+2.d0*r0+q0)
         w0(3) = 0.125d0*(1.d0+2.d0*r0-q0)

         X_o_dx = (X(1,s)+Xshift(1,l)-x_lower(1))/dx(1)
         ic_lower(1) = NINT(X_o_dx)+ilower1-2
         ic_upper(1) = ic_lower(1) + 3
         r1 = X_o_dx - ((ic_lower(1)+1-ilower1)+0.5d0)
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
c     Spread V onto u at the positions specified by X using the IB
c     4-point delta function using standard (double) precision
c     accumulation on the Cartesian grid.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_ib_4_spread2d(
     &     dx,x_lower,x_upper,depth,
     &     indices,Xshift,nindices,
     &     X,V,
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

      REAL Xshift(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1),x_lower(0:NDIM-1),x_upper(0:NDIM-1)
      REAL u(CELL2dVECG(ilower,iupper,nugc),0:depth-1)
      REAL X(0:NDIM-1,0:*)
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
      INTEGER d,k,l,s

      REAL X_o_dx,q0,q1,r0,r1
      REAL w0(0:3),w1(0:3),f(0:3)
      REAL w(0:3,0:3),wy
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
c     of X(s) within the cell and compute the interpolation weights.
c
         X_o_dx = (X(0,s)+Xshift(0,l)-x_lower(0))/dx(0)
         ic_lower(0) = NINT(X_o_dx)+ilower0-2
         ic_upper(0) = ic_lower(0) + 3
         r0 = X_o_dx - ((ic_lower(0)+1-ilower0)+0.5d0)
         q0 = sqrt(1.d0+4.d0*r0*(1.d0-r0))
         w0(0) = 0.125d0*(3.d0-2.d0*r0-q0)
         w0(1) = 0.125d0*(3.d0-2.d0*r0+q0)
         w0(2) = 0.125d0*(1.d0+2.d0*r0+q0)
         w0(3) = 0.125d0*(1.d0+2.d0*r0-q0)

         X_o_dx = (X(1,s)+Xshift(1,l)-x_lower(1))/dx(1)
         ic_lower(1) = NINT(X_o_dx)+ilower1-2
         ic_upper(1) = ic_lower(1) + 3
         r1 = X_o_dx - ((ic_lower(1)+1-ilower1)+0.5d0)
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
c     Interpolate u onto V at the positions specified by X using the IB
c     4-point delta function that has been broadened to have a support
c     of 8 meshwidths.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_ib_4_w8_interp2d(
     &     dx,x_lower,x_upper,depth,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u,
     &     indices,Xshift,nindices,
     &     X,V)
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

      REAL Xshift(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1),x_lower(0:NDIM-1),x_upper(0:NDIM-1)
      REAL u(CELL2dVECG(ilower,iupper,nugc),0:depth-1)
      REAL X(0:NDIM-1,0:*)
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

      REAL X_o_dx,q0,q1,r0,r1
      REAL w0(0:7),w1(0:7)
      REAL w(0:7,0:7),wy
c
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
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
c     of X(s) within the cell and compute the interpolation weights.
c
         X_o_dx = (X(0,s)+Xshift(0,l)-x_lower(0))/dx(0)
         ic_lower(0) = NINT(X_o_dx)+ilower0-4
         ic_upper(0) = ic_lower(0) + 7
         r0 = 0.5d0*(X_o_dx - ((ic_lower(0)+3-ilower0)+0.5d0))
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

         X_o_dx = (X(1,s)+Xshift(1,l)-x_lower(1))/dx(1)
         ic_lower(1) = NINT(X_o_dx)+ilower1-4
         ic_upper(1) = ic_lower(1) + 7
         r1 = 0.5d0*(X_o_dx - ((ic_lower(1)+3-ilower1)+0.5d0))
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
c     Spread V onto u at the positions specified by X using the IB
c     4-point delta function that has been broadened to have a support
c     of 8 meshwidths.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_ib_4_w8_spread2d(
     &     dx,x_lower,x_upper,depth,
     &     indices,Xshift,nindices,
     &     X,V,
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

      REAL Xshift(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1),x_lower(0:NDIM-1),x_upper(0:NDIM-1)
      REAL u(CELL2dVECG(ilower,iupper,nugc),0:depth-1)
      REAL X(0:NDIM-1,0:*)
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

      REAL X_o_dx,q0,q1,r0,r1
      REAL w0(0:7),w1(0:7)
      REAL w(0:7,0:7),wy
c
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
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
c     of X(s) within the cell and compute the interpolation weights.
c
         X_o_dx = (X(0,s)+Xshift(0,l)-x_lower(0))/dx(0)
         ic_lower(0) = NINT(X_o_dx)+ilower0-4
         ic_upper(0) = ic_lower(0) + 7
         r0 = 0.5d0*(X_o_dx - ((ic_lower(0)+3-ilower0)+0.5d0))
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

         X_o_dx = (X(1,s)+Xshift(1,l)-x_lower(1))/dx(1)
         ic_lower(1) = NINT(X_o_dx)+ilower1-4
         ic_upper(1) = ic_lower(1) + 7
         r1 = 0.5d0*(X_o_dx - ((ic_lower(1)+3-ilower1)+0.5d0))
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
c     Interpolate u onto V at the positions specified by X using the IB
c     6-point delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_ib_6_interp2d(
     &     dx,x_lower,x_upper,depth,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u,
     &     indices,Xshift,nindices,
     &     X,V)
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

      REAL Xshift(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1),x_lower(0:NDIM-1),x_upper(0:NDIM-1)
      REAL u(CELL2dVECG(ilower,iupper,nugc),0:depth-1)
      REAL X(0:NDIM-1,0:*)
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

      REAL X_o_dx,r,alpha,beta,gamma,discr,K
      REAL pm3,pm2,pm1,p,pp1,pp2
      REAL w0(0:5),w1(0:5)
      REAL w(0:5,0:5),wy

      PARAMETER (K = (59.d0/60.d0)*(1.d0-sqrt(1.d0-(3220.d0/3481.d0))))
c
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
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
c     of X(s) within the cell and compute the interpolation weights.
c
         X_o_dx = (X(0,s)+Xshift(0,l)-x_lower(0))/dx(0)
         ic_lower(0) = NINT(X_o_dx)+ilower0-3
         ic_upper(0) = ic_lower(0) + 5
         r = 1.d0 - X_o_dx + ((ic_lower(0)+2-ilower0)+0.5d0)

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

         X_o_dx = (X(1,s)+Xshift(1,l)-x_lower(1))/dx(1)
         ic_lower(1) = NINT(X_o_dx)+ilower1-3
         ic_upper(1) = ic_lower(1) + 5
         r = 1.d0 - X_o_dx + ((ic_lower(1)+2-ilower1)+0.5d0)

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
c     Spread V onto u at the positions specified by X using the IB
c     6-point delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_ib_6_spread2d(
     &     dx,x_lower,x_upper,depth,
     &     indices,Xshift,nindices,
     &     X,V,
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

      REAL Xshift(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1),x_lower(0:NDIM-1),x_upper(0:NDIM-1)
      REAL u(CELL2dVECG(ilower,iupper,nugc),0:depth-1)
      REAL X(0:NDIM-1,0:*)
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

      REAL X_o_dx,r,alpha,beta,gamma,discr,K
      REAL pm3,pm2,pm1,p,pp1,pp2
      REAL w0(0:5),w1(0:5)
      REAL w(0:5,0:5),wy

      PARAMETER (K = (59.d0/60.d0)*(1.d0-sqrt(1.d0-(3220.d0/3481.d0))))
c
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
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
c     of X(s) within the cell and compute the interpolation weights.
c
         X_o_dx = (X(0,s)+Xshift(0,l)-x_lower(0))/dx(0)
         ic_lower(0) = NINT(X_o_dx)+ilower0-3
         ic_upper(0) = ic_lower(0) + 5
         r = 1.d0 - X_o_dx + ((ic_lower(0)+2-ilower0)+0.5d0)

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

         X_o_dx = (X(1,s)+Xshift(1,l)-x_lower(1))/dx(1)
         ic_lower(1) = NINT(X_o_dx)+ilower1-3
         ic_upper(1) = ic_lower(1) + 5
         r = 1.d0 - X_o_dx + ((ic_lower(1)+2-ilower1)+0.5d0)

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
