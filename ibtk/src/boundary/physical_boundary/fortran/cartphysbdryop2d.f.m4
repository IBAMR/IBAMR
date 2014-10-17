c
c     Routines to set physical boundary condition values.
c
c     Created on 21 May 2007 by Boyce Griffith
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
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     For cell-centered values, we follow a similar approach as that
c     implemented in class SAMRAI::solv::CartesianRobinBcHelper.  Let
c     u_g denote the ghost cell value and let u_i denote the
c     mirror-image interior cell value, and let n be the number of cell
c     widths separating the ghost cell center and the interior cell
c     center.  We define
c
c     u_b = (u_g + u_i)/2
c     u_n = (u_g - u_i)/(n*h)
c
c     If
c
c     a*u_b + b*u_n = g
c
c     then
c
c     u_g = (-(a*n*h-2*b)/(a*n*h+2*b))*u_i + (2*n*h/(a*n*h+2*b))*g
c     = f_i*u_i + f_g*g
c
c     with
c
c     f_i = -(a*n*h-2*b)/(a*n*h+2*b)
c     f_g = 2*n*h/(a*n*h+2*b)
c
c     For side-centered values, we follow a similar approach.  In this
c     case, however, u_b can be a degree of freedom of the problem, so
c     that
c
c     u_g = u_i + (-a*n*h/b)*u_b + (n*h/b)*g
c     = f_i*u_i + f_b*u_b + f_g*g
c
c     with
c
c     f_i = 1
c     f_b = -a*n*h/b
c     f_g = n*h/b
c
c     For Dirichlet boundary conditions, b=0, and the foregoing
c     expressions are ill defined.  Consequently, in this case, we
c     eliminate u_b and simply set
c
c     u_b = g/a
c     u_g = 2*u_b - u_i
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Set cell centered boundary values using the supplied Robin
c     coefficients along the codimension 1 upper/lower x boundary.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ccrobinphysbdryop1x2d(
     &     U,U_gcw,
     &     acoef,bcoef,gcoef,
     &     location_index,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     blower1,bupper1,
     &     dx,
     &     adjoint_op)
c
      implicit none
c
c     Input.
c
      INTEGER U_gcw

      INTEGER location_index

      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1

      INTEGER blower1,bupper1

      REAL acoef(blower1:bupper1)
      REAL bcoef(blower1:bupper1)
      REAL gcoef(blower1:bupper1)

      REAL dx(0:NDIM-1)

      INTEGER adjoint_op
c
c     Input/Output.
c
      REAL U(CELL2d(ilower,iupper,U_gcw))
c
c     Local variables.
c
      INTEGER i,i_g,i_i
      INTEGER j
      INTEGER sgn
      REAL    a,b,g,h
      REAL    f_g,f_i,n,u_g,u_i
c
c     Initialize temporary variables to yield errors.
c
      u_g = 2.d0**15.d0
      u_i = 2.d0**15.d0
c
c     Set values along the upper/lower x side of the patch.
c
      if ( (location_index .eq. 0) .or.
     &     (location_index .eq. 1) ) then

         h = dx(location_index/NDIM)

         if (location_index .eq. 0) then
            sgn = -1
            i_g = ilower0-1     ! ghost    index
            i_i = ilower0       ! interior index
         else
            sgn = +1
            i_g = iupper0+1     ! ghost    index
            i_i = iupper0       ! interior index
         endif

         do j = blower1,bupper1
            a = acoef(j)
            b = bcoef(j)
            g = gcoef(j)
            do i = 0,U_gcw-1
               n = 1.d0+2.d0*i
               f_i = -(a*n*h-2.d0*b)/(a*n*h+2.d0*b)
               f_g = 2.d0*n*h/(a*n*h+2.d0*b)
               if (adjoint_op .eq. 1) then
                  u_g = U(i_g+sgn*i,j)
                  U(i_i-sgn*i,j) = U(i_i-sgn*i,j) + f_i*u_g
               else
                  u_i = U(i_i-sgn*i,j)
                  U(i_g+sgn*i,j) = f_i*u_i + f_g*g
               endif
            enddo
         enddo

      endif
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Set cell centered boundary values using the supplied Robin
c     coefficients along the codimension 1 upper/lower y boundary.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ccrobinphysbdryop1y2d(
     &     U,U_gcw,
     &     acoef,bcoef,gcoef,
     &     location_index,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     blower0,bupper0,
     &     dx,
     &     adjoint_op)
c
      implicit none
c
c     Input.
c
      INTEGER U_gcw

      INTEGER location_index

      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1

      INTEGER blower0,bupper0

      REAL acoef(blower0:bupper0)
      REAL bcoef(blower0:bupper0)
      REAL gcoef(blower0:bupper0)

      REAL dx(0:NDIM-1)

      INTEGER adjoint_op
c
c     Input/Output.
c
      REAL U(CELL2d(ilower,iupper,U_gcw))
c
c     Local variables.
c
      INTEGER i
      INTEGER j,j_g,j_i
      INTEGER sgn
      REAL    a,b,g,h
      REAL    f_g,f_i,n,u_g,u_i
c
c     Initialize temporary variables to yield errors.
c
      u_g = 2.d0**15.d0
      u_i = 2.d0**15.d0
c
c     Set values along the upper/lower y side of the patch.
c
      if ( (location_index .eq. 2) .or.
     &     (location_index .eq. 3) ) then

         h = dx(location_index/NDIM)

         if (location_index .eq. 2) then
            sgn = -1
            j_g = ilower1-1     ! ghost    index
            j_i = ilower1       ! interior index
         else
            sgn = +1
            j_g = iupper1+1     ! ghost    index
            j_i = iupper1       ! interior index
         endif

         do i = blower0,bupper0
            a = acoef(i)
            b = bcoef(i)
            g = gcoef(i)
            do j = 0,U_gcw-1
               n = 1.d0+2.d0*j
               f_i = -(a*n*h-2.d0*b)/(a*n*h+2.d0*b)
               f_g = 2.d0*n*h/(a*n*h+2.d0*b)
               if (adjoint_op .eq. 1) then
                  u_g = U(i,j_g+sgn*j)
                  U(i,j_i-sgn*j) = U(i,j_i-sgn*j) + f_i*u_g
               else
                  u_i = U(i,j_i-sgn*j)
                  U(i,j_g+sgn*j) = f_i*u_i + f_g*g
               endif
            enddo
         enddo

      endif
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Set cell centered boundary values along the codimension 2 boundary
c     by extrapolating values from the codimension 1 boundary.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ccrobinphysbdryop22d(
     &     U,U_gcw,
     &     location_index,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     blower0,bupper0,
     &     blower1,bupper1,
     &     adjoint_op)
c
      implicit none
c
c     Input.
c
      INTEGER U_gcw

      INTEGER location_index

      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1

      INTEGER blower0,bupper0
      INTEGER blower1,bupper1

      INTEGER adjoint_op
c
c     Input/Output.
c
      REAL U(CELL2d(ilower,iupper,U_gcw))
c
c     Local variables.
c
      REAL    U_g
      INTEGER i,i_bdry,i_mirror
      INTEGER j,j_bdry,j_mirror
      INTEGER sgn_x,sgn_y
c
c     Set the codimension 2 boundary values via linear extrapolation.
c
      if     (location_index .eq. 0) then

         i_bdry = ilower0       ! lower x, lower y
         j_bdry = ilower1

         sgn_x = -1
         sgn_y = -1

      elseif (location_index .eq. 1) then

         i_bdry = iupper0       ! upper x, lower y
         j_bdry = ilower1

         sgn_x = +1
         sgn_y = -1

      elseif (location_index .eq. 2) then

         i_bdry = ilower0       ! lower x, upper y
         j_bdry = iupper1

         sgn_x = -1
         sgn_y = +1

      else

         i_bdry = iupper0       ! upper x, upper y
         j_bdry = iupper1

         sgn_x = +1
         sgn_y = +1

      endif

      do j = blower1,bupper1
         j_mirror = j_bdry+(j_bdry-j+sgn_y)
         do i = blower0,bupper0
            i_mirror = i_bdry+(i_bdry-i+sgn_x)
            if (adjoint_op .eq. 1) then
               U_g = U(i,j)
               U(i_mirror,j_mirror) = U(i_mirror,j_mirror) + U_g
               U(i,j_bdry)          = U(i,j_bdry)          + U_g
               U(i_bdry,j)          = U(i_bdry,j)          + U_g
               U(i_mirror,j_bdry)   = U(i_mirror,j_bdry)   - U_g
               U(i_bdry,j_mirror)   = U(i_bdry,j_mirror)   - U_g
            else
               U(i,j) = U(i_mirror,j_mirror)
     &              + (U(i,j_bdry)-U(i_mirror,j_bdry))
     &              + (U(i_bdry,j)-U(i_bdry,j_mirror))
            endif
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Set side centered boundary values using the supplied Robin
c     coefficients along the codimension 1 upper/lower x boundary.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine scrobinphysbdryop1x2d(
     &     u0,u_gcw,
     &     acoef,bcoef,gcoef,
     &     location_index,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     blower1,bupper1,
     &     dx,
     &     adjoint_op)
c
      implicit none
c
c     Input.
c
      INTEGER u_gcw

      INTEGER location_index

      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1

      INTEGER blower1,bupper1

      REAL acoef(blower1:bupper1)
      REAL bcoef(blower1:bupper1)
      REAL gcoef(blower1:bupper1)

      REAL dx(0:NDIM-1)

      INTEGER adjoint_op
c
c     Input/Output.
c
      REAL u0(SIDE2d0(ilower,iupper,u_gcw))
c
c     Local variables.
c
      INTEGER i,i_b,i_i
      INTEGER j
      INTEGER sgn
      REAL    a,b,g,h,f_b,f_g,f_i,n,u_b,u_g,u_i
c
c     Initialize temporary variables to yield errors.
c
      u_b = 2.d0**15.d0
      u_g = 2.d0**15.d0
      u_i = 2.d0**15.d0
c
c     Set values along the upper/lower x side of the patch.
c
      if ( (location_index .eq. 0) .or.
     &     (location_index .eq. 1) ) then

         h = dx(location_index/NDIM)

         if (location_index .eq. 0) then
            sgn = -1
            i_b = ilower0       ! boundary index
            i_i = ilower0+1     ! interior index
         else
            sgn = +1
            i_b = iupper0+1     ! boundary index
            i_i = iupper0       ! interior index
         endif

         do j = blower1,bupper1
            a = acoef(j)
            b = bcoef(j)
            g = gcoef(j)
            if (abs(b) .lt. 1.d-12) then
c     Dirichlet boundary conditions
               u_b = g/a
               u0(i_b,j) = u_b
               do i = 1,u_gcw
                  f_i = -1.d0
                  f_b = 2.d0
                  if (adjoint_op .eq. 1) then
                     u_g = u0(i_b+sgn*i,j)
                     u0(i_b-sgn*i,j) = u0(i_b-sgn*i,j) + f_i*u_g
                     u0(i_b      ,j) = u0(i_b      ,j) + f_b*u_g
                  else
                     u_i = u0(i_b-sgn*i,j)
                     u0(i_b+sgn*i,j) = f_i*u_i + f_b*u_b
                  endif
               enddo
            else
c     Robin boundary conditions
               u_b = u0(i_b,j)
               do i = 1,u_gcw
                  n = 2.d0*i
                  f_i = 1.d0
                  f_b = -a*n*h/b
                  f_g = n*h/b
                  if (adjoint_op .eq. 1) then
                     u_g = u0(i_b+sgn*i,j)
                     u0(i_b-sgn*i,j) = u0(i_b-sgn*i,j) + f_i*u_g
                     u0(i_b      ,j) = u0(i_b      ,j) + f_b*u_g
                  else
                     u_i = u0(i_b-sgn*i,j)
                     u0(i_b+sgn*i,j) = f_i*u_i + f_b*u_b + f_g*g
                  endif
               enddo
            endif
         enddo

      endif
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Set side centered boundary values using the supplied Robin
c     coefficients along the codimension 1 upper/lower y boundary.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine scrobinphysbdryop1y2d(
     &     u1,u_gcw,
     &     acoef,bcoef,gcoef,
     &     location_index,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     blower0,bupper0,
     &     dx,
     &     adjoint_op)
c
      implicit none
c
c     Input.
c
      INTEGER u_gcw

      INTEGER location_index

      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1

      INTEGER blower0,bupper0

      REAL acoef(blower0:bupper0)
      REAL bcoef(blower0:bupper0)
      REAL gcoef(blower0:bupper0)

      REAL dx(0:NDIM-1)

      INTEGER adjoint_op
c
c     Input/Output.
c
      REAL u1(SIDE2d1(ilower,iupper,u_gcw))
c
c     Local variables.
c
      INTEGER i
      INTEGER j,j_b,j_i
      INTEGER sgn
      REAL    a,b,g,h,f_b,f_g,f_i,n,u_b,u_g,u_i
c
c     Initialize temporary variables to yield errors.
c
      u_b = 2.d0**15.d0
      u_g = 2.d0**15.d0
      u_i = 2.d0**15.d0
c
c     Set values along the upper/lower y side of the patch.
c
      if ( (location_index .eq. 2) .or.
     &     (location_index .eq. 3) ) then

         h = dx(location_index/NDIM)

         if (location_index .eq. 2) then
            sgn = -1
            j_b = ilower1       ! boundary index
            j_i = ilower1+1     ! interior index
         else
            sgn = +1
            j_b = iupper1+1     ! boundary index
            j_i = iupper1       ! interior index
         endif

         do i = blower0,bupper0
            a = acoef(i)
            b = bcoef(i)
            g = gcoef(i)
            if (abs(b) .lt. 1.d-12) then
c     Dirichlet boundary conditions
               u_b = g/a
               u1(i,j_b) = u_b
               do j = 1,u_gcw
                  f_i = -1.d0
                  f_b = 2.d0
                  if (adjoint_op .eq. 1) then
                     u_g = u1(i,j_b+sgn*j)
                     u1(i,j_b-sgn*j) = u1(i,j_b-sgn*j) + f_i*u_g
                     u1(i,j_b      ) = u1(i,j_b      ) + f_b*u_g
                  else
                     u_i = u1(i,j_b-sgn*j)
                     u1(i,j_b+sgn*j) = f_i*u_i + f_b*u_b
                  endif
               enddo
            else
c     Robin boundary conditions
               u_b = u1(i,j_b)
               do j = 1,u_gcw
                  n = 2.d0*j
                  f_i = 1.d0
                  f_b = -a*n*h/b
                  f_g = n*h/b
                  if (adjoint_op .eq. 1) then
                     u_g = u1(i,j_b+sgn*j)
                     u1(i,j_b-sgn*j) = u1(i,j_b-sgn*j) + f_i*u_g
                     u1(i,j_b      ) = u1(i,j_b      ) + f_b*u_g
                  else
                     u_i = u1(i,j_b-sgn*j)
                     u1(i,j_b+sgn*j) = f_i*u_i + f_b*u_b + f_g*g
                  endif
               enddo
            endif
         enddo

      endif
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Set side centered boundary values along the codimension 2 boundary
c     by extrapolating values from the codimension 1 boundary.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine scrobinphysbdryop22d(
     &     u0,u1,u_gcw,
     &     location_index,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     blower0,bupper0,
     &     blower1,bupper1,
     &     adjoint_op)
c
      implicit none
c
c     Input.
c
      INTEGER u_gcw

      INTEGER location_index

      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1

      INTEGER blower0,bupper0
      INTEGER blower1,bupper1

      INTEGER adjoint_op
c
c     Input/Output.
c
      REAL u0(SIDE2d0(ilower,iupper,u_gcw))
      REAL u1(SIDE2d1(ilower,iupper,u_gcw))
c
c     Local variables.
c
      INTEGER i,i_bdry,i_shift
      INTEGER j,j_bdry,j_shift
      REAL u_g,del
c
c     Initialize index variables to yield errors in most cases.
c
      i       = 2**15
      i_bdry  = 2**15
      i_shift = 2**15

      j       = 2**15
      j_bdry  = 2**15
      j_shift = 2**15
c
c     Set the codimension 2 boundary values via linear extrapolation.
c
      if     (location_index .eq. 0) then

         i_bdry = ilower0       ! lower x, lower y
         j_bdry = ilower1
         i_shift = +1
         j_shift = +1

      elseif (location_index .eq. 1) then

         i_bdry = iupper0       ! upper x, lower y
         j_bdry = ilower1
         i_shift = -1
         j_shift = +1

      elseif (location_index .eq. 2) then

         i_bdry = ilower0       ! lower x, upper y
         j_bdry = iupper1
         i_shift = +1
         j_shift = -1

      else

         i_bdry = iupper0       ! upper x, upper y
         j_bdry = iupper1
         i_shift = -1
         j_shift = -1

      endif

      do j = blower1,bupper1
         do i = blower0,bupper0+1
            del = dble(abs(j-j_bdry))
            if (adjoint_op .eq. 1) then
               u_g = u0(i,j)
               u0(i,j_bdry) = u0(i,j_bdry) + (1.d0+del)*u_g
               u0(i,j_bdry+j_shift) = u0(i,j_bdry+j_shift) - del*u_g
            else
               u0(i,j) = (1.d0+del)*u0(i,j_bdry)
     &              - del*u0(i,j_bdry+j_shift)
            endif
         enddo
      enddo

      do j = blower1,bupper1+1
         do i = blower0,bupper0
            del = dble(abs(i-i_bdry))
            if (adjoint_op .eq. 1) then
               u_g = u1(i,j)
               u1(i_bdry,j) = u1(i_bdry,j) + (1.d0+del)*u_g
               u1(i_bdry+i_shift,j) = u1(i_bdry+i_shift,j) - del*u_g
            else
               u1(i,j) = (1.d0+del)*u1(i_bdry,j)
     &              - del*u1(i_bdry+i_shift,j)
            endif
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
