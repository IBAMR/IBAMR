c
c     Routines to set physical boundary condition values.
c
c     Created on 21 May 2007 by Boyce Griffith
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
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     For cell centered values, we follow a similar approach as that
c     implemented in class SAMRAI::solv::CartesianRobinBcHelper.
c     Namely, with u_i denoting the interior cell, u_o denoting the
c     ghost cell, and u_b and u_n denoting the value and normal
c     derivative of u at the boundary,
c
c          u_b = (u_i + u_o)/2   and   u_n = (u_o - u_i)/h
c
c     Now, if
c
c          a*u_b + b*u_n = g
c
c     then
c
c          u_o = u_i*(-(a*h - 2*b)/(a*h + 2*b)) + g*(2*h/(a*h + 2*b))
c
c     The following routines evaluate this formula for use with
c     codimension 1 boundary boxes.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Set cell centered boundary values using the supplied Robin
c     coefficients along the codimension 1 upper/lower x patch boundary.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ccrobinphysbdryop1x3d(
     &     U,U_gcw,
     &     acoef,bcoef,gcoef,
     &     location_index,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
     &     blower1,bupper1,
     &     blower2,bupper2,
     &     dx)
c
      implicit none
c
c     Input.
c
      INTEGER U_gcw

      INTEGER location_index

      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2

      INTEGER blower1,bupper1
      INTEGER blower2,bupper2

      REAL acoef(blower1:bupper1,blower2:bupper2)
      REAL bcoef(blower1:bupper1,blower2:bupper2)
      REAL gcoef(blower1:bupper1,blower2:bupper2)

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL U(CELL3d(ilower,iupper,U_gcw))
c
c     Local variables.
c
      INTEGER i,i_g,i_i
      INTEGER j
      INTEGER k
      INTEGER sgn
      REAL    a,b,g,h,u_i
c
c     Correct the values along the upper/lower x side of the patch.
c
      if ( (location_index .eq. 0) .or.
     &     (location_index .eq. 1) ) then

         if (location_index .eq. 0) then
            sgn = -1
            i_g = ilower0-1     ! ghost    index
            i_i = ilower0       ! interior index
         else
            sgn = +1
            i_g = iupper0+1     ! ghost    index
            i_i = iupper0       ! interior index
         endif

         do k = blower2,bupper2
            do j = blower1,bupper1
               a = acoef(j,k)
               b = bcoef(j,k)
               g = gcoef(j,k)
               do i = 0,U_gcw-1
                  h = (1.d0+2.d0*dble(i))*dx(0)
                  u_i = U(i_i-sgn*i,j,k)
                  U(i_g+sgn*i,j,k) =
     &                 -(a*h*u_i - 2.d0*b*u_i - 2.d0*g*h)/(a*h + 2.d0*b)
               enddo
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
c     coefficients along the codimension 1 upper/lower y patch boundary.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ccrobinphysbdryop1y3d(
     &     U,U_gcw,
     &     acoef,bcoef,gcoef,
     &     location_index,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
     &     blower0,bupper0,
     &     blower2,bupper2,
     &     dx)
c
      implicit none
c
c     Input.
c
      INTEGER U_gcw

      INTEGER location_index

      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2

      INTEGER blower0,bupper0
      INTEGER blower2,bupper2

      REAL acoef(blower0:bupper0,blower2:bupper2)
      REAL bcoef(blower0:bupper0,blower2:bupper2)
      REAL gcoef(blower0:bupper0,blower2:bupper2)

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL U(CELL3d(ilower,iupper,U_gcw))
c
c     Local variables.
c
      INTEGER i
      INTEGER j,j_g,j_i
      INTEGER k
      INTEGER sgn
      REAL    a,b,g,h,u_i
c
c     Correct the values along the upper/lower y side of the patch.
c
      if ( (location_index .eq. 2) .or.
     &     (location_index .eq. 3) ) then

         if (location_index .eq. 2) then
            sgn = -1
            j_g = ilower1-1     ! ghost    index
            j_i = ilower1       ! interior index
         else
            sgn = +1
            j_g = iupper1+1     ! ghost    index
            j_i = iupper1       ! interior index
         endif

         do k = blower2,bupper2
            do i = blower0,bupper0
               a = acoef(i,k)
               b = bcoef(i,k)
               g = gcoef(i,k)
               do j = 0,U_gcw-1
                  h = (1.d0+2.d0*dble(j))*dx(1)
                  u_i = U(i,j_i-sgn*j,k)
                  U(i,j_g+sgn*j,k) =
     &                 -(a*h*u_i - 2.d0*b*u_i - 2.d0*g*h)/(a*h + 2.d0*b)
               enddo
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
c     coefficients along the codimension 1 upper/lower z patch boundary.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ccrobinphysbdryop1z3d(
     &     U,U_gcw,
     &     acoef,bcoef,gcoef,
     &     location_index,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
     &     blower0,bupper0,
     &     blower1,bupper1,
     &     dx)
c
      implicit none
c
c     Input.
c
      INTEGER U_gcw

      INTEGER location_index

      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2

      INTEGER blower0,bupper0
      INTEGER blower1,bupper1

      REAL acoef(blower0:bupper0,blower1:bupper1)
      REAL bcoef(blower0:bupper0,blower1:bupper1)
      REAL gcoef(blower0:bupper0,blower1:bupper1)

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL U(CELL3d(ilower,iupper,U_gcw))
c
c     Local variables.
c
      INTEGER i
      INTEGER j
      INTEGER k,k_g,k_i
      INTEGER sgn
      REAL    a,b,g,h,u_i
c
c     Correct the values along the upper/lower z side of the patch.
c
      if ( (location_index .eq. 4) .or.
     &     (location_index .eq. 5) ) then

         if (location_index .eq. 4) then
            sgn = -1
            k_g = ilower2-1     ! ghost    index
            k_i = ilower2       ! interior index
         else
            sgn = +1
            k_g = iupper2+1     ! ghost    index
            k_i = iupper2       ! interior index
         endif

         do j = blower1,bupper1
            do i = blower0,bupper0
               a = acoef(i,j)
               b = bcoef(i,j)
               g = gcoef(i,j)
               do k = 0,U_gcw-1
                  h = (1.d0+2.d0*dble(k))*dx(2)
                  u_i = U(i,j,k_i-sgn*k)
                  U(i,j,k_g+sgn*k) =
     &                 -(a*h*u_i - 2.d0*b*u_i - 2.d0*g*h)/(a*h + 2.d0*b)
               enddo
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
c     coefficients along the codimension 2 patch boundary.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ccrobinphysbdryop23d(
     &     U,U_gcw,
     &     location_index,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
     &     blower0,bupper0,
     &     blower1,bupper1,
     &     blower2,bupper2)
c
      implicit none
c
c     Input.
c
      INTEGER U_gcw

      INTEGER location_index

      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2

      INTEGER blower0,bupper0
      INTEGER blower1,bupper1
      INTEGER blower2,bupper2
c
c     Input/Output.
c
      REAL U(CELL3d(ilower,iupper,U_gcw))
c
c     Local variables.
c
      INTEGER i,i_bdry,i_mirror
      INTEGER j,j_bdry,j_mirror
      INTEGER k,k_bdry,k_mirror
      INTEGER sgn_x,sgn_y,sgn_z
c
c     Set the codimension 2 boundary values via linear extrapolation.
c

      i_bdry = 0
      j_bdry = 0
      k_bdry = 0

      if     (location_index .eq.  0) then

         j_bdry = ilower1       ! lower y, lower z
         k_bdry = ilower2

         sgn_x =  0
         sgn_y = -1
         sgn_z = -1

      elseif (location_index .eq.  1) then

         j_bdry = iupper1       ! upper y, lower z
         k_bdry = ilower2

         sgn_x =  0
         sgn_y = +1
         sgn_z = -1

      elseif (location_index .eq.  2) then

         j_bdry = ilower1       ! lower y, upper z
         k_bdry = iupper2

         sgn_x =  0
         sgn_y = -1
         sgn_z = +1

      elseif (location_index .eq.  3) then

         j_bdry = iupper1       ! upper y, upper z
         k_bdry = iupper2

         sgn_x =  0
         sgn_y = +1
         sgn_z = +1

      elseif (location_index .eq.  4) then

         i_bdry = ilower0       ! lower x, lower z
         k_bdry = ilower2

         sgn_x = -1
         sgn_y =  0
         sgn_z = -1

      elseif (location_index .eq.  5) then

         i_bdry = ilower0       ! lower x, upper z
         k_bdry = iupper2

         sgn_x = -1
         sgn_y =  0
         sgn_z = +1

      elseif (location_index .eq.  6) then

         i_bdry = iupper0       ! upper x, lower z
         k_bdry = ilower2

         sgn_x = +1
         sgn_y =  0
         sgn_z = -1

      elseif (location_index .eq.  7) then

         i_bdry = iupper0       ! upper x, upper z
         k_bdry = iupper2

         sgn_x = +1
         sgn_y =  0
         sgn_z = +1

      elseif (location_index .eq.  8) then

         i_bdry = ilower0       ! lower x, lower y
         j_bdry = ilower1

         sgn_x = -1
         sgn_y = -1
         sgn_z =  0

      elseif (location_index .eq.  9) then

         i_bdry = iupper0       ! upper x, lower y
         j_bdry = ilower1

         sgn_x = +1
         sgn_y = -1
         sgn_z =  0

      elseif (location_index .eq. 10) then

         i_bdry = ilower0       ! lower x, upper y
         j_bdry = iupper1

         sgn_x = -1
         sgn_y = +1
         sgn_z =  0

      else

         i_bdry = iupper0       ! upper x, upper y
         j_bdry = iupper1

         sgn_x = +1
         sgn_y = +1
         sgn_z =  0

      endif

      if     ((location_index.ge.0).and.(location_index.lt.4)) then
c
c     Set x-edge values.
c
         do k = blower2,bupper2
            k_mirror = k_bdry+(k_bdry-k+sgn_z)
            do j = blower1,bupper1
               j_mirror = j_bdry+(j_bdry-j+sgn_y)
               do i = blower0,bupper0
                  U(i,j,k) = U(i,j_mirror,k_mirror)
     &                 + (U(i,j_bdry,k)-U(i,j_bdry,k_mirror))
     &                 + (U(i,j,k_bdry)-U(i,j_mirror,k_bdry))
               enddo
            enddo
         enddo

      elseif ((location_index.ge.4).and.(location_index.lt.8)) then
c
c     Set y-edge values
c
         do k = blower2,bupper2
            k_mirror = k_bdry+(k_bdry-k+sgn_z)
            do j = blower1,bupper1
               do i = blower0,bupper0
                  i_mirror = i_bdry+(i_bdry-i+sgn_x)
                  U(i,j,k) = U(i_mirror,j,k_mirror)
     &                 + (U(i_bdry,j,k)-U(i_bdry,j,k_mirror))
     &                 + (U(i,j,k_bdry)-U(i_mirror,j,k_bdry))
               enddo
            enddo
         enddo

      elseif ((location_index.ge.8).and.(location_index.lt.12)) then
c
c     Set z-edge values
c
         do k = blower2,bupper2
            do j = blower1,bupper1
               j_mirror = j_bdry+(j_bdry-j+sgn_y)
               do i = blower0,bupper0
                  i_mirror = i_bdry+(i_bdry-i+sgn_x)
                  U(i,j,k) = U(i_mirror,j_mirror,k)
     &                 + (U(i_bdry,j,k)-U(i_bdry,j_mirror,k))
     &                 + (U(i,j_bdry,k)-U(i_mirror,j_bdry,k))
               enddo
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
c     coefficients along the codimension 3 patch boundary.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ccrobinphysbdryop33d(
     &     U,U_gcw,
     &     location_index,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
     &     blower0,bupper0,
     &     blower1,bupper1,
     &     blower2,bupper2)
c
      implicit none
c
c     Input.
c
      INTEGER U_gcw

      INTEGER location_index

      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2

      INTEGER blower0,bupper0
      INTEGER blower1,bupper1
      INTEGER blower2,bupper2
c
c     Input/Output.
c
      REAL U(CELL3d(ilower,iupper,U_gcw))
c
c     Local variables.
c
      INTEGER i,i_bdry,i_mirror
      INTEGER j,j_bdry,j_mirror
      INTEGER k,k_bdry,k_mirror
      INTEGER sgn_x,sgn_y,sgn_z
c
c     Set the codimension 3 boundary values via linear extrapolation.
c

      if     (location_index .eq.  0) then

         i_bdry = ilower0
         j_bdry = ilower1       ! lower x, lower y, lower z
         k_bdry = ilower2

         sgn_x = -1
         sgn_y = -1
         sgn_z = -1

      elseif (location_index .eq.  1) then

         i_bdry = iupper0
         j_bdry = ilower1       ! upper x, lower y, lower z
         k_bdry = ilower2

         sgn_x = +1
         sgn_y = -1
         sgn_z = -1

      elseif (location_index .eq.  2) then

         i_bdry = ilower0
         j_bdry = iupper1       ! lower x, upper y, lower z
         k_bdry = ilower2

         sgn_x = -1
         sgn_y = +1
         sgn_z = -1

      elseif (location_index .eq.  3) then

         i_bdry = iupper0
         j_bdry = iupper1       ! upper x, upper y, lower z
         k_bdry = ilower2

         sgn_x = +1
         sgn_y = +1
         sgn_z = -1

      elseif (location_index .eq.  4) then

         i_bdry = ilower0
         j_bdry = ilower1       ! lower x, lower y, upper z
         k_bdry = iupper2

         sgn_x = -1
         sgn_y = -1
         sgn_z = +1

      elseif (location_index .eq.  5) then

         i_bdry = iupper0
         j_bdry = ilower1       ! upper x, lower y, upper z
         k_bdry = iupper2

         sgn_x = +1
         sgn_y = -1
         sgn_z = +1

      elseif (location_index .eq.  6) then

         i_bdry = ilower0
         j_bdry = iupper1       ! lower x, upper y, upper z
         k_bdry = iupper2

         sgn_x = -1
         sgn_y = +1
         sgn_z = +1

      else

         i_bdry = iupper0
         j_bdry = iupper1       ! upper x, upper y, upper z
         k_bdry = iupper2

         sgn_x = +1
         sgn_y = +1
         sgn_z = +1

      endif

      do k = blower2,bupper2
         k_mirror = k_bdry+(k_bdry-k+sgn_z)
         do j = blower1,bupper1
            j_mirror = j_bdry+(j_bdry-j+sgn_y)
            do i = blower0,bupper0
               i_mirror = i_bdry+(i_bdry-i+sgn_x)
               U(i,j,k) = U(i_mirror,j_mirror,k_mirror)
     &              + (U(i,j_bdry,k_bdry)-U(i_mirror,j_bdry,k_bdry))
     &              + (U(i_bdry,j,k_bdry)-U(i_bdry,j_mirror,k_bdry))
     &              + (U(i_bdry,j_bdry,k)-U(i_bdry,j_bdry,k_mirror))
            enddo
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
      subroutine scrobinphysbdryop1x3d(
     &     u0,u_gcw,
     &     acoef,bcoef,gcoef,
     &     location_index,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
     &     blower1,bupper1,
     &     blower2,bupper2,
     &     dx)
c
      implicit none
c
c     Input.
c
      INTEGER u_gcw

      INTEGER location_index

      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2

      INTEGER blower1,bupper1
      INTEGER blower2,bupper2

      REAL acoef(blower1:bupper1,blower2:bupper2)
      REAL bcoef(blower1:bupper1,blower2:bupper2)
      REAL gcoef(blower1:bupper1,blower2:bupper2)

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL u0(SIDE3d0(ilower,iupper,u_gcw))
c
c     Local variables.
c
      INTEGER i,i_b,i_i
      INTEGER j
      INTEGER k
      INTEGER sgn
      REAL    a,b,g,u_b,u_i
c
c     Set values along the upper/lower x side of the patch.
c
      if ( (location_index .eq. 0) .or.
     &     (location_index .eq. 1) ) then

         if (location_index .eq. 0) then
            sgn = -1
            i_b = ilower0       ! boundary index
            i_i = ilower0+1     ! interior index
         else
            sgn = +1
            i_b = iupper0+1     ! boundary index
            i_i = iupper0       ! interior index
         endif

         do k = blower2,bupper2
            do j = blower1,bupper1
               a = acoef(j,k)
               b = bcoef(j,k)
               g = gcoef(j,k)

               if (abs(b) .lt. 1.d-12) then
c     Dirichlet boundary conditions
                  u_b = g/a
                  u0(i_b,j,k) = u_b
                  do i = 1,u_gcw
                     u_i = u0(i_b-sgn*i,j,k)
                     u0(i_b+sgn*i,j,k) = 2.d0*u_b-u_i
                  enddo
               else
c     Robin boundary conditions
                  u_b = u0(i_b,j,k)
                  do i = 1,u_gcw
                     u_i = u0(i_b-sgn*i,j,k)
                     u0(i_b+sgn*i,j,k) = u_i+
     &                    (2.d0*dble(i)*dx(0)*(g-a*u_b))/b
                  enddo
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
c     Set side centered boundary values using the supplied Robin
c     coefficients along the codimension 1 upper/lower y boundary.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine scrobinphysbdryop1y3d(
     &     u1,u_gcw,
     &     acoef,bcoef,gcoef,
     &     location_index,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
     &     blower0,bupper0,
     &     blower2,bupper2,
     &     dx)
c
      implicit none
c
c     Input.
c
      INTEGER u_gcw

      INTEGER location_index

      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2

      INTEGER blower0,bupper0
      INTEGER blower2,bupper2

      REAL acoef(blower0:bupper0,blower2:bupper2)
      REAL bcoef(blower0:bupper0,blower2:bupper2)
      REAL gcoef(blower0:bupper0,blower2:bupper2)

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL u1(SIDE3d1(ilower,iupper,u_gcw))
c
c     Local variables.
c
      INTEGER i
      INTEGER j,j_b,j_i
      INTEGER k
      INTEGER sgn
      REAL    a,b,g,u_b,u_i
c
c     Set values along the upper/lower y side of the patch.
c
      if ( (location_index .eq. 2) .or.
     &     (location_index .eq. 3) ) then

         if (location_index .eq. 2) then
            sgn = -1
            j_b = ilower1       ! boundary index
            j_i = ilower1+1     ! interior index
         else
            sgn = +1
            j_b = iupper1+1     ! boundary index
            j_i = iupper1       ! interior index
         endif

         do k = blower2,bupper2
            do i = blower0,bupper0
               a = acoef(i,k)
               b = bcoef(i,k)
               g = gcoef(i,k)

               if (abs(b) .lt. 1.d-12) then
c     Dirichlet boundary conditions
                  u_b = g/a
                  u1(i,j_b,k) = u_b
                  do j = 1,u_gcw
                     u_i = u1(i,j_b-sgn*j,k)
                     u1(i,j_b+sgn*j,k) = 2.d0*u_b-u_i
                  enddo
               else
c     Robin boundary conditions
                  u_b = u1(i,j_b,k)
                  do j = 1,u_gcw
                     u_i = u1(i,j_b-sgn*j,k)
                     u1(i,j_b+sgn*j,k) = u_i+
     &                    (2.d0*dble(j)*dx(1)*(g-a*u_b))/b
                  enddo
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
c     Set side centered boundary values using the supplied Robin
c     coefficients along the codimension 1 upper/lower z boundary.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine scrobinphysbdryop1z3d(
     &     u2,u_gcw,
     &     acoef,bcoef,gcoef,
     &     location_index,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
     &     blower0,bupper0,
     &     blower1,bupper1,
     &     dx)
c
      implicit none
c
c     Input.
c
      INTEGER u_gcw

      INTEGER location_index

      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2

      INTEGER blower0,bupper0
      INTEGER blower1,bupper1

      REAL acoef(blower0:bupper0,blower1:bupper1)
      REAL bcoef(blower0:bupper0,blower1:bupper1)
      REAL gcoef(blower0:bupper0,blower1:bupper1)

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL u2(SIDE3d2(ilower,iupper,u_gcw))
c
c     Local variables.
c
      INTEGER i
      INTEGER j
      INTEGER k,k_b,k_i
      INTEGER sgn
      REAL    a,b,g,u_b,u_i
c
c     Set values along the upper/lower z side of the patch.
c
      if ( (location_index .eq. 4) .or.
     &     (location_index .eq. 5) ) then

         if (location_index .eq. 4) then
            sgn = -1
            k_b = ilower2       ! boundary index
            k_i = ilower2+1     ! interior index
         else
            sgn = +1
            k_b = iupper2+1     ! boundary index
            k_i = iupper2       ! interior index
         endif

         do j = blower1,bupper1
            do i = blower0,bupper0
               a = acoef(i,j)
               b = bcoef(i,j)
               g = gcoef(i,j)

               if (abs(b) .lt. 1.d-12) then
c     Dirichlet boundary conditions
                  u_b = g/a
                  u2(i,j,k_b) = u_b
                  do k = 1,u_gcw
                     u_i = u2(i,j,k_b-sgn*k)
                     u2(i,j,k_b+sgn*k) = 2.d0*u_b-u_i
                  enddo
               else
c     Robin boundary conditions
                  u_b = u2(i,j,k_b)
                  do k = 1,u_gcw
                     u_i = u2(i,j,k_b-sgn*k)
                     u2(i,j,k_b+sgn*k) = u_i+
     &                    (2.d0*dble(k)*dx(2)*(g-a*u_b))/b
                  enddo
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
c     Set side centered boundary values along the codimension 2 boundary
c     by extrapolating values from the codimension 1 boundary.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine scrobinphysbdryop23d(
     &     u0,u1,u2,u_gcw,
     &     location_index,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
     &     blower0,bupper0,
     &     blower1,bupper1,
     &     blower2,bupper2)
c
      implicit none
c
c     Input.
c
      INTEGER u_gcw

      INTEGER location_index

      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2

      INTEGER blower0,bupper0
      INTEGER blower1,bupper1
      INTEGER blower2,bupper2
c
c     Input/Output.
c
      REAL u0(SIDE3d0(ilower,iupper,u_gcw))
      REAL u1(SIDE3d1(ilower,iupper,u_gcw))
      REAL u2(SIDE3d2(ilower,iupper,u_gcw))
c
c     Local variables.
c
      INTEGER i,i_bdry,i_shift
      INTEGER j,j_bdry,j_shift
      INTEGER k,k_bdry,k_shift
      REAL u_b,du
c
c     Initialize index variables to yield errors in most cases.
c
      i       = 2**15
      i_bdry  = 2**15
      i_shift = 2**15

      j       = 2**15
      j_bdry  = 2**15
      j_shift = 2**15

      k       = 2**15
      k_bdry  = 2**15
      k_shift = 2**15
c
c     Set the codimension 2 boundary values via linear extrapolation.
c
      if     (location_index .eq.  0) then

         j_bdry = ilower1       ! lower y, lower z
         k_bdry = ilower2

         j_shift = +1
         k_shift = +1

      elseif (location_index .eq.  1) then

         j_bdry = iupper1       ! upper y, lower z
         k_bdry = ilower2

         j_shift = -1
         k_shift = +1

      elseif (location_index .eq.  2) then

         j_bdry = ilower1       ! lower y, upper z
         k_bdry = iupper2

         j_shift = +1
         k_shift = -1

      elseif (location_index .eq.  3) then

         j_bdry = iupper1       ! upper y, upper z
         k_bdry = iupper2

         j_shift = -1
         k_shift = -1

      elseif (location_index .eq.  4) then

         i_bdry = ilower0       ! lower x, lower z
         k_bdry = ilower2

         i_shift = +1
         k_shift = +1

      elseif (location_index .eq.  5) then

         i_bdry = ilower0       ! lower x, upper z
         k_bdry = iupper2

         i_shift = +1
         k_shift = -1

      elseif (location_index .eq.  6) then

         i_bdry = iupper0       ! upper x, lower z
         k_bdry = ilower2

         i_shift = -1
         k_shift = +1

      elseif (location_index .eq.  7) then

         i_bdry = iupper0       ! upper x, upper z
         k_bdry = iupper2

         i_shift = -1
         k_shift = -1

      elseif (location_index .eq.  8) then

         i_bdry = ilower0       ! lower x, lower y
         j_bdry = ilower1

         i_shift = +1
         j_shift = +1

      elseif (location_index .eq.  9) then

         i_bdry = iupper0       ! upper x, lower y
         j_bdry = ilower1

         i_shift = -1
         j_shift = +1

      elseif (location_index .eq. 10) then

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

      if     ((location_index.ge.0).and.(location_index.lt.4)) then
c
c     Set x-edge values.
c
         do k = blower2,bupper2
            do j = blower1,bupper1+1
               do i = blower0,bupper0
                  u_b = u1(i,j,k_bdry)
                  du = u1(i,j,k_bdry)-u1(i,j,k_bdry+k_shift)
                  u1(i,j,k) = u_b + dble(abs(k-k_bdry))*du
               enddo
            enddo
         enddo

         do k = blower2,bupper2+1
            do j = blower1,bupper1
               do i = blower0,bupper0
                  u_b = u2(i,j_bdry,k)
                  du = u2(i,j_bdry,k)-u2(i,j_bdry+j_shift,k)
                  u2(i,j,k) = u_b + dble(abs(j-j_bdry))*du
               enddo
            enddo
         enddo

      elseif ((location_index.ge.4).and.(location_index.lt.8)) then
c
c     Set y-edge values
c
         do k = blower2,bupper2
            do j = blower1,bupper1
               do i = blower0,bupper0+1
                  u_b = u0(i,j,k_bdry)
                  du = u0(i,j,k_bdry)-u0(i,j,k_bdry+k_shift)
                  u0(i,j,k) = u_b + dble(abs(k-k_bdry))*du
               enddo
            enddo
         enddo

         do k = blower2,bupper2+1
            do j = blower1,bupper1
               do i = blower0,bupper0
                  u_b = u2(i_bdry,j,k)
                  du = u2(i_bdry,j,k)-u2(i_bdry+i_shift,j,k)
                  u2(i,j,k) = u_b + dble(abs(i-i_bdry))*du
               enddo
            enddo
         enddo

      elseif ((location_index.ge.8).and.(location_index.lt.12)) then
c
c     Set z-edge values
c
         do k = blower2,bupper2
            do j = blower1,bupper1
               do i = blower0,bupper0+1
                  u_b = u0(i,j_bdry,k)
                  du = u0(i,j_bdry,k)-u0(i,j_bdry+j_shift,k)
                  u0(i,j,k) = u_b + dble(abs(j-j_bdry))*du
               enddo
            enddo
         enddo

         do k = blower2,bupper2
            do j = blower1,bupper1+1
               do i = blower0,bupper0
                  u_b = u1(i_bdry,j,k)
                  du = u1(i_bdry,j,k)-u1(i_bdry+i_shift,j,k)
                  u1(i,j,k) = u_b + dble(abs(i-i_bdry))*du
               enddo
            enddo
         enddo

      endif
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Set side centered boundary values along the codimension 3 boundary
c     by extrapolating values from the codimension 1 and 2 boundary.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine scrobinphysbdryop33d(
     &     u0,u1,u2,u_gcw,
     &     location_index,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
     &     blower0,bupper0,
     &     blower1,bupper1,
     &     blower2,bupper2)
c
      implicit none
c
c     Input.
c
      INTEGER u_gcw

      INTEGER location_index

      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2

      INTEGER blower0,bupper0
      INTEGER blower1,bupper1
      INTEGER blower2,bupper2
c
c     Input/Output.
c
      REAL u0(SIDE3d0(ilower,iupper,u_gcw))
      REAL u1(SIDE3d1(ilower,iupper,u_gcw))
      REAL u2(SIDE3d2(ilower,iupper,u_gcw))
c
c     Local variables.
c
      INTEGER i,i_bdry,i_shift
      INTEGER j,j_bdry,j_shift
      INTEGER k,k_bdry,k_shift
      REAL u_b,du_i,du_j,du_k
c
c     Initialize index variables to yield errors in most cases.
c
      i       = 2**15
      i_bdry  = 2**15
      i_shift = 2**15

      j       = 2**15
      j_bdry  = 2**15
      j_shift = 2**15

      k       = 2**15
      k_bdry  = 2**15
      k_shift = 2**15

      du_i = 2.d0**15.d0
      du_j = 2.d0**15.d0
      du_k = 2.d0**15.d0
c
c     Set the codimension 3 boundary values via linear extrapolation.
c
      if     (location_index .eq.  0) then

         i_bdry = ilower0
         j_bdry = ilower1       ! lower x, lower y, lower z
         k_bdry = ilower2

         i_shift = +1
         j_shift = +1
         k_shift = +1

      elseif (location_index .eq.  1) then

         i_bdry = iupper0
         j_bdry = ilower1       ! upper x, lower y, lower z
         k_bdry = ilower2

         i_shift = -1
         j_shift = +1
         k_shift = +1

      elseif (location_index .eq.  2) then

         i_bdry = ilower0
         j_bdry = iupper1       ! lower x, upper y, lower z
         k_bdry = ilower2

         i_shift = +1
         j_shift = -1
         k_shift = +1

      elseif (location_index .eq.  3) then

         i_bdry = iupper0
         j_bdry = iupper1       ! upper x, upper y, lower z
         k_bdry = ilower2

         i_shift = -1
         j_shift = -1
         k_shift = +1

      elseif (location_index .eq.  4) then

         i_bdry = ilower0
         j_bdry = ilower1       ! lower x, lower y, upper z
         k_bdry = iupper2

         i_shift = +1
         j_shift = +1
         k_shift = -1

      elseif (location_index .eq.  5) then

         i_bdry = iupper0
         j_bdry = ilower1       ! upper x, lower y, upper z
         k_bdry = iupper2

         i_shift = -1
         j_shift = +1
         k_shift = -1

      elseif (location_index .eq.  6) then

         i_bdry = ilower0
         j_bdry = iupper1       ! lower x, upper y, upper z
         k_bdry = iupper2

         i_shift = +1
         j_shift = -1
         k_shift = -1

      else

         i_bdry = iupper0
         j_bdry = iupper1       ! upper x, upper y, upper z
         k_bdry = iupper2

         i_shift = -1
         j_shift = -1
         k_shift = -1

      endif

      do k = blower2,bupper2
         do j = blower1,bupper1
            do i = blower0,bupper0+1
               u_b = u0(i,j_bdry,k_bdry)
               du_j = u0(i,j_bdry,k_bdry)-u0(i,j_bdry+j_shift,k_bdry)
               du_k = u0(i,j_bdry,k_bdry)-u0(i,j_bdry,k_bdry+k_shift)
               u0(i,j,k) = u_b
     &              + dble(abs(j-j_bdry))*du_j
     &              + dble(abs(k-k_bdry))*du_k
            enddo
         enddo
      enddo

      do k = blower2,bupper2
         do j = blower1,bupper1+1
            do i = blower0,bupper0
               u_b = u1(i_bdry,j,k_bdry)
               du_i = u1(i_bdry,j,k_bdry)-u1(i_bdry+i_shift,j,k_bdry)
               du_k = u1(i_bdry,j,k_bdry)-u1(i_bdry,j,k_bdry+k_shift)
               u1(i,j,k) = u_b
     &              + dble(abs(i-i_bdry))*du_i
     &              + dble(abs(k-k_bdry))*du_k
            enddo
         enddo
      enddo

      do k = blower2,bupper2+1
         do j = blower1,bupper1
            do i = blower0,bupper0
               u_b = u2(i_bdry,j_bdry,k)
               du_j = u2(i_bdry,j_bdry,k)-u2(i_bdry+i_shift,j_bdry,k)
               du_k = u2(i_bdry,j_bdry,k)-u2(i_bdry,j_bdry+j_shift,k)
               u2(i,j,k) = u_b
     &              + dble(abs(i-i_bdry))*du_i
     &              + dble(abs(j-j_bdry))*du_j
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
