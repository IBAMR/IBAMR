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
define(NDIM,3)dnl
define(REAL,`double precision')dnl
define(INTEGER,`integer')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim3d.i)dnl
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Set inhomogeneous Robin boundary condition coefficients to enforce
c     incompressibility or no-stress at physical boundaries.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine navier_stokes_open_bc_coefs3d(
     &     U,U_gcw,
     &     acoef,bcoef,gcoef,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
     &     blower0,bupper0,
     &     blower1,bupper1,
     &     blower2,bupper2,
     &     location_index,
     &     comp_idx,
     &     dx)
c
      implicit none
c
c     Input.
c
      INTEGER U_gcw

      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2

      INTEGER blower0,bupper0
      INTEGER blower1,bupper1
      INTEGER blower2,bupper2

      INTEGER location_index,comp_idx

      REAL U(CELL3d(ilower,iupper,U_gcw),0:NDIM-1)

      REAL acoef(blower0:bupper0,
     &           blower1:bupper1,
     &           blower2:bupper2)
      REAL bcoef(blower0:bupper0,
     &           blower1:bupper1,
     &           blower2:bupper2)

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL gcoef(blower0:bupper0,
     &           blower1:bupper1,
     &           blower2:bupper2)
c
c     Local variables.
c
      INTEGER i,i_intr,i_bdry,i_s
      INTEGER j,j_intr,j_bdry,j_s
      INTEGER k,k_intr,k_bdry,k_s
      INTEGER bdry_normal_axis
      REAL F,sgn
c
c     Initialize index variables to yield errors in most cases.
c
      i      = 2**15
      i_intr = 2**15
      i_bdry = 2**15
      i_s    = 2**15

      j      = 2**15
      j_intr = 2**15
      j_bdry = 2**15
      j_s    = 2**15

      k      = 2**15
      k_intr = 2**15
      k_bdry = 2**15
      k_s    = 2**15
c
c     Set the boundary condition coefficients.
c
c     At "open" boundaries, modify the normal velocity boundary
c     conditions to enforce div u = 0, and modify the tangential
c     velocity boundary conditions to enforce zero stress.  This is done
c     by specifying a normal flux F at the boundary.
c

      bdry_normal_axis = location_index/2

      if ( (location_index .eq. 0) .or.
     &     (location_index .eq. 1) ) then

         if (location_index .eq. 0) then
            sgn = +1.d0
            i_intr = ilower0
            i_bdry = ilower0-1
            i_s    = ilower0
         else
            sgn = -1.d0
            i_intr = iupper0
            i_bdry = iupper0+1
            i_s    = iupper0+1
         endif

         if (comp_idx .eq. 0) then
c
c     Set F to enforce div u = 0.
c
            do k = ilower2,iupper2
               do j = ilower1,iupper1
                  if ( abs(acoef(i_s,j,k)       ) .lt. 1.0d-12 .or.
     &                 abs(bcoef(i_s,j,k) - 1.d0) .lt. 1.0d-12 ) then
                     F =  (0.25d0/dx(1))*(
     &                    + U(i_intr,j+1,k,1)
     &                    - U(i_intr,j-1,k,1)
     &                    + U(i_bdry,j+1,k,1)
     &                    - U(i_bdry,j-1,k,1)) +
     &                    (0.25d0/dx(2))*(
     &                    + U(i_intr,j,k+1,2)
     &                    - U(i_intr,j,k-1,2)
     &                    + U(i_bdry,j,k+1,2)
     &                    - U(i_bdry,j,k-1,2))
                     gcoef(i_s,j,k) = sgn*F
                  endif
               enddo
            enddo

         elseif (comp_idx .eq. 1) then
c
c     Set F to enforce t * sigma * n = 0.
c
            do k = ilower2,iupper2
               do j = ilower1,iupper1
                  if ( abs(acoef(i_s,j,k)       ) .lt. 1.0d-12 .or.
     &                 abs(bcoef(i_s,j,k) - 1.d0) .lt. 1.0d-12 ) then
                     F = (0.25d0/dx(1))*(
     &                    + U(i_intr,j+1,k,bdry_normal_axis)
     &                    - U(i_intr,j-1,k,bdry_normal_axis)
     &                    + U(i_bdry,j+1,k,bdry_normal_axis)
     &                    - U(i_bdry,j-1,k,bdry_normal_axis))
                     gcoef(i_s,j,k) = sgn*F
                  endif
               enddo
            enddo

         elseif (comp_idx .eq. 2) then
c
c     Set F to enforce t * sigma * n = 0.
c
            do k = ilower2,iupper2
               do j = ilower1,iupper1
                  if ( abs(acoef(i_s,j,k)       ) .lt. 1.0d-12 .or.
     &                 abs(bcoef(i_s,j,k) - 1.d0) .lt. 1.0d-12 ) then
                     F = (0.25d0/dx(2))*(
     &                    + U(i_intr,j,k+1,bdry_normal_axis)
     &                    - U(i_intr,j,k-1,bdry_normal_axis)
     &                    + U(i_bdry,j,k+1,bdry_normal_axis)
     &                    - U(i_bdry,j,k-1,bdry_normal_axis))
                     gcoef(i_s,j,k) = sgn*F
                  endif
               enddo
            enddo

         endif

      elseif ( (location_index .eq. 2) .or.
     &         (location_index .eq. 3) ) then

         if (location_index .eq. 2) then
            sgn = +1.d0
            j_intr = ilower1
            j_bdry = ilower1-1
            j_s    = ilower1
         else
            sgn = -1.d0
            j_intr = iupper1
            j_bdry = iupper1+1
            j_s    = iupper1+1
         endif

         if (comp_idx .eq. 0) then
c
c     Set F to enforce t * sigma * n = 0.
c
            do k = ilower2,iupper2
               do i = ilower0,iupper0
                  if ( abs(acoef(i,j_s,k)       ) .lt. 1.0d-12 .or.
     &                 abs(bcoef(i,j_s,k) - 1.d0) .lt. 1.0d-12 ) then
                     F = (0.25d0/dx(0))*(
     &                    + U(i+1,j_intr,k,bdry_normal_axis)
     &                    - U(i-1,j_intr,k,bdry_normal_axis)
     &                    + U(i+1,j_bdry,k,bdry_normal_axis)
     &                    - U(i-1,j_bdry,k,bdry_normal_axis))
                     gcoef(i,j_s,k) = sgn*F
                  endif
               enddo
            enddo

         elseif (comp_idx .eq. 1) then
c
c     Set F to enforce div u = 0.
c
            do k = ilower2,iupper2
               do i = ilower0,iupper0
                  if ( abs(acoef(i,j_s,k)       ) .lt. 1.0d-12 .or.
     &                 abs(bcoef(i,j_s,k) - 1.d0) .lt. 1.0d-12 ) then
                     F =  (0.25d0/dx(0))*(
     &                    + U(i+1,j_intr,k,0)
     &                    - U(i-1,j_intr,k,0)
     &                    + U(i+1,j_bdry,k,0)
     &                    - U(i-1,j_bdry,k,0)) +
     &                    (0.25d0/dx(2))*(
     &                    + U(i,j_intr,k+1,2)
     &                    - U(i,j_intr,k-1,2)
     &                    + U(i,j_bdry,k+1,2)
     &                    - U(i,j_bdry,k-1,2))
                     gcoef(i,j_s,k) = sgn*F
                  endif
               enddo
            enddo

         elseif (comp_idx .eq. 2) then
c
c     Set F to enforce t * sigma * n = 0.
c
            do k = ilower2,iupper2
               do i = ilower0,iupper0
                  if ( abs(acoef(i,j_s,k)       ) .lt. 1.0d-12 .or.
     &                 abs(bcoef(i,j_s,k) - 1.d0) .lt. 1.0d-12 ) then
                     F = (0.25d0/dx(2))*(
     &                    + U(i,j_intr,k+1,bdry_normal_axis)
     &                    - U(i,j_intr,k-1,bdry_normal_axis)
     &                    + U(i,j_bdry,k+1,bdry_normal_axis)
     &                    - U(i,j_bdry,k-1,bdry_normal_axis))
                     gcoef(i,j_s,k) = sgn*F
                  endif
               enddo
            enddo

         endif

      elseif ( (location_index .eq. 4) .or.
     &         (location_index .eq. 5) ) then

         if (location_index .eq. 4) then
            sgn = +1.d0
            k_intr = ilower2
            k_bdry = ilower2-1
            k_s    = ilower2
         else
            sgn = -1.d0
            k_intr = iupper2
            k_bdry = iupper2+1
            k_s    = iupper2+1
         endif

         if (comp_idx .eq. 0) then
c
c     Set F to enforce t * sigma * n = 0.
c
            do j = ilower1,iupper1
               do i = ilower0,iupper0
                  if ( abs(acoef(i,j,k_s)       ) .lt. 1.0d-12 .or.
     &                 abs(bcoef(i,j,k_s) - 1.d0) .lt. 1.0d-12 ) then
                     F = (0.25d0/dx(0))*(
     &                    + U(i+1,j,k_intr,bdry_normal_axis)
     &                    - U(i-1,j,k_intr,bdry_normal_axis)
     &                    + U(i+1,j,k_bdry,bdry_normal_axis)
     &                    - U(i-1,j,k_bdry,bdry_normal_axis))
                     gcoef(i,j,k_s) = sgn*F
                  endif
               enddo
            enddo

         elseif (comp_idx .eq. 1) then
c
c     Set F to enforce t * sigma * n = 0.
c
            do j = ilower1,iupper1
               do i = ilower0,iupper0
                  if ( abs(acoef(i,j,k_s)       ) .lt. 1.0d-12 .or.
     &                 abs(bcoef(i,j,k_s) - 1.d0) .lt. 1.0d-12 ) then
                     F = (0.25d0/dx(1))*(
     &                    + U(i,j+1,k_intr,bdry_normal_axis)
     &                    - U(i,j-1,k_intr,bdry_normal_axis)
     &                    + U(i,j+1,k_bdry,bdry_normal_axis)
     &                    - U(i,j-1,k_bdry,bdry_normal_axis))
                     gcoef(i,j,k_s) = sgn*F
                  endif
               enddo
            enddo

         elseif (comp_idx .eq. 2) then
c
c     Set F to enforce div u = 0.
c
            do j = ilower1,iupper1
               do i = ilower0,iupper0
                  if ( abs(acoef(i,j,k_s)       ) .lt. 1.0d-12 .or.
     &                 abs(bcoef(i,j,k_s) - 1.d0) .lt. 1.0d-12 ) then
                     F =  (0.25d0/dx(0))*(
     &                    + U(i+1,j,k_intr,0)
     &                    - U(i-1,j,k_intr,0)
     &                    + U(i+1,j,k_bdry,0)
     &                    - U(i-1,j,k_bdry,0)) +
     &                    (0.25d0/dx(1))*(
     &                    + U(i,j+1,k_intr,1)
     &                    - U(i,j-1,k_intr,1)
     &                    + U(i,j+1,k_bdry,1)
     &                    - U(i,j-1,k_bdry,1))
                     gcoef(i,j,k_s) = sgn*F
                  endif
               enddo
            enddo

         endif

      endif
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Set inhomogeneous Robin boundary condition coefficients to enforce
c     tangential boundary conditions at physical boundaries using lagged
c     values of phi.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine navier_stokes_tangential_bc_coefs3d(
     &     Phi,Phi_gcw,
     &     acoef,bcoef,gcoef,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
     &     blower0,bupper0,
     &     blower1,bupper1,
     &     blower2,bupper2,
     &     location_index,
     &     comp_idx,
     &     rho,dx,dt)
c
      implicit none
c
c     Input.
c
      INTEGER Phi_gcw

      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2

      INTEGER blower0,bupper0
      INTEGER blower1,bupper1
      INTEGER blower2,bupper2

      INTEGER location_index,comp_idx

      REAL Phi(CELL3d(ilower,iupper,Phi_gcw))

      REAL acoef(blower0:bupper0,
     &           blower1:bupper1,
     &           blower2:bupper2)
      REAL bcoef(blower0:bupper0,
     &           blower1:bupper1,
     &           blower2:bupper2)

      REAL rho
      REAL dx(0:NDIM-1)
      REAL dt
c
c     Input/Output.
c
      REAL gcoef(blower0:bupper0,
     &           blower1:bupper1,
     &           blower2:bupper2)
c
c     Local variables.
c
      INTEGER i,i_intr,i_bdry,i_s
      INTEGER j,j_intr,j_bdry,j_s
      INTEGER k,k_intr,k_bdry,k_s
      REAL grad_Phi,t_grad_grad_Phi_n,sgn
c
c     Initialize index variables to yield errors in most cases.
c
      i      = 2**15
      i_intr = 2**15
      i_bdry = 2**15
      i_s    = 2**15

      j      = 2**15
      j_intr = 2**15
      j_bdry = 2**15
      j_s    = 2**15

      k      = 2**15
      k_intr = 2**15
      k_bdry = 2**15
      k_s    = 2**15
c
c     Set the boundary condition coefficients.
c
      if ( (location_index .eq. 0) .or.
     &     (location_index .eq. 1) ) then

         if (location_index .eq. 0) then
            sgn = -1.d0
            i_intr = ilower0
            i_bdry = ilower0-1
            i_s    = ilower0
         else
            sgn = +1.d0
            i_intr = iupper0
            i_bdry = iupper0+1
            i_s    = iupper0+1
         endif

         do k = ilower2,iupper2
            do j = ilower1,iupper1
               if ( abs(acoef(i_s,j,k) - 1.d0) .lt. 1.0d-12 .or.
     &              abs(bcoef(i_s,j,k)       ) .lt. 1.0d-12 ) then

                  if (comp_idx .eq. 1) then
                     grad_Phi = (0.5d0/dx(1))*(
     &                    Phi(i_intr,j+1,k) - Phi(i_intr,j-1,k))
                     gcoef(i_s,j,k) = gcoef(i_s,j,k) + (dt/rho)*grad_Phi
                  elseif (comp_idx .eq. 2) then
                     grad_Phi = (0.5d0/dx(2))*(
     &                    Phi(i_intr,j,k+1) - Phi(i_intr,j,k-1))
                     gcoef(i_s,j,k) = gcoef(i_s,j,k) + (dt/rho)*grad_Phi
                  endif

               else

                  if (comp_idx .eq. 1) then
                     t_grad_grad_Phi_n = (0.5d0*sgn/dx(1))*(
     &                    +(Phi(i_bdry,j+1,k)-Phi(i_intr,j+1,k))/dx(0)
     &                    +(Phi(i_bdry,j-1,k)-Phi(i_intr,j-1,k))/dx(0))
                     gcoef(i_s,j,k) = gcoef(i_s,j,k) +
     &                    (dt/rho)*t_grad_grad_Phi_n
                  elseif (comp_idx .eq. 2) then
                     t_grad_grad_Phi_n = (0.5d0*sgn/dx(2))*(
     &                    +(Phi(i_bdry,j,k+1)-Phi(i_intr,j,k+1))/dx(0)
     &                    +(Phi(i_bdry,j,k-1)-Phi(i_intr,j,k-1))/dx(0))
                     gcoef(i_s,j,k) = gcoef(i_s,j,k) +
     &                    2.d0*(dt/rho)*t_grad_grad_Phi_n
                  endif

               endif
            enddo
         enddo

      elseif ( (location_index .eq. 2) .or.
     &         (location_index .eq. 3) ) then

         if (location_index .eq. 2) then
            sgn = -1.d0
            j_intr = ilower1
            j_bdry = ilower1-1
            j_s    = ilower1
         else
            sgn = +1.d0
            j_intr = iupper1
            j_bdry = iupper1+1
            j_s    = iupper1+1
         endif

         do k = ilower2,iupper2
            do i = ilower0,iupper0
               if ( abs(acoef(i,j_s,k) - 1.d0) .lt. 1.0d-12 .or.
     &              abs(bcoef(i,j_s,k)       ) .lt. 1.0d-12 ) then

                  if (comp_idx .eq. 0) then
                     grad_Phi = (0.5d0/dx(0))*(
     &                    Phi(i+1,j_intr,k) - Phi(i-1,j_intr,k))
                     gcoef(i,j_s,k) = gcoef(i,j_s,k) + (dt/rho)*grad_Phi
                  elseif (comp_idx .eq. 2) then
                     grad_Phi = (0.5d0/dx(2))*(
     &                    Phi(i,j_intr,k+1) - Phi(i,j_intr,k-1))
                     gcoef(i,j_s,k) = gcoef(i,j_s,k) + (dt/rho)*grad_Phi
                  endif

               else

                  if (comp_idx .eq. 0) then
                     t_grad_grad_Phi_n = (0.5d0*sgn/dx(0))*(
     &                    +(Phi(i+1,j_bdry,k)-Phi(i+1,j_intr,k))/dx(1)
     &                    +(Phi(i-1,j_bdry,k)-Phi(i-1,j_intr,k))/dx(1))
                     gcoef(i,j_s,k) = gcoef(i,j_s,k) +
     &                    (dt/rho)*t_grad_grad_Phi_n
                  elseif (comp_idx .eq. 2) then
                     t_grad_grad_Phi_n = (0.5d0*sgn/dx(2))*(
     &                    +(Phi(i,j_bdry,k+1)-Phi(i,j_intr,k+1))/dx(1)
     &                    +(Phi(i,j_bdry,k-1)-Phi(i,j_intr,k-1))/dx(1))
                     gcoef(i,j_s,k) = gcoef(i,j_s,k) +
     &                    2.d0*(dt/rho)*t_grad_grad_Phi_n
                  endif

               endif
            enddo
         enddo

      elseif ( (location_index .eq. 4) .or.
     &         (location_index .eq. 5) ) then

         if (location_index .eq. 4) then
            sgn = -1.d0
            k_intr = ilower2
            k_bdry = ilower2-1
            k_s    = ilower2
         else
            sgn = +1.d0
            k_intr = iupper2
            k_bdry = iupper2+1
            k_s    = iupper2+1
         endif

         do j = ilower1,iupper1
            do i = ilower0,iupper0
               if ( abs(acoef(i,j,k_s) - 1.d0) .lt. 1.0d-12 .or.
     &              abs(bcoef(i,j,k_s)       ) .lt. 1.0d-12 ) then

                  if (comp_idx .eq. 0) then
                     grad_Phi = (0.5d0/dx(0))*(
     &                    Phi(i+1,j,k_intr) - Phi(i-1,j,k_intr))
                     gcoef(i,j,k_s) = gcoef(i,j,k_s) + (dt/rho)*grad_Phi
                  elseif (comp_idx .eq. 1) then
                     grad_Phi = (0.5d0/dx(2))*(
     &                    Phi(i,j+1,k_intr) - Phi(i,j-1,k_intr))
                     gcoef(i,j,k_s) = gcoef(i,j,k_s) + (dt/rho)*grad_Phi
                  endif

               else

                  if (comp_idx .eq. 0) then
                     t_grad_grad_Phi_n = (0.5d0*sgn/dx(0))*(
     &                    +(Phi(i+1,j,k_bdry)-Phi(i+1,j,k_intr))/dx(2)
     &                    +(Phi(i-1,j,k_bdry)-Phi(i-1,j,k_intr))/dx(2))
                     gcoef(i,j,k_s) = gcoef(i,j,k_s) +
     &                    (dt/rho)*t_grad_grad_Phi_n
                  elseif (comp_idx .eq. 1) then
                     t_grad_grad_Phi_n = (0.5d0*sgn/dx(1))*(
     &                    +(Phi(i,j+1,k_bdry)-Phi(i,j+1,k_intr))/dx(2)
     &                    +(Phi(i,j-1,k_bdry)-Phi(i,j-1,k_intr))/dx(2))
                     gcoef(i,j,k_s) = gcoef(i,j,k_s) +
     &                    2.d0*(dt/rho)*t_grad_grad_Phi_n
                  endif

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
c     Set the homogeneous Robin boundary condition coefficients for the
c     projection Poisson equation.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine navier_stokes_homogeneous_projection_bc_coefs3d(
     &     acoef,bcoef,
     &     blower0,bupper0,
     &     blower1,bupper1,
     &     blower2,bupper2)
c
      implicit none
c
c     Input.
c
      INTEGER blower0,bupper0
      INTEGER blower1,bupper1
      INTEGER blower2,bupper2
c
c     Input/Output.
c
      REAL acoef(blower0:bupper0,
     &           blower1:bupper1,
     &           blower2:bupper2)
      REAL bcoef(blower0:bupper0,
     &           blower1:bupper1,
     &           blower2:bupper2)
c
c     Local variables.
c
      INTEGER i,j,k
c
c     "Flip" the values of the homogeneous Robin bc coefs.
c
      do k = blower2,bupper2
         do j = blower1,bupper1
            do i = blower0,bupper0
               acoef(i,j,k) = 1.d0 - acoef(i,j,k)
            enddo
         enddo
      enddo

      do k = blower2,bupper2
         do j = blower1,bupper1
            do i = blower0,bupper0
               bcoef(i,j,k) = 1.d0 - bcoef(i,j,k)
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Set the inhomogeneous Robin boundary condition coefficients for
c     the projection Poisson equation.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine navier_stokes_fc_inhomogeneous_projection_bc_coefs3d(
     &     u0,u1,u2,u_gcw,
     &     P,P_gcw,
     &     acoef,bcoef,gcoef,P_bdry,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
     &     blower0,bupper0,
     &     blower1,bupper1,
     &     blower2,bupper2,
     &     location_index,
     &     using_pressure_increment,
     &     rho,dt)
c
      implicit none
c
c     Input.
c
      INTEGER u_gcw,P_gcw

      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2

      INTEGER blower0,bupper0
      INTEGER blower1,bupper1
      INTEGER blower2,bupper2

      INTEGER location_index,using_pressure_increment

      REAL u0(FACE3d0(ilower,iupper,u_gcw))
      REAL u1(FACE3d1(ilower,iupper,u_gcw))
      REAL u2(FACE3d2(ilower,iupper,u_gcw))

      REAL P(CELL3d(ilower,iupper,P_gcw))

      REAL acoef(blower0:bupper0,
     &           blower1:bupper1,
     &           blower2:bupper2)
      REAL bcoef(blower0:bupper0,
     &           blower1:bupper1,
     &           blower2:bupper2)
      REAL P_bdry(blower0:bupper0,
     &            blower1:bupper1,
     &            blower2:bupper2)

      REAL rho,dt
c
c     Input/Output.
c
      REAL gcoef(blower0:bupper0,
     &           blower1:bupper1,
     &           blower2:bupper2)
c
c     Local variables.
c
      INTEGER i,i_intr,i_bdry,i_s
      INTEGER j,j_intr,j_bdry,j_s
      INTEGER k,k_intr,k_bdry,k_s
      REAL sgn
c
c     Initialize index variables to yield errors in most cases.
c
      i      = 2**15
      i_intr = 2**15
      i_bdry = 2**15
      i_s    = 2**15

      j      = 2**15
      j_intr = 2**15
      j_bdry = 2**15
      j_s    = 2**15

      k      = 2**15
      k_intr = 2**15
      k_bdry = 2**15
      k_s    = 2**15
c
c     Set the boundary condition coefficients.
c
      if ( (location_index .eq. 0) .or.
     &     (location_index .eq. 1) ) then

         if (location_index .eq. 0) then
            sgn = -1.d0
            i_intr = ilower0
            i_bdry = ilower0-1
            i_s    = ilower0
         else
            sgn = +1.d0
            i_intr = iupper0
            i_bdry = iupper0+1
            i_s    = iupper0+1
         endif

         do k = ilower2,iupper2
            do j = ilower1,iupper1
               if ( abs(acoef(i_s,j,k) - 1.d0) .lt. 1.0d-12 .or.
     &              abs(bcoef(i_s,j,k)       ) .lt. 1.0d-12 ) then
                  if (using_pressure_increment .eq. 1) then
                     gcoef(i_s,j,k) = P_bdry(i_s,j,k) -
     &                    0.5d0*(P(i_intr,j,k)+P(i_bdry,j,k))
                  else
                     gcoef(i_s,j,k) = P_bdry(i_s,j,k)
                  endif
               else
                  gcoef(i_s,j,k) =
     &                 sgn*(rho/dt)*(u0(i_s,j,k) - gcoef(i_s,j,k))
               endif
            enddo
         enddo

      elseif ( (location_index .eq. 2) .or.
     &         (location_index .eq. 3) ) then

         if (location_index .eq. 2) then
            sgn = -1.d0
            j_intr = ilower1
            j_bdry = ilower1-1
            j_s    = ilower1
         else
            sgn = +1.d0
            j_intr = iupper1
            j_bdry = iupper1+1
            j_s    = iupper1+1
         endif

         do k = ilower2,iupper2
            do i = ilower0,iupper0
               if ( abs(acoef(i,j_s,k) - 1.d0) .lt. 1.0d-12 .or.
     &              abs(bcoef(i,j_s,k)       ) .lt. 1.0d-12 ) then
                  if (using_pressure_increment .eq. 1) then
                     gcoef(i,j_s,k) = P_bdry(i,j_s,k) -
     &                    0.5d0*(P(i,j_intr,k)+P(i,j_bdry,k))
                  else
                     gcoef(i,j_s,k) = P_bdry(i,j_s,k)
                  endif
               else
                  gcoef(i,j_s,k) =
     &                 sgn*(rho/dt)*(u1(j_s,k,i) - gcoef(i,j_s,k))
               endif
            enddo
         enddo

      elseif ( (location_index .eq. 4) .or.
     &         (location_index .eq. 5) ) then

         if (location_index .eq. 4) then
            sgn = -1.d0
            k_intr = ilower2
            k_bdry = ilower2-1
            k_s    = ilower2
         else
            sgn = +1.d0
            k_intr = iupper2
            k_bdry = iupper2+1
            k_s    = iupper2+1
         endif

         do j = ilower1,iupper1
            do i = ilower0,iupper0
               if ( abs(acoef(i,j,k_s) - 1.d0) .lt. 1.0d-12 .or.
     &              abs(bcoef(i,j,k_s)       ) .lt. 1.0d-12 ) then
                  if (using_pressure_increment .eq. 1) then
                     gcoef(i,j,k_s) = P_bdry(i,j,k_s) -
     &                    0.5d0*(P(i,j,k_intr)+P(i,j,k_bdry))
                  else
                     gcoef(i,j,k_s) = P_bdry(i,j,k_s)
                  endif
               else
                  gcoef(i,j,k_s) =
     &                 sgn*(rho/dt)*(u2(k_s,i,j) - gcoef(i,j,k_s))
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
c     Set the inhomogeneous Robin boundary condition coefficients for
c     the projection Poisson equation.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine navier_stokes_sc_inhomogeneous_projection_bc_coefs3d(
     &     u0,u1,u2,u_gcw,
     &     P,P_gcw,
     &     acoef,bcoef,gcoef,P_bdry,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
     &     blower0,bupper0,
     &     blower1,bupper1,
     &     blower2,bupper2,
     &     location_index,
     &     using_pressure_increment,
     &     rho,dt)
c
      implicit none
c
c     Input.
c
      INTEGER u_gcw,P_gcw

      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2

      INTEGER blower0,bupper0
      INTEGER blower1,bupper1
      INTEGER blower2,bupper2

      INTEGER location_index,using_pressure_increment

      REAL u0(SIDE3d0(ilower,iupper,u_gcw))
      REAL u1(SIDE3d1(ilower,iupper,u_gcw))
      REAL u2(SIDE3d2(ilower,iupper,u_gcw))

      REAL P(CELL3d(ilower,iupper,P_gcw))

      REAL acoef(blower0:bupper0,
     &           blower1:bupper1,
     &           blower2:bupper2)
      REAL bcoef(blower0:bupper0,
     &           blower1:bupper1,
     &           blower2:bupper2)
      REAL P_bdry(blower0:bupper0,
     &            blower1:bupper1,
     &            blower2:bupper2)

      REAL rho,dt
c
c     Input/Output.
c
      REAL gcoef(blower0:bupper0,
     &           blower1:bupper1,
     &           blower2:bupper2)
c
c     Local variables.
c
      INTEGER i,i_intr,i_bdry,i_s
      INTEGER j,j_intr,j_bdry,j_s
      INTEGER k,k_intr,k_bdry,k_s
      REAL sgn
c
c     Initialize index variables to yield errors in most cases.
c
      i      = 2**15
      i_intr = 2**15
      i_bdry = 2**15
      i_s    = 2**15

      j      = 2**15
      j_intr = 2**15
      j_bdry = 2**15
      j_s    = 2**15

      k      = 2**15
      k_intr = 2**15
      k_bdry = 2**15
      k_s    = 2**15
c
c     Set the boundary condition coefficients.
c
      if ( (location_index .eq. 0) .or.
     &     (location_index .eq. 1) ) then

         if (location_index .eq. 0) then
            sgn = -1.d0
            i_intr = ilower0
            i_bdry = ilower0-1
            i_s    = ilower0
         else
            sgn = +1.d0
            i_intr = iupper0
            i_bdry = iupper0+1
            i_s    = iupper0+1
         endif

         do k = ilower2,iupper2
            do j = ilower1,iupper1
               if ( abs(acoef(i_s,j,k) - 1.d0) .lt. 1.0d-12 .or.
     &              abs(bcoef(i_s,j,k)       ) .lt. 1.0d-12 ) then
                  if (using_pressure_increment .eq. 1) then
                     gcoef(i_s,j,k) = P_bdry(i_s,j,k) -
     &                    0.5d0*(P(i_intr,j,k)+P(i_bdry,j,k))
                  else
                     gcoef(i_s,j,k) = P_bdry(i_s,j,k)
                  endif
               else
                  gcoef(i_s,j,k) =
     &                 sgn*(rho/dt)*(u0(i_s,j,k) - gcoef(i_s,j,k))
               endif
            enddo
         enddo

      elseif ( (location_index .eq. 2) .or.
     &         (location_index .eq. 3) ) then

         if (location_index .eq. 2) then
            sgn = -1.d0
            j_intr = ilower1
            j_bdry = ilower1-1
            j_s    = ilower1
         else
            sgn = +1.d0
            j_intr = iupper1
            j_bdry = iupper1+1
            j_s    = iupper1+1
         endif

         do k = ilower2,iupper2
            do i = ilower0,iupper0
               if ( abs(acoef(i,j_s,k) - 1.d0) .lt. 1.0d-12 .or.
     &              abs(bcoef(i,j_s,k)       ) .lt. 1.0d-12 ) then
                  if (using_pressure_increment .eq. 1) then
                     gcoef(i,j_s,k) = P_bdry(i,j_s,k) -
     &                    0.5d0*(P(i,j_intr,k)+P(i,j_bdry,k))
                  else
                     gcoef(i,j_s,k) = P_bdry(i,j_s,k)
                  endif
               else
                  gcoef(i,j_s,k) =
     &                 sgn*(rho/dt)*(u1(i,j_s,k) - gcoef(i,j_s,k))
               endif
            enddo
         enddo

      elseif ( (location_index .eq. 4) .or.
     &         (location_index .eq. 5) ) then

         if (location_index .eq. 4) then
            sgn = -1.d0
            k_intr = ilower2
            k_bdry = ilower2-1
            k_s    = ilower2
         else
            sgn = +1.d0
            k_intr = iupper2
            k_bdry = iupper2+1
            k_s    = iupper2+1
         endif

         do j = ilower1,iupper1
            do i = ilower0,iupper0
               if ( abs(acoef(i,j,k_s) - 1.d0) .lt. 1.0d-12 .or.
     &              abs(bcoef(i,j,k_s)       ) .lt. 1.0d-12 ) then
                  if (using_pressure_increment .eq. 1) then
                     gcoef(i,j,k_s) = P_bdry(i,j,k_s) -
     &                    0.5d0*(P(i,j,k_intr)+P(i,j,k_bdry))
                  else
                     gcoef(i,j,k_s) = P_bdry(i,j,k_s)
                  endif
               else
                  gcoef(i,j,k_s) =
     &                 sgn*(rho/dt)*(u2(i,j,k_s) - gcoef(i,j,k_s))
               endif
            enddo
         enddo

      endif
c
      return
      end
c
