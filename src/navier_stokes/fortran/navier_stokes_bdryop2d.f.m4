c
c     Routines to set physical boundary condition values.
c
c     Created on 26 Aug 2007
c             by Boyce Griffith (boyce@bigboy.nyconnect.com).
c
c     Last modified: <14.Dec.2007 18:51:03 griffith@box221.cims.nyu.edu>
c
define(NDIM,2)dnl
define(REAL,`double precision')dnl
define(INTEGER,`integer')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim2d.i)dnl
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Modify boundary conditions to enforce incompressibility or
c     no-stress at physical boundaries.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine navier_stokes_open_bc_coefs2d(
     &     U,U_gcw,
     &     acoef,bcoef,gcoef,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     blower0,bupper0,
     &     blower1,bupper1,
     &     location_index,
     &     bdry_normal_axis,
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

      INTEGER blower0,bupper0
      INTEGER blower1,bupper1

      INTEGER location_index,bdry_normal_axis,comp_idx

      REAL U(CELL2d(ilower,iupper,U_gcw),0:NDIM-1)

      REAL acoef(blower0:bupper0,blower1:bupper1)
      REAL bcoef(blower0:bupper0,blower1:bupper1)

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL gcoef(blower0:bupper0,blower1:bupper1)
c
c     Local variables.
c
      INTEGER i,i_intr,i_bdry,i_s
      INTEGER j,j_intr,j_bdry,j_s
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
c
c     Set the boundary condition coefficients.
c
c     At "open" boundaries, modify the normal velocity boundary
c     conditions to enforce div u = 0, and modify the tangential
c     velocity boundary conditions to enforce zero stress.  This is done
c     by specifying a normal flux F at the boundary.
c
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
            do j = ilower1,iupper1
               if ( abs(acoef(i_s,j)       ) .lt. 1.0d-12 .or.
     &              abs(bcoef(i_s,j) - 1.d0) .lt. 1.0d-12 ) then
                  F = (0.25d0/dx(1))*(
     &                 + U(i_intr,j+1,1)
     &                 - U(i_intr,j-1,1)
     &                 + U(i_bdry,j+1,1)
     &                 - U(i_bdry,j-1,1))
                  gcoef(i_s,j) = sgn*F
               endif
            enddo

         elseif (comp_idx .eq. 1) then
c
c     Set F to enforce t * sigma * n = 0.
c
            do j = ilower1,iupper1
               if ( abs(acoef(i_s,j)       ) .lt. 1.0d-12 .or.
     &              abs(bcoef(i_s,j) - 1.d0) .lt. 1.0d-12 ) then
                  F = (0.25d0/dx(1))*(
     &                 + U(i_intr,j+1,bdry_normal_axis)
     &                 - U(i_intr,j-1,bdry_normal_axis)
     &                 + U(i_bdry,j+1,bdry_normal_axis)
     &                 - U(i_bdry,j-1,bdry_normal_axis))
                  gcoef(i_s,j) = sgn*F
               endif
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
            do i = ilower0,iupper0
               if ( abs(acoef(i,j_s)       ) .lt. 1.0d-12 .or.
     &              abs(bcoef(i,j_s) - 1.d0) .lt. 1.0d-12 ) then
                  F = (0.25d0/dx(1))*(
     &                 + U(i+1,j_intr,bdry_normal_axis)
     &                 - U(i-1,j_intr,bdry_normal_axis)
     &                 + U(i+1,j_bdry,bdry_normal_axis)
     &                 - U(i-1,j_bdry,bdry_normal_axis))
                  gcoef(i,j_s) = sgn*F
               endif
            enddo

         elseif (comp_idx .eq. 1) then
c
c     Set F to enforce div u = 0.
c
            do i = ilower0,iupper0
               if ( abs(acoef(i,j_s)       ) .lt. 1.0d-12 .or.
     &              abs(bcoef(i,j_s) - 1.d0) .lt. 1.0d-12 ) then
                  F = (0.25d0/dx(0))*(
     &                 + U(i+1,j_intr,0)
     &                 - U(i-1,j_intr,0)
     &                 + U(i+1,j_bdry,0)
     &                 - U(i-1,j_bdry,0))
                  gcoef(i,j_s) = sgn*F
               endif
            enddo

         endif

      endif
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Modify boundary conditions to enforce tangential boundary
c     conditions at physical boundaries using lagged values of phi.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine navier_stokes_tangential_bc_coefs2d(
     &     Phi,Phi_gcw,
     &     acoef,bcoef,gcoef,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     blower0,bupper0,
     &     blower1,bupper1,
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

      INTEGER blower0,bupper0
      INTEGER blower1,bupper1

      INTEGER location_index,comp_idx

      REAL Phi(CELL2d(ilower,iupper,Phi_gcw))

      REAL acoef(blower0:bupper0,blower1:bupper1)
      REAL bcoef(blower0:bupper0,blower1:bupper1)

      REAL rho
      REAL dx(0:NDIM-1)
      REAL dt
c
c     Input/Output.
c
      REAL gcoef(blower0:bupper0,blower1:bupper1)
c
c     Local variables.
c
      INTEGER i,i_intr,i_bdry,i_s
      INTEGER j,j_intr,j_bdry,j_s
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

         do j = ilower1,iupper1
            if ( abs(acoef(i_s,j) - 1.d0) .lt. 1.0d-12 .or.
     &           abs(bcoef(i_s,j)       ) .lt. 1.0d-12 ) then

               if (comp_idx .eq. 1) then
                  grad_Phi = (0.5d0/dx(1))*(
     &                 Phi(i_intr,j+1) - Phi(i_intr,j-1))
                  gcoef(i_s,j) = gcoef(i_s,j) + (dt/rho)*grad_Phi
               endif

            else

               if (comp_idx .eq. 1) then
                  t_grad_grad_Phi_n = (0.5d0*sgn/dx(1))*(
     &                 + (Phi(i_bdry,j+1) - Phi(i_intr,j+1))/dx(0)
     &                 + (Phi(i_bdry,j-1) - Phi(i_intr,j-1))/dx(0))
                  gcoef(i_s,j) = gcoef(i_s,j) +
     &                 (dt/rho)*t_grad_grad_Phi_n
               endif

            endif
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

         do i = ilower0,iupper0
            if ( abs(acoef(i,j_s) - 1.d0) .lt. 1.0d-12 .or.
     &           abs(bcoef(i,j_s)       ) .lt. 1.0d-12 ) then

               if (comp_idx .eq. 0) then
                  grad_Phi = (0.5d0/dx(0))*(
     &                 Phi(i+1,j_intr) - Phi(i-1,j_intr))
                  gcoef(i,j_s) = gcoef(i,j_s) + (dt/rho)*grad_Phi
               endif

            else

               if (comp_idx .eq. 0) then
                  t_grad_grad_Phi_n = (0.5d0*sgn/dx(0))*(
     &                 + (Phi(i+1,j_bdry) - Phi(i+1,j_intr))/dx(1)
     &                 + (Phi(i-1,j_bdry) - Phi(i-1,j_intr))/dx(1))
                  gcoef(i,j_s) = gcoef(i,j_s) +
     &                 (dt/rho)*t_grad_grad_Phi_n
               endif

            endif
         enddo

      endif
c
      return
      end
c
