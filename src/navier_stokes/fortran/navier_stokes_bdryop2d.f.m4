c
c     Routines to set physical boundary condition values.
c
c     Created on 26 Aug 2007
c             by Boyce Griffith (boyce@bigboy.nyconnect.com).
c
c     Last modified: <26.Aug.2007 02:37:24 boyce@bigboy.nyconnect.com>
c
define(NDIM,2)dnl
define(REAL,`double precision')dnl
define(INTEGER,`integer')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim2d.i)dnl
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Set cell centered normal velocity boundary values by approximating
c     div(u) = 0 at the boundary.
c
c     WARNING: This function assumes that the grid spacing is equal in
c     all coordinate directions.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine navier_stokes_normal_velocity2d(
     &     U,U_gcw,
     &     V,V_gcw,
     &     location_index,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     blower0,bupper0,
     &     blower1,bupper1)
c
      implicit none
c
c     Input.
c
      INTEGER U_gcw,V_gcw

      INTEGER location_index

      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1

      INTEGER blower0,bupper0
      INTEGER blower1,bupper1
c
c     Input/Output.
c
      REAL U(CELL2d(ilower,iupper,U_gcw))
      REAL V(CELL2d(ilower,iupper,V_gcw))
c
c     Local variables.
c
      INTEGER i,j
c
c     Correct the values along the upper/lower x side of the patch.
c
      if (location_index .eq. 0) then
         i = ilower0            ! interior index
         do j = blower1,bupper1
            U(i-1,j) = U(i+1,j)
     &           +V(i,j+1)-V(i,j-1)
     &           +0.2d0*(+U(i+2,j  )-U(i  ,j  )
     &                   +V(i+1,j+1)-V(i+1,j-1))
         enddo
      elseif (location_index .eq. 1) then
         i = iupper0            ! interior index
         do j = blower1,bupper1
            U(i+1,j) = U(i-1,j)
     &           -V(i,j+1)+V(i,j-1)
     &           -0.2d0*(+U(i  ,j  )-U(i-2,j  )
     &                   +V(i-1,j+1)-V(i-1,j-1))
         enddo
c
c     Correct the values along the upper/lower y side of the patch.
c
      elseif (location_index .eq. 2) then
         j = ilower1            ! interior index
         do i = blower0,bupper0
            V(i,j-1) = V(i,j+1)
     &           +U(i+1,j)-U(i-1,j)
     &           +0.2d0*(+U(i+1,j+1)-U(i-1,j+1)
     &                   +V(i  ,j  )-V(i  ,j+2))
         enddo
      elseif (location_index .eq. 3) then
         j = iupper1            ! interior index
         do i = blower0,bupper0
            V(i,j+1) = V(i,j-1)
     &           -U(i+1,j)+U(i-1,j)
     &           -0.2d0*(+U(i+1,j-1)-U(i-1,j-1)
     &                   +V(i  ,j  )-V(i  ,j-2))
         enddo
      endif
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Set cell centered tangential velocity boundary values by
c     approximating t*[grad u + (grad u)^T]*n = 0 at the boundary.
c
c     WARNING: This function assumes that the grid spacing is equal in
c     all coordinate directions.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine navier_stokes_tangential_velocity2d(
     &     U,U_gcw,
     &     V,V_gcw,
     &     location_index,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     blower0,bupper0,
     &     blower1,bupper1)
c
      implicit none
c
c     Input.
c
      INTEGER U_gcw,V_gcw

      INTEGER location_index

      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1

      INTEGER blower0,bupper0
      INTEGER blower1,bupper1
c
c     Input/Output.
c
      REAL U(CELL2d(ilower,iupper,U_gcw))
      REAL V(CELL2d(ilower,iupper,V_gcw))
c
c     Local variables.
c
      INTEGER i,j
c
c     Correct the values along the upper/lower x side of the patch.
c
      if (location_index .eq. 0) then
         i = ilower0            ! interior index
         do j = blower1,bupper1
            V(i-1,j) = V(i+1,j)
     &           +U(i,j+1)-U(i,j-1)
     &           +0.2d0*(+V(i+2,j  )-V(i  ,j  )
     &                   +U(i+1,j+1)-U(i+1,j-1))
         enddo
      elseif (location_index .eq. 1) then
         i = iupper0            ! interior index
         do j = blower1,bupper1
            V(i+1,j) = V(i-1,j)
     &           -U(i,j+1)+U(i,j-1)
     &           +0.2d0*(-V(i  ,j  )+V(i-2,j  )
     &                   -U(i-1,j+1)+U(i-1,j-1))
         enddo
c
c     Correct the values along the upper/lower y side of the patch.
c
      elseif (location_index .eq. 2) then
         j = ilower1            ! interior index
         do i = blower0,bupper0
            U(i,j-1) = U(i,j+1)
     &           +V(i+1,j)-V(i-1,j)
     &           +0.5d0*(+U(i  ,j+2)-U(i  ,j  )
     &                   +V(i+1,j+1)-V(i-1,j+1))
         enddo
      elseif (location_index .eq. 3) then
         j = iupper1            ! interior index
         do i = blower0,bupper0
            U(i,j+1) = U(i,j-1)
     &           -V(i+1,j)+V(i-1,j)
     &           +0.2d0*(-U(i  ,j  )+U(i  ,j-2)
     &                   -V(i+1,j-1)+V(i-1,j-1))
         enddo
      endif
c
      return
      end
c
