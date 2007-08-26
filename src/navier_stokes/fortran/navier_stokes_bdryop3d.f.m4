c
c     Routines to set physical boundary condition values.
c
c     Created on 21 May 2007
c             by Boyce Griffith (griffith@box221.cims.nyu.edu).
c
c     Last modified: <26.Aug.2007 02:35:05 boyce@bigboy.nyconnect.com>
c
define(NDIM,3)dnl
define(REAL,`double precision')dnl
define(INTEGER,`integer')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim3d.i)dnl
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Set cell centered normal velocity boundary values by approximating
c     div(u) = 0 at the boundary.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine navier_stokes_normal_velocity3d(
     &     U,U_gcw,
     &     V,V_gcw,
     &     W,W_gcw,
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
      INTEGER U_gcw,V_gcw,W_gcw

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
      REAL V(CELL3d(ilower,iupper,V_gcw))
      REAL W(CELL3d(ilower,iupper,W_gcw))
c
c     Local variables.
c
      INTEGER i,j,k
c
c     Correct the values along the upper/lower x side of the patch.
c
      if (location_index .eq. 0) then
         i = ilower0            ! interior index
         do k = blower2,bupper2
            do j = blower1,bupper1
            enddo
         enddo
      elseif (location_index .eq. 1) then
         i = iupper0            ! interior index
         do k = blower2,bupper2
            do j = blower1,bupper1
            enddo
         enddo
c
c     Correct the values along the upper/lower y side of the patch.
c
      elseif (location_index .eq. 2) then
         j = ilower1            ! interior index
         do k = blower2,bupper2
            do i = blower0,bupper0
            enddo
         enddo
      elseif (location_index .eq. 3) then
         j = iupper1            ! interior index
         do k = blower2,bupper2
            do i = blower0,bupper0
            enddo
         enddo
c
c     Correct the values along the upper/lower z side of the patch.
c
      elseif (location_index .eq. 4) then
         k = ilower2            ! interior index
         do j = blower1,bupper1
            do i = blower0,bupper0
            enddo
         enddo
      elseif (location_index .eq. 5) then
         k = iupper2            ! interior index
         do j = blower1,bupper1
            do i = blower0,bupper0
            enddo
         enddo
      endif
c
      return
      end
c
