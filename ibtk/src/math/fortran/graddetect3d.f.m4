c ---------------------------------------------------------------------
c
c Copyright (c) 2011 - 2019 by the IBAMR developers
c All rights reserved.
c
c This file is part of IBAMR.
c
c IBAMR is free software and is distributed under the 3-clause BSD
c license. The full text of the license can be found in the file
c COPYRIGHT at the top level directory of IBAMR.
c
c ---------------------------------------------------------------------

define(NDIM,3)dnl
define(REAL,`double precision')dnl
define(INTEGER,`integer')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim3d.i)dnl
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Tag cells for refinement based on the gradients of u.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine graddetect3d(
     &     tags,tags_gcw,
     &     u,u_gcw,tol,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
     &     dx)
c
      implicit none
c
c     Input.
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2
      INTEGER tags_gcw,u_gcw

      REAL u(CELL3d(ilower,iupper,u_gcw))
      REAL tol

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      INTEGER tags(CELL3d(ilower,iupper,tags_gcw))
c
c     Local variables
c
      INTEGER i0,i1,i2
      LOGICAL tagcell
      REAL diag(0:NDIM-1),diag012
      REAL presm1,presp1,facejump
c
c     Tag cells with large gradients in u.
c
      diag(0) = sqrt(dx(2)**2+dx(1)**2)
      diag(1) = sqrt(dx(0)**2+dx(2)**2)
      diag(2) = sqrt(dx(0)**2+dx(1)**2)
      diag012 = sqrt(dx(0)**2+dx(1)**2+dx(2)**2)

      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0

               tagcell = .false.

               presm1 = u(i0-1,i1,i2)
               presp1 = u(i0+1,i1,i2)
               facejump = abs(u(i0,i1,i2)-presm1)
               facejump = max(facejump,abs(u(i0,i1,i2)-presp1))

               tagcell = ((facejump).gt.(tol*dx(0)))

               if (.not.tagcell) then
                  presm1 = u(i0,i1-1,i2)
                  presp1 = u(i0,i1+1,i2)
                  facejump = abs(u(i0,i1,i2)-presm1)
                  facejump = max(facejump,abs(u(i0,i1,i2)-presp1))
                  tagcell = ((facejump).gt.(tol*dx(1)))
               endif

               if (.not.tagcell) then
                  presm1 = u(i0,i1,i2-1)
                  presp1 = u(i0,i1,i2+1)
                  facejump = abs(u(i0,i1,i2)-presm1)
                  facejump = max(facejump,abs(u(i0,i1,i2)-presp1))
                  tagcell = ((facejump).gt.(tol*dx(2)))
               endif

c     2Dimensional diagonals

               if (.not.tagcell) then
                  presm1 = u(i0,i1-1,i2-1)
                  presp1 = u(i0,i1+1,i2+1)
                  facejump = abs(u(i0,i1,i2)-presm1)
                  facejump = max(facejump,abs(u(i0,i1,i2)-presp1))
                  tagcell = ((facejump).gt.(tol*diag(0)))
               endif

               if (.not.tagcell) then
                  presm1 = u(i0,i1+1,i2-1)
                  presp1 = u(i0,i1-1,i2+1)
                  facejump = abs(u(i0,i1,i2)-presm1)
                  facejump = max(facejump,abs(u(i0,i1,i2)-presp1))
                  tagcell = ((facejump).gt.(tol*diag(0)))
               endif

               if (.not.tagcell) then
                  presm1 = u(i0-1,i1,i2-1)
                  presp1 = u(i0+1,i1,i2+1)
                  facejump = abs(u(i0,i1,i2)-presm1)
                  facejump = max(facejump,abs(u(i0,i1,i2)-presp1))
                  tagcell = ((facejump).gt.(tol*diag(1)))
               endif

               if (.not.tagcell) then
                  presm1 = u(i0-1,i1,i2+1)
                  presp1 = u(i0+1,i1,i2-1)
                  facejump = abs(u(i0,i1,i2)-presm1)
                  facejump = max(facejump,abs(u(i0,i1,i2)-presp1))
                  tagcell = ((facejump).gt.(tol*diag(1)))
               endif

               if (.not.tagcell) then
                  presm1 = u(i0-1,i1-1,i2)
                  presp1 = u(i0+1,i1+1,i2)
                  facejump = abs(u(i0,i1,i2)-presm1)
                  facejump = max(facejump,abs(u(i0,i1,i2)-presp1))
                  tagcell = ((facejump).gt.(tol*diag(2)))
               endif

               if (.not.tagcell) then
                  presm1 = u(i0-1,i1+1,i2)
                  presp1 = u(i0+1,i1-1,i2)
                  facejump = abs(u(i0,i1,i2)-presm1)
                  facejump = max(facejump,abs(u(i0,i1,i2)-presp1))
                  tagcell = ((facejump).gt.(tol*diag(2)))
               endif

c     3Dimensional diagonals

               if (.not.tagcell) then
                  presm1 = u(i0-1,i1-1,i2-1)
                  presp1 = u(i0+1,i1+1,i2+1)
                  facejump = abs(u(i0,i1,i2)-presm1)
                  facejump = max(facejump,abs(u(i0,i1,i2)-presp1))
                  tagcell = ((facejump).gt.(tol*diag012))
               endif

               if (.not.tagcell) then
                  presm1 = u(i0-1,i1-1,i2+1)
                  presp1 = u(i0+1,i1+1,i2-1)
                  facejump = abs(u(i0,i1,i2)-presm1)
                  facejump = max(facejump,abs(u(i0,i1,i2)-presp1))
                  tagcell = ((facejump).gt.(tol*diag012))
               endif

               if (.not.tagcell) then
                  presm1 = u(i0-1,i1+1,i2-1)
                  presp1 = u(i0+1,i1-1,i2+1)
                  facejump = abs(u(i0,i1,i2)-presm1)
                  facejump = max(facejump,abs(u(i0,i1,i2)-presp1))
                  tagcell = ((facejump).gt.(tol*diag012))
               endif

               if (.not.tagcell) then
                  presm1 = u(i0-1,i1+1,i2+1)
                  presp1 = u(i0+1,i1-1,i2-1)
                  facejump = abs(u(i0,i1,i2)-presm1)
                  facejump = max(facejump,abs(u(i0,i1,i2)-presp1))
                  tagcell = ((facejump).gt.(tol*diag012))
               endif
c
c     Tag cells.
c
               if ( tagcell ) then
                  tags(i0,i1,i2) = 1
               endif

            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
