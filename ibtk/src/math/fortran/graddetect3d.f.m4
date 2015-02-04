c
c     Routines for detecting steep gradients.
c
c     Created on 27 Nov 2003 by Boyce Griffith
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
