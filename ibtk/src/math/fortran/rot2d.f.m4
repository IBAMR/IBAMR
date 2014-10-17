c
c     Routines to compute discrete rots on patches.
c
c     Created on 06 Aug 2013 by Amneet Bhalla
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
c     Computes W = rot U = [dU/dx1, -dU/dx0].
c
c     Uses centered differences to compute the side centered rot of a
c     node centered scalar field U.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ntosrot2d(
     &     W0,W1,W_gcw,
     &     U,U_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     dx)
c
      implicit none
c
c     Input.
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER W_gcw,U_gcw

      REAL U(NODE2d(ilower,iupper,U_gcw))

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL W0(SIDE2d0(ilower,iupper,W_gcw))
      REAL W1(SIDE2d1(ilower,iupper,W_gcw))
c
c     Local variables.
c
      INTEGER i0,i1
c
c     Compute the side centered rot of node centered U.
c

      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0+1
            W0(i0,i1) = (U(i0,i1+1)-U(i0,i1))/dx(1)
         enddo
      enddo

      do i1 = ilower1,iupper1+1
         do i0 = ilower0,iupper0
            W1(i0,i1) = -(U(i0+1,i1)-U(i0,i1))/dx(0)
         enddo
      enddo

c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes W = rot U = [dU/dy, -dU/dx]
c
c     Uses centered differences to compute the side centered rot of a
c     cell centered scalar field U.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ctosrot2d(
     &     W0,W1,W_gcw,
     &     U,U_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     dx)
c
      implicit none
c
c     Input.
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER W_gcw,U_gcw

      REAL U(CELL2d(ilower,iupper,U_gcw))

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL W0(SIDE2d0(ilower,iupper,W_gcw))
      REAL W1(SIDE2d1(ilower,iupper,W_gcw))
c
c     Local variables.
c
      INTEGER i0,i1
c
c     Compute the side centered rot of cell centered U.
c

      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0+1
            W0(i0,i1) = (U(i0-1,i1+1)-U(i0-1,i1-1) + 
     &                   U(i0,i1+1)-U(i0,i1-1))/(0.4d1*dx(1))
         enddo
      enddo

      do i1 = ilower1,iupper1+1
         do i0 = ilower0,iupper0
            W1(i0,i1) = -(U(i0+1,i1-1)-U(i0-1,i1-1) + 
     &                    U(i0+1,i1)-U(i0-1,i1))/(0.4d1*dx(0))
         enddo
      enddo

c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
