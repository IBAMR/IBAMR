c
c     Routines to interpolate cell data to and from face/side data on a
c     patch.
c
c     Created on 12 Jun 2003 by Boyce Griffith
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
c     Compute the face centered normal vector field (u0,u1) from the
c     cell centered vector field V using simple averaging.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ctofinterp2nd2d(
     &     u0,u1,u_gcw,
     &     V,V_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1)
c
      implicit none
c
c     Input.
c
      INTEGER u_gcw,V_gcw

      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1

      REAL V(CELL2d(ilower,iupper,V_gcw),0:NDIM-1)
c
c     Output.
c
      REAL u0(FACE2d0(ilower,iupper,u_gcw))
      REAL u1(FACE2d1(ilower,iupper,u_gcw))
c
c     Local variables.
c
      INTEGER i0,i1
c
c     Compute the face centered vector field (u0,u1) from the cell
c     centered vector field V.
c
      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0+1
            u0(i0,i1) = 0.5d0*(V(i0-1,i1,0)+V(i0,i1,0))
         enddo
      enddo
      do i0 = ilower0,iupper0
         do i1 = ilower1,iupper1+1
            u1(i1,i0) = 0.5d0*(V(i0,i1-1,1)+V(i0,i1,1))
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the face centered normal vector field (u0,u1) from the
c     cell centered scalar field V using simple averaging.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ctofcwiseinterp2nd2d(
     &     u0,u1,u_gcw,
     &     V,V_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1)
c
      implicit none
c
c     Input.
c
      INTEGER u_gcw,V_gcw

      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1

      REAL V(CELL2d(ilower,iupper,V_gcw))
c
c     Output.
c
      REAL u0(FACE2d0(ilower,iupper,u_gcw))
      REAL u1(FACE2d1(ilower,iupper,u_gcw))
c
c     Local variables.
c
      INTEGER i0,i1
c
c     Compute the face centered vector field (u0,u1) from the cell
c     centered scalar field V.
c
      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0+1
            u0(i0,i1) = 0.5d0*(V(i0-1,i1)+V(i0,i1))
         enddo
      enddo
      do i0 = ilower0,iupper0
         do i1 = ilower1,iupper1+1
            u1(i1,i0) = 0.5d0*(V(i0,i1-1)+V(i0,i1))
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the side centered normal vector field (u0,u1) from the
c     cell centered vector field V using simple averaging.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ctosinterp2nd2d(
     &     u0,u1,u_gcw,
     &     V,V_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1)
c
      implicit none
c
c     Input.
c
      INTEGER u_gcw,V_gcw

      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1

      REAL V(CELL2d(ilower,iupper,V_gcw),0:NDIM-1)
c
c     Output.
c
      REAL u0(SIDE2d0(ilower,iupper,u_gcw))
      REAL u1(SIDE2d1(ilower,iupper,u_gcw))
c
c     Local variables.
c
      INTEGER i0,i1
c
c     Compute the side centered vector field (u0,u1) from the cell
c     centered vector field V.
c
      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0+1
            u0(i0,i1) = 0.5d0*(V(i0-1,i1,0)+V(i0,i1,0))
         enddo
      enddo
      do i1 = ilower1,iupper1+1
         do i0 = ilower0,iupper0
            u1(i0,i1) = 0.5d0*(V(i0,i1-1,1)+V(i0,i1,1))
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the side centered normal vector field (u0,u1) from the
c     cell centered scalar field V using simple averaging.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ctoscwiseinterp2nd2d(
     &     u0,u1,u_gcw,
     &     V,V_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1)
c
      implicit none
c
c     Input.
c
      INTEGER u_gcw,V_gcw

      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1

      REAL V(CELL2d(ilower,iupper,V_gcw))
c
c     Output.
c
      REAL u0(SIDE2d0(ilower,iupper,u_gcw))
      REAL u1(SIDE2d1(ilower,iupper,u_gcw))
c
c     Local variables.
c
      INTEGER i0,i1
c
c     Compute the side centered vector field (u0,u1) from the cell
c     centered scalar field V.
c
      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0+1
            u0(i0,i1) = 0.5d0*(V(i0-1,i1)+V(i0,i1))
         enddo
      enddo
      do i1 = ilower1,iupper1+1
         do i0 = ilower0,iupper0
            u1(i0,i1) = 0.5d0*(V(i0,i1-1)+V(i0,i1))
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the cell centered vector field U from the face centered
c     normal vector field (v0,v1) using simple averaging.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ftocinterp2nd2d(
     &     U,U_gcw,
     &     v0,v1,v_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1)
c
      implicit none
c
c     Input.
c
      INTEGER U_gcw,v_gcw

      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1

      REAL v0(FACE2d0(ilower,iupper,v_gcw))
      REAL v1(FACE2d1(ilower,iupper,v_gcw))
c
c     Output.
c
      REAL U(CELL2d(ilower,iupper,U_gcw),0:NDIM-1)
c
c     Local variables.
c
      INTEGER i0,i1
c
c     Compute the cell centered vector field U from the face centered
c     vector field (v0,v1).
c
      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0
            U(i0,i1,0) = 0.5d0*(v0(i0,i1)+v0(i0+1,i1))
         enddo
      enddo
      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0
            U(i0,i1,1) = 0.5d0*(v1(i1,i0)+v1(i1+1,i0))
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the cell centered vector field U from the side centered
c     normal vector field (v0,v1) using simple averaging.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine stocinterp2nd2d(
     &     U,U_gcw,
     &     v0,v1,v_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1)
c
      implicit none
c
c     Input.
c
      INTEGER U_gcw,v_gcw

      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1

      REAL v0(SIDE2d0(ilower,iupper,v_gcw))
      REAL v1(SIDE2d1(ilower,iupper,v_gcw))
c
c     Output.
c
      REAL U(CELL2d(ilower,iupper,U_gcw),0:NDIM-1)
c
c     Local variables.
c
      INTEGER i0,i1
c
c     Compute the cell centered vector field U from the side centered
c     vector field (v0,v1).
c
      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0
            U(i0,i1,0) = 0.5d0*(v0(i0,i1)+v0(i0+1,i1))
         enddo
      enddo
      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0
            U(i0,i1,1) = 0.5d0*(v1(i0,i1)+v1(i0,i1+1))
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the cell centered vector field U from the side centered
c     normal vector field (v0,v1) using simple averaging.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine stocinterp2ndspecial2d(
     &     direction,
     &     U,U_gcw,
     &     alpha,
     &     v0,v1,v_gcw,
     &     beta,
     &     W,W_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1)
c
      implicit none
c
c     Input.
c
      INTEGER direction

      INTEGER U_gcw,v_gcw,W_gcw

      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1

      REAL alpha

      REAL v0(SIDE2d0(ilower,iupper,v_gcw))
      REAL v1(SIDE2d1(ilower,iupper,v_gcw))

      REAL beta

      REAL W(CELL2d(ilower,iupper,W_gcw))
c
c     Output.
c
      REAL U(CELL2d(ilower,iupper,U_gcw))
c
c     Local variables.
c
      INTEGER i0,i1
c
c     Compute the cell centered vector field U from the side centered
c     vector field (v0,v1).
c
      if ( direction.eq.0 ) then
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
               U(i0,i1) = 0.5d0*alpha*(v0(i0,i1)+v0(i0+1,i1))
     &              + beta*W(i0,i1)
            enddo
         enddo
      else
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
               U(i0,i1) = 0.5d0*alpha*(v1(i0,i1)+v1(i0,i1+1))
     &              + beta*W(i0,i1)
            enddo
         enddo
      endif
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
