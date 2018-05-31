c
c     Routines to interpolate cell data to and from face/side data on a
c     patch.
c
c     Created on 12 Jun 2003 by Boyce Griffith
c
c     Copyright (c) 2002-2017, Boyce Griffith
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
c     Compute the face centered normal vector field (u0,u1,u2) from the
c     cell centered vector field V using simple averaging.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ctofinterp2nd3d(
     &     u0,u1,u2,u_gcw,
     &     V,V_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2)
c
      implicit none
c
c     Input.
c
      INTEGER u_gcw,V_gcw

      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2

      REAL V(CELL3d(ilower,iupper,V_gcw),0:NDIM-1)
c
c     Output.
c
      REAL u0(FACE3d0(ilower,iupper,u_gcw))
      REAL u1(FACE3d1(ilower,iupper,u_gcw))
      REAL u2(FACE3d2(ilower,iupper,u_gcw))
c
c     Local variables.
c
      INTEGER i0,i1,i2
c
c     Compute the face centered vector field (u0,u1,u2) from the cell
c     centered vector field V.
c
      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0+1
               u0(i0,i1,i2) = 0.5d0*(V(i0-1,i1,i2,0)+V(i0,i1,i2,0))
            enddo
         enddo
      enddo
      do i0 = ilower0,iupper0
         do i2 = ilower2,iupper2
            do i1 = ilower1,iupper1+1
               u1(i1,i2,i0) = 0.5d0*(V(i0,i1-1,i2,1)+V(i0,i1,i2,1))
            enddo
         enddo
      enddo
      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0
            do i2 = ilower2,iupper2+1
               u2(i2,i0,i1) = 0.5d0*(V(i0,i1,i2-1,2)+V(i0,i1,i2,2))
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the face centered normal vector field (u0,u1,u2) from the
c     cell centered scalar field V using simple averaging.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ctofcwiseinterp2nd3d(
     &     u0,u1,u2,u_gcw,
     &     V,V_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2)
c
      implicit none
c
c     Input.
c
      INTEGER u_gcw,V_gcw

      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2

      REAL V(CELL3d(ilower,iupper,V_gcw))
c
c     Output.
c
      REAL u0(FACE3d0(ilower,iupper,u_gcw))
      REAL u1(FACE3d1(ilower,iupper,u_gcw))
      REAL u2(FACE3d2(ilower,iupper,u_gcw))
c
c     Local variables.
c
      INTEGER i0,i1,i2
c
c     Compute the face centered vector field (u0,u1,u2) from the cell
c     centered scalar field V.
c
      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0+1
               u0(i0,i1,i2) = 0.5d0*(V(i0-1,i1,i2)+V(i0,i1,i2))
            enddo
         enddo
      enddo
      do i0 = ilower0,iupper0
         do i2 = ilower2,iupper2
            do i1 = ilower1,iupper1+1
               u1(i1,i2,i0) = 0.5d0*(V(i0,i1-1,i2)+V(i0,i1,i2))
            enddo
         enddo
      enddo
      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0
            do i2 = ilower2,iupper2+1
               u2(i2,i0,i1) = 0.5d0*(V(i0,i1,i2-1)+V(i0,i1,i2))
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the side centered normal vector field (u0,u1,u2) from the
c     cell centered vector field V using simple averaging.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ctosinterp2nd3d(
     &     u0,u1,u2,u_gcw,
     &     V,V_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2)
c
      implicit none
c
c     Input.
c
      INTEGER u_gcw,V_gcw

      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2

      REAL V(CELL3d(ilower,iupper,V_gcw),0:NDIM-1)
c
c     Output.
c
      REAL u0(SIDE3d0(ilower,iupper,u_gcw))
      REAL u1(SIDE3d1(ilower,iupper,u_gcw))
      REAL u2(SIDE3d2(ilower,iupper,u_gcw))
c
c     Local variables.
c
      INTEGER i0,i1,i2
c
c     Compute the side centered vector field (u0,u1,u2) from the cell
c     centered vector field V.
c
      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0+1
               u0(i0,i1,i2) = 0.5d0*(V(i0-1,i1,i2,0)+V(i0,i1,i2,0))
            enddo
         enddo
      enddo
      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1+1
            do i0 = ilower0,iupper0
               u1(i0,i1,i2) = 0.5d0*(V(i0,i1-1,i2,1)+V(i0,i1,i2,1))
            enddo
         enddo
      enddo
      do i2 = ilower2,iupper2+1
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
               u2(i0,i1,i2) = 0.5d0*(V(i0,i1,i2-1,2)+V(i0,i1,i2,2))
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the side centered normal vector field (u0,u1,u2) from the
c     cell centered scalar field V using simple averaging.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ctoscwiseinterp2nd3d(
     &     u0,u1,u2,u_gcw,
     &     V,V_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2)
c
      implicit none
c
c     Input.
c
      INTEGER u_gcw,V_gcw

      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2

      REAL V(CELL3d(ilower,iupper,V_gcw))
c
c     Output.
c
      REAL u0(SIDE3d0(ilower,iupper,u_gcw))
      REAL u1(SIDE3d1(ilower,iupper,u_gcw))
      REAL u2(SIDE3d2(ilower,iupper,u_gcw))
c
c     Local variables.
c
      INTEGER i0,i1,i2
c
c     Compute the side centered vector field (u0,u1,u2) from the cell
c     centered scalar field V.
c
      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0+1
               u0(i0,i1,i2) = 0.5d0*(V(i0-1,i1,i2)+V(i0,i1,i2))
            enddo
         enddo
      enddo
      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1+1
            do i0 = ilower0,iupper0
               u1(i0,i1,i2) = 0.5d0*(V(i0,i1-1,i2)+V(i0,i1,i2))
            enddo
         enddo
      enddo
      do i2 = ilower2,iupper2+1
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
               u2(i0,i1,i2) = 0.5d0*(V(i0,i1,i2-1)+V(i0,i1,i2))
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the cell centered vector field U from the face centered
c     normal vector field (v0,v1,v2) using simple averaging.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ftocinterp2nd3d(
     &     U,U_gcw,
     &     v0,v1,v2,v_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2)
c
      implicit none
c
c     Input.
c
      INTEGER U_gcw,v_gcw

      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2

      REAL v0(FACE3d0(ilower,iupper,v_gcw))
      REAL v1(FACE3d1(ilower,iupper,v_gcw))
      REAL v2(FACE3d2(ilower,iupper,v_gcw))
c
c     Output.
c
      REAL U(CELL3d(ilower,iupper,U_gcw),0:NDIM-1)
c
c     Local variables.
c
      INTEGER i0,i1,i2
c
c     Compute the cell centered vector field U from the face centered
c     vector field (v0,v1,v2).
c
      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
               U(i0,i1,i2,0) = 0.5d0*(v0(i0,i1,i2)+v0(i0+1,i1,i2))
            enddo
         enddo
      enddo
      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
               U(i0,i1,i2,1) = 0.5d0*(v1(i1,i2,i0)+v1(i1+1,i2,i0))
            enddo
         enddo
      enddo
      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
               U(i0,i1,i2,2) = 0.5d0*(v2(i2,i0,i1)+v2(i2+1,i0,i1))
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the cell centered vector field U from the side centered
c     normal vector field (v0,v1,v2) using simple averaging.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine stocinterp2nd3d(
     &     U,U_gcw,
     &     v0,v1,v2,v_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2)
c
      implicit none
c
c     Input.
c
      INTEGER U_gcw,v_gcw

      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2

      REAL v0(SIDE3d0(ilower,iupper,v_gcw))
      REAL v1(SIDE3d1(ilower,iupper,v_gcw))
      REAL v2(SIDE3d2(ilower,iupper,v_gcw))
c
c     Output.
c
      REAL U(CELL3d(ilower,iupper,U_gcw),0:NDIM-1)
c
c     Local variables.
c
      INTEGER i0,i1,i2
c
c     Compute the cell centered vector field U from the side centered
c     vector field (v0,v1,v2).
c
      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
               U(i0,i1,i2,0) = 0.5d0*(v0(i0,i1,i2)+v0(i0+1,i1,i2))
            enddo
         enddo
      enddo
      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
               U(i0,i1,i2,1) = 0.5d0*(v1(i0,i1,i2)+v1(i0,i1+1,i2))
            enddo
         enddo
      enddo
      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
               U(i0,i1,i2,2) = 0.5d0*(v2(i0,i1,i2)+v2(i0,i1,i2+1))
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the cell centered vector field U from the side centered
c     normal vector field (v0,v1,v2) using simple averaging.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine stocinterp2ndspecial3d(
     &     direction,
     &     U,U_gcw,
     &     alpha,
     &     v0,v1,v2,v_gcw,
     &     beta,
     &     W,W_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2)
c
      implicit none
c
c     Input.
c
      INTEGER direction

      INTEGER U_gcw,v_gcw,W_gcw

      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2

      REAL alpha

      REAL v0(SIDE3d0(ilower,iupper,v_gcw))
      REAL v1(SIDE3d1(ilower,iupper,v_gcw))
      REAL v2(SIDE3d2(ilower,iupper,v_gcw))

      REAL beta

      REAL W(CELL3d(ilower,iupper,W_gcw))
c
c     Output.
c
      REAL U(CELL3d(ilower,iupper,U_gcw))
c
c     Local variables.
c
      INTEGER i0,i1,i2
c
c     Compute the cell centered vector field U from the side centered
c     vector field (v0,v1,v2).
c
      if ( direction.eq.0 ) then
         do i2 = ilower2,iupper2
            do i1 = ilower1,iupper1
               do i0 = ilower0,iupper0
                  U(i0,i1,i2)=0.5d0*alpha*(v0(i0,i1,i2)+v0(i0+1,i1,i2))
     &                 + beta*W(i0,i1,i2)
               enddo
            enddo
         enddo
      elseif ( direction.eq.1 ) then
         do i2 = ilower2,iupper2
            do i1 = ilower1,iupper1
               do i0 = ilower0,iupper0
                  U(i0,i1,i2)=0.5d0*alpha*(v1(i0,i1,i2)+v1(i0,i1+1,i2))
     &                 + beta*W(i0,i1,i2)
               enddo
            enddo
         enddo
      else
         do i2 = ilower2,iupper2
            do i1 = ilower1,iupper1
               do i0 = ilower0,iupper0
                  U(i0,i1,i2)=0.5d0*alpha*(v2(i0,i1,i2)+v2(i0,i1,i2+1))
     &                 + beta*W(i0,i1,i2)
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
c     Compute the cell centered field U from the edge centered
c     field (v0,v1,v2) using simple averaging.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine etocinterp3d(
     &     U,U_gcw,
     &     v0,v1,v2,v_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2)
c
      implicit none
c
c     Input.
c
      INTEGER U_gcw,v_gcw

      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2

      REAL v0(EDGE3d0(ilower,iupper,v_gcw))
      REAL v1(EDGE3d1(ilower,iupper,v_gcw))
      REAL v2(EDGE3d2(ilower,iupper,v_gcw))
c
c     Output.
c
      REAL U(CELL3d(ilower,iupper,U_gcw))
c
c     Local variables.
c
      INTEGER i0,i1,i2
      REAL avg0,avg1,avg2
c
c     Compute the cell centered field U from the edge centered
c     field (v0,v1,v2).
c
      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
              avg0 = v0(i0,i1,i2)+v0(i0,i1+1,i2)+
     &               v0(i0,i1,i2+1)+v0(i0,i1+1,i2+1)
              avg1 = v1(i0,i1,i2)+v1(i0+1,i1,i2)+
     &               v1(i0,i1,i2+1)+v1(i0+1,i1,i2+1)
              avg2 = v2(i0,i1,i2)+v2(i0+1,i1,i2)+
     &               v2(i0,i1+1,i2)+v2(i0+1,i1+1,i2)
              U(i0,i1,i2) = (avg0+avg1+avg2)/12.d0
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the edge centered vector field (u0,u1,u2) from the cell centered
c     field V using simple averaging.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ctoeinterp3d(
     &     u0,u1,u2,u_gcw,
     &     V,V_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
     &     U_ghost_interp)
     
c
      implicit none
c
c     Input.
c
      INTEGER u_gcw,V_gcw

      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2
      INTEGER U_ghost_interp

      REAL V(CELL3d(ilower,iupper,V_gcw))
c
c     Output.
c
      REAL u0(EDGE3d0(ilower,iupper,u_gcw))
      REAL u1(EDGE3d1(ilower,iupper,u_gcw))
      REAL u2(EDGE3d2(ilower,iupper,u_gcw))
c
c     Local variables.
c
      INTEGER i0,i1,i2
      INTEGER gcw_shift

c     Compute the edge centered interpolation of V
      gcw_shift = 0
      if (U_ghost_interp .eq. 1) then
         gcw_shift = U_gcw
      endif
      
      do i2 = ilower2-gcw_shift,iupper2+gcw_shift+1
         do i1 = ilower1-gcw_shift,iupper1+gcw_shift+1
            do i0 = ilower0-gcw_shift,iupper0+gcw_shift
               u0(i0,i1,i2) = 0.25d0*(V(i0,i1-1,i2) + V(i0,i1,i2-1)  
     &                      + V(i0,i1,i2) + V(i0,i1-1,i2-1))
     
            enddo
         enddo
      enddo 
  
      do i2 = ilower2-gcw_shift,iupper2+gcw_shift+1
         do i1 = ilower1-gcw_shift,iupper1+gcw_shift
            do i0 = ilower0-gcw_shift,iupper0+gcw_shift+1
               u1(i0,i1,i2) = 0.25d0*(V(i0,i1,i2) + V(i0,i1,i2-1)  
     &                      + V(i0-1,i1,i2) + V(i0-1,i1,i2-1))

            enddo
         enddo
      enddo

      do i2 = ilower2-gcw_shift,iupper2+gcw_shift
         do i1 = ilower1-gcw_shift,iupper1+gcw_shift+1
            do i0 = ilower0-gcw_shift,iupper0+gcw_shift+1
               u2(i0,i1,i2) = 0.25d0*(V(i0,i1,i2) + V(i0,i1-1,i2)  
     &                      + V(i0-1,i1,i2) + V(i0-1,i1-1,i2))

            enddo
         enddo
      enddo

c
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the side centered normal vector field (u0,u1,u2) from the
c     cell centered vector field V using harmonic averaging.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ctosharmonicinterp2nd3d(
     &     u0,u1,u2,u_gcw,
     &     V,V_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2)
c
      implicit none
c
c     Input.
c
      INTEGER u_gcw,V_gcw

      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2

      REAL V(CELL3d(ilower,iupper,V_gcw),0:NDIM-1)
c
c     Output.
c
      REAL u0(SIDE3d0(ilower,iupper,u_gcw))
      REAL u1(SIDE3d1(ilower,iupper,u_gcw))
      REAL u2(SIDE3d2(ilower,iupper,u_gcw))
c
c     Local variables.
c
      INTEGER i0,i1,i2
      REAL nmr,dmr
c
c     Compute the side centered vector field (u0,u1,u2) from the cell
c     centered vector field V.
c
      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0+1
               nmr = 2.d0*V(i0-1,i1,i2,0)*V(i0,i1,i2,0)
               dmr = V(i0-1,i1,i2,0)+V(i0,i1,i2,0)
               u0(i0,i1,i2) = nmr/dmr
            enddo
         enddo
      enddo
      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1+1
            do i0 = ilower0,iupper0
               nmr = 2.d0*V(i0,i1-1,i2,1)*V(i0,i1,i2,1)
               dmr = V(i0,i1-1,i2,1)+V(i0,i1,i2,1)
               u1(i0,i1,i2) = nmr/dmr
            enddo
         enddo
      enddo
      do i2 = ilower2,iupper2+1
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
               nmr = 2*V(i0,i1,i2-1,2)*V(i0,i1,i2,2)
               dmr = V(i0,i1,i2-1,2)+V(i0,i1,i2,2)
               u2(i0,i1,i2) = nmr/dmr
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the edge centered vector field (u0,u1,u2) from the cell centered
c     field V using harmonic averaging.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ctoeharmonicinterp3d(
     &     u0,u1,u2,u_gcw,
     &     V,V_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
     &     U_ghost_interp)
     
c
      implicit none
c
c     Input.
c
      INTEGER u_gcw,V_gcw

      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2
      INTEGER U_ghost_interp

      REAL V(CELL3d(ilower,iupper,V_gcw))
c
c     Output.
c
      REAL u0(EDGE3d0(ilower,iupper,u_gcw))
      REAL u1(EDGE3d1(ilower,iupper,u_gcw))
      REAL u2(EDGE3d2(ilower,iupper,u_gcw))
c
c     Local variables.
c
      INTEGER i0,i1,i2
      INTEGER gcw_shift
      REAL    nmr,dmr

c     Compute the edge centered interpolation of V
      gcw_shift = 0
      if (U_ghost_interp .eq. 1) then
         gcw_shift = U_gcw
      endif

      nmr = 4.d0
      
      do i2 = ilower2-gcw_shift,iupper2+gcw_shift+1
         do i1 = ilower1-gcw_shift,iupper1+gcw_shift+1
            do i0 = ilower0-gcw_shift,iupper0+gcw_shift
               dmr = 1.d0/V(i0,i1-1,i2) + 1.d0/V(i0,i1,i2-1)  
     &               + 1.d0/V(i0,i1,i2) + 1.d0/V(i0,i1-1,i2-1)
               u0(i0,i1,i2) = nmr/dmr
     
            enddo
         enddo
      enddo 
  
      do i2 = ilower2-gcw_shift,iupper2+gcw_shift+1
         do i1 = ilower1-gcw_shift,iupper1+gcw_shift
            do i0 = ilower0-gcw_shift,iupper0+gcw_shift+1
               dmr = 1.d0/V(i0,i1,i2) + 1.d0/V(i0,i1,i2-1)  
     &              + 1.d0/V(i0-1,i1,i2) + 1.d0/V(i0-1,i1,i2-1)
               u1(i0,i1,i2) = nmr/dmr

            enddo
         enddo
      enddo

      do i2 = ilower2-gcw_shift,iupper2+gcw_shift
         do i1 = ilower1-gcw_shift,iupper1+gcw_shift+1
            do i0 = ilower0-gcw_shift,iupper0+gcw_shift+1
               dmr = 1.d0/V(i0,i1,i2) + 1.d0/V(i0,i1-1,i2)  
     &               + 1.d0/V(i0-1,i1,i2) + 1.d0/V(i0-1,i1-1,i2)
               u2(i0,i1,i2) = nmr/dmr

            enddo
         enddo
      enddo

c
      return
      end