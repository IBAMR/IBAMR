c ---------------------------------------------------------------------
c
c Copyright (c) 2011 - 2023 by the IBAMR developers
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
c     Compute the node centered vector field U from the face centered
c     normal vector field (v0,v1,v2) using simple averaging.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ftoninterp2nd3d(
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
      REAL U(NODE3d(ilower,iupper,U_gcw),0:NDIM-1)
c
c     Local variables.
c
      INTEGER i0,i1,i2
c
c     Compute the node centered vector field U from the face centered
c     vector field (v0,v1,v2).
c
      do i2 = ilower2,iupper2+1
         do i1 = ilower1,iupper1+1
            do i0 = ilower0,iupper0+1
               U(i0,i1,i2,0) = 0.25d0*(v0(i0,i1-1,i2-1)+
     &                                 v0(i0,i1  ,i2-1)+
     &                                 v0(i0,i1-1,i2  )+
     &                                 v0(i0,i1  ,i2  ))
            enddo
         enddo
      enddo
      do i2 = ilower2,iupper2+1
         do i1 = ilower1,iupper1+1
            do i0 = ilower0,iupper0+1
               U(i0,i1,i2,1) = 0.25d0*(v1(i1,i2-1,i0-1)+
     &                                 v1(i1,i2-1,i0  )+
     &                                 v1(i1,i2  ,i0-1)+
     &                                 v1(i1,i2  ,i0  ))
            enddo
         enddo
      enddo
      do i2 = ilower2,iupper2+1
         do i1 = ilower1,iupper1+1
            do i0 = ilower0,iupper0+1
               U(i0,i1,i2,2) = 0.25d0*(v2(i2,i0-1,i1-1)+
     &                                 v2(i2,i0  ,i1-1)+
     &                                 v2(i2,i0-1,i1  )+
     &                                 v2(i2,i0  ,i1  ))
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
c     Compute the node centered vector field U from the side centered
c     normal vector field (v0,v1,v2) using simple averaging.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine stoninterp2nd3d(
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
      REAL U(NODE3d(ilower,iupper,U_gcw),0:NDIM-1)
c
c     Local variables.
c
      INTEGER i0,i1,i2
c
c     Compute the node centered vector field U from the side centered
c     vector field (v0,v1,v2).
c
      do i2 = ilower2,iupper2+1
         do i1 = ilower1,iupper1+1
            do i0 = ilower0,iupper0+1
               U(i0,i1,i2,0) = 0.25d0*(v0(i0,i1-1,i2-1)+
     &                                 v0(i0,i1  ,i2-1)+
     &                                 v0(i0,i1-1,i2  )+
     &                                 v0(i0,i1  ,i2  ))
            enddo
         enddo
      enddo
      do i2 = ilower2,iupper2+1
         do i1 = ilower1,iupper1+1
            do i0 = ilower0,iupper0+1
               U(i0,i1,i2,1) = 0.25d0*(v1(i0-1,i1,i2-1)+
     &                                 v1(i0  ,i1,i2-1)+
     &                                 v1(i0-1,i1,i2  )+
     &                                 v1(i0  ,i1,i2  ))
            enddo
         enddo
      enddo
      do i2 = ilower2,iupper2+1
         do i1 = ilower1,iupper1+1
            do i0 = ilower0,iupper0+1
               U(i0,i1,i2,2) = 0.25d0*(v2(i0-1,i1-1,i2)+
     &                                 v2(i0  ,i1-1,i2)+
     &                                 v2(i0-1,i1  ,i2)+
     &                                 v2(i0  ,i1  ,i2))
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the node centered field U from the cell centered
c     field V using simple averaging.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ctoninterp3d(
     &     U,U_gcw,
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
      INTEGER U_gcw,V_gcw
      INTEGER U_ghost_interp

      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2

      REAL V(CELL3d(ilower,iupper,V_gcw))

c
c     Output.
c
      REAL U(NODE3d(ilower,iupper,U_gcw))
c
c     Local variables.
c
      INTEGER i0,i1,i2
      INTEGER gcw_shift
c
c     Compute the node centered scalar field U from the cell centered
c     scalar field V.
c
      gcw_shift = 0
      if (U_ghost_interp .eq. 1) then
         gcw_shift = U_gcw
      endif

      do i2 = ilower2-gcw_shift,iupper2+gcw_shift+1
         do i1 = ilower1-gcw_shift,iupper1+gcw_shift+1
            do i0 = ilower0-gcw_shift,iupper0+gcw_shift+1
               U(i0,i1,i2) = 0.125d0*(V(i0,i1,i2)+V(i0-1,i1,i2)
     &              +V(i0,i1-1,i2)+V(i0-1,i1-1,i2)+V(i0,i1,i2-1)
     &              +V(i0-1,i1,i2-1)+V(i0,i1-1,i2-1)+V(i0-1,i1-1,i2-1))
            enddo
         enddo
      enddo
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
c
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
c
c     Compute the edge centered interpolation of V
c
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
c
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
c     Function.
c
      REAL h_avg2
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
               u0(i0,i1,i2) = h_avg2(V(i0-1,i1,i2,0),V(i0,i1,i2,0))
            enddo
         enddo
      enddo
      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1+1
            do i0 = ilower0,iupper0
               u1(i0,i1,i2) = h_avg2(V(i0,i1-1,i2,1),V(i0,i1,i2,1))
            enddo
         enddo
      enddo
      do i2 = ilower2,iupper2+1
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
               u2(i0,i1,i2) = h_avg2(V(i0,i1,i2-1,2),V(i0,i1,i2,2))
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
c     cell centered scalar field V using harmonic averaging.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ctoscwiseharmonicinterp2nd3d(
     &     u0,u1,u2,u_gcw,
     &     V,V_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2)
c
      implicit none
c
c     Function.
c
      REAL h_avg2
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
               u0(i0,i1,i2) = h_avg2(V(i0-1,i1,i2),V(i0,i1,i2))
            enddo
         enddo
      enddo
      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1+1
            do i0 = ilower0,iupper0
               u1(i0,i1,i2) = h_avg2(V(i0,i1-1,i2),V(i0,i1,i2))
            enddo
         enddo
      enddo
      do i2 = ilower2,iupper2+1
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
               u2(i0,i1,i2) = h_avg2(V(i0,i1,i2-1),V(i0,i1,i2))
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
c     Function.
c
      REAL h_avg4
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
               u0(i0,i1,i2) = h_avg4(V(i0,i1-1,i2), V(i0,i1,i2-1),
     &                          V(i0,i1,i2), V(i0,i1-1,i2-1))

            enddo
         enddo
      enddo

      do i2 = ilower2-gcw_shift,iupper2+gcw_shift+1
         do i1 = ilower1-gcw_shift,iupper1+gcw_shift
            do i0 = ilower0-gcw_shift,iupper0+gcw_shift+1
               u1(i0,i1,i2) = h_avg4(V(i0,i1,i2),V(i0,i1,i2-1),
     &                          V(i0-1,i1,i2),V(i0-1,i1,i2-1))

            enddo
         enddo
      enddo

      do i2 = ilower2-gcw_shift,iupper2+gcw_shift
         do i1 = ilower1-gcw_shift,iupper1+gcw_shift+1
            do i0 = ilower0-gcw_shift,iupper0+gcw_shift+1
               u2(i0,i1,i2) = h_avg4(V(i0,i1,i2), V(i0,i1-1,i2),
     &                          V(i0-1,i1,i2),V(i0-1,i1-1,i2))

            enddo
         enddo
      enddo

c
      return
      end
