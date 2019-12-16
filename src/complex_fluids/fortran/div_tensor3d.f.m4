c ---------------------------------------------------------------------
c
c Copyright (c) 2019 - 2019 by the IBAMR developers
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c	Computes d = alpha div s
c
c	where d is vector valued side centered
c       and s is symmetric tensor valued cell centered
c       using centered differencing
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine div_tensor_c_to_s_3d(dx, d_data_0, d_data_1, d_data_2,
     &        d_gcw, s_data, s_gcw, ilower0, iupper0,
     &        ilower1,  iupper1, ilower2, iupper2, alpha)
      implicit none
c     INPUTS
      INTEGER ilower0,  iupper0
      INTEGER iupper1,  ilower1
      INTEGER iupper2, ilower2
      INTEGER s_gcw,  d_gcw

      REAL alpha
c     RETURNS
      REAL d_data_0(SIDE3d0(ilower,iupper,d_gcw))
      REAL d_data_1(SIDE3d1(ilower,iupper,d_gcw))
      REAL d_data_2(SIDE3d2(ilower,iupper,d_gcw))
c     TAU DATA
      REAL s_data(CELL3d(ilower,iupper,s_gcw),0:5)
      REAL dx(0:2)

      INTEGER i0, i1, i2
      REAL scale0_x, scale0_y, scale0_z
      REAL scale1_x, scale1_y, scale1_z
      REAL scale2_x, scale2_y, scale2_z

      scale0_x = alpha/dx(0)
      scale0_y = alpha/(dx(1)*4.d0)
      scale0_z = alpha/(dx(2)*4.d0)
      scale1_y = alpha/dx(1)
      scale1_x = alpha/(dx(0)*4.d0)
      scale1_z = alpha/(dx(2)*4.d0)
      scale2_z = alpha/dx(2)
      scale2_x = alpha/(dx(0)*4.d0)
      scale2_y = alpha/(dx(1)*4.d0)

      do i2 = ilower2, iupper2
        do i1 = ilower1, iupper1
          do i0 = ilower0, (iupper0+1)
            d_data_0(i0,i1,i2) =
     &       scale0_x*(s_data(i0, i1, i2, 0)-s_data(i0-1, i1, i2, 0))+
     &       scale0_y*(s_data(i0-1, i1+1, i2,5)+s_data(i0, i1+1, i2,5)
     &              -s_data(i0-1, i1-1, i2,5)-s_data(i0, i1-1, i2,5))+
     &       scale0_z*(s_data(i0-1, i1, i2+1,4)+s_data(i0, i1, i2+1,4)
     &              -s_data(i0-1,i1,i2-1,4) - s_data(i0,i1,i2-1,4))
          enddo
        enddo
      enddo
      do i2 = ilower2, iupper2
        do i1 = ilower1, (iupper1+1)
          do i0 = ilower0, iupper0
            d_data_1(i0,i1,i2) =
     &       scale1_y*(s_data(i0, i1, i2,1)-s_data(i0, i1-1, i2,1)) + 
     &       scale1_x*(s_data(i0+1, i1, i2,5)+s_data(i0+1, i1-1, i2,5)
     &              -s_data(i0-1, i1-1, i2,5)-s_data(i0-1, i1, i2,5))+
     &       scale1_z*(s_data(i0, i1,i2+1,3)+s_data(i0,i1-1,i2+1,3)
     &              -s_data(i0,i1,i2-1,3)-s_data(i0,i1-1,i2-1,3))
          enddo
        enddo
      enddo
      do i2 = ilower2, (iupper2+1)
        do i1 = ilower1, iupper1
          do i0 = ilower0, iupper0
            d_data_2(i0,i1,i2) = 
     &       scale2_x*(s_data(i0+1,i1,i2,4)+s_data(i0+1,i1,i2-1,4)
     &                 -s_data(i0-1,i1,i2,4)-s_data(i0-1,i1,i2-1,4))+
     &       scale2_y*(s_data(i0,i1+1,i2,3)+s_data(i0,i1+1,i2-1,3)
     &                 -s_data(i0,i1-1,i2,3)-s_data(i0,i1-1,i2-1,3))+
     &       scale2_z*(s_data(i0,i1,i2,2)-s_data(i0,i1,i2-1,2))
          enddo
        enddo
      enddo
      end subroutine
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       Computes d = alpha div s
c
c       where d is vector valued cell centered
c       and s is symmetric tensor valued cell centered
c       using centered differencing
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine div_tensor_c_to_c_3d(dx, d_data,
     &        d_gcw, s_data, s_gcw, ilower0, iupper0,
     &        ilower1,  iupper1, ilower2, iupper2, alpha)

      implicit none
c     INPUTS
      INTEGER ilower0,  iupper0
      INTEGER iupper1,  ilower1
      INTEGER iupper2, ilower2
      INTEGER s_gcw,  d_gcw

      REAL alpha
c     RETURNS
      REAL d_data(CELL3d(ilower,iupper,d_gcw),0:2)
c     TAU DATA
      REAL s_data(CELL3d(ilower,iupper,s_gcw),0:5)
      REAL dx(0:2)

      INTEGER i0, i1, i2
      REAL scale_x, scale_y, scale_z

      scale_x = alpha/(2.d0*dx(0))
      scale_y = alpha/(2.d0*dx(1))
      scale_z = alpha/(2.d0*dx(2))

      do i2 = ilower2, iupper2
        do i1 = ilower1, iupper1
          do i0 = ilower0, iupper0
            d_data(i0,i1,i2,0) =
     &        scale_x*(s_data(i0+1,i1,i2,0)-s_data(i0-1,i1,i2,0))+
     &        scale_y*(s_data(i0,i1+1,i2,5)-s_data(i0,i1-1,i2,5))+
     &        scale_z*(s_data(i0,i1,i2+1,4)-s_data(i0,i1,i2-1,4))
            d_data(i0,i1,i2,1) =
     &        scale_y*(s_data(i0,i1+1,i2,1)-s_data(i0,i1-1,i2,1))+
     &        scale_x*(s_data(i0+1,i1,i2,5)-s_data(i0-1,i1,i2,5))+
     &        scale_z*(s_data(i0,i1,i2+1,3)-s_data(i0,i1,i2-1,3))
            d_data(i0,i1,i2,2) =
     &        scale_x*(s_data(i0+1,i1,i2,4)-s_data(i0-1,i1,i2,4))+
     &        scale_y*(s_data(i0,i1+1,i2,3)-s_data(i0,i1-1,i2,3))+
     &        scale_z*(s_data(i0,i1,i2+1,2)-s_data(i0,i1,i2-1,2))
          enddo
        enddo
      enddo
      end subroutine
