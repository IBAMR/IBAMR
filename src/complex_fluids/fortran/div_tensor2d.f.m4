c Filename: div_tensor2d.f.m4
c Created on 23 Oct 2019 by Aaron Barrett
c
c Copyright (c) 2002-2019, Boyce Griffith
c All rights reserved.
c
c Redistribution and use in source and binary forms, with or without
c modification, are permitted provided that the following conditions are met:
c
c    * Redistributions of source code must retain the above copyright notice,
c      this list of conditions and the following disclaimer.
c
c    * Redistributions in binary form must reproduce the above copyright
c      notice, this list of conditions and the following disclaimer in the
c      documentation and/or other materials provided with the distribution.
c
c    * Neither the name of The University of North Carolina nor the names of
c      its contributors may be used to endorse or promote products derived from
c      this software without specific prior written permission.
c
c THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
c AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
c IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
c ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
c LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
c CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
c SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
c INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
c CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
c ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
c POSSIBILITY OF SUCH DAMAGE.

define(NDIM,2)dnl
define(REAL,`double precision')dnl
define(INTEGER,`integer')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim2d.i)dnl
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       Computes d = alpha div s
c
c       where d is vector valued side centered
c       and s is tensor valued cell centered
c       using centered differencing
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine div_tensor_c_to_s_2d(dx, d_data_0, d_data_1, d_gcw,
     &        s_data, s_gcw,  ilower0,
     &        iupper0, ilower1,  iupper1, alpha)
      implicit none
c     INPUTS
      INTEGER ilower0,  iupper0
      INTEGER iupper1,  ilower1
      INTEGER s_gcw,  d_gcw

      REAL alpha
c     RETURNS
      REAL d_data_0(SIDE2d0(ilower,iupper,d_gcw))
      REAL d_data_1(SIDE2d1(ilower,iupper,d_gcw))
c     TAU DATA
      REAL s_data(CELL2d(ilower,iupper,s_gcw),0:2)
      REAL dx(0:1)

      INTEGER i0, i1
      REAL scale0_x, scale0_y
      REAL scale1_x, scale1_y

      scale0_x = alpha/dx(0)
      scale0_y = alpha/(dx(1)*4.d0)
      scale1_y = alpha/dx(1)
      scale1_x = alpha/(dx(0)*4.d0)

      do i1 = ilower1, iupper1
        do i0 = ilower0, iupper0+1
          d_data_0(i0,i1) =
     &      scale0_x*(s_data(i0, i1,0) - s_data(i0-1, i1,0)) +
     &      scale0_y*(s_data(i0-1, i1+1,2)+s_data(i0, i1+1,2)
     &                -s_data(i0-1, i1-1,2)-s_data(i0, i1-1,2))
        enddo
      enddo
      do i1 = ilower1, (iupper1+1)
        do i0 = ilower0, iupper0
          d_data_1(i0,i1) =
     &      scale1_y*(s_data(i0, i1,1) - s_data(i0, i1-1,1)) +
     &      scale1_x*(s_data(i0+1, i1,2) + s_data(i0+1, i1-1,2)
     &                -s_data(i0-1, i1-1,2) - s_data(i0-1, i1,2))
        enddo
      enddo
      end subroutine
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       Computes d = alpha div s
c
c       where d is vector valued cell centered
c       and s is tensor valued cell centered
c       using centered differencing
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine div_tensor_c_to_c_2d(dx, d_data, d_gcw,
     &        s_data, s_gcw,  ilower0,
     &        iupper0, ilower1,  iupper1, alpha)
      implicit none
c     INPUTS
      INTEGER ilower0,  iupper0
      INTEGER iupper1,  ilower1
      INTEGER s_gcw,  d_gcw
      REAL alpha
c     RETURNS
      REAL d_data(CELL2d(ilower,iupper,d_gcw),0:1)
c     TAU DATA
      REAL s_data(CELL2d(ilower,iupper,s_gcw),0:2)
      REAL dx(0:1)
      INTEGER i0, i1
      REAL scale_x, scale_y

      scale_x = alpha/(2.d0*dx(0))
      scale_y = alpha/(2.d0*dx(1))

      do i1 = ilower1, iupper1
        do i0 = ilower0, iupper0
          d_data(i0,i1,0) =
     &     scale_x*(s_data(i0+1, i1,0) - s_data(i0-1, i1,0)) +
     &     scale_y*(s_data(i0, i1+1,2) - s_data(i0, i1-1,2))
          d_data(i0,i1,1) =
     &     scale_y*(s_data(i0, i1+1,1) - s_data(i0, i1-1,1)) +
     &     scale_x*(s_data(i0+1, i1,2) - s_data(i0-1, i1,2))
        enddo
      enddo
      end subroutine
