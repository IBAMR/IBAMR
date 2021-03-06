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
c     Computes F = alpha div grad U.
c
c     Uses the seven point stencil to compute the discrete Laplacian of
c     a variable U.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine laplace3d(
     &     F,F_gcw,
     &     alpha,
     &     U,U_gcw,
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
      INTEGER F_gcw,U_gcw

      REAL alpha

      REAL U(ilower0-U_gcw:iupper0+U_gcw,
     &       ilower1-U_gcw:iupper1+U_gcw,
     &       ilower2-U_gcw:iupper2+U_gcw)

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL F(ilower0-F_gcw:iupper0+F_gcw,
     &       ilower1-F_gcw:iupper1+F_gcw,
     &       ilower2-F_gcw:iupper2+F_gcw)
c
c     Local variables.
c
      INTEGER i0,i1,i2
      REAL    fac0,fac1,fac2
c
c     Compute the discrete Laplacian of U.
c
      fac0 = alpha/(dx(0)*dx(0))
      fac1 = alpha/(dx(1)*dx(1))
      fac2 = alpha/(dx(2)*dx(2))

      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
               F(i0,i1,i2) =
     &              fac0*(U(i0-1,i1,i2)+U(i0+1,i1,i2)-2.d0*U(i0,i1,i2))+
     &              fac1*(U(i0,i1-1,i2)+U(i0,i1+1,i2)-2.d0*U(i0,i1,i2))+
     &              fac2*(U(i0,i1,i2-1)+U(i0,i1,i2+1)-2.d0*U(i0,i1,i2))
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes F = alpha div grad U + beta V.
c
c     Uses the seven point stencil to compute the discrete Laplacian of
c     a variable U.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine laplaceadd3d(
     &     F,F_gcw,
     &     alpha,
     &     U,U_gcw,
     &     beta,
     &     V,V_gcw,
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
      INTEGER F_gcw,U_gcw,V_gcw

      REAL alpha

      REAL U(ilower0-U_gcw:iupper0+U_gcw,
     &       ilower1-U_gcw:iupper1+U_gcw,
     &       ilower2-U_gcw:iupper2+U_gcw)

      REAL beta

      REAL V(ilower0-V_gcw:iupper0+V_gcw,
     &       ilower1-V_gcw:iupper1+V_gcw,
     &       ilower2-V_gcw:iupper2+V_gcw)

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL F(ilower0-F_gcw:iupper0+F_gcw,
     &       ilower1-F_gcw:iupper1+F_gcw,
     &       ilower2-F_gcw:iupper2+F_gcw)
c
c     Local variables.
c
      INTEGER i0,i1,i2
      REAL    fac0,fac1,fac2
c
c     Compute the discrete Laplacian of U.
c
      fac0 = alpha/(dx(0)*dx(0))
      fac1 = alpha/(dx(1)*dx(1))
      fac2 = alpha/(dx(2)*dx(2))

      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
               F(i0,i1,i2) =
     &              fac0*(U(i0-1,i1,i2)+U(i0+1,i1,i2)-2.d0*U(i0,i1,i2))+
     &              fac1*(U(i0,i1-1,i2)+U(i0,i1+1,i2)-2.d0*U(i0,i1,i2))+
     &              fac2*(U(i0,i1,i2-1)+U(i0,i1,i2+1)-2.d0*U(i0,i1,i2))+
     &              beta* V(i0,i1,i2)
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes F = alpha div grad U + beta U.
c
c     Uses the seven point stencil to compute the damped discrete
c     Laplacian of a variable U.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine dampedlaplace3d(
     &     F,F_gcw,
     &     alpha,beta,
     &     U,U_gcw,
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
      INTEGER F_gcw,U_gcw

      REAL alpha,beta

      REAL U(ilower0-U_gcw:iupper0+U_gcw,
     &       ilower1-U_gcw:iupper1+U_gcw,
     &       ilower2-U_gcw:iupper2+U_gcw)

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL F(ilower0-F_gcw:iupper0+F_gcw,
     &       ilower1-F_gcw:iupper1+F_gcw,
     &       ilower2-F_gcw:iupper2+F_gcw)
c
c     Local variables.
c
      INTEGER i0,i1,i2
      REAL    fac0,fac1,fac2
c
c     Compute the discrete Laplacian of U.
c
      fac0 = alpha/(dx(0)*dx(0))
      fac1 = alpha/(dx(1)*dx(1))
      fac2 = alpha/(dx(2)*dx(2))

      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
               F(i0,i1,i2) =
     &              fac0*(U(i0-1,i1,i2)+U(i0+1,i1,i2)-2.d0*U(i0,i1,i2))+
     &              fac1*(U(i0,i1-1,i2)+U(i0,i1+1,i2)-2.d0*U(i0,i1,i2))+
     &              fac2*(U(i0,i1,i2-1)+U(i0,i1,i2+1)-2.d0*U(i0,i1,i2))+
     &              beta* U(i0,i1,i2)
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes F = alpha div grad U + beta U + gamma V.
c
c     Uses the seven point stencil to compute the damped discrete
c     Laplacian of a variable U.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine dampedlaplaceadd3d(
     &     F,F_gcw,
     &     alpha,beta,
     &     U,U_gcw,
     &     gamma,
     &     V,V_gcw,
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
      INTEGER F_gcw,U_gcw,V_gcw

      REAL alpha,beta

      REAL U(ilower0-U_gcw:iupper0+U_gcw,
     &       ilower1-U_gcw:iupper1+U_gcw,
     &       ilower2-U_gcw:iupper2+U_gcw)

      REAL gamma

      REAL V(ilower0-V_gcw:iupper0+V_gcw,
     &       ilower1-V_gcw:iupper1+V_gcw,
     &       ilower2-V_gcw:iupper2+V_gcw)

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL F(ilower0-F_gcw:iupper0+F_gcw,
     &       ilower1-F_gcw:iupper1+F_gcw,
     &       ilower2-F_gcw:iupper2+F_gcw)
c
c     Local variables.
c
      INTEGER i0,i1,i2
      REAL    fac0,fac1,fac2
c
c     Compute the discrete Laplacian of U.
c
      fac0 = alpha/(dx(0)*dx(0))
      fac1 = alpha/(dx(1)*dx(1))
      fac2 = alpha/(dx(2)*dx(2))

      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
               F(i0,i1,i2) =
     &              fac0*(U(i0-1,i1,i2)+U(i0+1,i1,i2)-2.d0*U(i0,i1,i2))+
     &              fac1*(U(i0,i1-1,i2)+U(i0,i1+1,i2)-2.d0*U(i0,i1,i2))+
     &              fac2*(U(i0,i1,i2-1)+U(i0,i1,i2+1)-2.d0*U(i0,i1,i2))+
     &              beta* U(i0,i1,i2)                                  +
     &              gamma*V(i0,i1,i2)
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes F = div alpha grad U.
c
c     Uses a seven point stencil to compute the cell centered
c     anisotropic Laplacian of a cell centered variable U, where the
c     grid aligned anisotropic diffusion coefficient alpha is specified
c     on cell faces.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ctocanisoflaplace3d(
     &     F,F_gcw,
     &     alpha0,alpha1,alpha2,alpha_gcw,
     &     U,U_gcw,
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
      INTEGER F_gcw,alpha_gcw,U_gcw

      REAL alpha0(FACE3d0(ilower,iupper,alpha_gcw))
      REAL alpha1(FACE3d1(ilower,iupper,alpha_gcw))
      REAL alpha2(FACE3d2(ilower,iupper,alpha_gcw))

      REAL U(CELL3d(ilower,iupper,U_gcw))

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL F(CELL3d(ilower,iupper,F_gcw))
c
c     Local variables.
c
      INTEGER i0,i1,i2,k
      REAL    f0(0:1),f1(0:1),f2(0:1)
c
c     Compute the cell centered anisotropic Laplacian of U.
c
      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
               do k = 0,1
                  f0(k) = alpha0(i0+k,i1,i2)*
     &                 (U(i0+k,i1,i2)-U(i0-1+k,i1,i2))/dx(0)
                  f1(k) = alpha1(i1+k,i2,i0)*
     &                 (U(i0,i1+k,i2)-U(i0,i1-1+k,i2))/dx(1)
                  f2(k) = alpha2(i2+k,i0,i1)*
     &                 (U(i0,i1,i2+k)-U(i0,i1,i2-1+k))/dx(2)
               enddo

               F(i0,i1,i2) =
     &              (f0(1)-f0(0))/dx(0) +
     &              (f1(1)-f1(0))/dx(1) +
     &              (f2(1)-f2(0))/dx(2)
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes F = div alpha grad U + beta V.
c
c     Uses a seven point stencil to compute the cell centered
c     anisotropic Laplacian of a cell centered variable U, where the
c     grid aligned anisotropic diffusion coefficient alpha is specified
c     on cell faces.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ctocanisoflaplaceadd3d(
     &     F,F_gcw,
     &     alpha0,alpha1,alpha2,alpha_gcw,
     &     U,U_gcw,
     &     beta,
     &     V,V_gcw,
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
      INTEGER F_gcw,alpha_gcw,U_gcw,V_gcw

      REAL alpha0(FACE3d0(ilower,iupper,alpha_gcw))
      REAL alpha1(FACE3d1(ilower,iupper,alpha_gcw))
      REAL alpha2(FACE3d2(ilower,iupper,alpha_gcw))

      REAL U(CELL3d(ilower,iupper,U_gcw))

      REAL beta

      REAL V(CELL3d(ilower,iupper,V_gcw))

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL F(CELL3d(ilower,iupper,F_gcw))
c
c     Local variables.
c
      INTEGER i0,i1,i2,k
      REAL    f0(0:1),f1(0:1),f2(0:1)
c
c     Compute the cell centered anisotropic Laplacian of U.
c
      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
               do k = 0,1
                  f0(k) = alpha0(i0+k,i1,i2)*
     &                 (U(i0+k,i1,i2)-U(i0-1+k,i1,i2))/dx(0)
                  f1(k) = alpha1(i1+k,i2,i0)*
     &                 (U(i0,i1+k,i2)-U(i0,i1-1+k,i2))/dx(1)
                  f2(k) = alpha2(i2+k,i0,i1)*
     &                 (U(i0,i1,i2+k)-U(i0,i1,i2-1+k))/dx(2)
               enddo

               F(i0,i1,i2) = beta*V(i0,i1,i2) +
     &              (f0(1)-f0(0))/dx(0) +
     &              (f1(1)-f1(0))/dx(1) +
     &              (f2(1)-f2(0))/dx(2)
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes F = div alpha grad U + beta U.
c
c     Uses a seven point stencil to compute the damped cell centered
c     anisotropic Laplacian of a cell centered variable U, where the
c     grid aligned anisotropic diffusion coefficient alpha is specified
c     on cell faces.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ctocanisofdampedlaplace3d(
     &     F,F_gcw,
     &     alpha0,alpha1,alpha2,alpha_gcw,
     &     beta,
     &     U,U_gcw,
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
      INTEGER F_gcw,alpha_gcw,U_gcw

      REAL alpha0(FACE3d0(ilower,iupper,alpha_gcw))
      REAL alpha1(FACE3d1(ilower,iupper,alpha_gcw))
      REAL alpha2(FACE3d2(ilower,iupper,alpha_gcw))

      REAL beta

      REAL U(CELL3d(ilower,iupper,U_gcw))

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL F(CELL3d(ilower,iupper,F_gcw))
c
c     Local variables.
c
      INTEGER i0,i1,i2,k
      REAL    f0(0:1),f1(0:1),f2(0:1)
c
c     Compute the cell centered anisotropic Laplacian of U.
c
      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
               do k = 0,1
                  f0(k) = alpha0(i0+k,i1,i2)*
     &                 (U(i0+k,i1,i2)-U(i0-1+k,i1,i2))/dx(0)
                  f1(k) = alpha1(i1+k,i2,i0)*
     &                 (U(i0,i1+k,i2)-U(i0,i1-1+k,i2))/dx(1)
                  f2(k) = alpha2(i2+k,i0,i1)*
     &                 (U(i0,i1,i2+k)-U(i0,i1,i2-1+k))/dx(2)
               enddo

               F(i0,i1,i2) = beta*U(i0,i1,i2) +
     &              (f0(1)-f0(0))/dx(0) +
     &              (f1(1)-f1(0))/dx(1) +
     &              (f2(1)-f2(0))/dx(2)
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes F = div alpha grad U + beta U + gamma V.
c
c     Uses a seven point stencil to compute the damped cell centered
c     anisotropic Laplacian of a cell centered variable U, where the
c     grid aligned anisotropic diffusion coefficient alpha is specified
c     on cell faces.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ctocanisofdampedlaplaceadd3d(
     &     F,F_gcw,
     &     alpha0,alpha1,alpha2,alpha_gcw,
     &     beta,
     &     U,U_gcw,
     &     gamma,
     &     V,V_gcw,
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
      INTEGER F_gcw,alpha_gcw,U_gcw,V_gcw

      REAL alpha0(FACE3d0(ilower,iupper,alpha_gcw))
      REAL alpha1(FACE3d1(ilower,iupper,alpha_gcw))
      REAL alpha2(FACE3d2(ilower,iupper,alpha_gcw))

      REAL beta

      REAL U(CELL3d(ilower,iupper,U_gcw))

      REAL gamma

      REAL V(CELL3d(ilower,iupper,V_gcw))

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL F(CELL3d(ilower,iupper,F_gcw))
c
c     Local variables.
c
      INTEGER i0,i1,i2,k
      REAL    f0(0:1),f1(0:1),f2(0:1)
c
c     Compute the cell centered anisotropic Laplacian of U.
c
      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
               do k = 0,1
                  f0(k) = alpha0(i0+k,i1,i2)*
     &                 (U(i0+k,i1,i2)-U(i0-1+k,i1,i2))/dx(0)
                  f1(k) = alpha1(i1+k,i2,i0)*
     &                 (U(i0,i1+k,i2)-U(i0,i1-1+k,i2))/dx(1)
                  f2(k) = alpha2(i2+k,i0,i1)*
     &                 (U(i0,i1,i2+k)-U(i0,i1,i2-1+k))/dx(2)
               enddo

               F(i0,i1,i2) = beta*U(i0,i1,i2) + gamma*V(i0,i1,i2) +
     &              (f0(1)-f0(0))/dx(0) +
     &              (f1(1)-f1(0))/dx(1) +
     &              (f2(1)-f2(0))/dx(2)
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes F = div alpha grad U.
c
c     Uses a seven point stencil to compute the cell centered
c     anisotropic Laplacian of a cell centered variable U, where the
c     grid aligned anisotropic diffusion coefficient alpha is specified
c     on cell sides.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ctocanisoslaplace3d(
     &     F,F_gcw,
     &     alpha0,alpha1,alpha2,alpha_gcw,
     &     U,U_gcw,
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
      INTEGER F_gcw,alpha_gcw,U_gcw

      REAL alpha0(SIDE3d0(ilower,iupper,alpha_gcw))
      REAL alpha1(SIDE3d1(ilower,iupper,alpha_gcw))
      REAL alpha2(SIDE3d2(ilower,iupper,alpha_gcw))

      REAL U(CELL3d(ilower,iupper,U_gcw))

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL F(CELL3d(ilower,iupper,F_gcw))
c
c     Local variables.
c
      INTEGER i0,i1,i2,k
      REAL    f0(0:1),f1(0:1),f2(0:1)
c
c     Compute the cell centered anisotropic Laplacian of U.
c
      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
               do k = 0,1
                  f0(k) = alpha0(i0+k,i1,i2)*
     &                 (U(i0+k,i1,i2)-U(i0-1+k,i1,i2))/dx(0)
                  f1(k) = alpha1(i0,i1+k,i2)*
     &                 (U(i0,i1+k,i2)-U(i0,i1-1+k,i2))/dx(1)
                  f2(k) = alpha2(i0,i1,i2+k)*
     &                 (U(i0,i1,i2+k)-U(i0,i1,i2-1+k))/dx(2)
               enddo

               F(i0,i1,i2) =
     &              (f0(1)-f0(0))/dx(0) +
     &              (f1(1)-f1(0))/dx(1) +
     &              (f2(1)-f2(0))/dx(2)
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes F = div alpha grad U + beta V.
c
c     Uses a seven point stencil to compute the cell centered
c     anisotropic Laplacian of a cell centered variable U, where the
c     grid aligned anisotropic diffusion coefficient alpha is specified
c     on cell sides.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ctocanisoslaplaceadd3d(
     &     F,F_gcw,
     &     alpha0,alpha1,alpha2,alpha_gcw,
     &     U,U_gcw,
     &     beta,
     &     V,V_gcw,
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
      INTEGER F_gcw,alpha_gcw,U_gcw,V_gcw

      REAL alpha0(SIDE3d0(ilower,iupper,alpha_gcw))
      REAL alpha1(SIDE3d1(ilower,iupper,alpha_gcw))
      REAL alpha2(SIDE3d2(ilower,iupper,alpha_gcw))

      REAL U(CELL3d(ilower,iupper,U_gcw))

      REAL beta

      REAL V(CELL3d(ilower,iupper,V_gcw))

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL F(CELL3d(ilower,iupper,F_gcw))
c
c     Local variables.
c
      INTEGER i0,i1,i2,k
      REAL    f0(0:1),f1(0:1),f2(0:1)
c
c     Compute the cell centered anisotropic Laplacian of U.
c
      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
               do k = 0,1
                  f0(k) = alpha0(i0+k,i1,i2)*
     &                 (U(i0+k,i1,i2)-U(i0-1+k,i1,i2))/dx(0)
                  f1(k) = alpha1(i0,i1+k,i2)*
     &                 (U(i0,i1+k,i2)-U(i0,i1-1+k,i2))/dx(1)
                  f2(k) = alpha2(i0,i1,i2+k)*
     &                 (U(i0,i1,i2+k)-U(i0,i1,i2-1+k))/dx(2)
               enddo

               F(i0,i1,i2) = beta*V(i0,i1,i2) +
     &              (f0(1)-f0(0))/dx(0) +
     &              (f1(1)-f1(0))/dx(1) +
     &              (f2(1)-f2(0))/dx(2)
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes F = div alpha grad U + beta U.
c
c     Uses a seven point stencil to compute the damped cell centered
c     anisotropic Laplacian of a cell centered variable U, where the
c     grid aligned anisotropic diffusion coefficient alpha is specified
c     on cell sides.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ctocanisosdampedlaplace3d(
     &     F,F_gcw,
     &     alpha0,alpha1,alpha2,alpha_gcw,
     &     beta,
     &     U,U_gcw,
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
      INTEGER F_gcw,alpha_gcw,U_gcw

      REAL alpha0(SIDE3d0(ilower,iupper,alpha_gcw))
      REAL alpha1(SIDE3d1(ilower,iupper,alpha_gcw))
      REAL alpha2(SIDE3d2(ilower,iupper,alpha_gcw))

      REAL beta

      REAL U(CELL3d(ilower,iupper,U_gcw))

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL F(CELL3d(ilower,iupper,F_gcw))
c
c     Local variables.
c
      INTEGER i0,i1,i2,k
      REAL    f0(0:1),f1(0:1),f2(0:1)
c
c     Compute the cell centered anisotropic Laplacian of U.
c
      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
               do k = 0,1
                  f0(k) = alpha0(i0+k,i1,i2)*
     &                 (U(i0+k,i1,i2)-U(i0-1+k,i1,i2))/dx(0)
                  f1(k) = alpha1(i0,i1+k,i2)*
     &                 (U(i0,i1+k,i2)-U(i0,i1-1+k,i2))/dx(1)
                  f2(k) = alpha2(i0,i1,i2+k)*
     &                 (U(i0,i1,i2+k)-U(i0,i1,i2-1+k))/dx(2)
               enddo

               F(i0,i1,i2) = beta*U(i0,i1,i2) +
     &              (f0(1)-f0(0))/dx(0) +
     &              (f1(1)-f1(0))/dx(1) +
     &              (f2(1)-f2(0))/dx(2)
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes F = div alpha grad U + beta U + gamma V.
c
c     Uses a seven point stencil to compute the damped cell centered
c     anisotropic Laplacian of a cell centered variable U, where the
c     grid aligned anisotropic diffusion coefficient alpha is specified
c     on cell sides.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ctocanisosdampedlaplaceadd3d(
     &     F,F_gcw,
     &     alpha0,alpha1,alpha2,alpha_gcw,
     &     beta,
     &     U,U_gcw,
     &     gamma,
     &     V,V_gcw,
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
      INTEGER F_gcw,alpha_gcw,U_gcw,V_gcw

      REAL alpha0(SIDE3d0(ilower,iupper,alpha_gcw))
      REAL alpha1(SIDE3d1(ilower,iupper,alpha_gcw))
      REAL alpha2(SIDE3d2(ilower,iupper,alpha_gcw))

      REAL beta

      REAL U(CELL3d(ilower,iupper,U_gcw))

      REAL gamma

      REAL V(CELL3d(ilower,iupper,V_gcw))

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL F(CELL3d(ilower,iupper,F_gcw))
c
c     Local variables.
c
      INTEGER i0,i1,i2,k
      REAL    f0(0:1),f1(0:1),f2(0:1)
c
c     Compute the cell centered anisotropic Laplacian of U.
c
      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
               do k = 0,1
                  f0(k) = alpha0(i0+k,i1,i2)*
     &                 (U(i0+k,i1,i2)-U(i0-1+k,i1,i2))/dx(0)
                  f1(k) = alpha1(i0,i1+k,i2)*
     &                 (U(i0,i1+k,i2)-U(i0,i1-1+k,i2))/dx(1)
                  f2(k) = alpha2(i0,i1,i2+k)*
     &                 (U(i0,i1,i2+k)-U(i0,i1,i2-1+k))/dx(2)
               enddo

               F(i0,i1,i2) = beta*U(i0,i1,i2) + gamma*V(i0,i1,i2) +
     &              (f0(1)-f0(0))/dx(0) +
     &              (f1(1)-f1(0))/dx(1) +
     &              (f2(1)-f2(0))/dx(2)
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes F = div alpha grad U.
c
c     Uses a nineteen point stencil to compute the cell centered
c     anisotropic Laplacian of a cell centered variable U, where the
c     general (non-grid aligned) anisotropic diffusion coefficient alpha
c     is specified on cell faces.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ctocgenanisoflaplace3d(
     &     F,F_gcw,
     &     alpha0,alpha1,alpha2,alpha_gcw,
     &     U,U_gcw,
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
      INTEGER F_gcw,alpha_gcw,U_gcw

      REAL alpha0(FACE3d0(ilower,iupper,alpha_gcw),0:NDIM-1)
      REAL alpha1(FACE3d1(ilower,iupper,alpha_gcw),0:NDIM-1)
      REAL alpha2(FACE3d2(ilower,iupper,alpha_gcw),0:NDIM-1)

      REAL U(CELL3d(ilower,iupper,U_gcw))

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL F(CELL3d(ilower,iupper,F_gcw))
c
c     Local variables.
c
      INTEGER i0,i1,i2,k
      REAL    f0(0:1),f1(0:1),f2(0:1)
c
c     Compute the cell centered anisotropic Laplacian of U.
c
      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
               do k = 0,1
                  f0(k) =
     &                 alpha0(i0+k,i1,i2,0)*
     &                 ( U(i0+k,i1,i2)-U(i0-1+k,i1,i2) )      /dx(0) +
     &
     &                 alpha0(i0+k,i1,i2,1)*0.125d0*
     &                 ( U(i0  +k,i1+1,i2)-U(i0  +k,i1-1,i2)
     &                 + U(i0-1+k,i1+1,i2)-U(i0-1+k,i1-1,i2) )/dx(1) +
     &
     &                 alpha0(i0+k,i1,i2,2)*0.125d0*
     &                 ( U(i0  +k,i1,i2+1)-U(i0  +k,i1,i2-1)
     &                 + U(i0-1+k,i1,i2+1)-U(i0-1+k,i1,i2-1) )/dx(2)

                  f1(k) =
     &                 alpha1(i1+k,i2,i0,0)*0.125d0*
     &                 ( U(i0+1,i1  +k,i2)-U(i0-1,i1  +k,i2)
     &                 + U(i0+1,i1-1+k,i2)-U(i0-1,i1-1+k,i2) )/dx(0) +
     &
     &                 alpha1(i1+k,i2,i0,1)*
     &                 ( U(i0,i1+k,i2)-U(i0,i1-1+k,i2) )      /dx(1) +
     &
     &                 alpha1(i1+k,i2,i0,2)*0.125d0*
     &                 ( U(i0,i1  +k,i2+1)-U(i0,i1  +k,i2-1)
     &                 + U(i0,i1-1+k,i2+1)-U(i0,i1-1+k,i2-1) )/dx(2)

                  f2(k) =
     &                 alpha2(i2+k,i0,i1,0)*0.125d0*
     &                 ( U(i0+1,i1,i2  +k)-U(i0-1,i1,i2  +k)
     &                 + U(i0+1,i1,i2-1+k)-U(i0-1,i1,i2-1+k) )/dx(0) +
     &
     &                 alpha2(i2+k,i0,i1,1)*0.125d0*
     &                 ( U(i0,i1+1,i2  +k)-U(i0,i1-1,i2  +k)
     &                 + U(i0,i1+1,i2-1+k)-U(i0,i1-1,i2-1+k) )/dx(1) +
     &
     &                 alpha2(i2+k,i0,i1,2)*
     &                 ( U(i0,i1,i2+k)-U(i0,i1,i2-1+k) )      /dx(2)
               enddo

               F(i0,i1,i2) =
     &              (f0(1)-f0(0))/dx(0) +
     &              (f1(1)-f1(0))/dx(1) +
     &              (f2(1)-f2(0))/dx(2)
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes F = div alpha grad U + beta V.
c
c     Uses a nineteen point stencil to compute the cell centered
c     anisotropic Laplacian of a cell centered variable U, where the
c     general (non-grid aligned) anisotropic diffusion coefficient alpha
c     is specified on cell faces.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ctocgenanisoflaplaceadd3d(
     &     F,F_gcw,
     &     alpha0,alpha1,alpha2,alpha_gcw,
     &     U,U_gcw,
     &     beta,
     &     V,V_gcw,
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
      INTEGER F_gcw,alpha_gcw,U_gcw,V_gcw

      REAL alpha0(FACE3d0(ilower,iupper,alpha_gcw),0:NDIM-1)
      REAL alpha1(FACE3d1(ilower,iupper,alpha_gcw),0:NDIM-1)
      REAL alpha2(FACE3d2(ilower,iupper,alpha_gcw),0:NDIM-1)

      REAL U(CELL3d(ilower,iupper,U_gcw))

      REAL beta

      REAL V(CELL3d(ilower,iupper,V_gcw))

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL F(CELL3d(ilower,iupper,F_gcw))
c
c     Local variables.
c
      INTEGER i0,i1,i2,k
      REAL    f0(0:1),f1(0:1),f2(0:1)
c
c     Compute the cell centered anisotropic Laplacian of U.
c
      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
               do k = 0,1
                  f0(k) =
     &                 alpha0(i0+k,i1,i2,0)*
     &                 ( U(i0+k,i1,i2)-U(i0-1+k,i1,i2) )      /dx(0) +
     &
     &                 alpha0(i0+k,i1,i2,1)*0.125d0*
     &                 ( U(i0  +k,i1+1,i2)-U(i0  +k,i1-1,i2)
     &                 + U(i0-1+k,i1+1,i2)-U(i0-1+k,i1-1,i2) )/dx(1) +
     &
     &                 alpha0(i0+k,i1,i2,2)*0.125d0*
     &                 ( U(i0  +k,i1,i2+1)-U(i0  +k,i1,i2-1)
     &                 + U(i0-1+k,i1,i2+1)-U(i0-1+k,i1,i2-1) )/dx(2)

                  f1(k) =
     &                 alpha1(i1+k,i2,i0,0)*0.125d0*
     &                 ( U(i0+1,i1  +k,i2)-U(i0-1,i1  +k,i2)
     &                 + U(i0+1,i1-1+k,i2)-U(i0-1,i1-1+k,i2) )/dx(0) +
     &
     &                 alpha1(i1+k,i2,i0,1)*
     &                 ( U(i0,i1+k,i2)-U(i0,i1-1+k,i2) )      /dx(1) +
     &
     &                 alpha1(i1+k,i2,i0,2)*0.125d0*
     &                 ( U(i0,i1  +k,i2+1)-U(i0,i1  +k,i2-1)
     &                 + U(i0,i1-1+k,i2+1)-U(i0,i1-1+k,i2-1) )/dx(2)

                  f2(k) =
     &                 alpha2(i2+k,i0,i1,0)*0.125d0*
     &                 ( U(i0+1,i1,i2  +k)-U(i0-1,i1,i2  +k)
     &                 + U(i0+1,i1,i2-1+k)-U(i0-1,i1,i2-1+k) )/dx(0) +
     &
     &                 alpha2(i2+k,i0,i1,1)*0.125d0*
     &                 ( U(i0,i1+1,i2  +k)-U(i0,i1-1,i2  +k)
     &                 + U(i0,i1+1,i2-1+k)-U(i0,i1-1,i2-1+k) )/dx(1) +
     &
     &                 alpha2(i2+k,i0,i1,2)*
     &                 ( U(i0,i1,i2+k)-U(i0,i1,i2-1+k) )      /dx(2)
               enddo

               F(i0,i1,i2) = beta*V(i0,i1,i2) +
     &              (f0(1)-f0(0))/dx(0) +
     &              (f1(1)-f1(0))/dx(1) +
     &              (f2(1)-f2(0))/dx(2)
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes F = div alpha grad U + beta U.
c
c     Uses a nineteen point stencil to compute the damped cell centered
c     anisotropic Laplacian of a cell centered variable U, where the
c     general (non-grid aligned) anisotropic diffusion coefficient alpha
c     is specified on cell faces.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ctocgenanisofdampedlaplace3d(
     &     F,F_gcw,
     &     alpha0,alpha1,alpha2,alpha_gcw,
     &     beta,
     &     U,U_gcw,
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
      INTEGER F_gcw,alpha_gcw,U_gcw

      REAL alpha0(FACE3d0(ilower,iupper,alpha_gcw),0:NDIM-1)
      REAL alpha1(FACE3d1(ilower,iupper,alpha_gcw),0:NDIM-1)
      REAL alpha2(FACE3d2(ilower,iupper,alpha_gcw),0:NDIM-1)

      REAL beta

      REAL U(CELL3d(ilower,iupper,U_gcw))

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL F(CELL3d(ilower,iupper,F_gcw))
c
c     Local variables.
c
      INTEGER i0,i1,i2,k
      REAL    f0(0:1),f1(0:1),f2(0:1)
c
c     Compute the cell centered anisotropic Laplacian of U.
c
      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
               do k = 0,1
                  f0(k) =
     &                 alpha0(i0+k,i1,i2,0)*
     &                 ( U(i0+k,i1,i2)-U(i0-1+k,i1,i2) )      /dx(0) +
     &
     &                 alpha0(i0+k,i1,i2,1)*0.125d0*
     &                 ( U(i0  +k,i1+1,i2)-U(i0  +k,i1-1,i2)
     &                 + U(i0-1+k,i1+1,i2)-U(i0-1+k,i1-1,i2) )/dx(1) +
     &
     &                 alpha0(i0+k,i1,i2,2)*0.125d0*
     &                 ( U(i0  +k,i1,i2+1)-U(i0  +k,i1,i2-1)
     &                 + U(i0-1+k,i1,i2+1)-U(i0-1+k,i1,i2-1) )/dx(2)

                  f1(k) =
     &                 alpha1(i1+k,i2,i0,0)*0.125d0*
     &                 ( U(i0+1,i1  +k,i2)-U(i0-1,i1  +k,i2)
     &                 + U(i0+1,i1-1+k,i2)-U(i0-1,i1-1+k,i2) )/dx(0) +
     &
     &                 alpha1(i1+k,i2,i0,1)*
     &                 ( U(i0,i1+k,i2)-U(i0,i1-1+k,i2) )      /dx(1) +
     &
     &                 alpha1(i1+k,i2,i0,2)*0.125d0*
     &                 ( U(i0,i1  +k,i2+1)-U(i0,i1  +k,i2-1)
     &                 + U(i0,i1-1+k,i2+1)-U(i0,i1-1+k,i2-1) )/dx(2)

                  f2(k) =
     &                 alpha2(i2+k,i0,i1,0)*0.125d0*
     &                 ( U(i0+1,i1,i2  +k)-U(i0-1,i1,i2  +k)
     &                 + U(i0+1,i1,i2-1+k)-U(i0-1,i1,i2-1+k) )/dx(0) +
     &
     &                 alpha2(i2+k,i0,i1,1)*0.125d0*
     &                 ( U(i0,i1+1,i2  +k)-U(i0,i1-1,i2  +k)
     &                 + U(i0,i1+1,i2-1+k)-U(i0,i1-1,i2-1+k) )/dx(1) +
     &
     &                 alpha2(i2+k,i0,i1,2)*
     &                 ( U(i0,i1,i2+k)-U(i0,i1,i2-1+k) )      /dx(2)
               enddo

               F(i0,i1,i2) = beta*U(i0,i1,i2) +
     &              (f0(1)-f0(0))/dx(0) +
     &              (f1(1)-f1(0))/dx(1) +
     &              (f2(1)-f2(0))/dx(2)
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes F = div alpha grad U + beta U + gamma V.
c
c     Uses a nineteen point stencil to compute the damped cell centered
c     anisotropic Laplacian of a cell centered variable U, where the
c     general (non-grid aligned) anisotropic diffusion coefficient alpha
c     is specified on cell faces.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ctocgenanisofdampedlaplaceadd3d(
     &     F,F_gcw,
     &     alpha0,alpha1,alpha2,alpha_gcw,
     &     beta,
     &     U,U_gcw,
     &     gamma,
     &     V,V_gcw,
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
      INTEGER F_gcw,alpha_gcw,U_gcw,V_gcw

      REAL alpha0(FACE3d0(ilower,iupper,alpha_gcw),0:NDIM-1)
      REAL alpha1(FACE3d1(ilower,iupper,alpha_gcw),0:NDIM-1)
      REAL alpha2(FACE3d2(ilower,iupper,alpha_gcw),0:NDIM-1)

      REAL beta

      REAL U(CELL3d(ilower,iupper,U_gcw))

      REAL gamma

      REAL V(CELL3d(ilower,iupper,V_gcw))

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL F(CELL3d(ilower,iupper,F_gcw))
c
c     Local variables.
c
      INTEGER i0,i1,i2,k
      REAL    f0(0:1),f1(0:1),f2(0:1)
c
c     Compute the cell centered anisotropic Laplacian of U.
c
      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
               do k = 0,1
                  f0(k) =
     &                 alpha0(i0+k,i1,i2,0)*
     &                 ( U(i0+k,i1,i2)-U(i0-1+k,i1,i2) )      /dx(0) +
     &
     &                 alpha0(i0+k,i1,i2,1)*0.125d0*
     &                 ( U(i0  +k,i1+1,i2)-U(i0  +k,i1-1,i2)
     &                 + U(i0-1+k,i1+1,i2)-U(i0-1+k,i1-1,i2) )/dx(1) +
     &
     &                 alpha0(i0+k,i1,i2,2)*0.125d0*
     &                 ( U(i0  +k,i1,i2+1)-U(i0  +k,i1,i2-1)
     &                 + U(i0-1+k,i1,i2+1)-U(i0-1+k,i1,i2-1) )/dx(2)

                  f1(k) =
     &                 alpha1(i1+k,i2,i0,0)*0.125d0*
     &                 ( U(i0+1,i1  +k,i2)-U(i0-1,i1  +k,i2)
     &                 + U(i0+1,i1-1+k,i2)-U(i0-1,i1-1+k,i2) )/dx(0) +
     &
     &                 alpha1(i1+k,i2,i0,1)*
     &                 ( U(i0,i1+k,i2)-U(i0,i1-1+k,i2) )      /dx(1) +
     &
     &                 alpha1(i1+k,i2,i0,2)*0.125d0*
     &                 ( U(i0,i1  +k,i2+1)-U(i0,i1  +k,i2-1)
     &                 + U(i0,i1-1+k,i2+1)-U(i0,i1-1+k,i2-1) )/dx(2)

                  f2(k) =
     &                 alpha2(i2+k,i0,i1,0)*0.125d0*
     &                 ( U(i0+1,i1,i2  +k)-U(i0-1,i1,i2  +k)
     &                 + U(i0+1,i1,i2-1+k)-U(i0-1,i1,i2-1+k) )/dx(0) +
     &
     &                 alpha2(i2+k,i0,i1,1)*0.125d0*
     &                 ( U(i0,i1+1,i2  +k)-U(i0,i1-1,i2  +k)
     &                 + U(i0,i1+1,i2-1+k)-U(i0,i1-1,i2-1+k) )/dx(1) +
     &
     &                 alpha2(i2+k,i0,i1,2)*
     &                 ( U(i0,i1,i2+k)-U(i0,i1,i2-1+k) )      /dx(2)
               enddo

               F(i0,i1,i2) = beta*U(i0,i1,i2) + gamma*V(i0,i1,i2) +
     &              (f0(1)-f0(0))/dx(0) +
     &              (f1(1)-f1(0))/dx(1) +
     &              (f2(1)-f2(0))/dx(2)
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes F = div alpha grad U.
c
c     Uses a nineteen point stencil to compute the cell centered
c     anisotropic Laplacian of a cell centered variable U, where the
c     general (non-grid aligned) anisotropic diffusion coefficient alpha
c     is specified on cell sides.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ctocgenanisoslaplace3d(
     &     F,F_gcw,
     &     alpha0,alpha1,alpha2,alpha_gcw,
     &     U,U_gcw,
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
      INTEGER F_gcw,alpha_gcw,U_gcw

      REAL alpha0(SIDE3d0(ilower,iupper,alpha_gcw),0:NDIM-1)
      REAL alpha1(SIDE3d1(ilower,iupper,alpha_gcw),0:NDIM-1)
      REAL alpha2(SIDE3d2(ilower,iupper,alpha_gcw),0:NDIM-1)

      REAL U(CELL3d(ilower,iupper,U_gcw))

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL F(CELL3d(ilower,iupper,F_gcw))
c
c     Local variables.
c
      INTEGER i0,i1,i2,k
      REAL    f0(0:1),f1(0:1),f2(0:1)
c
c     Compute the cell centered anisotropic Laplacian of U.
c
      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
               do k = 0,1
                  f0(k) =
     &                 alpha0(i0+k,i1,i2,0)*
     &                 ( U(i0+k,i1,i2)-U(i0-1+k,i1,i2) )      /dx(0) +
     &
     &                 alpha0(i0+k,i1,i2,1)*0.125d0*
     &                 ( U(i0  +k,i1+1,i2)-U(i0  +k,i1-1,i2)
     &                 + U(i0-1+k,i1+1,i2)-U(i0-1+k,i1-1,i2) )/dx(1) +
     &
     &                 alpha0(i0+k,i1,i2,2)*0.125d0*
     &                 ( U(i0  +k,i1,i2+1)-U(i0  +k,i1,i2-1)
     &                 + U(i0-1+k,i1,i2+1)-U(i0-1+k,i1,i2-1) )/dx(2)

                  f1(k) =
     &                 alpha1(i0,i1+k,i2,0)*0.125d0*
     &                 ( U(i0+1,i1  +k,i2)-U(i0-1,i1  +k,i2)
     &                 + U(i0+1,i1-1+k,i2)-U(i0-1,i1-1+k,i2) )/dx(0) +
     &
     &                 alpha1(i0,i1+k,i2,1)*
     &                 ( U(i0,i1+k,i2)-U(i0,i1-1+k,i2) )      /dx(1) +
     &
     &                 alpha1(i0,i1+k,i2,2)*0.125d0*
     &                 ( U(i0,i1  +k,i2+1)-U(i0,i1  +k,i2-1)
     &                 + U(i0,i1-1+k,i2+1)-U(i0,i1-1+k,i2-1) )/dx(2)

                  f2(k) =
     &                 alpha2(i0,i1,i2+k,0)*0.125d0*
     &                 ( U(i0+1,i1,i2  +k)-U(i0-1,i1,i2  +k)
     &                 + U(i0+1,i1,i2-1+k)-U(i0-1,i1,i2-1+k) )/dx(0) +
     &
     &                 alpha2(i0,i1,i2+k,1)*0.125d0*
     &                 ( U(i0,i1+1,i2  +k)-U(i0,i1-1,i2  +k)
     &                 + U(i0,i1+1,i2-1+k)-U(i0,i1-1,i2-1+k) )/dx(1) +
     &
     &                 alpha2(i0,i1,i2+k,2)*
     &                 ( U(i0,i1,i2+k)-U(i0,i1,i2-1+k) )      /dx(2)
               enddo

               F(i0,i1,i2) =
     &              (f0(1)-f0(0))/dx(0) +
     &              (f1(1)-f1(0))/dx(1) +
     &              (f2(1)-f2(0))/dx(2)
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes F = div alpha grad U + beta V.
c
c     Uses a nineteen point stencil to compute the cell centered
c     anisotropic Laplacian of a cell centered variable U, where the
c     general (non-grid aligned) anisotropic diffusion coefficient alpha
c     is specified on cell sides.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ctocgenanisoslaplaceadd3d(
     &     F,F_gcw,
     &     alpha0,alpha1,alpha2,alpha_gcw,
     &     U,U_gcw,
     &     beta,
     &     V,V_gcw,
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
      INTEGER F_gcw,alpha_gcw,U_gcw,V_gcw

      REAL alpha0(SIDE3d0(ilower,iupper,alpha_gcw),0:NDIM-1)
      REAL alpha1(SIDE3d1(ilower,iupper,alpha_gcw),0:NDIM-1)
      REAL alpha2(SIDE3d2(ilower,iupper,alpha_gcw),0:NDIM-1)

      REAL U(CELL3d(ilower,iupper,U_gcw))

      REAL beta

      REAL V(CELL3d(ilower,iupper,V_gcw))

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL F(CELL3d(ilower,iupper,F_gcw))
c
c     Local variables.
c
      INTEGER i0,i1,i2,k
      REAL    f0(0:1),f1(0:1),f2(0:1)
c
c     Compute the cell centered anisotropic Laplacian of U.
c
      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
               do k = 0,1
                  f0(k) =
     &                 alpha0(i0+k,i1,i2,0)*
     &                 ( U(i0+k,i1,i2)-U(i0-1+k,i1,i2) )      /dx(0) +
     &
     &                 alpha0(i0+k,i1,i2,1)*0.125d0*
     &                 ( U(i0  +k,i1+1,i2)-U(i0  +k,i1-1,i2)
     &                 + U(i0-1+k,i1+1,i2)-U(i0-1+k,i1-1,i2) )/dx(1) +
     &
     &                 alpha0(i0+k,i1,i2,2)*0.125d0*
     &                 ( U(i0  +k,i1,i2+1)-U(i0  +k,i1,i2-1)
     &                 + U(i0-1+k,i1,i2+1)-U(i0-1+k,i1,i2-1) )/dx(2)

                  f1(k) =
     &                 alpha1(i0,i1+k,i2,0)*0.125d0*
     &                 ( U(i0+1,i1  +k,i2)-U(i0-1,i1  +k,i2)
     &                 + U(i0+1,i1-1+k,i2)-U(i0-1,i1-1+k,i2) )/dx(0) +
     &
     &                 alpha1(i0,i1+k,i2,1)*
     &                 ( U(i0,i1+k,i2)-U(i0,i1-1+k,i2) )      /dx(1) +
     &
     &                 alpha1(i0,i1+k,i2,2)*0.125d0*
     &                 ( U(i0,i1  +k,i2+1)-U(i0,i1  +k,i2-1)
     &                 + U(i0,i1-1+k,i2+1)-U(i0,i1-1+k,i2-1) )/dx(2)

                  f2(k) =
     &                 alpha2(i0,i1,i2+k,0)*0.125d0*
     &                 ( U(i0+1,i1,i2  +k)-U(i0-1,i1,i2  +k)
     &                 + U(i0+1,i1,i2-1+k)-U(i0-1,i1,i2-1+k) )/dx(0) +
     &
     &                 alpha2(i0,i1,i2+k,1)*0.125d0*
     &                 ( U(i0,i1+1,i2  +k)-U(i0,i1-1,i2  +k)
     &                 + U(i0,i1+1,i2-1+k)-U(i0,i1-1,i2-1+k) )/dx(1) +
     &
     &                 alpha2(i0,i1,i2+k,2)*
     &                 ( U(i0,i1,i2+k)-U(i0,i1,i2-1+k) )      /dx(2)
               enddo

               F(i0,i1,i2) = beta*V(i0,i1,i2) +
     &              (f0(1)-f0(0))/dx(0) +
     &              (f1(1)-f1(0))/dx(1) +
     &              (f2(1)-f2(0))/dx(2)
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes F = div alpha grad U + beta U.
c
c     Uses a nineteen point stencil to compute the damped cell centered
c     anisotropic Laplacian of a cell centered variable U, where the
c     general (non-grid aligned) anisotropic diffusion coefficient alpha
c     is specified on cell sides.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ctocgenanisosdampedlaplace3d(
     &     F,F_gcw,
     &     alpha0,alpha1,alpha2,alpha_gcw,
     &     beta,
     &     U,U_gcw,
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
      INTEGER F_gcw,alpha_gcw,U_gcw

      REAL alpha0(SIDE3d0(ilower,iupper,alpha_gcw),0:NDIM-1)
      REAL alpha1(SIDE3d1(ilower,iupper,alpha_gcw),0:NDIM-1)
      REAL alpha2(SIDE3d2(ilower,iupper,alpha_gcw),0:NDIM-1)

      REAL beta

      REAL U(CELL3d(ilower,iupper,U_gcw))

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL F(CELL3d(ilower,iupper,F_gcw))
c
c     Local variables.
c
      INTEGER i0,i1,i2,k
      REAL    f0(0:1),f1(0:1),f2(0:1)
c
c     Compute the cell centered anisotropic Laplacian of U.
c
      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
               do k = 0,1
                  f0(k) =
     &                 alpha0(i0+k,i1,i2,0)*
     &                 ( U(i0+k,i1,i2)-U(i0-1+k,i1,i2) )      /dx(0) +
     &
     &                 alpha0(i0+k,i1,i2,1)*0.125d0*
     &                 ( U(i0  +k,i1+1,i2)-U(i0  +k,i1-1,i2)
     &                 + U(i0-1+k,i1+1,i2)-U(i0-1+k,i1-1,i2) )/dx(1) +
     &
     &                 alpha0(i0+k,i1,i2,2)*0.125d0*
     &                 ( U(i0  +k,i1,i2+1)-U(i0  +k,i1,i2-1)
     &                 + U(i0-1+k,i1,i2+1)-U(i0-1+k,i1,i2-1) )/dx(2)

                  f1(k) =
     &                 alpha1(i0,i1+k,i2,0)*0.125d0*
     &                 ( U(i0+1,i1  +k,i2)-U(i0-1,i1  +k,i2)
     &                 + U(i0+1,i1-1+k,i2)-U(i0-1,i1-1+k,i2) )/dx(0) +
     &
     &                 alpha1(i0,i1+k,i2,1)*
     &                 ( U(i0,i1+k,i2)-U(i0,i1-1+k,i2) )      /dx(1) +
     &
     &                 alpha1(i0,i1+k,i2,2)*0.125d0*
     &                 ( U(i0,i1  +k,i2+1)-U(i0,i1  +k,i2-1)
     &                 + U(i0,i1-1+k,i2+1)-U(i0,i1-1+k,i2-1) )/dx(2)

                  f2(k) =
     &                 alpha2(i0,i1,i2+k,0)*0.125d0*
     &                 ( U(i0+1,i1,i2  +k)-U(i0-1,i1,i2  +k)
     &                 + U(i0+1,i1,i2-1+k)-U(i0-1,i1,i2-1+k) )/dx(0) +
     &
     &                 alpha2(i0,i1,i2+k,1)*0.125d0*
     &                 ( U(i0,i1+1,i2  +k)-U(i0,i1-1,i2  +k)
     &                 + U(i0,i1+1,i2-1+k)-U(i0,i1-1,i2-1+k) )/dx(1) +
     &
     &                 alpha2(i0,i1,i2+k,2)*
     &                 ( U(i0,i1,i2+k)-U(i0,i1,i2-1+k) )      /dx(2)
               enddo

               F(i0,i1,i2) = beta*U(i0,i1,i2) +
     &              (f0(1)-f0(0))/dx(0) +
     &              (f1(1)-f1(0))/dx(1) +
     &              (f2(1)-f2(0))/dx(2)
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes F = div alpha grad U + beta U + gamma V.
c
c     Uses a nineteen point stencil to compute the damped cell centered
c     anisotropic Laplacian of a cell centered variable U, where the
c     general (non-grid aligned) anisotropic diffusion coefficient alpha
c     is specified on cell sides.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ctocgenanisosdampedlaplaceadd3d(
     &     F,F_gcw,
     &     alpha0,alpha1,alpha2,alpha_gcw,
     &     beta,
     &     U,U_gcw,
     &     gamma,
     &     V,V_gcw,
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
      INTEGER F_gcw,alpha_gcw,U_gcw,V_gcw

      REAL alpha0(SIDE3d0(ilower,iupper,alpha_gcw),0:NDIM-1)
      REAL alpha1(SIDE3d1(ilower,iupper,alpha_gcw),0:NDIM-1)
      REAL alpha2(SIDE3d2(ilower,iupper,alpha_gcw),0:NDIM-1)

      REAL beta

      REAL U(CELL3d(ilower,iupper,U_gcw))

      REAL gamma

      REAL V(CELL3d(ilower,iupper,V_gcw))

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL F(CELL3d(ilower,iupper,F_gcw))
c
c     Local variables.
c
      INTEGER i0,i1,i2,k
      REAL    f0(0:1),f1(0:1),f2(0:1)
c
c     Compute the cell centered anisotropic Laplacian of U.
c
      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
               do k = 0,1
                  f0(k) =
     &                 alpha0(i0+k,i1,i2,0)*
     &                 ( U(i0+k,i1,i2)-U(i0-1+k,i1,i2) )      /dx(0) +
     &
     &                 alpha0(i0+k,i1,i2,1)*0.125d0*
     &                 ( U(i0  +k,i1+1,i2)-U(i0  +k,i1-1,i2)
     &                 + U(i0-1+k,i1+1,i2)-U(i0-1+k,i1-1,i2) )/dx(1) +
     &
     &                 alpha0(i0+k,i1,i2,2)*0.125d0*
     &                 ( U(i0  +k,i1,i2+1)-U(i0  +k,i1,i2-1)
     &                 + U(i0-1+k,i1,i2+1)-U(i0-1+k,i1,i2-1) )/dx(2)

                  f1(k) =
     &                 alpha1(i0,i1+k,i2,0)*0.125d0*
     &                 ( U(i0+1,i1  +k,i2)-U(i0-1,i1  +k,i2)
     &                 + U(i0+1,i1-1+k,i2)-U(i0-1,i1-1+k,i2) )/dx(0) +
     &
     &                 alpha1(i0,i1+k,i2,1)*
     &                 ( U(i0,i1+k,i2)-U(i0,i1-1+k,i2) )      /dx(1) +
     &
     &                 alpha1(i0,i1+k,i2,2)*0.125d0*
     &                 ( U(i0,i1  +k,i2+1)-U(i0,i1  +k,i2-1)
     &                 + U(i0,i1-1+k,i2+1)-U(i0,i1-1+k,i2-1) )/dx(2)

                  f2(k) =
     &                 alpha2(i0,i1,i2+k,0)*0.125d0*
     &                 ( U(i0+1,i1,i2  +k)-U(i0-1,i1,i2  +k)
     &                 + U(i0+1,i1,i2-1+k)-U(i0-1,i1,i2-1+k) )/dx(0) +
     &
     &                 alpha2(i0,i1,i2+k,1)*0.125d0*
     &                 ( U(i0,i1+1,i2  +k)-U(i0,i1-1,i2  +k)
     &                 + U(i0,i1+1,i2-1+k)-U(i0,i1-1,i2-1+k) )/dx(1) +
     &
     &                 alpha2(i0,i1,i2+k,2)*
     &                 ( U(i0,i1,i2+k)-U(i0,i1,i2-1+k) )      /dx(2)
               enddo

               F(i0,i1,i2) = beta*U(i0,i1,i2) + gamma*V(i0,i1,i2) +
     &              (f0(1)-f0(0))/dx(0) +
     &              (f1(1)-f1(0))/dx(1) +
     &              (f2(1)-f2(0))/dx(2)
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
