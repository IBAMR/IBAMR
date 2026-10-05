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
c     Computes D = alpha div grad U + beta U + gamma W on the cells of
c     a patch, evaluating each flux once.
c
c     Uses the seven point stencil.  The side centered fluxes alpha grad U
c     are computed one row at a time into the scratch arrays gx, gy, and
c     gz.  The terms beta U and gamma W are omitted when beta and gamma
c     are zero.  The fluxes on the sides of the patch box that lie on its
c     boundary are returned in FL0, FU0, FL1, FU1, FL2, FU2 (lower and
c     upper sides in each coordinate direction).
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine fluxlaplace3d(
     &     D,D_gcw,
     &     alpha,beta,
     &     U,U_gcw,
     &     gamma,
     &     W,W_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
     &     dx,
     &     FL0,FU0,FL1,FU1,FL2,FU2,
     &     gx,gy,gz)
c
      implicit none
c
c     Input.
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2
      INTEGER D_gcw,U_gcw,W_gcw

      REAL alpha,beta

      REAL U(CELL3d(ilower,iupper,U_gcw))

      REAL gamma

      REAL W(CELL3d(ilower,iupper,W_gcw))

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL D(CELL3d(ilower,iupper,D_gcw))

      REAL FL0(ilower1:iupper1,ilower2:iupper2)
      REAL FU0(ilower1:iupper1,ilower2:iupper2)
      REAL FL1(ilower0:iupper0,ilower2:iupper2)
      REAL FU1(ilower0:iupper0,ilower2:iupper2)
      REAL FL2(ilower0:iupper0,ilower1:iupper1)
      REAL FU2(ilower0:iupper0,ilower1:iupper1)

      REAL gx(ilower0:iupper0+1)
      REAL gy(ilower0:iupper0,0:1)
      REAL gz(ilower0:iupper0,ilower1:iupper1,0:1)
c
c     Local variables.
c
      INTEGER i0,i1,i2,jl,jh,kl,kh
      REAL    fac0,fac1,fac2,dfac0,dfac1,dfac2
c
c     Compute D.
c
      fac0 = alpha/dx(0)
      fac1 = alpha/dx(1)
      fac2 = alpha/dx(2)
      dfac0 = 1.d0/dx(0)
      dfac1 = 1.d0/dx(1)
      dfac2 = 1.d0/dx(2)

      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0
            gz(i0,i1,0) = fac2*(U(i0,i1,ilower2)-U(i0,i1,ilower2-1))
            FL2(i0,i1) = gz(i0,i1,0)
         enddo
      enddo
      kh = 0
      do i2 = ilower2,iupper2
         kl = mod(i2-ilower2,2)
         kh = 1-kl
         do i0 = ilower0,iupper0
            gy(i0,0) = fac1*(U(i0,ilower1,i2)-U(i0,ilower1-1,i2))
            FL1(i0,i2) = gy(i0,0)
         enddo
         jh = 0
         do i1 = ilower1,iupper1
            jl = mod(i1-ilower1,2)
            jh = 1-jl
            do i0 = ilower0,iupper0+1
               gx(i0) = fac0*(U(i0,i1,i2)-U(i0-1,i1,i2))
            enddo
            do i0 = ilower0,iupper0
               gy(i0,jh) = fac1*(U(i0,i1+1,i2)-U(i0,i1,i2))
            enddo
            do i0 = ilower0,iupper0
               gz(i0,i1,kh) = fac2*(U(i0,i1,i2+1)-U(i0,i1,i2))
            enddo
            FL0(i1,i2) = gx(ilower0)
            FU0(i1,i2) = gx(iupper0+1)
            if (beta.eq.0.d0) then
               if (gamma.eq.0.d0) then
                  do i0 = ilower0,iupper0
                     D(i0,i1,i2) =
     &                    dfac0*(gx(i0+1)-gx(i0)) +
     &                    dfac1*(gy(i0,jh)-gy(i0,jl)) +
     &                    dfac2*(gz(i0,i1,kh)-gz(i0,i1,kl))
                  enddo
               else
                  do i0 = ilower0,iupper0
                     D(i0,i1,i2) =
     &                    dfac0*(gx(i0+1)-gx(i0)) +
     &                    dfac1*(gy(i0,jh)-gy(i0,jl)) +
     &                    dfac2*(gz(i0,i1,kh)-gz(i0,i1,kl)) +
     &                    gamma*W(i0,i1,i2)
                  enddo
               endif
            else
               if (gamma.eq.0.d0) then
                  do i0 = ilower0,iupper0
                     D(i0,i1,i2) =
     &                    dfac0*(gx(i0+1)-gx(i0)) +
     &                    dfac1*(gy(i0,jh)-gy(i0,jl)) +
     &                    dfac2*(gz(i0,i1,kh)-gz(i0,i1,kl)) +
     &                    beta*U(i0,i1,i2)
                  enddo
               else
                  do i0 = ilower0,iupper0
                     D(i0,i1,i2) =
     &                    (dfac0*(gx(i0+1)-gx(i0)) +
     &                     dfac1*(gy(i0,jh)-gy(i0,jl)) +
     &                     dfac2*(gz(i0,i1,kh)-gz(i0,i1,kl)) +
     &                     beta*U(i0,i1,i2)) +
     &                    gamma*W(i0,i1,i2)
                  enddo
               endif
            endif
         enddo
         do i0 = ilower0,iupper0
            FU1(i0,i2) = gy(i0,jh)
         enddo
      enddo
      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0
            FU2(i0,i1) = gz(i0,i1,kh)
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
c
c     Computes F = alpha div grad U + beta U + gamma grad P + delta V.
c
c     Uses the seven point stencil to compute the damped discrete
c     Laplacian of a side centered variable U=(U0,U1,U2) and centered
c     differences to compute the partial gradient of a cell centered
c     variable P on each side.  The sum is evaluated as
c     (alpha div grad U + beta U) + (gamma grad P + delta V).
c     The term beta U is omitted when beta is zero.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
define(BETA_TERM,`beta*U$1(i0,i1,i2) +
     &         ')dnl
define(LAPLACE_GRAD_LOOP,`dnl
         do i2 = ilower2,iupper2$4
            do i1 = ilower1,iupper1$3
               do i0 = ilower0,iupper0$2
                  F$1(i0,i1,i2) =
     &         fac0*(U$1(i0-1,i1,i2)+U$1(i0+1,i1,i2)-2.d0*U$1(i0,i1,i2))+
     &         fac1*(U$1(i0,i1-1,i2)+U$1(i0,i1+1,i2)-2.d0*U$1(i0,i1,i2))+
     &         fac2*(U$1(i0,i1,i2-1)+U$1(i0,i1,i2+1)-2.d0*U$1(i0,i1,i2))+
     &         ifelse($6,1,`BETA_TERM($1)')dnl
(gfac$1*(P(i0,i1,i2)-P($5)) + delta*V$1(i0,i1,i2))
               enddo
            enddo
         enddo
')dnl
define(LAPLACE_GRAD,`dnl
      if (beta.eq.0.d0) then
LAPLACE_GRAD_LOOP($1,`$2',`$3',`$4',`$5',0)dnl
      else
LAPLACE_GRAD_LOOP($1,`$2',`$3',`$4',`$5',1)dnl
      endif
')dnl
      subroutine stoslaplacegrad3d(
     &     F0,F1,F2,F_gcw,
     &     alpha,beta,
     &     U0,U1,U2,U_gcw,
     &     gamma,
     &     P,P_gcw,
     &     delta,
     &     V0,V1,V2,V_gcw,
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
      INTEGER F_gcw,U_gcw,P_gcw,V_gcw

      REAL alpha,beta

      REAL U0(SIDE3d0(ilower,iupper,U_gcw))
      REAL U1(SIDE3d1(ilower,iupper,U_gcw))
      REAL U2(SIDE3d2(ilower,iupper,U_gcw))

      REAL gamma

      REAL P(CELL3d(ilower,iupper,P_gcw))

      REAL delta

      REAL V0(SIDE3d0(ilower,iupper,V_gcw))
      REAL V1(SIDE3d1(ilower,iupper,V_gcw))
      REAL V2(SIDE3d2(ilower,iupper,V_gcw))

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL F0(SIDE3d0(ilower,iupper,F_gcw))
      REAL F1(SIDE3d1(ilower,iupper,F_gcw))
      REAL F2(SIDE3d2(ilower,iupper,F_gcw))
c
c     Local variables.
c
      INTEGER i0,i1,i2
      REAL    fac0,fac1,fac2,gfac0,gfac1,gfac2
c
c     Compute F.
c
      fac0 = alpha/(dx(0)*dx(0))
      fac1 = alpha/(dx(1)*dx(1))
      fac2 = alpha/(dx(2)*dx(2))
      gfac0 = gamma/dx(0)
      gfac1 = gamma/dx(1)
      gfac2 = gamma/dx(2)

LAPLACE_GRAD(0,`+1',,,`i0-1,i1,i2')dnl
LAPLACE_GRAD(1,,`+1',,`i0,i1-1,i2')dnl
LAPLACE_GRAD(2,,,`+1',`i0,i1,i2-1')dnl
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
