c ---------------------------------------------------------------------
c c
c  File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/patchdata/fortran/pdat_m4arrdim3d.i $
c  Package:     SAMRAI patchdata
c  Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision: 1917 $
c  Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
c  Description: m4 include file for dimensioning 3d arrays in FORTRAN routines.
c

c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    Reference Map Technique (3D)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c------------------------------ div_t_3d ---------------------------------
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine div_tmain3d(
     &     f11, f12, f13, f21, f22, f23, f31, f32, f33, f_gcw, 
     &     GGS, Txx_s, Txy_s, Txz_s, Tyx_s, Tyy_s,  
     &     Tyz_s, Tzx_s, Tzy_s, Tzz_s, ts_gcw,
     &     phi, phi_gcw, H_phi, H_phi_gcw,
     &     xi0, xi1, xi2, xi_gcw,
     &     ilower0, iupper0, ilower1, iupper1, ilower2, iupper2,
     &     dx, div_t_gcw,
     &     div_t1, div_t2, div_t3)

      implicit none

      integer ilower0, iupper0
      integer ilower1, iupper1
      integer ilower2, iupper2

      double precision dx(0:3-1)

      integer f_gcw, ts_gcw, div_t_gcw
      integer phi_gcw, H_phi_gcw, xi_gcw

      integer i0, i1, i2
      double precision fac0, fac1, fac2
      double precision GGS

      double precision  f11(ilower0-f_gcw:iupper0+f_gcw,
     &          ilower1-f_gcw:iupper1+f_gcw,
     &          ilower2-f_gcw:iupper2+f_gcw)
      double precision  f12(ilower0-f_gcw:iupper0+f_gcw,
     &          ilower1-f_gcw:iupper1+f_gcw,
     &          ilower2-f_gcw:iupper2+f_gcw)
      double precision  f13(ilower0-f_gcw:iupper0+f_gcw,
     &          ilower1-f_gcw:iupper1+f_gcw,
     &          ilower2-f_gcw:iupper2+f_gcw)
      double precision  f21(ilower0-f_gcw:iupper0+f_gcw,
     &          ilower1-f_gcw:iupper1+f_gcw,
     &          ilower2-f_gcw:iupper2+f_gcw)
      double precision  f22(ilower0-f_gcw:iupper0+f_gcw,
     &          ilower1-f_gcw:iupper1+f_gcw,
     &          ilower2-f_gcw:iupper2+f_gcw)
      double precision  f23(ilower0-f_gcw:iupper0+f_gcw,
     &          ilower1-f_gcw:iupper1+f_gcw,
     &          ilower2-f_gcw:iupper2+f_gcw)
      double precision  f31(ilower0-f_gcw:iupper0+f_gcw,
     &          ilower1-f_gcw:iupper1+f_gcw,
     &          ilower2-f_gcw:iupper2+f_gcw)
      double precision  f32(ilower0-f_gcw:iupper0+f_gcw,
     &          ilower1-f_gcw:iupper1+f_gcw,
     &          ilower2-f_gcw:iupper2+f_gcw)
      double precision  f33(ilower0-f_gcw:iupper0+f_gcw,
     &          ilower1-f_gcw:iupper1+f_gcw,
     &          ilower2-f_gcw:iupper2+f_gcw)

      double precision Txx_s(ilower0-ts_gcw:iupper0+ts_gcw,
     &          ilower1-ts_gcw:iupper1+ts_gcw,
     &          ilower2-ts_gcw:iupper2+ts_gcw)
      double precision Txy_s(ilower0-ts_gcw:iupper0+ts_gcw,
     &          ilower1-ts_gcw:iupper1+ts_gcw,
     &          ilower2-ts_gcw:iupper2+ts_gcw)
      double precision Txz_s(ilower0-ts_gcw:iupper0+ts_gcw,
     &          ilower1-ts_gcw:iupper1+ts_gcw,
     &          ilower2-ts_gcw:iupper2+ts_gcw)
      double precision Tyx_s(ilower0-ts_gcw:iupper0+ts_gcw,
     &          ilower1-ts_gcw:iupper1+ts_gcw,
     &          ilower2-ts_gcw:iupper2+ts_gcw)
      double precision Tyy_s(ilower0-ts_gcw:iupper0+ts_gcw,
     &          ilower1-ts_gcw:iupper1+ts_gcw,
     &          ilower2-ts_gcw:iupper2+ts_gcw)
      double precision Tyz_s(ilower0-ts_gcw:iupper0+ts_gcw,
     &          ilower1-ts_gcw:iupper1+ts_gcw,
     &          ilower2-ts_gcw:iupper2+ts_gcw)
      double precision Tzx_s(ilower0-ts_gcw:iupper0+ts_gcw,
     &          ilower1-ts_gcw:iupper1+ts_gcw,
     &          ilower2-ts_gcw:iupper2+ts_gcw)
      double precision Tzy_s(ilower0-ts_gcw:iupper0+ts_gcw,
     &          ilower1-ts_gcw:iupper1+ts_gcw,
     &          ilower2-ts_gcw:iupper2+ts_gcw)
      double precision Tzz_s(ilower0-ts_gcw:iupper0+ts_gcw,
     &          ilower1-ts_gcw:iupper1+ts_gcw,
     &          ilower2-ts_gcw:iupper2+ts_gcw)

      double precision phi  (ilower0-phi_gcw:iupper0+phi_gcw,
     &          ilower1-phi_gcw:iupper1+phi_gcw,
     &          ilower2-phi_gcw:iupper2+phi_gcw)
      double precision H_phi(ilower0-H_phi_gcw:iupper0+H_phi_gcw,
     &          ilower1-H_phi_gcw:iupper1+H_phi_gcw,
     &          ilower2-H_phi_gcw:iupper2+H_phi_gcw)

      double precision xi0(ilower0-xi_gcw:iupper0+xi_gcw,
     &          ilower1-xi_gcw:iupper1+xi_gcw,
     &          ilower2-xi_gcw:iupper2+xi_gcw)
      double precision xi1(ilower0-xi_gcw:iupper0+xi_gcw,
     &          ilower1-xi_gcw:iupper1+xi_gcw,
     &          ilower2-xi_gcw:iupper2+xi_gcw)
      double precision xi2(ilower0-xi_gcw:iupper0+xi_gcw,
     &          ilower1-xi_gcw:iupper1+xi_gcw,
     &          ilower2-xi_gcw:iupper2+xi_gcw)

      double precision div_t1(ilower0-div_t_gcw:iupper0+1+div_t_gcw,
     &          ilower1-div_t_gcw:iupper1+div_t_gcw,
     &          ilower2-div_t_gcw:iupper2+div_t_gcw)
      double precision div_t2(ilower0-div_t_gcw:iupper0+div_t_gcw,
     &          ilower1-div_t_gcw:iupper1+1+div_t_gcw,
     &          ilower2-div_t_gcw:iupper2+div_t_gcw)
      double precision div_t3(ilower0-div_t_gcw:iupper0+div_t_gcw,
     &          ilower1-div_t_gcw:iupper1+div_t_gcw,
     &          ilower2-div_t_gcw:iupper2+1+div_t_gcw)

      fac0 = 1.0d0/dx(0)
      fac1 = 1.0d0/dx(1)
      fac2 = 1.0d0/dx(2)

      call solid_stress_tensor_3d(
     &     Txx_s, Txy_s, Txz_s, Tyx_s, Tyy_s, 
     &     Tyz_s, Tzx_s, Tzy_s, Tzz_s,
     &     ilower0, iupper0, ilower1, iupper1, 
     &     ilower2, iupper2, ts_gcw,
     &     phi, phi_gcw, xi0, xi1, xi2, xi_gcw, 
     &     H_phi, H_phi_gcw, dx,
     &     f11, f12, f13, f21, f22, f23, 
     &     f31, f32, f33, f_gcw, GGS)

c     div_t1 on x-faces: d/dx Txx + d/dy Tyx + d/dz Tzx
      do i2 = ilower2, iupper2
        do i1 = ilower1, iupper1
          do i0 = ilower0, iupper0+1
            div_t1(i0,i1,i2) =
     &          (Txx_s(i0,i1,i2) - Txx_s(i0-1,i1,i2))*fac0
     &        + 0.25d0 * (
     &          (Tyx_s(i0,i1+1,i2) + Tyx_s(i0-1,i1+1,i2))
     &        - (Tyx_s(i0,i1-1,i2) + Tyx_s(i0-1,i1-1,i2))
     &        ) * fac1
     &        + 0.25d0 * (
     &          (Tzx_s(i0,i1,i2+1) + Tzx_s(i0-1,i1,i2+1))
     &        - (Tzx_s(i0,i1,i2-1) + Tzx_s(i0-1,i1,i2-1))
     &        ) * fac2
          end do
        end do
      end do

c     div_t2 on y-faces: d/dx Txy + d/dy Tyy + d/dz Tzy
      do i2 = ilower2, iupper2
        do i1 = ilower1, iupper1+1
          do i0 = ilower0, iupper0
            div_t2(i0,i1,i2) =
     &          (Tyy_s(i0,i1,i2) - Tyy_s(i0,i1-1,i2))*fac1
     &        + 0.25d0 * (
     &          (Txy_s(i0+1,i1,i2) + Txy_s(i0+1,i1-1,i2))
     &        - (Txy_s(i0-1,i1,i2) + Txy_s(i0-1,i1-1,i2))
     &        ) * fac0
     &        + 0.25d0 * (
     &          (Tzy_s(i0,i1,i2+1) + Tzy_s(i0,i1-1,i2+1))
     &        - (Tzy_s(i0,i1,i2-1) + Tzy_s(i0,i1-1,i2-1))
     &        ) * fac2
          end do
        end do
      end do

c     div_t3 on z-faces: d/dx Txz + d/dy Tyz + d/dz Tzz
      do i2 = ilower2, iupper2+1
        do i1 = ilower1, iupper1
          do i0 = ilower0, iupper0
            div_t3(i0,i1,i2) =
     &          (Tzz_s(i0,i1,i2) - Tzz_s(i0,i1,i2-1))*fac2
     &        + 0.25d0 * (
     &          (Txz_s(i0+1,i1,i2) + Txz_s(i0+1,i1,i2-1))
     &        - (Txz_s(i0-1,i1,i2) + Txz_s(i0-1,i1,i2-1))
     &        ) * fac0
     &        + 0.25d0 * (
     &          (Tyz_s(i0,i1+1,i2) + Tyz_s(i0,i1+1,i2-1))
     &        - (Tyz_s(i0,i1-1,i2) + Tyz_s(i0,i1-1,i2-1))
     &        ) * fac1
          end do
        end do
      end do


c       ==== DEBUG: constant div_t forcing ====
      do i2 = ilower2, iupper2
      do i1 = ilower1, iupper1
        do i0 = ilower0, iupper0+1
          div_t1(i0,i1, i2) = 0.d0
        enddo 
      enddo
      end do

      do i2 = ilower2, iupper2
      do i1 = ilower1, iupper1+1
        do i0 = ilower0, iupper0
          div_t2(i0,i1, i2) = 0.d0
        enddo
      enddo
      end do

       do i2 = ilower2+1, iupper2
      do i1 = ilower1, iupper1
        do i0 = ilower0, iupper0
          div_t3(i0,i1, i2) = 0.d0
        enddo
      enddo
      end do



      return
      end

c-------------------- solid_stress_tensor_3d --------------------
      subroutine solid_stress_tensor_3d(
     &     Txx_s, Txy_s, Txz_s, Tyx_s, Tyy_s, Tyz_s, 
     &     Tzx_s, Tzy_s, Tzz_s,
     &     ilower0, iupper0, ilower1, iupper1, ilower2, 
     &     iupper2, ts_gcw,
     &     phi, phi_gcw, xi0, xi1, xi2, xi_gcw, H_phi, 
     &     H_phi_gcw, dx,
     &     f11, f12, f13, f21, f22, f23, f31, 
     &     f32, f33, f_gcw, GGS)

      implicit none

      integer ilower0, iupper0
      integer ilower1, iupper1
      integer ilower2, iupper2

      double precision dx(0:3-1)

      integer ts_gcw, phi_gcw, xi_gcw, H_phi_gcw, f_gcw
      integer i0, i1, i2

      double precision GGS
      double precision B11,B12,B13,B21,B22,B23,B31,B32,B33
      double precision trB,iso_term

      double precision phi(ilower0-phi_gcw:iupper0+phi_gcw,
     &          ilower1-phi_gcw:iupper1+phi_gcw,
     &          ilower2-phi_gcw:iupper2+phi_gcw)
      double precision xi0(ilower0-xi_gcw:iupper0+xi_gcw,
     &          ilower1-xi_gcw:iupper1+xi_gcw,
     &          ilower2-xi_gcw:iupper2+xi_gcw)
      double precision xi1(ilower0-xi_gcw:iupper0+xi_gcw,
     &          ilower1-xi_gcw:iupper1+xi_gcw,
     &          ilower2-xi_gcw:iupper2+xi_gcw)
      double precision xi2(ilower0-xi_gcw:iupper0+xi_gcw,
     &          ilower1-xi_gcw:iupper1+xi_gcw,
     &          ilower2-xi_gcw:iupper2+xi_gcw)

      double precision Txx_s(ilower0-ts_gcw:iupper0+ts_gcw,
     &          ilower1-ts_gcw:iupper1+ts_gcw,
     &          ilower2-ts_gcw:iupper2+ts_gcw)
      double precision Txy_s(ilower0-ts_gcw:iupper0+ts_gcw,
     &          ilower1-ts_gcw:iupper1+ts_gcw,
     &          ilower2-ts_gcw:iupper2+ts_gcw)
      double precision Txz_s(ilower0-ts_gcw:iupper0+ts_gcw,
     &          ilower1-ts_gcw:iupper1+ts_gcw,
     &          ilower2-ts_gcw:iupper2+ts_gcw)
      double precision Tyx_s(ilower0-ts_gcw:iupper0+ts_gcw,
     &          ilower1-ts_gcw:iupper1+ts_gcw,
     &          ilower2-ts_gcw:iupper2+ts_gcw)
      double precision Tyy_s(ilower0-ts_gcw:iupper0+ts_gcw,
     &          ilower1-ts_gcw:iupper1+ts_gcw,
     &          ilower2-ts_gcw:iupper2+ts_gcw)
      double precision Tyz_s(ilower0-ts_gcw:iupper0+ts_gcw,
     &          ilower1-ts_gcw:iupper1+ts_gcw,
     &          ilower2-ts_gcw:iupper2+ts_gcw)
      double precision Tzx_s(ilower0-ts_gcw:iupper0+ts_gcw,
     &          ilower1-ts_gcw:iupper1+ts_gcw,
     &          ilower2-ts_gcw:iupper2+ts_gcw)
      double precision Tzy_s(ilower0-ts_gcw:iupper0+ts_gcw,
     &          ilower1-ts_gcw:iupper1+ts_gcw,
     &          ilower2-ts_gcw:iupper2+ts_gcw)
      double precision Tzz_s(ilower0-ts_gcw:iupper0+ts_gcw,
     &          ilower1-ts_gcw:iupper1+ts_gcw,
     &          ilower2-ts_gcw:iupper2+ts_gcw)

      double precision H_phi(ilower0-H_phi_gcw:iupper0+H_phi_gcw,
     &          ilower1-H_phi_gcw:iupper1+H_phi_gcw,
     &          ilower2-H_phi_gcw:iupper2+H_phi_gcw)

      double precision f11(ilower0-f_gcw:iupper0+f_gcw,
     &          ilower1-f_gcw:iupper1+f_gcw,
     &          ilower2-f_gcw:iupper2+f_gcw)
      double precision f12(ilower0-f_gcw:iupper0+f_gcw,
     &          ilower1-f_gcw:iupper1+f_gcw,
     &          ilower2-f_gcw:iupper2+f_gcw)
      double precision f13(ilower0-f_gcw:iupper0+f_gcw,
     &          ilower1-f_gcw:iupper1+f_gcw,
     &          ilower2-f_gcw:iupper2+f_gcw)
      double precision f21(ilower0-f_gcw:iupper0+f_gcw,
     &          ilower1-f_gcw:iupper1+f_gcw,
     &          ilower2-f_gcw:iupper2+f_gcw)
      double precision f22(ilower0-f_gcw:iupper0+f_gcw,
     &          ilower1-f_gcw:iupper1+f_gcw,
     &          ilower2-f_gcw:iupper2+f_gcw)
      double precision f23(ilower0-f_gcw:iupper0+f_gcw,
     &          ilower1-f_gcw:iupper1+f_gcw,
     &          ilower2-f_gcw:iupper2+f_gcw)
      double precision f31(ilower0-f_gcw:iupper0+f_gcw,
     &          ilower1-f_gcw:iupper1+f_gcw,
     &          ilower2-f_gcw:iupper2+f_gcw)
      double precision f32(ilower0-f_gcw:iupper0+f_gcw,
     &          ilower1-f_gcw:iupper1+f_gcw,
     &          ilower2-f_gcw:iupper2+f_gcw)
      double precision f33(ilower0-f_gcw:iupper0+f_gcw,
     &          ilower1-f_gcw:iupper1+f_gcw,
     &          ilower2-f_gcw:iupper2+f_gcw)

     
c     Compute smoothed Heaviside
      call smoothed_heaviside_3d(phi, H_phi, dx,
     &  ilower0,iupper0, ilower1,iupper1, ilower2,iupper2,
     &  phi_gcw, H_phi_gcw)

c     Compute deformation gradient tensor F
      call deformation_gradient_tensor_F_3d(
     &  xi0, xi1, xi2,
     &  f11, f12, f13, f21, f22, f23, f31, f32, f33,
     &  xi_gcw, f_gcw, dx,
     &  ilower0,iupper0, ilower1,iupper1, ilower2,iupper2)

c     Compute Neo-Hookean deviatoric stress: GGS*(B - (trB/3)I) masked by (1-H)
      do i2 = ilower2, iupper2
        do i1 = ilower1, iupper1
          do i0 = ilower0, iupper0

        B11 = f11(i0,i1,i2)**2 + f12(i0,i1,i2)**2 + f13(i0,i1,i2)**2
        B12 = f11(i0,i1,i2)*f21(i0,i1,i2) + f12(i0,i1,i2)*f22(i0,i1,i2)
     &          + f13(i0,i1,i2)*f23(i0,i1,i2)
        B13 = f11(i0,i1,i2)*f31(i0,i1,i2) + f12(i0,i1,i2)*f32(i0,i1,i2)
     &          + f13(i0,i1,i2)*f33(i0,i1,i2)

        B21 = f21(i0,i1,i2)*f11(i0,i1,i2) + f22(i0,i1,i2)*f12(i0,i1,i2)
     &          + f23(i0,i1,i2)*f13(i0,i1,i2)
        B22 = f21(i0,i1,i2)**2 + f22(i0,i1,i2)**2 + f23(i0,i1,i2)**2
        B23 = f21(i0,i1,i2)*f31(i0,i1,i2) + f22(i0,i1,i2)*f32(i0,i1,i2)
     &          + f23(i0,i1,i2)*f33(i0,i1,i2)

        B31 = f31(i0,i1,i2)*f11(i0,i1,i2) + f32(i0,i1,i2)*f12(i0,i1,i2)
     &          + f33(i0,i1,i2)*f13(i0,i1,i2)
        B32 = f31(i0,i1,i2)*f21(i0,i1,i2) + f32(i0,i1,i2)*f22(i0,i1,i2)
     &          + f33(i0,i1,i2)*f23(i0,i1,i2)
        B33 = f31(i0,i1,i2)**2 + f32(i0,i1,i2)**2 + f33(i0,i1,i2)**2

            trB = B11 + B22 + B33
            iso_term = trB/3.0d0

        Txx_s(i0,i1,i2) = GGS*(B11 - iso_term) 
     &      * (1.0d0 - H_phi(i0,i1,i2))

        Txy_s(i0,i1,i2) = GGS* B12            
     &      * (1.0d0 - H_phi(i0,i1,i2))

        Txz_s(i0,i1,i2) = GGS* B13            
     &      * (1.0d0 - H_phi(i0,i1,i2))

        Tyx_s(i0,i1,i2) = GGS* B21            
     &      * (1.0d0 - H_phi(i0,i1,i2))

        Tyy_s(i0,i1,i2) = GGS*(B22 - iso_term) 
     &      * (1.0d0 - H_phi(i0,i1,i2))

        Tyz_s(i0,i1,i2) = GGS* B23            
     &      * (1.0d0 - H_phi(i0,i1,i2))

        Tzx_s(i0,i1,i2) = GGS* B31            
     &      * (1.0d0 - H_phi(i0,i1,i2))

        Tzy_s(i0,i1,i2) = GGS* B32            
     &      * (1.0d0 - H_phi(i0,i1,i2))

        Tzz_s(i0,i1,i2) = GGS*(B33 - iso_term) 
     &      * (1.0d0 - H_phi(i0,i1,i2))

          end do
        end do
      end do

      return
      end

c-------------------------- smoothed_heaviside_3d -------------------------
      subroutine smoothed_heaviside_3d(phi, H_phi, dx,
     &  ilower0,iupper0, ilower1,iupper1, ilower2,iupper2,
     &  phi_gcw, H_phi_gcw)

      implicit none

      integer ilower0, iupper0, ilower1, iupper1, ilower2, iupper2
      integer phi_gcw, H_phi_gcw

      double precision phi  (ilower0-phi_gcw:iupper0+phi_gcw,
     &          ilower1-phi_gcw:iupper1+phi_gcw,
     &          ilower2-phi_gcw:iupper2+phi_gcw)
      double precision H_phi(ilower0-H_phi_gcw:iupper0+H_phi_gcw,
     &          ilower1-H_phi_gcw:iupper1+H_phi_gcw,
     &          ilower2-H_phi_gcw:iupper2+H_phi_gcw)

      integer i0, i1, i2
      double precision dx(0:3-1)
      double precision epsilon, phi_abs
      double precision, parameter :: pi = 3.14159265358979323846

      epsilon = 1.5d0 * max(dx(0), max(dx(1), dx(2)))

      do i2 = ilower2, iupper2
        do i1 = ilower1, iupper1
          do i0 = ilower0, iupper0
            phi_abs = abs(phi(i0,i1,i2))
            if (phi(i0,i1,i2) .le. -epsilon) then
              H_phi(i0,i1,i2) = 0.0d0
            else if (phi_abs .lt. epsilon) then
        H_phi(i0,i1,i2) = 0.5d0 * ( 1.0d0 + phi(i0,i1,i2)/epsilon
     &           + sin(pi*phi(i0,i1,i2)/epsilon)/pi )
            else
              H_phi(i0,i1,i2) = 1.0d0
            end if
          end do
        end do
      end do

      return
      end

c-------------------- deformation_gradient_tensor_F_3d --------------------
c  Here: compute grad(xi) with centered differences and invert to get F = inv(grad(xi))
      subroutine deformation_gradient_tensor_F_3d(
     &     xi0, xi1, xi2,
     &     f11, f12, f13, f21, f22, f23, f31, f32, f33,
     &     xi_gcw, f_gcw, dx,
     &     ilower0,iupper0, ilower1,iupper1, ilower2,iupper2)

      implicit none

      integer xi_gcw, f_gcw
      integer ilower0, iupper0, ilower1, iupper1, ilower2, iupper2
      double precision dx(0:3-1)

      double precision xi0(ilower0-xi_gcw:iupper0+xi_gcw,
     &          ilower1-xi_gcw:iupper1+xi_gcw,
     &          ilower2-xi_gcw:iupper2+xi_gcw)
      double precision xi1(ilower0-xi_gcw:iupper0+xi_gcw,
     &          ilower1-xi_gcw:iupper1+xi_gcw,
     &          ilower2-xi_gcw:iupper2+xi_gcw)
      double precision xi2(ilower0-xi_gcw:iupper0+xi_gcw,
     &          ilower1-xi_gcw:iupper1+xi_gcw,
     &          ilower2-xi_gcw:iupper2+xi_gcw)

      double precision f11(ilower0-f_gcw:iupper0+f_gcw,
     &          ilower1-f_gcw:iupper1+f_gcw,
     &          ilower2-f_gcw:iupper2+f_gcw)
      double precision f12(ilower0-f_gcw:iupper0+f_gcw,
     &          ilower1-f_gcw:iupper1+f_gcw,
     &          ilower2-f_gcw:iupper2+f_gcw)
      double precision f13(ilower0-f_gcw:iupper0+f_gcw,
     &          ilower1-f_gcw:iupper1+f_gcw,
     &          ilower2-f_gcw:iupper2+f_gcw)
      double precision f21(ilower0-f_gcw:iupper0+f_gcw,
     &          ilower1-f_gcw:iupper1+f_gcw,
     &          ilower2-f_gcw:iupper2+f_gcw)
      double precision f22(ilower0-f_gcw:iupper0+f_gcw,
     &          ilower1-f_gcw:iupper1+f_gcw,
     &          ilower2-f_gcw:iupper2+f_gcw)
      double precision f23(ilower0-f_gcw:iupper0+f_gcw,
     &          ilower1-f_gcw:iupper1+f_gcw,
     &          ilower2-f_gcw:iupper2+f_gcw)
      double precision f31(ilower0-f_gcw:iupper0+f_gcw,
     &          ilower1-f_gcw:iupper1+f_gcw,
     &          ilower2-f_gcw:iupper2+f_gcw)
      double precision f32(ilower0-f_gcw:iupper0+f_gcw,
     &          ilower1-f_gcw:iupper1+f_gcw,
     &          ilower2-f_gcw:iupper2+f_gcw)
      double precision f33(ilower0-f_gcw:iupper0+f_gcw,
     &          ilower1-f_gcw:iupper1+f_gcw,
     &          ilower2-f_gcw:iupper2+f_gcw)

      integer i0, i1, i2
      double precision fac0, fac1, fac2

c     grad(xi) = A, then F = inv(A)
      double precision a11,a12,a13,a21,a22,a23,a31,a32,a33
      double precision detA
      double precision inv11,inv12,inv13,inv21,inv22,inv23
      double precision inv31,inv32,inv33

      fac0 = 0.5d0/dx(0)
      fac1 = 0.5d0/dx(1)
      fac2 = 0.5d0/dx(2)

      do i2 = ilower2, iupper2
        do i1 = ilower1, iupper1
          do i0 = ilower0, iupper0

c           A = grad(xi)
            a11 = (xi0(i0+1,i1,i2) - xi0(i0-1,i1,i2))*fac0
            a12 = (xi0(i0,i1+1,i2) - xi0(i0,i1-1,i2))*fac1
            a13 = (xi0(i0,i1,i2+1) - xi0(i0,i1,i2-1))*fac2

            a21 = (xi1(i0+1,i1,i2) - xi1(i0-1,i1,i2))*fac0
            a22 = (xi1(i0,i1+1,i2) - xi1(i0,i1-1,i2))*fac1
            a23 = (xi1(i0,i1,i2+1) - xi1(i0,i1,i2-1))*fac2

            a31 = (xi2(i0+1,i1,i2) - xi2(i0-1,i1,i2))*fac0
            a32 = (xi2(i0,i1+1,i2) - xi2(i0,i1-1,i2))*fac1
            a33 = (xi2(i0,i1,i2+1) - xi2(i0,i1,i2-1))*fac2

            detA =
     &          a11*(a22*a33 - a23*a32)
     &        - a12*(a21*a33 - a23*a31)
     &        + a13*(a21*a32 - a22*a31)

c           Inverse of A (with a tiny safeguard)
            if (abs(detA) .lt. 1.0d-14) then
              inv11 = 1.0d0; inv12 = 0.0d0; inv13 = 0.0d0
              inv21 = 0.0d0; inv22 = 1.0d0; inv23 = 0.0d0
              inv31 = 0.0d0; inv32 = 0.0d0; inv33 = 1.0d0
            else
              inv11 =  (a22*a33 - a23*a32)/detA
              inv12 = -(a12*a33 - a13*a32)/detA
              inv13 =  (a12*a23 - a13*a22)/detA

              inv21 = -(a21*a33 - a23*a31)/detA
              inv22 =  (a11*a33 - a13*a31)/detA
              inv23 = -(a11*a23 - a13*a21)/detA

              inv31 =  (a21*a32 - a22*a31)/detA
              inv32 = -(a11*a32 - a12*a31)/detA
              inv33 =  (a11*a22 - a12*a21)/detA
            end if

c           F = inv(grad(xi))
            f11(i0,i1,i2) = inv11
            f12(i0,i1,i2) = inv12
            f13(i0,i1,i2) = inv13

            f21(i0,i1,i2) = inv21
            f22(i0,i1,i2) = inv22
            f23(i0,i1,i2) = inv23

            f31(i0,i1,i2) = inv31
            f32(i0,i1,i2) = inv32
            f33(i0,i1,i2) = inv33

          end do
        end do
      end do

      return
      end
