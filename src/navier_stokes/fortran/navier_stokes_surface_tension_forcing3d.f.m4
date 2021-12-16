c ---------------------------------------------------------------------
c
c Copyright (c) 2017 - 2019 by the IBAMR developers
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
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Mollify indicator function using IB_4 kernel
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine mollify_ib_4_3d(
     &     V,V_gcw,
     &     U,U_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2)
c
      implicit none
include(TOP_SRCDIR/src/fortran/const.i)dnl
c
c     Input.
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2
      INTEGER V_gcw,U_gcw
    
c
c     Input/Output.
c
      REAL U(CELL3d(ilower,iupper,U_gcw))
      REAL V(CELL3d(ilower,iupper,V_gcw))
      REAL w0(-1:1),w1(-1:1),w2(-1:1),wy,wz
      REAL w(-1:1,-1:1,-1:1)

c
c     Local variables.
c
      INTEGER k0,k1,k2
      INTEGER i0,i1,i2
      

c     Compute 1D weights.
      w0(-1) = fourth; w1(-1) = fourth; w2(-1) = fourth
      w0(0) = half; w1(0) = half; w2(0) = half
      w0(1) = fourth; w1(1) = fourth; w2(1) = fourth

c     Compute the tensor product weight
      do k2 = -1,1
         wz = w2(k2)
         do k1 = -1,1
           wy = w1(k1)
           do k0 = -1,1
              w(k0,k1,k2) = w0(k0)*wy*wz
           enddo
        enddo  
      enddo
    
c     Mollify U to V.
      do i2 = ilower2,iupper2
        do i1 = ilower1,iupper1
          do i0 = ilower0,iupper0

              V(i0,i1,i2) = 0.d0
              do k2 = -1,1
                do k1 = -1,1
                  do k0 = -1,1
                    V(i0,i1,i2) = V(i0,i1,i2) + 
     &                   U(i0+k0,i1+k1,i2+k2)*w(k0,k1,k2)
                  enddo
                enddo
              enddo
 
          enddo
        enddo
      enddo
      

      return
      end

c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute gradient of the indicator function to estimate 
c     interface normal
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine sc_normal_3d(
     &     N00,N01,N02,
     &     N10,N11,N12,
     &     N20,N21,N22,
     &     N_gcw,
     &     U,U_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
     &     dx)
c
      implicit none
include(TOP_SRCDIR/src/fortran/const.i)dnl
c
c     Input.
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2
      INTEGER N_gcw,U_gcw
    
c
c     Input/Output.
c
      REAL N00(SIDE3d0(ilower,iupper,N_gcw))
      REAL N01(SIDE3d0(ilower,iupper,N_gcw))
      REAL N02(SIDE3d0(ilower,iupper,N_gcw))
      REAL N10(SIDE3d1(ilower,iupper,N_gcw))
      REAL N11(SIDE3d1(ilower,iupper,N_gcw))
      REAL N12(SIDE3d1(ilower,iupper,N_gcw))
      REAL N20(SIDE3d2(ilower,iupper,N_gcw))
      REAL N21(SIDE3d2(ilower,iupper,N_gcw))
      REAL N22(SIDE3d2(ilower,iupper,N_gcw))
      REAL U(CELL3d(ilower,iupper,U_gcw))
     
      REAL dx(0:NDIM-1)

c
c     Local variables.
c
      INTEGER i0,i1,i2
      REAL fac0,fac1,fac2
      
      fac0 = 1.d0/(dx(0))
      fac1 = 1.d0/(dx(1))
      fac2 = 1.d0/(dx(2))      

c
c     Find face normal gradients first and then interpolate 
c     face tangential gradients

c     Do N00.
      do i2 = ilower2 - 2, iupper2 + 2
        do i1 = ilower1 - 2, iupper1 + 2
          do i0 = ilower0 - 1, iupper0 + 2
             
            N00(i0,i1,i2) = fac0*(U(i0,i1,i2) - U(i0-1,i1,i2))

          enddo
        enddo
      enddo

      
c     Do N11.
      do i2 = ilower2 - 2, iupper2 + 2
        do i1 = ilower1 - 1, iupper1 + 2
          do i0 = ilower0 - 2, iupper0 + 2
             
              N11(i0,i1,i2) = fac1*(U(i0,i1,i2) - U(i0,i1-1,i2))

          enddo
        enddo
      enddo

c     Do N22.
      do i2 = ilower2 - 1, iupper2 + 2
        do i1 = ilower1 - 2, iupper1 + 2
          do i0 = ilower0 - 2, iupper0 + 2
             
              N22(i0,i1,i2) = fac2*(U(i0,i1,i2) - U(i0,i1,i2-1))

          enddo
        enddo
      enddo

c     Interpolate N11 to N01
      do i2 = ilower2 - 1, iupper2 + 1
        do i1 = ilower1 - 1, iupper1 + 1
          do i0 = ilower0 - 1, iupper0 + 2
             
              N01(i0,i1,i2) = fourth*(N11(i0-1,i1,i2) + N11(i0,i1,i2) + 
     &                         N11(i0-1,i1+1,i2) + N11(i0,i1+1,i2)) 
          enddo
        enddo
      enddo

c     Interpolate N22 to N02
      do i2 = ilower2 - 1, iupper2 + 1
        do i1 = ilower1 - 1, iupper1 + 1
          do i0 = ilower0 - 1, iupper0 + 2
             
              N02(i0,i1,i2) = fourth*(N22(i0-1,i1,i2) + N22(i0,i1,i2) + 
     &                         N22(i0-1,i1,i2+1) + N22(i0,i1,i2+1)) 
          enddo
        enddo
      enddo

c     Interpolate N00 to N10
      do i2 = ilower2 - 1, iupper2 + 1
        do i1 = ilower1 - 1, iupper1 + 2
          do i0 = ilower0 - 1, iupper0 + 1

              N10(i0,i1,i2) = fourth*(N00(i0,i1,i2) + N00(i0+1,i1,i2) +
     &                         N00(i0,i1-1,i2) + N00(i0+1,i1-1,i2))

          enddo
        enddo
      enddo

c     Interpolate N22 to N12
      do i2 = ilower2 - 1, iupper2 + 1
        do i1 = ilower1 - 1, iupper1 + 2
          do i0 = ilower0 - 1, iupper0 + 1

              N12(i0,i1,i2) = fourth*(N22(i0,i1,i2) + N22(i0,i1,i2+1) +
     &                         N22(i0,i1-1,i2) + N22(i0,i1-1,i2+1))

          enddo
        enddo
      enddo

c     Interpolate N00 to N20
      do i2 = ilower2 - 1, iupper2 + 2
        do i1 = ilower1 - 1, iupper1 + 1
          do i0 = ilower0 - 1, iupper0 + 1

              N20(i0,i1,i2) = fourth*(N00(i0,i1,i2) + N00(i0+1,i1,i2) +
     &                         N00(i0,i1,i2-1) + N00(i0+1,i1,i2-1))

          enddo
        enddo
      enddo

c     Interpolate N11 to N21
      do i2 = ilower2 - 1, iupper2 + 2
        do i1 = ilower1 - 1, iupper1 + 1
          do i0 = ilower0 - 1, iupper0 + 1
             
              N21(i0,i1,i2) = fourth*(N11(i0,i1,i2-1) + N11(i0,i1,i2) + 
     &                         N11(i0,i1+1,i2-1) + N11(i0,i1+1,i2)) 
          enddo
        enddo
      enddo

      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute cell center curvature of the interface.
c
c           K = - div (n/|n|)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine cc_curvature_3d(
     &     K,K_gcw,
     &     N00,N01,N02,
     &     N10,N11,N12,
     &     N20,N21,N22,
     &     N_gcw,
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
      INTEGER K_gcw,N_gcw
    
c
c     Input/Output.
c
      REAL N00(SIDE3d0(ilower,iupper,N_gcw))
      REAL N01(SIDE3d0(ilower,iupper,N_gcw))
      REAL N02(SIDE3d0(ilower,iupper,N_gcw))
      REAL N10(SIDE3d1(ilower,iupper,N_gcw))
      REAL N11(SIDE3d1(ilower,iupper,N_gcw))
      REAL N12(SIDE3d1(ilower,iupper,N_gcw))
      REAL N20(SIDE3d2(ilower,iupper,N_gcw))
      REAL N21(SIDE3d2(ilower,iupper,N_gcw))
      REAL N22(SIDE3d2(ilower,iupper,N_gcw))
      REAL K(CELL3d(ilower,iupper,K_gcw))
     
      REAL dx(0:NDIM-1)

c
c     Local variables.
c
      INTEGER i0,i1,i2
      REAL fac0,fac1,fac2
      REAL norm_grad_upper,norm_grad_lower,eps
 
      fac0 = 1.d0/(dx(0))
      fac1 = 1.d0/(dx(1))   
      fac2 = 1.d0/(dx(2))   
      eps = 1.d-10

c
c     Compute curvature K = -div (n/|n|) 
c     
      do i2 = ilower2 - 1, iupper2 + 1
        do i1 = ilower1 - 1, iupper1 + 1
          do i0 = ilower0 - 1, iupper0 + 1
            
c           compute -d/dx term.

            norm_grad_upper = sqrt(N00(i0+1,i1,i2)**2 +
     &                             N01(i0+1,i1,i2)**2 +
     &                             N02(i0+1,i1,i2)**2)
            norm_grad_lower = sqrt(N00(i0,i1,i2)**2 +
     &                             N01(i0,i1,i2)**2 +
     &                             N02(i0,i1,i2)**2)
            
            if (norm_grad_upper .gt. eps) then
                norm_grad_upper = 1.d0/norm_grad_upper
            else
                norm_grad_upper = 0.d0
            endif

            if (norm_grad_lower .gt. eps) then
                norm_grad_lower = 1.d0/norm_grad_lower
            else
                norm_grad_lower = 0.d0
            endif

c           Compute -div
            K(i0,i1,i2) = fac0*(N00(i0,i1,i2)*norm_grad_lower - 
     &              N00(i0+1,i1,i2)*norm_grad_upper)
            

c           compute -d/dy term.

            norm_grad_upper = sqrt(N10(i0,i1+1,i2)**2 +
     &                             N11(i0,i1+1,i2)**2 +
     &                             N12(i0,i1+1,i2)**2)
            norm_grad_lower = sqrt(N10(i0,i1,i2)**2 + 
     &                             N11(i0,i1,i2)**2 +
     &                             N12(i0,i1,i2)**2)

            if (norm_grad_upper .gt. eps) then
                norm_grad_upper = 1.d0/norm_grad_upper
            else
                norm_grad_upper = 0.d0
            endif
            
            if (norm_grad_lower .gt. eps) then
                norm_grad_lower = 1.d0/norm_grad_lower
            else
                norm_grad_lower = 0.d0
            endif    
              
c           Compute -div
            K(i0,i1,i2) = K(i0,i1,i2) + fac1*(N11(i0,i1,i2)
     &         * norm_grad_lower - N11(i0,i1+1,i2)*norm_grad_upper)

c           compute -d/dz term.

            norm_grad_upper = sqrt(N20(i0,i1,i2+1)**2 +
     &                             N21(i0,i1,i2+1)**2 +
     &                             N22(i0,i1,i2+1)**2)
            norm_grad_lower = sqrt(N20(i0,i1,i2)**2 + 
     &                             N21(i0,i1,i2)**2 +
     &                             N22(i0,i1,i2)**2)

            if (norm_grad_upper .gt. eps) then
                norm_grad_upper = 1.d0/norm_grad_upper
            else
                norm_grad_upper = 0.d0
            endif
            
            if (norm_grad_lower .gt. eps) then
                norm_grad_lower = 1.d0/norm_grad_lower
            else
                norm_grad_lower = 0.d0
            endif    
              
c           Compute -div
            K(i0,i1,i2) = K(i0,i1,i2) + fac2*(N22(i0,i1,i2)
     &         * norm_grad_lower - N22(i0,i1,i2+1)*norm_grad_upper)           

          enddo
        enddo
      enddo
      
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute surface tension forcing
c
c           F = sigma * K * grad(f)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine sc_surface_tension_force_3d(
     &     F0,F1,F2,F_gcw,
     &     K,K_gcw,
     &     N00,N11,N22,N_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
     &     sigma)
c
      implicit none
c
c     Input.
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2
      INTEGER F_gcw,K_gcw,N_gcw
    
c
c     Input/Output.
c
      REAL F0(SIDE3d0(ilower,iupper,F_gcw))
      REAL F1(SIDE3d1(ilower,iupper,F_gcw))
      REAL F2(SIDE3d2(ilower,iupper,F_gcw))

      REAL N00(SIDE3d0(ilower,iupper,N_gcw))
      REAL N11(SIDE3d1(ilower,iupper,N_gcw))
      REAL N22(SIDE3d2(ilower,iupper,N_gcw))
      REAL K(CELL3d(ilower,iupper,K_gcw))
     
      REAL sigma

c
c     Local variables.
c
      INTEGER i0,i1,i2
      REAL kappa
  
c
c     Compute F0  = sigma * K_x * N_x 
c     
      do i2 = ilower2, iupper2
        do i1 = ilower1, iupper1
          do i0 = ilower0, iupper0 + 1
            
              kappa = 0.5d0*(K(i0-1,i1,i2)+K(i0,i1,i2))
              F0(i0,i1,i2) = sigma*kappa*N00(i0,i1,i2)
            
          enddo
        enddo
      enddo

c
c     Compute F1  = sigma * K_y * N_y
c
      do i2 = ilower2, iupper2
        do i1 = ilower1, iupper1 + 1
          do i0 = ilower0, iupper0

              kappa = 0.5d0*(K(i0,i1-1,i2)+K(i0,i1,i2))
              F1(i0,i1,i2) = sigma*kappa*N11(i0,i1,i2)
                    
          enddo
        enddo
      enddo

c
c     Compute F2  = sigma * K_z * N_z
c
      do i2 = ilower2, iupper2 + 1
        do i1 = ilower1, iupper1
          do i0 = ilower0, iupper0

              kappa = 0.5d0*(K(i0,i1,i2-1)+K(i0,i1,i2))
              F2(i0,i1,i2) = sigma*kappa*N22(i0,i1,i2)
                    
          enddo
        enddo
      enddo
      
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute Marangoni force
c
c           F = grad T |grad C| - (grad T dot grad \phi) grad C
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine sc_marangoni_force_3d(
     &     F0,F1,F2,F_gcw,
     &     gradT00,gradT01,gradT02,
     &     gradT10,gradT11,gradT12,
     &     gradT20,gradT21,gradT22,
     &     gradT_gcw,
     &     N00,N01,N02,
     &     N10,N11,N12,
     &     N20,N21,N22,
     &     N_gcw,
     &     gradC00,gradC01,gradC02,
     &     gradC10,gradC11,gradC12,
     &     gradC20,gradC21,gradC22,
     &     gradC_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
     &     marangoni_coefficient)
c
      implicit none
c
c     Input.
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2
      INTEGER F_gcw,gradT_gcw,N_gcw,gradC_gcw

c
c     Input/Output.
c
      REAL F0(SIDE3d0(ilower,iupper,F_gcw))
      REAL F1(SIDE3d1(ilower,iupper,F_gcw))
      REAL F2(SIDE3d2(ilower,iupper,F_gcw))
      REAL gradT00(SIDE3d0(ilower,iupper,gradT_gcw))
      REAL gradT01(SIDE3d0(ilower,iupper,gradT_gcw))
      REAL gradT02(SIDE3d0(ilower,iupper,gradT_gcw))
      REAL gradT10(SIDE3d1(ilower,iupper,gradT_gcw))
      REAL gradT11(SIDE3d1(ilower,iupper,gradT_gcw))
      REAL gradT12(SIDE3d1(ilower,iupper,gradT_gcw))
      REAL gradT20(SIDE3d2(ilower,iupper,gradT_gcw))
      REAL gradT21(SIDE3d2(ilower,iupper,gradT_gcw))
      REAL gradT22(SIDE3d2(ilower,iupper,gradT_gcw))
      REAL N00(SIDE3d0(ilower,iupper,N_gcw))
      REAL N01(SIDE3d0(ilower,iupper,N_gcw))
      REAL N02(SIDE3d0(ilower,iupper,N_gcw))
      REAL N10(SIDE3d1(ilower,iupper,N_gcw))
      REAL N11(SIDE3d1(ilower,iupper,N_gcw))
      REAL N12(SIDE3d1(ilower,iupper,N_gcw))
      REAL N20(SIDE3d2(ilower,iupper,N_gcw))
      REAL N21(SIDE3d2(ilower,iupper,N_gcw))
      REAL N22(SIDE3d2(ilower,iupper,N_gcw))
      REAL gradC00(SIDE3d0(ilower,iupper,gradC_gcw))
      REAL gradC01(SIDE3d0(ilower,iupper,gradC_gcw))
      REAL gradC02(SIDE3d0(ilower,iupper,gradC_gcw))
      REAL gradC10(SIDE3d1(ilower,iupper,gradC_gcw))
      REAL gradC11(SIDE3d1(ilower,iupper,gradC_gcw))
      REAL gradC12(SIDE3d1(ilower,iupper,gradC_gcw))
      REAL gradC20(SIDE3d2(ilower,iupper,gradC_gcw))
      REAL gradC21(SIDE3d2(ilower,iupper,gradC_gcw))
      REAL gradC22(SIDE3d2(ilower,iupper,gradC_gcw))

      REAL marangoni_coefficient

c
c     Local variables.
c
      INTEGER i0,i1,i2
      REAL gradC_mag
      REAL gradT_dot_gradphi

c
      do i2 = ilower2, iupper2
        do i1 = ilower1, iupper1
          do i0 = ilower0, iupper0 + 1
            gradC_mag = sqrt(gradC00(i0,i1,i2)**2+gradC01(i0,i1,i2)**2
     &                  +gradC02(i0,i1,i2)**2)
            gradT_dot_gradphi = gradT00(i0,i1,i2)*N00(i0,i1,i2) +
     &                          gradT01(i0,i1,i2)*N01(i0,i1,i2) +
     &                          gradT02(i0,i1,i2)*N02(i0,i1,i2)
            F0(i0,i1,i2)=F0(i0,i1,i2)+marangoni_coefficient
     &                  *(gradT00(i0,i1,i2)*gradC_mag-gradT_dot_gradphi
     &                  *gradC00(i0,i1,i2))
          enddo
        enddo
      enddo
c
      do i2 = ilower2, iupper2
        do i1 = ilower1, iupper1+1
          do i0 = ilower0, iupper0
            gradC_mag = sqrt(gradC10(i0,i1,i2)**2+gradC11(i0,i1,i2)**2
     &                + gradC12(i0,i1,i2)**2)
            gradT_dot_gradphi = gradT10(i0,i1,i2)*N10(i0,i1,i2) +
     &                          gradT11(i0,i1,i2)*N11(i0,i1,i2) +
     &                          gradT12(i0,i1,i2)*N12(i0,i1,i2)
            F1(i0,i1,i2)=F1(i0,i1,i2)+marangoni_coefficient
     &                  *(gradT11(i0,i1,i2)*gradC_mag-gradT_dot_gradphi
     &                  *gradC11(i0,i1,i2))
          enddo
        enddo
      enddo
c
      do i2 = ilower2, iupper2+1
        do i1 = ilower1, iupper1
          do i0 = ilower0, iupper0
            gradC_mag = sqrt(gradC20(i0,i1,i2)**2+gradC21(i0,i1,i2)**2
     &                + gradC22(i0,i1,i2)**2)
            gradT_dot_gradphi = gradT20(i0,i1,i2)*N20(i0,i1,i2) +
     &                          gradT21(i0,i1,i2)*N21(i0,i1,i2) +
     &                          gradT22(i0,i1,i2)*N22(i0,i1,i2)
            F2(i0,i1,i2)=F2(i0,i1,i2)+marangoni_coefficient
     &                  *(gradT22(i0,i1,i2)*gradC_mag-gradT_dot_gradphi
     &                  *gradC22(i0,i1,i2))
          enddo
        enddo
      enddo
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute temperature dependence surface tension
c
c     F = marangoni_coef (T - T0) kappa grad C
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine sc_temperature_surface_tension_force_3d(
     &     F0,F1,F2,F_gcw,
     &     T,T_gcw,
     &     K,K_gcw,
     &     gradC00,gradC11,gradC22,
     &     gradC_gcw,
     &     ref_temperature,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
     &     marangoni_coefficient)
c
      implicit none
c
c     Input.
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2
      INTEGER F_gcw,T_gcw,K_gcw,gradC_gcw
c
c     Input/Output.
c
      REAL F0(SIDE3d0(ilower,iupper,F_gcw))
      REAL F1(SIDE3d1(ilower,iupper,F_gcw))
      REAL F2(SIDE3d2(ilower,iupper,F_gcw))
      REAL T(CELL3d(ilower,iupper,T_gcw))
      REAL K(CELL3d(ilower,iupper,K_gcw))
      REAL gradC00(SIDE3d0(ilower,iupper,gradC_gcw))
      REAL gradC11(SIDE3d1(ilower,iupper,gradC_gcw))
      REAL gradC22(SIDE3d2(ilower,iupper,gradC_gcw))
c
      REAL ref_temperature
      REAL marangoni_coefficient
c
c     Local variables.
c
      INTEGER i0,i1,i2
      REAL kappa, Temperature
c
      do i2 = ilower2, iupper2
        do i1 = ilower1, iupper1
         do i0 = ilower0, iupper0 + 1
            kappa = 0.5d0*(K(i0-1,i1,i2)+K(i0,i1,i2))
            Temperature = 0.5d0*(T(i0-1,i1,i2)+T(i0,i1,i2))
            F0(i0,i1,i2) = F0(i0,i1,i2)+marangoni_coefficient*
     &                   (Temperature - ref_temperature)*kappa
     &                   *gradC00(i0,i1,i2)
         enddo
       enddo
      enddo
c
      do i2 = ilower2, iupper2
        do i1 = ilower1, iupper1 + 1
         do i0 = ilower0, iupper0
            kappa = 0.5d0*(K(i0,i1-1,i2)+K(i0,i1,i2))
            Temperature = 0.5d0*(T(i0,i1-1,i2)+T(i0,i1,i2))
            F1(i0,i1,i2) = F1(i0,i1,i2)+marangoni_coefficient*
     &                   (Temperature - ref_temperature)*kappa
     &                    *gradC11(i0,i1,i2)
         enddo
       enddo
      enddo
c
      do i2 = ilower2, iupper2 + 1
         do i1 = ilower1, iupper1
          do i0 = ilower0, iupper0
             kappa = 0.5d0*(K(i0,i1,i2-1)+K(i0,i1,i2))
             Temperature = 0.5d0*(T(i0,i1,i2-1)+T(i0,i1,i2))
             F2(i0,i1,i2) = F2(i0,i1,i2)+marangoni_coefficient*
     &                   (Temperature - ref_temperature)*kappa
     &                    *gradC22(i0,i1,i2)
          enddo
        enddo
      enddo
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute Marangoni force
c     marangoni_coef is a function of T based on the work of Sahoo 1988
c     F = grad T |grad C| - (grad T dot grad \phi) grad C
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine sc_variable_coef_marangoni_force_3d(
     &     F0,F1,F2,F_gcw,
     &     gradT00,gradT01,gradT02,
     &     gradT10,gradT11,gradT12,
     &     gradT20,gradT21,gradT22,
     &     gradT_gcw,
     &     N00,N01,N02,
     &     N10,N11,N12,
     &     N20,N21,N22,
     &     N_gcw,
     &     gradC00,gradC01,gradC02,
     &     gradC10,gradC11,gradC12,
     &     gradC20,gradC21,gradC22,
     &     gradC_gcw,
     &     T,T_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
     &     marangoni_coefficient)
c
      implicit none
c
c     Input.
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2
      INTEGER F_gcw,gradT_gcw,N_gcw,gradC_gcw,T_gcw

c
c     Input/Output.
c
      REAL F0(SIDE3d0(ilower,iupper,F_gcw))
      REAL F1(SIDE3d1(ilower,iupper,F_gcw))
      REAL F2(SIDE3d2(ilower,iupper,F_gcw))
      REAL gradT00(SIDE3d0(ilower,iupper,gradT_gcw))
      REAL gradT01(SIDE3d0(ilower,iupper,gradT_gcw))
      REAL gradT02(SIDE3d0(ilower,iupper,gradT_gcw))
      REAL gradT10(SIDE3d1(ilower,iupper,gradT_gcw))
      REAL gradT11(SIDE3d1(ilower,iupper,gradT_gcw))
      REAL gradT12(SIDE3d1(ilower,iupper,gradT_gcw))
      REAL gradT20(SIDE3d2(ilower,iupper,gradT_gcw))
      REAL gradT21(SIDE3d2(ilower,iupper,gradT_gcw))
      REAL gradT22(SIDE3d2(ilower,iupper,gradT_gcw))
      REAL N00(SIDE3d0(ilower,iupper,N_gcw))
      REAL N01(SIDE3d0(ilower,iupper,N_gcw))
      REAL N02(SIDE3d0(ilower,iupper,N_gcw))
      REAL N10(SIDE3d1(ilower,iupper,N_gcw))
      REAL N11(SIDE3d1(ilower,iupper,N_gcw))
      REAL N12(SIDE3d1(ilower,iupper,N_gcw))
      REAL N20(SIDE3d2(ilower,iupper,N_gcw))
      REAL N21(SIDE3d2(ilower,iupper,N_gcw))
      REAL N22(SIDE3d2(ilower,iupper,N_gcw))
      REAL gradC00(SIDE3d0(ilower,iupper,gradC_gcw))
      REAL gradC01(SIDE3d0(ilower,iupper,gradC_gcw))
      REAL gradC02(SIDE3d0(ilower,iupper,gradC_gcw))
      REAL gradC10(SIDE3d1(ilower,iupper,gradC_gcw))
      REAL gradC11(SIDE3d1(ilower,iupper,gradC_gcw))
      REAL gradC12(SIDE3d1(ilower,iupper,gradC_gcw))
      REAL gradC20(SIDE3d2(ilower,iupper,gradC_gcw))
      REAL gradC21(SIDE3d2(ilower,iupper,gradC_gcw))
      REAL gradC22(SIDE3d2(ilower,iupper,gradC_gcw))
      REAL T(CELL3d(ilower,iupper,T_gcw))

      REAL marangoni_coefficient

c
c     Local variables.
c
      INTEGER i0,i1,i2
      REAL gradC_mag
      REAL gradT_dot_gradphi
      REAL R,ai,gamma_s,ki,deltaH0,deltaHi,K,dsigma_dT0
      REAL T_sc
      R=8314
      ai=0.015
      gamma_s=1.3E-8
      ki=3.18E-3
      deltaH0=-1.66E8
      deltaHi = 0.5d0
      dsigma_dT0 = -5E-4

c
      do i2 = ilower2, iupper2
        do i1 = ilower1, iupper1
          do i0 = ilower0, iupper0 + 1
            T_sc = 0.5d0*(T(i0-1,i1,i2)+T(i0,i1,i2))
            K = ki*EXP(-deltaH0/(R*T_sc))
            marangoni_coefficient = dsigma_dT0-R*gamma_s*log(1+K*ai)
     &       - K*ai*gamma_s*(deltaH0 - deltaHi)/(T_sc*(1+K*ai))
            gradC_mag = sqrt(gradC00(i0,i1,i2)**2+gradC01(i0,i1,i2)**2
     &                  +gradC02(i0,i1,i2)**2)
            gradT_dot_gradphi = gradT00(i0,i1,i2)*N00(i0,i1,i2) +
     &                          gradT01(i0,i1,i2)*N01(i0,i1,i2) +
     &                          gradT02(i0,i1,i2)*N02(i0,i1,i2)
            F0(i0,i1,i2)=F0(i0,i1,i2)+marangoni_coefficient
     &                  *(gradT00(i0,i1,i2)*gradC_mag-gradT_dot_gradphi
     &                  *gradC00(i0,i1,i2))
          enddo
        enddo
      enddo
c
      do i2 = ilower2, iupper2
        do i1 = ilower1, iupper1+1
          do i0 = ilower0, iupper0
            T_sc = 0.5d0*(T(i0,i1-1,i2)+T(i0,i1,i2))
            K = ki*EXP(-deltaH0/(R*T_sc))
            marangoni_coefficient = dsigma_dT0-R*gamma_s*log(1+K*ai)
     &       - K*ai*gamma_s*(deltaH0 - deltaHi)/(T_sc*(1+K*ai))
            gradC_mag = sqrt(gradC10(i0,i1,i2)**2+gradC11(i0,i1,i2)**2
     &                + gradC12(i0,i1,i2)**2)
            gradT_dot_gradphi = gradT10(i0,i1,i2)*N10(i0,i1,i2) +
     &                          gradT11(i0,i1,i2)*N11(i0,i1,i2) +
     &                          gradT12(i0,i1,i2)*N12(i0,i1,i2)
            F1(i0,i1,i2)=F1(i0,i1,i2)+marangoni_coefficient
     &                  *(gradT11(i0,i1,i2)*gradC_mag-gradT_dot_gradphi
     &                  *gradC11(i0,i1,i2))
          enddo
        enddo
      enddo
c
      do i2 = ilower2, iupper2+1
        do i1 = ilower1, iupper1
          do i0 = ilower0, iupper0
            T_sc = 0.5d0*(T(i0,i1,i2-1)+T(i0,i1,i2))
            K = ki*EXP(-deltaH0/(R*T_sc))
            marangoni_coefficient = dsigma_dT0-R*gamma_s*log(1+K*ai)
     &       - K*ai*gamma_s*(deltaH0 - deltaHi)/(T_sc*(1+K*ai))
            gradC_mag = sqrt(gradC20(i0,i1,i2)**2+gradC21(i0,i1,i2)**2
     &                + gradC22(i0,i1,i2)**2)
            gradT_dot_gradphi = gradT20(i0,i1,i2)*N20(i0,i1,i2) +
     &                          gradT21(i0,i1,i2)*N21(i0,i1,i2) +
     &                          gradT22(i0,i1,i2)*N22(i0,i1,i2)
            F2(i0,i1,i2)=F2(i0,i1,i2)+marangoni_coefficient
     &                  *(gradT22(i0,i1,i2)*gradC_mag-gradT_dot_gradphi
     &                  *gradC22(i0,i1,i2))
          enddo
        enddo
      enddo
      return
      end



