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

define(NDIM,2)dnl
define(REAL,`double precision')dnl
define(INTEGER,`integer')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim2d.i)dnl

c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Mollify indicator function using IB_4 kernel
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine mollify_ib_4_2d(
     &     V,V_gcw,
     &     U,U_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1)
c
      implicit none
include(TOP_SRCDIR/src/fortran/const.i)dnl
c
c     Input.
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER V_gcw,U_gcw
    
c
c     Input/Output.
c
      REAL U(CELL2d(ilower,iupper,U_gcw))
      REAL V(CELL2d(ilower,iupper,V_gcw))
      REAL w0(-1:1),w1(-1:1),wy
      REAL w(-1:1,-1:1)

c
c     Local variables.
c
      INTEGER k0,k1
      INTEGER i0,i1
      

c     Compute 1D weights.
      w0(-1) = fourth; w1(-1) = fourth
      w0(0) = half; w1(0) = half
      w0(1) = fourth; w1(1) = fourth

c     Compute the tensor product weight
      do k1 = -1,1
         wy = w1(k1)
         do k0 = -1,1
            w(k0,k1) = w0(k0)*wy
         enddo
      enddo  
    
c     Mollify U to V.
      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0

            V(i0,i1) = 0.d0
            do k1 = -1,1
               do k0 = -1,1
               V(i0,i1) = V(i0,i1) + U(i0+k0,i1+k1)*w(k0,k1)
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
      subroutine sc_normal_2d(
     &     N00,N01,
     &     N10,N11,
     &     N_gcw,
     &     U,U_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     dx)
c
      implicit none
include(TOP_SRCDIR/src/fortran/const.i)dnl
c
c     Input.
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER N_gcw,U_gcw
    
c
c     Input/Output.
c
      REAL N00(SIDE2d0(ilower,iupper,N_gcw))
      REAL N01(SIDE2d0(ilower,iupper,N_gcw))
      REAL N10(SIDE2d1(ilower,iupper,N_gcw))
      REAL N11(SIDE2d1(ilower,iupper,N_gcw))
      REAL U(CELL2d(ilower,iupper,U_gcw))
     
      REAL dx(0:NDIM-1)

c
c     Local variables.
c
      INTEGER i0,i1
      REAL fac0,fac1
      
      fac0 = 1.d0/(dx(0))
      fac1 = 1.d0/(dx(1))      

c
c     Find face normal gradients first and then interpolate 
c     face tangential gradients

c     Do N00.
      do i1 = ilower1 - 2, iupper1 + 2
         do i0 = ilower0 - 1, iupper0 + 2
             
            N00(i0,i1) = fac0*(U(i0,i1) - U(i0-1,i1))

         enddo
      enddo
      
c     Do N11.
      do i1 = ilower1 - 1, iupper1 + 2
         do i0 = ilower0 - 2, iupper0 + 2
             
            N11(i0,i1) = fac1*(U(i0,i1) - U(i0,i1-1))

         enddo
      enddo

c     Interpolate N11 to N01
      do i1 = ilower1 - 1, iupper1 + 1
         do i0 = ilower0 - 1, iupper0 + 2
             
            N01(i0,i1) = fourth*(N11(i0-1,i1) + N11(i0,i1) + 
     &                    N11(i0-1,i1+1) + N11(i0,i1+1)) 

         enddo
      enddo

c     Interpolate N00 to N10
      do i1 = ilower1 - 1, iupper1 + 2
         do i0 = ilower0 - 1, iupper0 + 1

            N10(i0,i1) = fourth*(N00(i0,i1) + N00(i0+1,i1) +
     &                   N00(i0,i1-1) + N00(i0+1,i1-1))

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
      subroutine cc_curvature_2d(
     &     K,K_gcw,
     &     N00,N01,
     &     N10,N11,
     &     N_gcw,
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
      INTEGER K_gcw,N_gcw
    
c
c     Input/Output.
c
      REAL N00(SIDE2d0(ilower,iupper,N_gcw))
      REAL N01(SIDE2d0(ilower,iupper,N_gcw))
      REAL N10(SIDE2d1(ilower,iupper,N_gcw))
      REAL N11(SIDE2d1(ilower,iupper,N_gcw))
      REAL K(CELL2d(ilower,iupper,K_gcw))
     
      REAL dx(0:NDIM-1)

c
c     Local variables.
c
      INTEGER i0,i1
      REAL fac0,fac1
      REAL norm_grad_upper,norm_grad_lower,eps
 
      fac0 = 1.d0/(dx(0))
      fac1 = 1.d0/(dx(1))      
      eps = 1.d-10

c
c     Compute curvature K = -div (n/|n|) 
c     
      do i1 = ilower1 - 1, iupper1 + 1
         do i0 = ilower0 - 1, iupper0 + 1
            
c           compute -d/dx term.

            norm_grad_upper = sqrt(N00(i0+1,i1)**2+N01(i0+1,i1)**2)
            norm_grad_lower = sqrt(N00(i0,i1)**2+N01(i0,i1)**2)
            
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
            K(i0,i1) = fac0*(N00(i0,i1)*norm_grad_lower - 
     &              N00(i0+1,i1)*norm_grad_upper)
            

c           compute -d/dy term.

            norm_grad_upper = sqrt(N10(i0,i1+1)**2+N11(i0,i1+1)**2)
            norm_grad_lower = sqrt(N10(i0,i1)**2+N11(i0,i1)**2)

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
            K(i0,i1) = K(i0,i1) + fac1*(N11(i0,i1)*norm_grad_lower -
     &              N11(i0,i1+1)*norm_grad_upper)          

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
      subroutine sc_surface_tension_force_2d(
     &     F0,F1,F_gcw,
     &     K,K_gcw,
     &     N00,N11,N_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     sigma)
c
      implicit none
c
c     Input.
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER F_gcw,K_gcw,N_gcw
    
c
c     Input/Output.
c
      REAL F0(SIDE2d0(ilower,iupper,F_gcw))
      REAL F1(SIDE2d1(ilower,iupper,F_gcw))
      REAL N00(SIDE2d0(ilower,iupper,N_gcw))
      REAL N11(SIDE2d1(ilower,iupper,N_gcw))
      REAL K(CELL2d(ilower,iupper,K_gcw))
     
      REAL sigma

c
c     Local variables.
c
      INTEGER i0,i1
      REAL kappa
  
c
c     Compute F0  = sigma * K_x * N_x 
c     
      do i1 = ilower1, iupper1
         do i0 = ilower0, iupper0 + 1
            
            kappa = 0.5d0*(K(i0-1,i1)+K(i0,i1))
            F0(i0,i1) = sigma*kappa*N00(i0,i1)
            
         enddo
      enddo

c
c     Compute F1  = sigma * K_y * N_y
c
      do i1 = ilower1, iupper1 + 1
         do i0 = ilower0, iupper0

            kappa = 0.5d0*(K(i0,i1-1)+K(i0,i1))
            F1(i0,i1) = sigma*kappa*N11(i0,i1)
                    
         enddo
      enddo
      
      return
      end




