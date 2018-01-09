c
c     Routines to compute continuum surface tension forces.
c
c     Created on 25 Dec 2017 by Nishant Nangia
c
c     Copyright (c) 2002-2017, Amneet Bhalla
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




