c ---------------------------------------------------------------------
c
c Copyright (c) 2017 - 2023 by the IBAMR developers
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
      do i1 = ilower1 - N_gcw, iupper1 + N_gcw
         do i0 = ilower0 + 1 - N_gcw, iupper0 + N_gcw
             
            N00(i0,i1) = fac0*(U(i0,i1) - U(i0-1,i1))

         enddo
      enddo
      
c     Do N11.
      do i1 = ilower1 + 1 - N_gcw, iupper1 + N_gcw
         do i0 = ilower0 - N_gcw, iupper0 + N_gcw
             
            N11(i0,i1) = fac1*(U(i0,i1) - U(i0,i1-1))

         enddo
      enddo

c     Interpolate N11 to N01
      do i1 = ilower1 + 1 - N_gcw, iupper1 + N_gcw - 1
         do i0 = ilower0 + 1 - N_gcw, iupper0 + N_gcw
             
            N01(i0,i1) = fourth*(N11(i0-1,i1) + N11(i0,i1) + 
     &                    N11(i0-1,i1+1) + N11(i0,i1+1)) 

         enddo
      enddo

c     Interpolate N00 to N10
      do i1 = ilower1 + 1 - N_gcw, iupper1 + N_gcw
         do i0 = ilower0 + 1 - N_gcw, iupper0 + N_gcw - 1

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
c           F = K * grad(f)
c
c     It is the caller's responsibility to multiply F by the (possibly
c     spatially-dependent) sigma coefficient.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine sc_surface_tension_force_2d(
     &     F0,F1,F_gcw,
     &     K,K_gcw,
     &     N00,N11,N_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1)
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

c
c     Local variables.
c
      INTEGER i0,i1
      REAL kappa
  
c
c     Compute F0  = K_x * N_x
c     
      do i1 = ilower1, iupper1
         do i0 = ilower0, iupper0 + 1
            
            kappa = 0.5d0*(K(i0-1,i1)+K(i0,i1))
            F0(i0,i1) = kappa*N00(i0,i1)
            
         enddo
      enddo

c
c     Compute F1  = K_y * N_y
c
      do i1 = ilower1, iupper1 + 1
         do i0 = ilower0, iupper0

            kappa = 0.5d0*(K(i0,i1-1)+K(i0,i1))
            F1(i0,i1) = kappa*N11(i0,i1)
                    
         enddo
      enddo
      
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute Marangoni force
c
c     F = grad T |grad C| - (grad T dot grad \phi)grad C
c
c     It is the caller's responsibility to multiply F by the (possibly
c     spatially-dependent) Marangoni coefficient.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine sc_marangoni_force_2d(
     &     F0,F1,F_gcw,
     &     gradT00,gradT01,
     &     gradT10,gradT11,
     &     gradT_gcw,
     &     N00,N01,
     &     N10,N11,
     &     N_gcw,
     &     gradC00,gradC01,
     &     gradC10,gradC11,
     &     gradC_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1)
c
      implicit none
c
c     Input.
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER F_gcw,gradT_gcw,N_gcw,gradC_gcw
c
c     Input/Output.
c
      REAL F0(SIDE2d0(ilower,iupper,F_gcw))
      REAL F1(SIDE2d1(ilower,iupper,F_gcw))
      REAL gradT00(SIDE2d0(ilower,iupper,gradT_gcw))
      REAL gradT01(SIDE2d0(ilower,iupper,gradT_gcw))
      REAL gradT10(SIDE2d1(ilower,iupper,gradT_gcw))
      REAL gradT11(SIDE2d1(ilower,iupper,gradT_gcw))
      REAL N00(SIDE2d0(ilower,iupper,N_gcw))
      REAL N01(SIDE2d0(ilower,iupper,N_gcw))
      REAL N10(SIDE2d1(ilower,iupper,N_gcw))
      REAL N11(SIDE2d1(ilower,iupper,N_gcw))
      REAL gradC00(SIDE2d0(ilower,iupper,gradC_gcw))
      REAL gradC01(SIDE2d0(ilower,iupper,gradC_gcw))
      REAL gradC10(SIDE2d1(ilower,iupper,gradC_gcw))
      REAL gradC11(SIDE2d1(ilower,iupper,gradC_gcw))
c
c     Local variables.
c
      INTEGER i0,i1
      REAL gradC_mag
      REAL gradT_dot_gradphi
      REAL norm_grad,eps
      eps = 1.d-10
c
      do i1 = ilower1, iupper1
        do i0 = ilower0, iupper0 + 1
          norm_grad = sqrt(N00(i0,i1)**2 + N01(i0,i1)**2)
          if (norm_grad .gt. eps) then
              norm_grad = 1.d0/norm_grad
          else
              norm_grad = 0.d0
          endif
          N00(i0,i1) = N00(i0,i1)*norm_grad
          N01(i0,i1) = N01(i0,i1)*norm_grad

          gradC_mag = sqrt(gradC00(i0,i1)**2 + gradC01(i0,i1)**2)
          gradT_dot_gradphi = gradT00(i0,i1)*N00(i0,i1) +
     &                          gradT01(i0,i1)*N01(i0,i1)
          F0(i0,i1) = gradT00(i0,i1)
     &                *gradC_mag-gradT_dot_gradphi*gradC00(i0,i1)
        enddo
      enddo
c
      do i1 = ilower1, iupper1 + 1
        do i0 = ilower0, iupper0
          norm_grad = sqrt(N10(i0,i1)**2 + N11(i0,i1)**2)
          if (norm_grad .gt. eps) then
              norm_grad = 1.d0/norm_grad
          else
              norm_grad = 0.d0
          endif
          N10(i0,i1) = N10(i0,i1)*norm_grad
          N11(i0,i1) = N11(i0,i1)*norm_grad

          gradC_mag = sqrt(gradC10(i0,i1)**2 + gradC11(i0,i1)**2)
          gradT_dot_gradphi = gradT10(i0,i1)*N10(i0,i1) +
     &                          gradT11(i0,i1)*N11(i0,i1)
          F1(i0,i1) = gradT11(i0,i1)
     &                *gradC_mag-gradT_dot_gradphi*gradC11(i0,i1)
        enddo
      enddo
      return
      end

c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Height-function curvature for a 2-D Cartesian grid.
c
c     A fixed 3 x 7 or 7 x 3 stencil is used.
c
c     K_valid = 1 when height-function curvature is available.
c     K_valid = 0 otherwise.
c
c     Curvature convention:
c
c              kappa = -div(grad(C)/|grad(C)|)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine hf_curvature_2d(
     &     K,K_valid,K_gcw,
     &     C,C_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     dx,hf_radius,vof_tol,gradient_tol)

      implicit none

c
c     Patch information.
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER K_gcw,C_gcw
      INTEGER hf_radius

c
c     Data arrays.
c
      REAL K(CELL2d(ilower,iupper,K_gcw))
      REAL K_valid(CELL2d(ilower,iupper,K_gcw))
      REAL C(CELL2d(ilower,iupper,C_gcw))

c
c     Grid spacing and tolerances.
c
      REAL dx(0:NDIM-1)
      REAL vof_tol
      REAL gradient_tol

c
c     Local variables.
c
      INTEGER i0,i1
      INTEGER m,kk
      INTEGER all_ok
      INTEGER ok

      REAL alpha
      REAL gx,gy,gmag
      REAL h(-1:1)
      REAL clow,chigh
      REAL hd,hdd
      REAL denom

c
c     This implementation is specifically the standard 3 x 7
c     height function.
c
      
      if (hf_radius .ne. 3) then
         write(*,*) 'hf_curvature_2d requires hf_radius = 3'
         stop
      endif

      do i1 = ilower1-1,iupper1+1
         do i0 = ilower0-1,iupper0+1

            K(i0,i1) = 0.d0
            K_valid(i0,i1) = 0.d0
            alpha = C(i0,i1)

            if (alpha .gt. vof_tol .and.
     &          alpha .lt. 1.d0-vof_tol) then

c
c              Estimate the interface normal.
c
               gx = (C(i0+1,i1)-C(i0-1,i1))
     &              /(2.d0*dx(0))

               gy = (C(i0,i1+1)-C(i0,i1-1))
     &              /(2.d0*dx(1))

               gmag = sqrt(gx*gx+gy*gy)

               if (gmag .gt. gradient_tol) then

                  if (abs(gy) .ge. abs(gx)) then

                     all_ok = 1

                     do m = -1,1

                        clow = max(0.d0,
     &                       min(1.d0,C(i0+m,i1-3)))

                        chigh = max(0.d0,
     &                        min(1.d0,C(i0+m,i1+3)))

                        ok = 0

c
c                       grad(C) points from C=0 toward C=1.
c
                        if (gy .gt. 0.d0) then

                           if (clow .le. vof_tol .and.
     &                         chigh .ge.
     &                         1.d0-vof_tol) then
                              ok = 1
                           endif

                        else

c
c                          gy < 0: opposite phase ordering.
c
                           if (clow .ge.
     &                         1.d0-vof_tol .and.
     &                         chigh .le. vof_tol) then
                              ok = 1
                           endif

                        endif

                        if (ok .eq. 0) all_ok = 0

c
c                       Construct the physical height.
c
                        h(m) = 0.d0

                        do kk = -3,3

                           h(m) = h(m)
     &                          + dx(1)*max(0.d0,
     &                            min(1.d0,
     &                            C(i0+m,i1+kk)))

                        enddo

                     enddo

c
c                    All three columns must have valid pure ends.
c
                     if (all_ok .eq. 1) then

c
c                       h varies in the x direction.
c
                        hd = (h(1)-h(-1))
     &                       /(2.d0*dx(0))

                        hdd = (h(1)-2.d0*h(0)+h(-1))
     &                        /(dx(0)*dx(0))

                        denom = (1.d0+hd*hd)**1.5d0

                        K(i0,i1) = -hdd/denom
                        K_valid(i0,i1) = 1.d0

                     endif

c
c                 |gx| > |gy| now check for the x cells
c
                  else 

                     all_ok = 1

                     do m = -1,1

c
c                       Fixed left and right ends.
c
                        clow = max(0.d0,
     &                       min(1.d0,C(i0-3,i1+m)))

                        chigh = max(0.d0,
     &                        min(1.d0,C(i0+3,i1+m)))

                        ok = 0

c
c                       If gx > 0, C increases from left to right.
c
                        if (gx .gt. 0.d0) then

                           if (clow .le. vof_tol .and.
     &                         chigh .ge.
     &                         1.d0-vof_tol) then
                              ok = 1
                           endif

                        else

c
c                          gx < 0: opposite phase ordering.
c
                           if (clow .ge.
     &                         1.d0-vof_tol .and.
     &                         chigh .le. vof_tol) then
                              ok = 1
                           endif

                        endif

                        if (ok .eq. 0) all_ok = 0

c
c                       Construct physical height by integrating
c                       over exactly seven cells in x.
c
                        h(m) = 0.d0

                        do kk = -3,3

                           h(m) = h(m)
     &                          + dx(0)*max(0.d0,
     &                            min(1.d0,
     &                            C(i0+kk,i1+m)))

                        enddo

                     enddo

                     if (all_ok .eq. 1) then

c
c                       h varies in the y direction.
c
                        hd = (h(1)-h(-1))
     &                       /(2.d0*dx(1))

                        hdd = (h(1)-2.d0*h(0)+h(-1))
     &                        /(dx(1)*dx(1))

                        denom = (1.d0+hd*hd)**1.5d0

                        K(i0,i1) = -hdd/denom
                        K_valid(i0,i1) = 1.d0

                     endif

                  endif

               endif

            endif

         enddo
      enddo

      return
      end

      subroutine sc_surface_tension_force_vof_2d(
     &     F0,F1,F_gcw,
     &     K,K_valid,K_gcw,
     &     N00,N11,N_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1)

      implicit none

      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER F_gcw,K_gcw,N_gcw

      REAL F0(SIDE2d0(ilower,iupper,F_gcw))
      REAL F1(SIDE2d1(ilower,iupper,F_gcw))

      REAL K(CELL2d(ilower,iupper,K_gcw))
      REAL K_valid(CELL2d(ilower,iupper,K_gcw))

      REAL N00(SIDE2d0(ilower,iupper,N_gcw))
      REAL N11(SIDE2d1(ilower,iupper,N_gcw))

      INTEGER i0,i1
      INTEGER vl,vr
      REAL kappa

c
c     X faces.
c
      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0+1

            vl = 0
            vr = 0

            if (K_valid(i0-1,i1) .gt. 0.5d0) vl = 1
            if (K_valid(i0,i1)   .gt. 0.5d0) vr = 1

            if (vl+vr .eq. 2) then

               kappa =
     &              0.5d0*(K(i0-1,i1)+K(i0,i1))

            else if (vl .eq. 1) then

               kappa = K(i0-1,i1)

            else if (vr .eq. 1) then

               kappa = K(i0,i1)

            else

               kappa = 0.d0
c we can add normal gradient calculation here 

            endif

            F0(i0,i1) = kappa*N00(i0,i1)

         enddo
      enddo

c
c     Y faces.
c
      do i1 = ilower1,iupper1+1
         do i0 = ilower0,iupper0

            vl = 0
            vr = 0

            if (K_valid(i0,i1-1) .gt. 0.5d0) vl = 1
            if (K_valid(i0,i1)   .gt. 0.5d0) vr = 1

            if (vl+vr .eq. 2) then

               kappa =
     &              0.5d0*(K(i0,i1-1)+K(i0,i1))

            else if (vl .eq. 1) then

               kappa = K(i0,i1-1)

            else if (vr .eq. 1) then

               kappa = K(i0,i1)

            else

               kappa = 0.d0 
c we can add normal gradient calculation here 

            endif

            F1(i0,i1) = kappa*N11(i0,i1)

         enddo
      enddo

      return
      end