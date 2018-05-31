c
c     Patch smoothers for use in multigrid preconditioners.
c
c     Created on 25 Aug 2010 by Boyce Griffith
c     Modified on 23 Aug 2017 by Amneet Bhalla
c
c     Copyright (c) 2002-2017, Boyce Griffith
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
define(NDIM,2)dnl
define(REAL,`double precision')dnl
define(INTEGER,`integer')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim2d.i)dnl
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Perform a single Gauss-Seidel sweep for F = alpha div grad U +
c     beta U.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine gssmooth2d(
     &     U,U_gcw,
     &     alpha,beta,
     &     F,F_gcw,
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
      INTEGER U_gcw,F_gcw

      REAL alpha,beta

      REAL F(ilower0-F_gcw:iupper0+F_gcw,
     &       ilower1-F_gcw:iupper1+F_gcw)

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL U(ilower0-U_gcw:iupper0+U_gcw,
     &       ilower1-U_gcw:iupper1+U_gcw)
c
c     Local variables.
c
      INTEGER i0,i1
      REAL    fac0,fac1,fac
c
c     Perform a single Gauss-Seidel sweep.
c
      fac0 = alpha/(dx(0)*dx(0))
      fac1 = alpha/(dx(1)*dx(1))
      fac = 0.5d0/(fac0+fac1-0.5d0*beta)

      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0
            U(i0,i1) = fac*(
     &           fac0*(U(i0-1,i1)+U(i0+1,i1)) +
     &           fac1*(U(i0,i1-1)+U(i0,i1+1)) -
     &           F(i0,i1))
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Perform a single "red" or "black" Gauss-Seidel sweep for F = alpha
c     div grad U + beta U.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine rbgssmooth2d(
     &     U,U_gcw,
     &     alpha,beta,
     &     F,F_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     dx,
     &     red_or_black)
c
      implicit none
c
c     Input.
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER U_gcw,F_gcw
      INTEGER red_or_black

      REAL alpha,beta

      REAL F(ilower0-F_gcw:iupper0+F_gcw,
     &       ilower1-F_gcw:iupper1+F_gcw)

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL U(ilower0-U_gcw:iupper0+U_gcw,
     &       ilower1-U_gcw:iupper1+U_gcw)
c
c     Local variables.
c
      INTEGER i0,i1
      REAL    fac0,fac1,fac
c
c     Perform a single "red" or "black" Gauss-Seidel sweep.
c
      red_or_black = mod(red_or_black,2) ! "red" = 0, "black" = 1

      fac0 = alpha/(dx(0)*dx(0))
      fac1 = alpha/(dx(1)*dx(1))
      fac = 0.5d0/(fac0+fac1-0.5d0*beta)

      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0
            if ( mod(i0+i1,2) .eq. red_or_black ) then
               U(i0,i1) = fac*(
     &              fac0*(U(i0-1,i1)+U(i0+1,i1)) +
     &              fac1*(U(i0,i1-1)+U(i0,i1+1)) -
     &              F(i0,i1))
            endif
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Perform a single Gauss-Seidel sweep for F = alpha div grad U +
c     beta U with masking of certain degrees of freedom.
c
c     NOTE: The solution U is unmodified at masked degrees of freedom.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine gssmoothmask2d(
     &     U,U_gcw,
     &     alpha,beta,
     &     F,F_gcw,
     &     mask,mask_gcw,
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
      INTEGER U_gcw,F_gcw,mask_gcw

      REAL alpha,beta

      REAL F(ilower0-F_gcw:iupper0+F_gcw,
     &       ilower1-F_gcw:iupper1+F_gcw)

      INTEGER mask(ilower0-mask_gcw:iupper0+mask_gcw,
     &             ilower1-mask_gcw:iupper1+mask_gcw)

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL U(ilower0-U_gcw:iupper0+U_gcw,
     &       ilower1-U_gcw:iupper1+U_gcw)
c
c     Local variables.
c
      INTEGER i0,i1
      REAL    fac0,fac1,fac
c
c     Perform a single Gauss-Seidel sweep.
c
      fac0 = alpha/(dx(0)*dx(0))
      fac1 = alpha/(dx(1)*dx(1))
      fac = 0.5d0/(fac0+fac1-0.5d0*beta)

      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0
            if (mask(i0,i1) .eq. 0) then
               U(i0,i1) = fac*(
     &              fac0*(U(i0-1,i1)+U(i0+1,i1)) +
     &              fac1*(U(i0,i1-1)+U(i0,i1+1)) -
     &              F(i0,i1))
            endif
         enddo
      enddo
c
      return
      end
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Perform a single "red" or "black" Gauss-Seidel sweep for F = alpha
c     div grad U + beta U with masking of certain degrees of freedom.
c
c     NOTE: The solution U is unmodified at masked degrees of freedom.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine rbgssmoothmask2d(
     &     U,U_gcw,
     &     alpha,beta,
     &     F,F_gcw,
     &     mask,mask_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     dx,
     &     red_or_black)
c
      implicit none
c
c     Input.
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER U_gcw,F_gcw,mask_gcw
      INTEGER red_or_black

      REAL alpha,beta

      REAL F(ilower0-F_gcw:iupper0+F_gcw,
     &       ilower1-F_gcw:iupper1+F_gcw)

      INTEGER mask(ilower0-mask_gcw:iupper0+mask_gcw,
     &             ilower1-mask_gcw:iupper1+mask_gcw)

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL U(ilower0-U_gcw:iupper0+U_gcw,
     &       ilower1-U_gcw:iupper1+U_gcw)
c
c     Local variables.
c
      INTEGER i0,i1
      REAL    fac0,fac1,fac
c
c     Perform a single "red" or "black" Gauss-Seidel sweep.
c
      red_or_black = mod(red_or_black,2) ! "red" = 0, "black" = 1

      fac0 = alpha/(dx(0)*dx(0))
      fac1 = alpha/(dx(1)*dx(1))
      fac = 0.5d0/(fac0+fac1-0.5d0*beta)

      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0
            if ( (mod(i0+i1,2) .eq. red_or_black) .and.
     &           (mask(i0,i1) .eq. 0) ) then
               U(i0,i1) = fac*(
     &              fac0*(U(i0-1,i1)+U(i0+1,i1)) +
     &              fac1*(U(i0,i1-1)+U(i0,i1+1)) -
     &              F(i0,i1))
            endif
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c              Variable coefficient patch smoothers

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Perform a single Gauss-Seidel sweep for F = div alpha grad U +
c     beta U.
c
c     The smoother is written for cell-centered U and side-centered
c     alpha = (alpha0,alpha1)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine vccellgssmooth2d(
     &     U,U_gcw,
     &     alpha0,alpha1,alpha_gcw,
     &     beta,
     &     F,F_gcw,
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
      INTEGER U_gcw,F_gcw,alpha_gcw

      REAL beta

      REAL F(ilower0-F_gcw:iupper0+F_gcw,
     &       ilower1-F_gcw:iupper1+F_gcw)

      REAL alpha0(SIDE2d0(ilower,iupper,alpha_gcw))
      REAL alpha1(SIDE2d1(ilower,iupper,alpha_gcw))

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL U(ilower0-U_gcw:iupper0+U_gcw,
     &       ilower1-U_gcw:iupper1+U_gcw)
c
c     Local variables.
c
      INTEGER i0,i1
      REAL    hx,hy
      REAL    facu0,facl0
      REAL    facu1,facl1
      REAL    fac
c
c     Perform a single Gauss-Seidel sweep.
c
      hx = dx(0)
      hy = dx(1)

      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0
            facu0 = alpha0(i0+1,i1)/(hx*hx)
            facl0 = alpha0(i0,i1)/(hx*hx)
            facu1 = alpha1(i0,i1+1)/(hy*hy)
            facl1 = alpha1(i0,i1)/(hy*hy)
            fac   = 1.d0/(facu0+facl0+facu1+facl1-beta)
            U(i0,i1) = fac*(
     &           facu0*U(i0+1,i1) +
     &           facl0*U(i0-1,i1) +
     &           facu1*U(i0,i1+1) +
     &           facl1*U(i0,i1-1) -
     &           F(i0,i1))
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Perform a single "red" or "black" Gauss-Seidel sweep for F =
c     div alpha grad U + beta U.
c
c     The smoother is written for cell-centered U and side-centered
c     alpha = (alpha0,alpha1)
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine vccellrbgssmooth2d(
     &     U,U_gcw,
     &     alpha0,alpha1,alpha_gcw,
     &     beta,
     &     F,F_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     dx,
     &     red_or_black)
c
      implicit none
c
c     Input.
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER U_gcw,F_gcw,alpha_gcw
      INTEGER red_or_black
      REAL beta

      REAL F(ilower0-F_gcw:iupper0+F_gcw,
     &       ilower1-F_gcw:iupper1+F_gcw)

      REAL alpha0(SIDE2d0(ilower,iupper,alpha_gcw))
      REAL alpha1(SIDE2d1(ilower,iupper,alpha_gcw))

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL U(ilower0-U_gcw:iupper0+U_gcw,
     &       ilower1-U_gcw:iupper1+U_gcw)
c
c     Local variables.
c
      INTEGER i0,i1
      REAL    hx,hy
      REAL    facu0,facl0
      REAL    facu1,facl1
      REAL    fac
c
c     Perform a single "red" or "black" Gauss-Seidel sweep.
c
      red_or_black = mod(red_or_black,2) ! "red" = 0, "black" = 1

      hx = dx(0)
      hy = dx(1)

      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0
            if ( mod(i0+i1,2) .eq. red_or_black ) then
               facu0 = alpha0(i0+1,i1)/(hx*hx)
               facl0 = alpha0(i0,i1)/(hx*hx)
               facu1 = alpha1(i0,i1+1)/(hy*hy)
               facl1 = alpha1(i0,i1)/(hy*hy)
               fac   = 1.d0/(facu0+facl0+facu1+facl1-beta)
               U(i0,i1) = fac*(
     &             facu0*U(i0+1,i1) +
     &             facl0*U(i0-1,i1) +
     &             facu1*U(i0,i1+1) +
     &             facl1*U(i0,i1-1) -
     &             F(i0,i1))
            endif
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  Perform a single Gauss-Seidel sweep for 
c     (f0,f1) = alpha div mu (grad (u0,u1) + grad (u0, u1)^T) + beta c (u0,u1).
c
c  The smoother is written for side-centered vector fields (u0, u1) and (f0, f1)
c  with node-centered coefficient mu and side-centered coefficient (c0,c1)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine vcgssmooth2d(
     &     u0,u1,u_gcw,
     &     f0,f1,f_gcw,
     &     c0,c1,c_gcw,
     &     mu,mu_gcw,
     &     alpha, beta,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     dx,
     &     var_c,
     &     use_harmonic_interp)
c
      implicit none
c
c     Functions.
c
      REAL a_avg, h_avg

c
c     Input.
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER u_gcw,f_gcw,c_gcw,mu_gcw
      INTEGER var_c,use_harmonic_interp

      REAL alpha,beta

      REAL mu(NODE2d(ilower,iupper,mu_gcw))
      
      REAL f0(SIDE2d0(ilower,iupper,f_gcw))
      REAL f1(SIDE2d1(ilower,iupper,f_gcw))

      REAL c0(SIDE2d0(ilower,iupper,c_gcw))
      REAL c1(SIDE2d1(ilower,iupper,c_gcw))

      REAL dx(0:NDIM-1)

c
c     Input/Output.
c

      REAL u0(SIDE2d0(ilower,iupper,u_gcw))
      REAL u1(SIDE2d1(ilower,iupper,u_gcw))
c
c     Local variables.
c
      INTEGER i0,i1
      REAL fac0,fac1,fac,nmr,dnr,mu_lower,mu_upper,c
c
c     Perform a single Gauss-Seidel sweep.
c
      fac0 = 1.d0/(dx(0))
      fac1 = 1.d0/(dx(1))

      fac = 2.d0*fac0**2.d0
      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0+1

         c = beta
         if (var_c .eq. 1) then
            c = c0(i0,i1)*beta
         endif       
         
         if (use_harmonic_interp .eq. 1) then
            mu_upper = h_avg(mu(i0,i1),mu(i0+1,i1),
     &                       mu(i0,i1+1),mu(i0+1,i1+1))
            mu_lower = h_avg(mu(i0,i1),mu(i0-1,i1),
     &                       mu(i0,i1+1),mu(i0-1,i1+1))
         else
            mu_upper = a_avg(mu(i0,i1),mu(i0+1,i1),
     &                       mu(i0,i1+1),mu(i0+1,i1+1))
            mu_lower = a_avg(mu(i0,i1),mu(i0-1,i1),
     &                       mu(i0,i1+1),mu(i0-1,i1+1))
         endif
         
         dnr =  alpha*(fac*(mu_upper + mu_lower) + 
     &       fac1**2.d0*(mu(i0,i1+1) + mu(i0,i1))) - c

         nmr = -f0(i0,i1) + alpha*(fac*(
     &     mu_upper*u0(i0+1,i1) + mu_lower*u0(i0-1,i1))+  
     &     fac1**2.d0*(mu(i0,i1+1)*u0(i0,i1+1) + mu(i0,i1)*u0(i0,i1-1))+   
     &     fac0*fac1*(mu(i0,i1+1)*(u1(i0,i1+1)-u1(i0-1,i1+1))-
     &        mu(i0,i1)*(u1(i0,i1)-u1(i0-1,i1))))               

         u0(i0,i1) = nmr/dnr
         enddo
      enddo

      fac = 2.d0*fac1**2.d0
      do i1 = ilower1,iupper1+1
         do i0 = ilower0,iupper0
  
         c = beta
         if (var_c .eq. 1) then
            c = c1(i0,i1)*beta
         endif

         if (use_harmonic_interp .eq. 1) then
            mu_upper = h_avg(mu(i0,i1),mu(i0+1,i1),
     &                       mu(i0,i1+1),mu(i0+1,i1+1))
            mu_lower = h_avg(mu(i0,i1),mu(i0+1,i1),
     &                       mu(i0,i1-1),mu(i0+1,i1-1))
         else
            mu_upper = a_avg(mu(i0,i1),mu(i0+1,i1),
     &                       mu(i0,i1+1),mu(i0+1,i1+1))
            mu_lower = a_avg(mu(i0,i1),mu(i0+1,i1),
     &                       mu(i0,i1-1),mu(i0+1,i1-1))
         endif
        
         dnr = alpha*(fac*(mu_upper+ mu_lower)+
     &      fac0**2.d0*(mu(i0+1,i1) + mu(i0,i1))) - c

         nmr = -f1(i0,i1) + alpha*(fac*(
     &    mu_upper*u1(i0,i1+1) + mu_lower*u1(i0,i1-1))+
     &    fac0**2.d0*(mu(i0+1,i1)*u1(i0+1,i1) + mu(i0,i1)*u1(i0-1,i1))+
     &    fac0*fac1*(mu(i0+1,i1)*(u0(i0+1,i1)-u0(i0+1,i1-1))-
     &       mu(i0,i1)*(u0(i0,i1)-u0(i0,i1-1))))  

         u1(i0,i1) = nmr/dnr

         enddo
      enddo

c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  Perform a single "red" or "black" Gauss-Seidel sweep for 
c     (f0,f1) = alpha div mu (grad (u0,u1) + grad (u0, u1)^T) + beta c (u0,u1).
c
c  The smoother is written for side-centered vector fields (u0, u1) and (f0, f1)
c  with node-centered coefficient mu and side-centered coefficient (c0,c1)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine vcrbgssmooth2d(
     &     u0,u1,u_gcw,
     &     f0,f1,f_gcw,
     &     c0,c1,c_gcw,
     &     mu,mu_gcw,
     &     alpha,beta,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     dx,
     &     var_c,
     &     use_harmonic_interp,
     &     red_or_black)
c
      implicit none
c
c     Functions.
c
      REAL a_avg, h_avg
c
c     Input.
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER u_gcw,f_gcw,c_gcw,mu_gcw
      INTEGER var_c,use_harmonic_interp
      INTEGER red_or_black

      REAL alpha,beta

      REAL mu(NODE2d(ilower,iupper,mu_gcw))
      
      REAL f0(SIDE2d0(ilower,iupper,f_gcw))
      REAL f1(SIDE2d1(ilower,iupper,f_gcw))

      REAL c0(SIDE2d0(ilower,iupper,c_gcw))
      REAL c1(SIDE2d1(ilower,iupper,c_gcw))

      REAL dx(0:NDIM-1)

c
c     Input/Output.
c

      REAL u0(SIDE2d0(ilower,iupper,u_gcw))
      REAL u1(SIDE2d1(ilower,iupper,u_gcw))
c
c     Local variables.
c
      INTEGER i0,i1
      REAL fac0,fac1,fac,nmr,dnr,mu_lower,mu_upper,c
c
c     Perform a single "red" or "black"  Gauss-Seidel sweep.
c
      red_or_black = mod(red_or_black,2) ! "red" = 0, "black" = 1
  
      fac0 = 1.d0/(dx(0))
      fac1 = 1.d0/(dx(1))

      fac = 2.d0*fac0**2.d0
      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0+1
            if (mod(i0+i1,2) .eq. red_or_black) then
         
            c = beta
            if (var_c .eq. 1) then
               c = c0(i0,i1)*beta
            endif       
         
            if (use_harmonic_interp .eq. 1) then
                mu_upper = h_avg(mu(i0,i1),mu(i0+1,i1),
     &                           mu(i0,i1+1),mu(i0+1,i1+1))
                mu_lower = h_avg(mu(i0,i1),mu(i0-1,i1),
     &                       mu(i0,i1+1),mu(i0-1,i1+1))
            else
                mu_upper = a_avg(mu(i0,i1),mu(i0+1,i1),
     &                           mu(i0,i1+1),mu(i0+1,i1+1))
                mu_lower = a_avg(mu(i0,i1),mu(i0-1,i1),
     &                       mu(i0,i1+1),mu(i0-1,i1+1))
            endif
            
            dnr =  alpha*(fac*(mu_upper + mu_lower) + 
     &         fac1**2.d0*(mu(i0,i1+1) + mu(i0,i1))) - c

            nmr = -f0(i0,i1) + alpha*(fac*(
     &         mu_upper*u0(i0+1,i1) + mu_lower*u0(i0-1,i1))+  
     &         fac1**2.d0*(mu(i0,i1+1)*u0(i0,i1+1)+ 
     &            mu(i0,i1)*u0(i0,i1-1))+   
     &         fac0*fac1*(mu(i0,i1+1)*(u1(i0,i1+1)-
     &            u1(i0-1,i1+1))-
     &            mu(i0,i1)*(u1(i0,i1)-u1(i0-1,i1))))               

            u0(i0,i1) = nmr/dnr
            
          endif
         enddo
      enddo

      fac = 2.d0*fac1**2.d0
      do i1 = ilower1,iupper1+1
         do i0 = ilower0,iupper0
            if (mod(i0+i1,2) .eq. red_or_black) then  
  
            c = beta
            if (var_c .eq. 1) then
               c = c1(i0,i1)*beta
            endif

            if (use_harmonic_interp .eq. 1) then
                mu_upper = h_avg(mu(i0,i1),mu(i0+1,i1),
     &                           mu(i0,i1+1),mu(i0+1,i1+1))
                mu_lower = h_avg(mu(i0,i1),mu(i0+1,i1),
     &                           mu(i0,i1-1),mu(i0+1,i1-1))
            else
                mu_upper = a_avg(mu(i0,i1),mu(i0+1,i1),
     &                           mu(i0,i1+1),mu(i0+1,i1+1))
                mu_lower = a_avg(mu(i0,i1),mu(i0+1,i1),
     &                           mu(i0,i1-1),mu(i0+1,i1-1))
            endif
        
            dnr = alpha*(fac*(mu_upper+ mu_lower)+
     &         fac0**2.d0*(mu(i0+1,i1) + mu(i0,i1))) - c

            nmr = -f1(i0,i1) + alpha*(fac*(
     &         mu_upper*u1(i0,i1+1) + mu_lower*u1(i0,i1-1))+
     &         fac0**2.d0*(mu(i0+1,i1)*u1(i0+1,i1)+ 
     &            mu(i0,i1)*u1(i0-1,i1))+
     &         fac0*fac1*(mu(i0+1,i1)*(u0(i0+1,i1)
     &            -u0(i0+1,i1-1))-
     &         mu(i0,i1)*(u0(i0,i1)-u0(i0,i1-1))))  

            u1(i0,i1) = nmr/dnr
       
           endif
         enddo
      enddo

c
      return
      end
c

c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  Perform a single Gauss-Seidel sweep for 
c     (f0,f1) = alpha div mu (grad (u0,u1) + grad (u0, u1)^T) + beta c (u0,u1),
c  with masking of certain degrees of freedom.
c
c     NOTE: The solution (u0,u1) is unmodified at masked degrees of freedom.
c
c  The smoother is written for side-centered vector fields (u0, u1) and (f0, f1)
c  with node-centered coefficient mu and side-centered coefficient (c0,c1)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine vcgssmoothmask2d(
     &     u0,u1,u_gcw,
     &     f0,f1,f_gcw,
     &     mask0,mask1,mask_gcw,
     &     c0,c1,c_gcw,
     &     mu,mu_gcw,
     &     alpha, beta,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     dx,
     &     var_c,
     &     use_harmonic_interp)
c
      implicit none
c
c     Functions.
c
      REAL a_avg, h_avg
c
c     Input.
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER u_gcw,f_gcw,c_gcw,mu_gcw,mask_gcw
      INTEGER var_c,use_harmonic_interp

      REAL alpha,beta

      REAL mu(NODE2d(ilower,iupper,mu_gcw))
      
      REAL f0(SIDE2d0(ilower,iupper,f_gcw))
      REAL f1(SIDE2d1(ilower,iupper,f_gcw))

      INTEGER mask0(SIDE2d0(ilower,iupper,mask_gcw))
      INTEGER mask1(SIDE2d1(ilower,iupper,mask_gcw))
      
      REAL c0(SIDE2d0(ilower,iupper,c_gcw))
      REAL c1(SIDE2d1(ilower,iupper,c_gcw))

      REAL dx(0:NDIM-1)

c
c     Input/Output.
c

      REAL u0(SIDE2d0(ilower,iupper,u_gcw))
      REAL u1(SIDE2d1(ilower,iupper,u_gcw))
c
c     Local variables.
c
      INTEGER i0,i1
      REAL fac0,fac1,fac,nmr,dnr,mu_lower,mu_upper,c
c
c     Perform a single Gauss-Seidel sweep.
c
      fac0 = 1.d0/(dx(0))
      fac1 = 1.d0/(dx(1))

      fac = 2.d0*fac0**2.d0
      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0+1
            if (mask0(i0,i1) .eq. 0) then
         
            c = beta
            if (var_c .eq. 1) then
               c = c0(i0,i1)*beta
            endif       
         
            if (use_harmonic_interp .eq. 1) then
                mu_upper = h_avg(mu(i0,i1),mu(i0+1,i1),
     &                           mu(i0,i1+1),mu(i0+1,i1+1))
                mu_lower = h_avg(mu(i0,i1),mu(i0-1,i1),
     &                           mu(i0,i1+1),mu(i0-1,i1+1))
            else
                mu_upper = a_avg(mu(i0,i1),mu(i0+1,i1),
     &                           mu(i0,i1+1),mu(i0+1,i1+1))
                mu_lower = a_avg(mu(i0,i1),mu(i0-1,i1),
     &                           mu(i0,i1+1),mu(i0-1,i1+1))
            endif
            
            dnr =  alpha*(fac*(mu_upper + mu_lower) + 
     &         fac1**2.d0*(mu(i0,i1+1) + mu(i0,i1))) - c

            nmr = -f0(i0,i1) + alpha*(fac*(
     &         mu_upper*u0(i0+1,i1) + mu_lower*u0(i0-1,i1))+  
     &         fac1**2.d0*(mu(i0,i1+1)*u0(i0,i1+1)+ 
     &            mu(i0,i1)*u0(i0,i1-1))+   
     &         fac0*fac1*(mu(i0,i1+1)*(u1(i0,i1+1)-
     &            u1(i0-1,i1+1))-
     &            mu(i0,i1)*(u1(i0,i1)-u1(i0-1,i1))))               

            u0(i0,i1) = nmr/dnr

            endif
         enddo
      enddo

      fac = 2.d0*fac1**2.d0
      do i1 = ilower1,iupper1+1
         do i0 = ilower0,iupper0
            if (mask1(i0,i1) .eq. 0) then
         
            c = beta
            if (var_c .eq. 1) then
               c = c1(i0,i1)*beta
            endif

            if (use_harmonic_interp .eq. 1) then
                mu_upper = h_avg(mu(i0,i1),mu(i0+1,i1),
     &                           mu(i0,i1+1),mu(i0+1,i1+1))
                mu_lower = h_avg(mu(i0,i1),mu(i0+1,i1),
     &                           mu(i0,i1-1),mu(i0+1,i1-1))
            else
                mu_upper = a_avg(mu(i0,i1),mu(i0+1,i1),
     &                           mu(i0,i1+1),mu(i0+1,i1+1))
                mu_lower = a_avg(mu(i0,i1),mu(i0+1,i1),
     &                           mu(i0,i1-1),mu(i0+1,i1-1))
            endif
        
            dnr = alpha*(fac*(mu_upper+ mu_lower)+
     &         fac0**2.d0*(mu(i0+1,i1) + mu(i0,i1))) - c

            nmr = -f1(i0,i1) + alpha*(fac*(
     &         mu_upper*u1(i0,i1+1) + mu_lower*u1(i0,i1-1))+
     &         fac0**2.d0*(mu(i0+1,i1)*u1(i0+1,i1)+ 
     &            mu(i0,i1)*u1(i0-1,i1))+
     &         fac0*fac1*(mu(i0+1,i1)*(u0(i0+1,i1)-
     &            u0(i0+1,i1-1))-
     &         mu(i0,i1)*(u0(i0,i1)-u0(i0,i1-1))))  

            u1(i0,i1) = nmr/dnr
 
            endif
         enddo
      enddo

c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  Perform a single "red" or "black" Gauss-Seidel sweep for 
c     (f0,f1) = alpha div mu (grad (u0,u1) + grad (u0, u1)^T) + beta c (u0,u1),
c  with masking of certain degrees of freedom.
c
c     NOTE: The solution (u0,u1) is unmodified at masked degrees of freedom.
c
c  The smoother is written for side-centered vector fields (u0, u1) and (f0, f1)
c  with node-centered coefficient mu and side-centered coefficient (c0,c1)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine vcrbgssmoothmask2d(
     &     u0,u1,u_gcw,
     &     f0,f1,f_gcw,
     &     mask0,mask1,mask_gcw,
     &     c0,c1,c_gcw,
     &     mu,mu_gcw,
     &     alpha,beta,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     dx,
     &     var_c,
     &     use_harmonic_interp,
     &     red_or_black)
c
      implicit none
c
c     Functions.
c
      REAL a_avg, h_avg
c
c     Input.
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER u_gcw,f_gcw,c_gcw,mu_gcw,mask_gcw
      INTEGER var_c,use_harmonic_interp
      INTEGER red_or_black

      REAL alpha,beta

      REAL mu(NODE2d(ilower,iupper,mu_gcw))
      
      REAL f0(SIDE2d0(ilower,iupper,f_gcw))
      REAL f1(SIDE2d1(ilower,iupper,f_gcw))

      INTEGER mask0(SIDE2d0(ilower,iupper,mask_gcw))
      INTEGER mask1(SIDE2d1(ilower,iupper,mask_gcw))
      
      REAL c0(SIDE2d0(ilower,iupper,c_gcw))
      REAL c1(SIDE2d1(ilower,iupper,c_gcw))

      REAL dx(0:NDIM-1)

c
c     Input/Output.
c

      REAL u0(SIDE2d0(ilower,iupper,u_gcw))
      REAL u1(SIDE2d1(ilower,iupper,u_gcw))
c
c     Local variables.
c
      INTEGER i0,i1
      REAL fac0,fac1,fac,nmr,dnr,mu_lower,mu_upper,c
c
c     Perform a single "red" or "black"  Gauss-Seidel sweep.
c
      red_or_black = mod(red_or_black,2) ! "red" = 0, "black" = 1
  
      fac0 = 1.d0/(dx(0))
      fac1 = 1.d0/(dx(1))

      fac = 2.d0*fac0**2.d0
      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0+1
            if ( (mod(i0+i1,2) .eq. red_or_black) .and. 
     &           (mask0(i0,i1) .eq. 0) ) then
         
            c = beta      
            if (var_c .eq. 1) then
               c = c0(i0,i1)*beta
            endif       
         
            if (use_harmonic_interp .eq. 1) then
                mu_upper = h_avg(mu(i0,i1),mu(i0+1,i1),
     &                           mu(i0,i1+1),mu(i0+1,i1+1))
                mu_lower = h_avg(mu(i0,i1),mu(i0-1,i1),
     &                           mu(i0,i1+1),mu(i0-1,i1+1))
            else
                mu_upper = a_avg(mu(i0,i1),mu(i0+1,i1),
     &                           mu(i0,i1+1),mu(i0+1,i1+1))
                mu_lower = a_avg(mu(i0,i1),mu(i0-1,i1),
     &                       mu(i0,i1+1),mu(i0-1,i1+1))
            endif
            
            dnr =  alpha*(fac*(mu_upper + mu_lower) + 
     &         fac1**2.d0*(mu(i0,i1+1) + mu(i0,i1))) - c

            nmr = -f0(i0,i1) + alpha*(fac*(
     &         mu_upper*u0(i0+1,i1) + mu_lower*u0(i0-1,i1))+  
     &         fac1**2.d0*(mu(i0,i1+1)*u0(i0,i1+1)+ 
     &            mu(i0,i1)*u0(i0,i1-1))+   
     &         fac0*fac1*(mu(i0,i1+1)*(u1(i0,i1+1)-
     &            u1(i0-1,i1+1))-
     &            mu(i0,i1)*(u1(i0,i1)-u1(i0-1,i1))))               

            u0(i0,i1) = nmr/dnr
            
          endif
         enddo
      enddo

      fac = 2.d0*fac1**2.d0
      do i1 = ilower1,iupper1+1
         do i0 = ilower0,iupper0
            if ( (mod(i0+i1,2) .eq. red_or_black) .and. 
     &           (mask1(i0,i1) .eq. 0) ) then  
  
            c = beta
            if (var_c .eq. 1) then
               c = c1(i0,i1)*beta
            endif

            if (use_harmonic_interp .eq. 1) then
                mu_upper = h_avg(mu(i0,i1),mu(i0+1,i1),
     &                           mu(i0,i1+1),mu(i0+1,i1+1))
                mu_lower = h_avg(mu(i0,i1),mu(i0+1,i1),
     &                           mu(i0,i1-1),mu(i0+1,i1-1))
            else
                mu_upper = a_avg(mu(i0,i1),mu(i0+1,i1),
     &                           mu(i0,i1+1),mu(i0+1,i1+1))
                mu_lower = a_avg(mu(i0,i1),mu(i0+1,i1),
     &                           mu(i0,i1-1),mu(i0+1,i1-1))
            endif
        
            dnr = alpha*(fac*(mu_upper+ mu_lower)+
     &         fac0**2.d0*(mu(i0+1,i1) + mu(i0,i1))) - c

            nmr = -f1(i0,i1) + alpha*(fac*(
     &         mu_upper*u1(i0,i1+1) + mu_lower*u1(i0,i1-1))+
     &         fac0**2.d0*(mu(i0+1,i1)*u1(i0+1,i1)+ 
     &            mu(i0,i1)*u1(i0-1,i1))+
     &         fac0*fac1*(mu(i0+1,i1)*(u0(i0+1,i1)
     &            -u0(i0+1,i1-1))-
     &         mu(i0,i1)*(u0(i0,i1)-u0(i0,i1-1))))  

            u1(i0,i1) = nmr/dnr
       
           endif
         enddo
      enddo

c
      return
      end
c

c
