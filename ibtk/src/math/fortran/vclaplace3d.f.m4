c
c     Routines to compute quantities related to variable coefficient
c     generalized Laplace operators.
c
c     Created on 14 Sep 2017 by Nishant Nangia
c
c     Copyright (c) 2002-2014, Boyce Griffith
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
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the arithmetic average of twelve inputs
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      REAL function a_avg(
     &     a0,a1,a2,a3,
     &     a4,a5,a6,a7,
     &     a8,a9,a10,a11)
      implicit none
      REAL a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11
      a_avg = (1.d0/12.d0)*(a0+a1+a2+a3+a4+a5+a6+a7+a8+a9+a10+a11)
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the harmonic average of twelve inputs
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      REAL function h_avg(
     &     a0,a1,a2,a3,
     &     a4,a5,a6,a7,
     &     a8,a9,a10,a11)
      implicit none
      REAL a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11
      REAL dmr
      dmr = 1.d0/a0+1.d0/a1+1.d0/a2+1.d0/a3+1.d0/a4+1.d0/a5+
     &      1.d0/a6+1.d0/a7+1.d0/a8+1.d0/a9+1.d0/a10+1.d0/a11
      h_avg = 12.d0/dmr
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes (f0,f1,f2) = alpha div mu (grad (u0,u1,u1) + grad(u0,u1,u2)^T) 
c     + beta (u0,u1,f2) + gamma (v0,v1,v2)
c
c     Computes the side-centered variable coefficient generalized
c     Laplacian, with edge-centered coefficient mu and side-centered
c     vector fields (u0,u1,u2) and (v0,v1,v2).
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine stosvclaplace3d(
     &     f0,f1,f2,f_gcw,
     &     alpha,beta,
     &     mu0,mu1,mu2,mu_gcw,
     &     rho0,rho1,rho2,rho_gcw,
     &     u0,u1,u2,u_gcw,
     &     gamma,
     &     v0,v1,v2,v_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
     &     dx,
     &     var_rho,
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
      INTEGER ilower2,iupper2
      INTEGER f_gcw,mu_gcw,rho_gcw,u_gcw,v_gcw,var_rho
      INTEGER use_harmonic_interp

      REAL alpha,beta,gamma

      REAL mu0(EDGE3d0(ilower,iupper,mu_gcw))
      REAL mu1(EDGE3d1(ilower,iupper,mu_gcw))
      REAL mu2(EDGE3d2(ilower,iupper,mu_gcw))

      REAL rho0(SIDE3d0(ilower,iupper,rho_gcw))
      REAL rho1(SIDE3d1(ilower,iupper,rho_gcw))
      REAL rho2(SIDE3d2(ilower,iupper,rho_gcw))

      REAL u0(SIDE3d0(ilower,iupper,u_gcw))
      REAL u1(SIDE3d1(ilower,iupper,u_gcw))
      REAL u2(SIDE3d2(ilower,iupper,u_gcw))

      REAL v0(SIDE3d0(ilower,iupper,v_gcw))
      REAL v1(SIDE3d1(ilower,iupper,v_gcw))
      REAL v2(SIDE3d2(ilower,iupper,v_gcw))
      
      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL f0(SIDE3d0(ilower,iupper,f_gcw))
      REAL f1(SIDE3d1(ilower,iupper,f_gcw))
      REAL f2(SIDE3d2(ilower,iupper,f_gcw))
c
c     Local variables.
c
      INTEGER i0,i1,i2
      REAL fac0,fac1,fac2,mu_lower,mu_upper,rho
c
c     Compute the discrete divergence of mu (grad (u0,u1,u2) + grad (u0,u1,u2)^T).
c
      fac0 = 1.d0/(dx(0))
      fac1 = 1.d0/(dx(1))
      fac2 = 1.d0/(dx(2))

      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0+1

               rho = beta
               if (var_rho .eq. 1) then
                  rho = rho0(i0,i1,i2)*beta
               endif

               if (use_harmonic_interp .eq. 1) then
                   mu_upper = h_avg(mu0(i0,i1,i2),mu0(i0,i1+1,i2),
     &                              mu0(i0,i1,i2+1),mu0(i0,i1+1,i2+1),
     &                              mu1(i0,i1,i2),mu1(i0+1,i1,i2),
     &                              mu1(i0,i1,i2+1),mu1(i0+1,i1,i2+1),
     &                              mu2(i0,i1,i2),mu2(i0+1,i1,i2),
     &                              mu2(i0,i1+1,i2),mu2(i0+1,i1+1,i2))

                   mu_lower = h_avg(mu0(i0-1,i1,i2),mu0(i0-1,i1+1,i2),
     &                            mu0(i0-1,i1,i2+1),mu0(i0-1,i1+1,i2+1),
     &                            mu1(i0-1,i1,i2),mu1(i0,i1,i2),
     &                            mu1(i0-1,i1,i2+1),mu1(i0,i1,i2+1),
     &                            mu2(i0-1,i1,i2),mu2(i0,i1,i2),
     &                            mu2(i0-1,i1+1,i2),mu2(i0,i1+1,i2))
               else
                    mu_upper = a_avg(mu0(i0,i1,i2),mu0(i0,i1+1,i2),
     &                              mu0(i0,i1,i2+1),mu0(i0,i1+1,i2+1),
     &                              mu1(i0,i1,i2),mu1(i0+1,i1,i2),
     &                              mu1(i0,i1,i2+1),mu1(i0+1,i1,i2+1),
     &                              mu2(i0,i1,i2),mu2(i0+1,i1,i2),
     &                              mu2(i0,i1+1,i2),mu2(i0+1,i1+1,i2))

                    mu_lower = a_avg(mu0(i0-1,i1,i2),mu0(i0-1,i1+1,i2),
     &                            mu0(i0-1,i1,i2+1),mu0(i0-1,i1+1,i2+1),
     &                            mu1(i0-1,i1,i2),mu1(i0,i1,i2),
     &                            mu1(i0-1,i1,i2+1),mu1(i0,i1,i2+1),
     &                            mu2(i0-1,i1,i2),mu2(i0,i1,i2),
     &                            mu2(i0-1,i1+1,i2),mu2(i0,i1+1,i2))
               endif

               f0(i0,i1,i2) = alpha*(
     &              2.d0*fac0**2.d0*(
     &              mu_upper*(u0(i0+1,i1,i2)-u0(i0,i1,i2))-
     &              mu_lower*(u0(i0,i1,i2)-u0(i0-1,i1,i2)))+
     &              fac1**2.d0*(
     &              mu2(i0,i1+1,i2)*(u0(i0,i1+1,i2)-u0(i0,i1,i2))-
     &              mu2(i0,i1,i2)*(u0(i0,i1,i2)-u0(i0,i1-1,i2)))+
     &              fac0*fac1*(
     &              mu2(i0,i1+1,i2)*(u1(i0,i1+1,i2)-u1(i0-1,i1+1,i2))-
     &              mu2(i0,i1,i2)*(u1(i0,i1,i2)-u1(i0-1,i1,i2)))+
     &              fac2**2.d0*(
     &              mu1(i0,i1,i2+1)*(u0(i0,i1,i2+1)-u0(i0,i1,i2))-
     &              mu1(i0,i1,i2)*(u0(i0,i1,i2)-u0(i0,i1,i2-1)))+
     &              fac0*fac2*(
     &              mu1(i0,i1,i2+1)*(u2(i0,i1,i2+1)-u2(i0-1,i1,i2+1))-
     &              mu1(i0,i1,i2)*(u2(i0,i1,i2)-u2(i0-1,i1,i2))))+
     &              rho*u0(i0,i1,i2) + gamma*v0(i0,i1,i2)
            enddo
         enddo
      enddo

      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1+1
            do i0 = ilower0,iupper0

               rho = beta
               if (var_rho .eq. 1) then
                  rho = rho1(i0,i1,i2)*beta
               endif

               if (use_harmonic_interp .eq. 1) then
                  mu_upper = h_avg(mu0(i0,i1,i2),mu0(i0,i1+1,i2),
     &                             mu0(i0,i1,i2+1),mu0(i0,i1+1,i2+1),
     &                             mu1(i0,i1,i2),mu1(i0+1,i1,i2),
     &                             mu1(i0,i1,i2+1),mu1(i0+1,i1,i2+1),
     &                             mu2(i0,i1,i2),mu2(i0+1,i1,i2),
     &                             mu2(i0,i1+1,i2),mu2(i0+1,i1+1,i2))
                  mu_lower = h_avg(mu0(i0,i1-1,i2),mu0(i0,i1,i2),
     &                             mu0(i0,i1-1,i2+1),mu0(i0,i1,i2+1),
     &                             mu1(i0,i1-1,i2),mu1(i0+1,i1-1,i2),
     &                            mu1(i0,i1-1,i2+1),mu1(i0+1,i1-1,i2+1),
     &                            mu2(i0,i1-1,i2),mu2(i0+1,i1-1,i2),
     &                            mu2(i0,i1,i2),mu2(i0+1,i1,i2))
               else
                  mu_upper = a_avg(mu0(i0,i1,i2),mu0(i0,i1+1,i2),
     &                             mu0(i0,i1,i2+1),mu0(i0,i1+1,i2+1),
     &                             mu1(i0,i1,i2),mu1(i0+1,i1,i2),
     &                             mu1(i0,i1,i2+1),mu1(i0+1,i1,i2+1),
     &                             mu2(i0,i1,i2),mu2(i0+1,i1,i2),
     &                             mu2(i0,i1+1,i2),mu2(i0+1,i1+1,i2))
                  mu_lower = a_avg(mu0(i0,i1-1,i2),mu0(i0,i1,i2),
     &                             mu0(i0,i1-1,i2+1),mu0(i0,i1,i2+1),
     &                             mu1(i0,i1-1,i2),mu1(i0+1,i1-1,i2),
     &                            mu1(i0,i1-1,i2+1),mu1(i0+1,i1-1,i2+1),
     &                            mu2(i0,i1-1,i2),mu2(i0+1,i1-1,i2),
     &                            mu2(i0,i1,i2),mu2(i0+1,i1,i2))
                endif       

               f1(i0,i1,i2) = alpha*(
     &              2.d0*fac1**2.d0*(
     &              mu_upper*(u1(i0,i1+1,i2)-u1(i0,i1,i2))-
     &              mu_lower*(u1(i0,i1,i2)-u1(i0,i1-1,i2)))+
     &              fac0**2.d0*(
     &              mu2(i0+1,i1,i2)*(u1(i0+1,i1,i2)-u1(i0,i1,i2))-
     &              mu2(i0,i1,i2)*(u1(i0,i1,i2)-u1(i0-1,i1,i2)))+
     &              fac0*fac1*(
     &              mu2(i0+1,i1,i2)*(u0(i0+1,i1,i2)-u0(i0+1,i1-1,i2))-
     &              mu2(i0,i1,i2)*(u0(i0,i1,i2)-u0(i0,i1-1,i2)))+
     &              fac2**2.d0*(
     &              mu0(i0,i1,i2+1)*(u1(i0,i1,i2+1)-u1(i0,i1,i2))-
     &              mu0(i0,i1,i2)*(u1(i0,i1,i2)-u1(i0,i1,i2-1)))+
     &              fac1*fac2*(
     &              mu0(i0,i1,i2+1)*(u2(i0,i1,i2+1)-u2(i0,i1-1,i2+1))-
     &              mu0(i0,i1,i2)*(u2(i0,i1,i2)-u2(i0,i1-1,i2))))+
     &              rho*u1(i0,i1,i2) + gamma*v1(i0,i1,i2)
            enddo
         enddo
      enddo

      do i2 = ilower2,iupper2+1
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0

               rho = beta
               if (var_rho .eq. 1) then
                  rho = rho2(i0,i1,i2)*beta
               endif

               if (use_harmonic_interp .eq. 1) then
                  mu_upper = h_avg(mu0(i0,i1,i2),mu0(i0,i1+1,i2),
     &                             mu0(i0,i1,i2+1),mu0(i0,i1+1,i2+1),
     &                             mu1(i0,i1,i2),mu1(i0+1,i1,i2),
     &                             mu1(i0,i1,i2+1),mu1(i0+1,i1,i2+1),
     &                             mu2(i0,i1,i2),mu2(i0+1,i1,i2),
     &                             mu2(i0,i1+1,i2),mu2(i0+1,i1+1,i2))

                  mu_lower = h_avg(mu0(i0,i1,i2-1),mu0(i0,i1+1,i2-1),
     &                             mu0(i0,i1,i2),mu0(i0,i1+1,i2),
     &                             mu1(i0,i1,i2-1),mu1(i0+1,i1,i2-1),
     &                             mu1(i0,i1,i2),mu1(i0+1,i1,i2),
     &                             mu2(i0,i1,i2-1),mu2(i0+1,i1,i2-1),
     &                            mu2(i0,i1+1,i2-1),mu2(i0+1,i1+1,i2-1))
               else
                  mu_upper = a_avg(mu0(i0,i1,i2),mu0(i0,i1+1,i2),
     &                             mu0(i0,i1,i2+1),mu0(i0,i1+1,i2+1),
     &                             mu1(i0,i1,i2),mu1(i0+1,i1,i2),
     &                             mu1(i0,i1,i2+1),mu1(i0+1,i1,i2+1),
     &                             mu2(i0,i1,i2),mu2(i0+1,i1,i2),
     &                             mu2(i0,i1+1,i2),mu2(i0+1,i1+1,i2))

                  mu_lower = a_avg(mu0(i0,i1,i2-1),mu0(i0,i1+1,i2-1),
     &                             mu0(i0,i1,i2),mu0(i0,i1+1,i2),
     &                             mu1(i0,i1,i2-1),mu1(i0+1,i1,i2-1),
     &                             mu1(i0,i1,i2),mu1(i0+1,i1,i2),
     &                             mu2(i0,i1,i2-1),mu2(i0+1,i1,i2-1),
     &                            mu2(i0,i1+1,i2-1),mu2(i0+1,i1+1,i2-1))
               endif

               f2(i0,i1,i2) = alpha*(
     &              2.d0*fac2**2.d0*(
     &              mu_upper*(u2(i0,i1,i2+1)-u2(i0,i1,i2))-
     &              mu_lower*(u2(i0,i1,i2)-u2(i0,i1,i2-1)))+
     &              fac1**2.d0*(
     &              mu0(i0,i1+1,i2)*(u2(i0,i1+1,i2)-u2(i0,i1,i2))-
     &              mu0(i0,i1,i2)*(u2(i0,i1,i2)-u2(i0,i1-1,i2)))+
     &              fac1*fac2*(
     &              mu0(i0,i1+1,i2)*(u1(i0,i1+1,i2)-u1(i0,i1+1,i2-1))-
     &              mu0(i0,i1,i2)*(u1(i0,i1,i2)-u1(i0,i1,i2-1)))+
     &              fac0**2.d0*(
     &              mu1(i0+1,i1,i2)*(u2(i0+1,i1,i2)-u2(i0,i1,i2))-
     &              mu1(i0,i1,i2)*(u2(i0,i1,i2)-u2(i0-1,i1,i2)))+
     &              fac0*fac2*(
     &              mu1(i0+1,i1,i2)*(u0(i0+1,i1,i2)-u0(i0+1,i1,i2-1))-
     &              mu1(i0,i1,i2)*(u0(i0,i1,i2)-u0(i0,i1,i2-1))))+
     &              rho*u2(i0,i1,i2) + gamma*v2(i0,i1,i2)
            enddo
         enddo
      enddo
c
      return
      end
c

