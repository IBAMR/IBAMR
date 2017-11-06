c
c     Routines to compute quantities related to variable coefficient
c     generalized Laplace operators.
c
c     Created on 27 May 2010 by Thomas Fai
c     Modified on 12 Feb 2016 by Nishant Nangia
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
c     Compute the arithmetic average of four inputs
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      REAL function a_avg(a0,a1,a2,a3)
      implicit none
      REAL a0,a1,a2,a3
      a_avg = 0.25d0*(a0+a1+a2+a3)
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the harmonic average of four inputs
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      REAL function h_avg(a0,a1,a2,a3)
      implicit none
      REAL a0,a1,a2,a3
      REAL dmr
      dmr = 1.d0/a0+1.d0/a1+1.d0/a2+1.d0/a3
      h_avg = 4.d0/dmr
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes (f0,f1) = alpha div mu (grad (u0,u1) + grad(u0,u1)^T) + beta (u0,u1) +
c     gamma (v0,v1).
c
c     Computes the side-centered variable coefficient generalized
c     Laplacian, with node-centered coefficient mu and side-centered
c     vector fields (u0,u1) and (v0,v1).
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine stosvclaplace2d(
     &     f0,f1,f_gcw,
     &     alpha,beta,
     &     mu,mu_gcw,
     &     rho0,rho1,rho_gcw,
     &     u0,u1,u_gcw,
     &     gamma,
     &     v0,v1,v_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
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
      INTEGER f_gcw,mu_gcw,rho_gcw,u_gcw,v_gcw,var_rho
      INTEGER use_harmonic_interp

      REAL alpha,beta,gamma

      REAL mu(NODE2d(ilower,iupper,mu_gcw))

      REAL rho0(SIDE2d0(ilower,iupper,rho_gcw))
      REAL rho1(SIDE2d1(ilower,iupper,rho_gcw))
      
      REAL u0(SIDE2d0(ilower,iupper,u_gcw))
      REAL u1(SIDE2d1(ilower,iupper,u_gcw))

      REAL v0(SIDE2d0(ilower,iupper,v_gcw))
      REAL v1(SIDE2d1(ilower,iupper,v_gcw))

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL f0(SIDE2d0(ilower,iupper,f_gcw))
      REAL f1(SIDE2d1(ilower,iupper,f_gcw))
c
c     Local variables.
c
      INTEGER i0,i1
      REAL    fac0,fac1,rho
      REAL    mu_lower,mu_upper
c
c     Compute the discrete divergence of mu (grad (u0,u1) + grad (u0,u1)^T).
c
      fac0 = 1.d0/(dx(0))
      fac1 = 1.d0/(dx(1))

      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0+1

            rho = beta
            if (var_rho .eq. 1) then
               rho = rho0(i0,i1)*beta
            endif

            if (use_harmonic_interp .eq. 1) then
                mu_upper = h_avg(mu(i0,i1),mu(i0+1,i1),
     &                          mu(i0,i1+1),mu(i0+1,i1+1))
                mu_lower = h_avg(mu(i0,i1),mu(i0-1,i1),
     &                       mu(i0,i1+1),mu(i0-1,i1+1))
            else
                mu_upper = a_avg(mu(i0,i1),mu(i0+1,i1),
     &                          mu(i0,i1+1),mu(i0+1,i1+1))
                mu_lower = a_avg(mu(i0,i1),mu(i0-1,i1),
     &                       mu(i0,i1+1),mu(i0-1,i1+1))
            endif

            f0(i0,i1) = alpha*(
     &           2.d0*fac0**2.d0*(
     &           mu_upper*
     &           (u0(i0+1,i1)-u0(i0,i1))-
     &           mu_lower*
     &           (u0(i0,i1)-u0(i0-1,i1)))+
     &           fac1**2.d0*(mu(i0,i1+1)*(u0(i0,i1+1)-u0(i0,i1))-
     &           mu(i0,i1)*(u0(i0,i1)-u0(i0,i1-1)))+
     &           fac0*fac1*(mu(i0,i1+1)*(u1(i0,i1+1)-u1(i0-1,i1+1))-
     &           mu(i0,i1)*(u1(i0,i1)-u1(i0-1,i1)))) +
     &           rho*u0(i0,i1) + gamma*v0(i0,i1)
         enddo
      enddo

      do i1 = ilower1,iupper1+1
         do i0 = ilower0,iupper0

            rho = beta
            if (var_rho .eq. 1) then
               rho = rho1(i0,i1)*beta
            endif

            if (use_harmonic_interp .eq. 1) then
               mu_upper = h_avg(mu(i0,i1),mu(i0+1,i1),
     &                          mu(i0,i1+1),mu(i0+1,i1+1))

               mu_lower = h_avg(mu(i0,i1),mu(i0+1,i1),
     &                          mu(i0,i1-1),mu(i0+1,i1-1))
            else
               mu_upper = a_avg(mu(i0,i1),mu(i0+1,i1),
     &                          mu(i0,i1+1),mu(i0+1,i1+1))

               mu_lower = a_avg(mu(i0,i1),mu(i0+1,i1),
     &                          mu(i0,i1-1),mu(i0+1,i1-1))
            endif


            f1(i0,i1) = alpha*(
     &           2.d0*fac1**2.d0*(
     &           mu_upper*
     &           (u1(i0,i1+1)-u1(i0,i1))-
     &           mu_lower*
     &           (u1(i0,i1)-u1(i0,i1-1)))+
     &           fac0**2.d0*(mu(i0+1,i1)*(u1(i0+1,i1)-u1(i0,i1))-
     &           mu(i0,i1)*(u1(i0,i1)-u1(i0-1,i1)))+
     &           fac0*fac1*(mu(i0+1,i1)*(u0(i0+1,i1)-u0(i0+1,i1-1))-
     &           mu(i0,i1)*(u0(i0,i1)-u0(i0,i1-1)))) +
     &           rho*u1(i0,i1) + gamma*v1(i0,i1)
         enddo
      enddo
c
      return
      end
c
