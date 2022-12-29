c ---------------------------------------------------------------------
c
c Copyright (c) 2011 - 2022 by the IBAMR developers
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
      REAL a_avg4, h_avg4
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
                mu_upper = h_avg4(mu(i0,i1),mu(i0+1,i1),
     &                          mu(i0,i1+1),mu(i0+1,i1+1))
                mu_lower = h_avg4(mu(i0,i1),mu(i0-1,i1),
     &                       mu(i0,i1+1),mu(i0-1,i1+1))
            else
                mu_upper = a_avg4(mu(i0,i1),mu(i0+1,i1),
     &                          mu(i0,i1+1),mu(i0+1,i1+1))
                mu_lower = a_avg4(mu(i0,i1),mu(i0-1,i1),
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
               mu_upper = h_avg4(mu(i0,i1),mu(i0+1,i1),
     &                          mu(i0,i1+1),mu(i0+1,i1+1))

               mu_lower = h_avg4(mu(i0,i1),mu(i0+1,i1),
     &                          mu(i0,i1-1),mu(i0+1,i1-1))
            else
               mu_upper = a_avg4(mu(i0,i1),mu(i0+1,i1),
     &                          mu(i0,i1+1),mu(i0+1,i1+1))

               mu_lower = a_avg4(mu(i0,i1),mu(i0+1,i1),
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
