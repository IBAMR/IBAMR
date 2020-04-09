c ---------------------------------------------------------------------
c
c Copyright (c) 2017 - 2017 by the IBAMR developers
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

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes P = mut*(grad U + grad U^T):grad_U
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine production_2d(
     &     P,P_gc,
     &     mut,mut_gc,
     &     U0,U1,U_gc0,U_gc1,
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
      INTEGER U_gc0,U_gc1
      INTEGER mut_gc

      REAL mut(CELL2d(ilower,iupper,mut_gc))
      REAL U0(SIDE2d0VECG(ilower,iupper,U_gc))
      REAL U1(SIDE2d1VECG(ilower,iupper,U_gc))
      REAL dx(0:NDIM-1)
c
c     Output.
      INTEGER P_gc
      REAL P(CELL2d(ilower,iupper,P_gc))
c
c     Local variables.
c
      INTEGER i0,i1
      REAL    fac0,fac1,fac2,fac3
      REAL del_u_del_x,del_v_del_y,del_u_del_y,del_v_del_x
c
c     Compute the Production term of k and omega equation.
c
      fac0 = 1.d0/dx(0)
      fac1 = 1.d0/dx(1)
      fac2 = 1.d0/(4.d0*dx(1))
      fac3 = 1.d0/(4.d0*dx(0))
c
      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0
            del_u_del_x = fac0*(U0(i0+1,i1)-U0(i0,i1))
            del_v_del_y = fac1*(U1(i0,i1+1)-U1(i0,i1))
            del_u_del_y = fac2*((U0(i0+1,i1+1)-U0(i0+1,i1-1))
     &               +(U0(i0,i1+1)-U0(i0,i1-1)))
            del_v_del_x = fac3*((U1(i0+1,i1+1)-U1(i0-1,i1+1))
     &               +(U1(i0+1,i1)-U1(i0-1,i1)))
            P(i0,i1) = mut(i0,i1)*(2.d0*(del_u_del_x**2.d0+
     &      del_v_del_y**2.d0)+(del_u_del_y+del_v_del_x)**2.d0)
         enddo
      enddo
c
      return
      end
cccccc
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Calculates blending function F1
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
      subroutine sst_blending_fcn_2d(
     &     F1,F1_gc,
     &     k,k_gc,
     &     w,w_gc,
     &     mu,mu_gc,
     &     rho,rho_gc,
     &     distance, distance_gc,
     &     beta_star,
     &     sigma_w2,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     dx)
c
c
      implicit none
c
c     Input
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER k_gc, w_gc, rho_gc, mu_gc, distance_gc
c
      REAL beta_star, sigma_w2
c
      REAL k(CELL2d(ilower,iupper,k_gc))
      REAL w(CELL2d(ilower,iupper,w_gc))
      REAL mu(CELL2d(ilower,iupper,mu_gc))
      REAL rho(CELL2d(ilower,iupper,rho_gc))
      REAL distance(CELL2d(ilower,iupper,distance_gc))
      REAL dx(0:NDIM-1)
c
c     Output.
c
      INTEGER F1_gc
      REAL F1(CELL2d(ilower,iupper,F1_gc))
c
c     Local variables.
      INTEGER i0,i1
      REAL grad_k_dot_w
      REAL cross_diffusion,first_term
      REAL second_term,gamma
      REAL fac0,fac1
c
      fac0 = 1.d0/(4.d0*dx(0)*dx(0))
      fac1 = 1.d0/(4.d0*dx(1)*dx(1))
c
      do i1 = ilower1,iupper1
        do i0 = ilower0,iupper0
            grad_k_dot_w = (fac0*(k(i0+1,i1)-k(i0-1,i1))*(w(i0+1,i1)
     &                     -w(i0-1,i1)))+ (fac1*(k(i0,i1+1)-k(i0,i1-1))
     &                     *(w(i0,i1+1)-w(i0,i1-1)))
            cross_diffusion = max((2.d0*rho(i0,i1)*sigma_w2
     &                        *grad_k_dot_w/w(i0,i1)), 1E-10)
            second_term = 4.d0*rho(i0,i1)*sigma_w2*k(i0,i1)
     &                    /(cross_diffusion*(distance(i0,i1)**2.d0))
            first_term = max((sqrt(k(i0,i1))/(beta_star*w(i0,i1)
     &                   *distance(i0,i1))),(500.d0*mu(i0,i1)
     &                   /(rho(i0,i1)*w(i0,i1)*(distance(i0,i1)
     &                   **2.d0))))
            gamma = min(first_term, second_term)
            F1(i0,i1) = tanh(gamma**4.d0)
        enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Calculates F2
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
      subroutine sst_f2_fcn_2d(
     &     F2,F2_gc,
     &     k,k_gc,
     &     w,w_gc,
     &     mu,mu_gc,
     &     rho,rho_gc,
     &     distance, distance_gc,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     beta_star)
c
c
      implicit none
c
c     Input
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER k_gc, w_gc, rho_gc, mu_gc, distance_gc
c
      REAL beta_star
c
      REAL k(CELL2d(ilower,iupper,k_gc))
      REAL w(CELL2d(ilower,iupper,w_gc))
      REAL mu(CELL2d(ilower,iupper,mu_gc))
      REAL rho(CELL2d(ilower,iupper,rho_gc))
      REAL distance(CELL2d(ilower,iupper,distance_gc))
c
c     Output.
c
      INTEGER F2_gc
      REAL F2(CELL2d(ilower,iupper,F2_gc))
c
c    Local variables.
      INTEGER i0,i1
      REAL gamma

      do i1 = ilower1,iupper1
        do i0 = ilower0,iupper0
            gamma = max((2.d0*sqrt(k(i0,i1))/(beta_star*w(i0,i1)
     &                   *distance(i0,i1))),(500.d0*mu(i0,i1)
     &                   /(rho(i0,i1)*w(i0,i1)*(distance(i0,i1)
     &                   **2.d0))))
            F2(i0,i1) = tanh(gamma**2.d0)
        enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Calculates turbulent viscosity, mu_t, for SST k-w model
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
      subroutine sst_mu_t_fcn_2d(
     &     mu_t,mu_t_gc,
     &     k,k_gc,
     &     w,w_gc,
     &     rho,rho_gc,
     &     F2, F2_gc,
     &     P, P_gc,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     a1)
c
c
      implicit none
c
c     Input
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER k_gc, w_gc, rho_gc, P_gc, F2_gc
c
      REAL a1
c
      REAL k(CELL2d(ilower,iupper,k_gc))
      REAL w(CELL2d(ilower,iupper,w_gc))
      REAL P(CELL2d(ilower,iupper,P_gc))
      REAL rho(CELL2d(ilower,iupper,rho_gc))
      REAL F2(CELL2d(ilower,iupper,F2_gc))
c
c     Output.
c
      INTEGER mu_t_gc
      REAL mu_t(CELL2d(ilower,iupper,mu_t_gc))
c
c    Local variables.
      INTEGER i0,i1
      REAL strain_rate_mag

      do i1 = ilower1,iupper1
        do i0 = ilower0,iupper0
            strain_rate_mag = sqrt(P(i0,i1)/(2.d0*mu_t(i0,i1)))
            mu_t(i0,i1) = a1*rho(i0,i1)*k(i0,i1)
     &                  / max(a1*w(i0,i1), sqrt(2.d0)
     &                  * strain_rate_mag*F2(i0,i1))
        enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes G_b = - (1/(rho*sigma_t)*mut*grad_rho . g
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
      subroutine sst_k_eqn_buoyancy_2d(
     &     P,P_gc,
     &     mut,mut_gc,
     &     rho,rho_gc,
     &     gx,
     &     gy,
     &     sigma_t,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     dx)
c
c
      implicit none
c
c     Input
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER rho_gc
      INTEGER mut_gc
c
      REAL gx,gy,sigma_t
c
      REAL mut(CELL2d(ilower,iupper,mut_gc))
      REAL rho(CELL2d(ilower,iupper,rho_gc))
      REAL dx(0:NDIM-1)
c
c     Output.
c
      INTEGER P_gc
      REAL P(CELL2d(ilower,iupper,P_gc))
c
c    Local variables.
c
      INTEGER i0,i1
      REAL    fac0,fac1

      fac0 = 1.d0/(2.d0*dx(0))
      fac1 = 1.d0/(2.d0*dx(1))
c
c
      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0
           P(i0,i1)=P(i0,i1)-((mut(i0,i1)/(rho(i0,i1)*sigma_t))
     &              *((fac0*gx*(rho(i0+1,i1)-rho(i0-1,i1)))
     &              +(fac1*gy*(rho(i0,i1+1)-rho(i0,i1-1)))))
         enddo
      enddo
c
      return
      end
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes P_k = min(G, 10*beta_star*k*w)
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
      subroutine sst_k_eqn_production_2d(
     & k_f,k_f_gc,
     & P,P_gc,
     & k,k_gc,
     & w,w_gc,
     & ilower0,iupper0,
     & ilower1,iupper1,
     & beta_star)
c
      implicit none
c
c     Input
      INTEGER ilower0, iupper0
      INTEGER ilower1, iupper1
      INTEGER k_gc, w_gc, P_gc

c
      REAL beta_star
c
      REAL k(CELL2d(ilower,iupper,k_gc))
      REAL w(CELL2d(ilower,iupper,w_gc))
      REAL P(CELL2d(ilower,iupper,P_gc))
c
c     output
c
c
      INTEGER k_f_gc
      REAL k_f(CELL2d(ilower,iupper,k_f_gc))
c
c
c    Local variables.
      INTEGER i0,i1
c
      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0
          k_f(i0,i1) = min(P(i0,i1),(10.d0*beta_star*k(i0,i1)*w(i0,i1)))
        enddo
      enddo
c
      return
      end
c
c Fortran routines for omega equation source terms
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes alpha*G/nu_t
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
      subroutine sst_w_eqn_production_2d(
     & w_f,w_f_gc,
     & P,P_gc,
     & mut,mut_gc,
     & rho,rho_gc,
     & F1,F1_gc,
     & ilower0,iupper0,
     & ilower1,iupper1,
     & alpha_1,alpha_2)
c
      implicit none
c
c     Input
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER F1_gc,rho_gc,mut_gc,P_gc
c
      REAL alpha_1,alpha_2
c
      REAL F1(CELL2d(ilower,iupper,F1_gc))
      REAL rho(CELL2d(ilower,iupper,rho_gc))
      REAL mut(CELL2d(ilower,iupper,mut_gc))
      REAL P(CELL2d(ilower,iupper,P_gc))
c
c     output
c
      INTEGER w_f_gc
      REAL w_f(CELL2d(ilower,iupper,w_f_gc))
c
c
c    Local variables.
      INTEGER i0,i1
      REAL alpha
      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0
            alpha = F1(i0,i1)*alpha_1 + (1.d0-F1(i0,i1))*alpha_2
            w_f(i0,i1) = P(i0,i1)*alpha*rho(i0,i1)/mut(i0,i1)
         enddo
      enddo
c
      return
      end
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes 2*(1-F1)*rho*(sigma_w2/omega)*grad k \dot \grad w
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
      subroutine sst_w_eqn_crossdiffusion_2d(
     & w_f,w_f_gc,
     & rho,rho_gc,
     & F1,F1_gc,
     & k,k_gc,
     & w,w_gc,
     & sigma_w2,
     & ilower0,iupper0,
     & ilower1,iupper1,
     & dx)
c
      implicit none
c
c     input variables
      REAL sigma_w2
      INTEGER F1_gc,rho_gc,k_gc,w_gc
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
c
      REAL F1(CELL2d(ilower,iupper,F1_gc))
      REAL rho(CELL2d(ilower,iupper,rho_gc))
      REAL k(CELL2d(ilower,iupper,k_gc))
      REAL w(CELL2d(ilower,iupper,w_gc))
      REAL dx(0:NDIM-1)

c     output variables
c
      INTEGER w_f_gc
      REAL w_f(CELL2d(ilower,iupper,w_f_gc))
c
c    Local variables.
      INTEGER i0,i1
      REAL grad_k_grad_w
      REAL fac0,fac1
c
      fac0 = 1.d0/(4.d0*dx(0)*dx(0))
      fac1 = 1.d0/(4.d0*dx(1)*dx(1))
c
c
      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0
         grad_k_grad_w = (fac0*((k(i0+1,i1)-k(i0-1,i1))*(w(i0+1,i1)
     &                     -w(i0-1,i1))))+(fac1*((k(i0,i1+1)-k(i0,i1-1))
     &                     *(w(i0,i1+1)-w(i0,i1-1))))
         w_f(i0,i1) = w_f(i0,i1) + ((2.d0*(1.d0-F1(i0,i1))*rho(i0,i1)
     &                *sigma_w2*grad_k_grad_w)/w(i0,i1))
         enddo
      enddo
c
      return
      end
