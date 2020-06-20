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
define(NDIM,3)dnl
define(REAL,`double precision')dnl
define(INTEGER,`integer')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim3d.i)dnl

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes P = mut*grad_U:(grad U + grad U^T)
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine production_3d(
     &     P,P_gc,
     &     mut,mut_gc,
     &     U0,U1,U2,U_gc0,U_gc1,U_gc2,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2, iupper2,
     &     dx)
c
      implicit none
c
c     Input.
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2
      INTEGER U_gc0,U_gc1,U_gc2
      INTEGER mut_gc

      REAL mut(CELL3d(ilower,iupper,mut_gc))

      REAL U0(SIDE3d0VECG(ilower,iupper,U_gc))
      REAL U1(SIDE3d1VECG(ilower,iupper,U_gc))
      REAL U2(SIDE3d2VECG(ilower,iupper,U_gc))
      REAL dx(0:NDIM-1)
c
c     Output.
      INTEGER P_gc

      REAL P(CELL3d(ilower,iupper,P_gc))
c
c     Local variables.
c
      INTEGER i0,i1,i2
      REAL    fac0,fac1,fac2,fac3,fac4,fac5
      REAL del_u_del_x,del_v_del_y,del_w_del_z,del_u_del_y
      REAL del_v_del_x,del_u_del_z,del_v_del_z,del_w_del_y
      REAL del_w_del_x
c
c     Compute the Production term of k and omega equation.
c
      fac0 = 1.d0/dx(0)
      fac1 = 1.d0/dx(1)
      fac2 = 1.d0/dx(2)
      fac3 = 1.d0/(4.d0*dx(1))
      fac4 = 1.d0/(4.d0*dx(0))
      fac5 = 1.d0/(4.d0*dx(2))

      do i2 = ilower2,iupper2
        do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0
            del_u_del_x = fac0*(U0(i0+1,i1,i2)-U0(i0,i1,i2))
            del_v_del_y = fac1*(U1(i0,i1+1,i2)-U1(i0,i1,i2))
            del_w_del_z = fac2*(U2(i0,i1,i2+1)-U2(i0,i1,i2))
            del_u_del_y = fac3*((U0(i0+1,i1+1,i2)-U0(i0+1,i1-1,i2))
     &               +(U0(i0,i1+1,i2)-U0(i0,i1-1,i2)))
            del_v_del_x = fac4*((U1(i0+1,i1+1,i2)-U1(i0-1,i1+1,i2))
     &               +(U1(i0+1,i1,i2)-U1(i0-1,i1,i2)))
            del_w_del_x = fac4*((U2(i0+1,i1,i2+1)-U2(i0-1,i1,i2+1))
     &               +(U2(i0+1,i1,i2)-U2(i0-1,i1,i2)))
            del_u_del_z = fac5*((U0(i0+1,i1,i2+1)-U0(i0+1,i1,i2-1))
     &               +(U0(i0,i1,i2+1)-U0(i0,i1,i2-1)))
            del_v_del_z = fac5*((U1(i0,i1+1,i2+1)-U1(i0,i1+1,i2-1))
     &               +(U1(i0,i1,i2+1)-U1(i0,i1,i2-1)))
            del_w_del_y = fac3*((U2(i0,i1+1,i2+1)-U2(i0,i1-1,i2+1))
     &               +(U2(i0,i1+1,i2)-U2(i0,i1-1,i2)))
c
            P(i0,i1,i2) = mut(i0,i1,i2)*((2.d0*(del_u_del_x**2.d0+
     &      del_v_del_y**2.d0+del_w_del_z**2.d0))+(del_u_del_y
     &       +del_v_del_x)**2.d0+(del_u_del_z+del_w_del_x)**2.d0
     &       +(del_v_del_z+del_w_del_y)**2.d0)
         enddo
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
      subroutine sst_blending_fcn_3d(
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
     &     ilower2,iupper2,
     &     dx)
c
c
      implicit none
c
c     Input
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2
      INTEGER k_gc, w_gc, rho_gc, mu_gc, distance_gc
c
      REAL beta_star, sigma_w2
c
      REAL k(CELL3d(ilower,iupper,k_gc))
      REAL w(CELL3d(ilower,iupper,w_gc))
      REAL mu(CELL3d(ilower,iupper,mu_gc))
      REAL rho(CELL3d(ilower,iupper,rho_gc))
      REAL distance(CELL3d(ilower,iupper,distance_gc))
      REAL dx(0:NDIM-1)
c
c     Output.
c
      INTEGER F1_gc
      REAL F1(CELL3d(ilower,iupper,F1_gc))
c
c    Local variables.
      INTEGER i0,i1,i2
      REAL grad_k_dot_w
      REAL cross_diffusion,first_term
      REAL second_term,gamma
      REAL fac0,fac1,fac2
c
      fac0 = 1.d0/(4.d0*dx(0)*dx(0))
      fac1 = 1.d0/(4.d0*dx(1)*dx(1))
      fac2 = 1.d0/(4.d0*dx(2)*dx(2))
c
      do i2 = ilower2,iupper2
        do i1 = ilower1,iupper1
          do i0 = ilower0,iupper0
            grad_k_dot_w = (fac0*(k(i0+1,i1,i2)-k(i0-1,i1,i2))*
     &      (w(i0+1,i1,i2)-w(i0-1,i1,i2)))+ (fac1*(k(i0,i1+1,i2)-
     &      k(i0,i1-1,i2))*(w(i0,i1+1,i2)-w(i0,i1-1,i2)))
     &      +(fac2*(k(i0,i1,i2+1)-k(i0,i1,i2-1))*(w(i0,i1,i2+1)
     &      -w(i0,i1,i2-1)))
     &
            cross_diffusion = max((2.d0*rho(i0,i1,i2)*sigma_w2
     &                        *grad_k_dot_w/w(i0,i1,i2)), 1E-10)
            second_term = 4.d0*rho(i0,i1,i2)*sigma_w2*k(i0,i1,i2)
     &                    /(cross_diffusion*(distance(i0,i1,i2)**2.d0))
            first_term = max((sqrt(k(i0,i1,i2))/(beta_star*w(i0,i1,i2)
     &                   *distance(i0,i1,i2))),(500.d0*mu(i0,i1,i2)
     &                   /(rho(i0,i1,i2)*w(i0,i1,i2)*(distance(i0,i1,i2)
     &                   **2.d0))))
            gamma = min(first_term, second_term)
            F1(i0,i1,i2) = tanh(gamma**4.d0)
          enddo
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
      subroutine sst_f2_fcn_3d(
     &     F2,F2_gc,
     &     k,k_gc,
     &     w,w_gc,
     &     mu,mu_gc,
     &     rho,rho_gc,
     &     distance, distance_gc,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
     &     beta_star)
c
c
      implicit none
c
c     Input
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2
      INTEGER k_gc, w_gc, rho_gc, mu_gc, distance_gc
c
      REAL beta_star
c
      REAL k(CELL3d(ilower,iupper,k_gc))
      REAL w(CELL3d(ilower,iupper,w_gc))
      REAL mu(CELL3d(ilower,iupper,mu_gc))
      REAL rho(CELL3d(ilower,iupper,rho_gc))
      REAL distance(CELL3d(ilower,iupper,distance_gc))
c
c     Output.
c
      INTEGER F2_gc
      REAL F2(CELL3d(ilower,iupper,F2_gc))
c
c    Local variables.
      INTEGER i0,i1,i2
      REAL gamma
      do i2 = ilower2,iupper2
        do i1 = ilower1,iupper1
          do i0 = ilower0,iupper0
            gamma = max((2.d0*sqrt(k(i0,i1,i2))/(beta_star*w(i0,i1,i2)
     &                   *distance(i0,i1,i2))),(500.d0*mu(i0,i1,i2)
     &                   /(rho(i0,i1,i2)*w(i0,i1,i2)*(distance(i0,i1,i2)
     &                   **2.d0))))
            F2(i0,i1,i2) = tanh(gamma**2.d0)
          enddo
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
      subroutine sst_mu_t_fcn_3d(
     &     mu_t,mu_t_gc,
     &     k,k_gc,
     &     w,w_gc,
     &     rho,rho_gc,
     &     F2, F2_gc,
     &     P, P_gc,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
     &     a1)
c
c
      implicit none
c
c     Input
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2
      INTEGER k_gc, w_gc, rho_gc, P_gc, F2_gc
c
      REAL a1
c
      REAL k(CELL3d(ilower,iupper,k_gc))
      REAL w(CELL3d(ilower,iupper,w_gc))
      REAL P(CELL3d(ilower,iupper,P_gc))
      REAL rho(CELL3d(ilower,iupper,rho_gc))
      REAL F2(CELL3d(ilower,iupper,F2_gc))
c
c     Output.
c
      INTEGER mu_t_gc
      REAL mu_t(CELL3d(ilower,iupper,mu_t_gc))
c
c    Local variables.
      INTEGER i0,i1,i2
      REAL strain_rate_mag
      do i2 = ilower2,iupper2
        do i1 = ilower1,iupper1
          do i0 = ilower0,iupper0
            strain_rate_mag = sqrt(P(i0,i1,i2)/(2.d0*mu_t(i0,i1,i2)))
            mu_t(i0,i1,i2) = a1*rho(i0,i1,i2)*k(i0,i1,i2)
     &                    /max(a1*w(i0,i1,i2), sqrt(2.d0)
     &                    *strain_rate_mag*F2(i0,i1,i2))
          enddo
        enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes G_b = - (1/(rho* sigma_t)*mut*grad_rho . g
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
      subroutine sst_k_eqn_buoyancy_3d(
     &     P,P_gc,
     &     mut,mut_gc,
     &     rho,rho_gc,
     &     gx,
     &     gy,
     &     gz,
     &     sigma_t,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
     &     dx)
c
c
      implicit none
c
c     Input
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2
      INTEGER rho_gc
      INTEGER mut_gc
c
      REAL gx,gy,gz,sigma_t

      REAL mut(CELL3d(ilower,iupper,mut_gc))
      REAL rho(CELL3d(ilower,iupper,rho_gc))
      REAL dx(0:NDIM-1)
c
c     Output.
c
      INTEGER P_gc
      REAL P(CELL3d(ilower,iupper,P_gc))
c
c    Local variables.
c
      INTEGER i0,i1,i2
      REAL    fac0,fac1,fac2

      fac0 = 1.d0/dx(0)
      fac1 = 1.d0/dx(1)
      fac2 = 1.d0/dx(2)

      do i2 = ilower2,iupper2
        do i1 = ilower1,iupper1
          do i0 = ilower0,iupper0
            P(i0,i1,i2)=P(i0,i1,i2)-((mut(i0,i1,i2)/(2.d0*rho(i0,i1,i2)
     &           *sigma_t))*((fac0*gx*(rho(i0+1,i1,i2)-rho(i0-1,i1,i2)))
     &              +(fac1*gy*(rho(i0,i1+1,i2)-rho(i0,i1-1,i2)))
     &              +(fac2*gz*(rho(i0,i1,i2+1)-rho(i0,i1,i2-1)))))
          enddo
        enddo
      enddo

      return
      end

c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes P_k = min (G, 10*beta_star*k*w)
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
      subroutine sst_k_eqn_production_3d(
     & k_f,k_f_gc,
     & P,P_gc,
     & k,k_gc,
     & w,w_gc,
     & ilower0,iupper0,
     & ilower1,iupper1,
     & ilower2,iupper2,
     & beta_star)

      implicit none

c     Input
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2
      INTEGER k_gc,w_gc

      REAL beta_star

      REAL k(CELL3d(ilower,iupper,k_gc))
      REAL w(CELL3d(ilower,iupper,w_gc))

c     output

      INTEGER P_gc
      REAL P(CELL3d(ilower,iupper,P_gc))

      INTEGER k_f_gc
      REAL k_f(CELL3d(ilower,iupper,k_f_gc))


c    Local variables.
      INTEGER i0,i1,i2
      do i2 = ilower2,iupper2
        do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0
          k_f(i0,i1,i2) = Min(P(i0,i1,i2),(10.d0*beta_star*k(i0,i1,i2)
     &       *w(i0,i1,i2)))
          enddo
        enddo
      enddo

      return
      end

c Fortran routines for omega equation source terms
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes alpha*G/nu_t
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
      subroutine sst_w_eqn_production_3d(
     & w_f,w_f_gc,
     & P,P_gc,
     & mut,mut_gc,
     & rho,rho_gc,
     & F1,F1_gc,
     & ilower0,iupper0,
     & ilower1,iupper1,
     & ilower2,iupper2,
     & alpha_1,alpha_2)

      implicit none

c     Input
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2
      INTEGER F1_gc,rho_gc,mut_gc

      REAL alpha_1,alpha_2

      REAL F1(CELL3d(ilower,iupper,F1_gc))
      REAL rho(CELL3d(ilower,iupper,rho_gc))
      REAL mut(CELL3d(ilower,iupper,mut_gc))

c     output

      INTEGER P_gc
      REAL P(CELL3d(ilower,iupper,P_gc))

      INTEGER w_f_gc
      REAL w_f(CELL3d(ilower,iupper,w_f_gc))


c    Local variables.
      INTEGER i0,i1,i2
      REAL alpha
      do i2 = ilower2,iupper2
        do i1 = ilower1,iupper1
          do i0 = ilower0,iupper0
            alpha=(F1(i0,i1,i2)*alpha_1)+((1.d0-F1(i0,i1,i2))*alpha_2)
          w_f(i0,i1,i2)=(alpha*rho(i0,i1,i2)/mut(i0,i1,i2))*P(i0,i1,i2)
          enddo
        enddo
      enddo

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
      subroutine sst_w_eqn_crossdiffusion_3d(
     & w_f,w_f_gc,
     & rho,rho_gc,
     & F1,F1_gc,
     & k,k_gc,
     & w,w_gc,
     & sigma_w2,
     & ilower0,iupper0,
     & ilower1,iupper1,
     & ilower2,iupper2,
     & dx)

      implicit none

c     input variables
      REAL sigma_w2
      INTEGER F1_gc,rho_gc,k_gc,w_gc
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2

      REAL F1(CELL3d(ilower,iupper,F1_gc))
      REAL rho(CELL3d(ilower,iupper,rho_gc))
      REAL k(CELL3d(ilower,iupper,k_gc))
      REAL w(CELL3d(ilower,iupper,w_gc))
      REAL dx(0:NDIM-1)

c     output variables

      INTEGER w_f_gc
      REAL w_f(CELL3d(ilower,iupper,w_f_gc))

c    Local variables.
      INTEGER i0,i1,i2
      REAL grad_k_grad_w
      REAL fac0,fac1,fac2

      fac0 = 1.d0/(4.d0*dx(0)*dx(0))
      fac1 = 1.d0/(4.d0*dx(1)*dx(1))
      fac2 = 1.d0/(4.d0*dx(2)*dx(2))

      do i2 = ilower2,iupper2
        do i1 = ilower1,iupper1
          do i0 = ilower0,iupper0
         grad_k_grad_w = (fac0*((k(i0+1,i1,i2)-k(i0-1,i1,i2))
     &          *(w(i0+1,i1,i2)-w(i0-1,i1,i2))))
     &    +(fac1*((k(i0,i1+1,i2)-k(i0,i1-1,i2))
     &                     *(w(i0,i1+1,i2)-w(i0,i1-1,i2))))
     &    +(fac2*((k(i0,i1,i2+1)-k(i0,i1,i2-1))
     &                     *(w(i0,i1,i2+1)-w(i0,i1,i2-1))))
          w_f(i0,i1,i2) = w_f(i0,i1,i2)+((2.d0*(1.d0-F1(i0,i1,i2))
     &             *rho(i0,i1,i2)*sigma_w2*grad_k_grad_w)/w(i0,i1,i2))
          enddo
        enddo
      enddo

      return
      end
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute tau_w = rho*u_tau*u_star. This is not implemented.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
      subroutine wall_shear_stress_3d(
     & tau_w,tau_w_gc,
     & U0,U_gc0,
     & U1,U_gc1,
     & U2,U_gc2,
     & k,k_gc,
     & rho,rho_gc,
     & mu,mu_gc,
     & U_tau0,U_tau_gc0,
     & U_tau1,U_tau_gc1,
     & U_tau2,U_tau_gc2,
     & yplus0,yplus_gc0,
     & yplus1,yplus_gc1,
     & yplus2,yplus_gc2,
     & kappa,
     & beta_star,
     & B,
     & wall_location_index,
     & ilower0,iupper0,
     & ilower1,iupper1,
     & ilower2,iupper2,
     & dx)
c
      implicit none
c
c     input variables
      REAL kappa, beta_star, B
      INTEGER U_gc0, U_gc1, U_gc2
      INTEGER k_gc, rho_gc, mu_gc, wall_location_index
      INTEGER yplus_gc0, yplus_gc1, yplus_gc2
      INTEGER U_tau_gc0, U_tau_gc1, U_tau_gc2
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2
c
      REAL U0(SIDE3d0VECG(ilower,iupper,U_gc))
      REAL U1(SIDE3d1VECG(ilower,iupper,U_gc))
      REAL U2(SIDE3d2VECG(ilower,iupper,U_gc))

      REAL U_tau0(SIDE3d0VECG(ilower,iupper,U_tau_gc))
      REAL U_tau1(SIDE3d1VECG(ilower,iupper,U_tau_gc))
      REAL U_tau2(SIDE3d2VECG(ilower,iupper,U_tau_gc))

      REAL yplus0(SIDE3d0VECG(ilower,iupper,yplus_gc))
      REAL yplus1(SIDE3d1VECG(ilower,iupper,yplus_gc))
      REAL yplus2(SIDE3d2VECG(ilower,iupper,yplus_gc))

      REAL k(CELL3d(ilower,iupper,k_gc))
      REAL rho(CELL3d(ilower,iupper,rho_gc))
      REAL mu(CELL3d(ilower,iupper,mu_gc))
      REAL dx(0:NDIM-1)

c     output variables
c
      INTEGER tau_w_gc
      REAL tau_w(EDGE3d0(ilower,iupper,tau_w_gc))

c     Local variables.
      INTEGER i0,i1
      REAL U_tau_vis, U_tau_log, U_tau
      REAL U_star_vis, U_star_log, U_star
      REAL mu_nc, rho_nc
      REAL U_mag, y_plus
c
c
c      do i2 = ilower2, iupper2+1
c        do i1 = ilower1,iupper1+1
c          do i0 = ilower0,iupper0+1
c          if (((wall_location_index .eq. 0) .and. (i0 .eq. ilower0)).or.
c     &  ((wall_location_index .eq. 1) .and.
c     & (i0 .eq. (number_of_indices(0) + 1))) .or.
c     &  ((wall_location_index .eq. 2) .and.  (i1 .eq. 0))
c     & .or. ((wall_location_index .eq. 3) .and.  (i1 .eq.
c     & (number_of_indices(1) + 1))) .or. ((wall_location_index .eq. 4)
c     &  .and. (number_of_indices(2) = ilower2)) .or.
c     &  ((wall_location_index .eq. 5) .and.  (number_of_indices(2) =
c     &  iupper2+1))) then
c            if ((wall_location_index .eq. 0) .or.
c     &            ((wall_location_index .eq. 1))) then
c               U_mag = sqrt(U1(i0,i1)*U1(i0,i1))
c             else
c               U_mag = sqrt(U0(i0,i1)*U0(i0,i1))
c             endif
c             mu_nc = 0.25d0*(mu(i0, i1) + mu(i0-1,i1) + mu(i0,i1-1)
c       &                   + mu(i0-1, i1-1))
c             rho_nc = 0.25d0*(rho(i0, i1) + rho(i0-1,i1) + rho(i0,i1-1)
c       &                   + rho(i0-1, i1-1))
c
c            U_tau_vis = sqrt(mu_nc*U_mag/(rho_nc*d(i0,i1)))
c            U_star_vis = sqrt(mu_nc*U_mag/(rho_nc*d(i0,i1)))
c             U_star_log = beta_star**0.25d0 * k(i0,i1)**0.5d0
c
c       This step needs to be checked surely because I use U_star
c        in yplus
c
c             y_plus = rho(i0,i1)*U_star_log*d(i0,i1)/mu(i0,i1)
c             U_tau_log = U_mag / (log(y_plus)/kappa + B)
c             U_tau = (U_tau_vis**4.d0 + U_tau_log**4.d0)**0.25d0
c             U_star = (U_star_vis**4.d0 + U_star_log**4.d0)**0.25d0
c             tau_w(i0,i1) = rho_nc*U_tau*U_star
c
c             P_k(i0, i1) = tau_w(i0,i1)*U_tau/(kappa*d(i0,i1))
c           endif
c         enddo
c        enddo
c      enddo
c
      return
      end
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute P_k = tau_w*U_tau/(kappa*U_mag). Not implemented.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
      subroutine near_wall_production_3d(
     & P_k, P_k_gc,
     & tau_w,tau_w_gc,
     & U_tau0,U_tau_gc0,
     & U_tau1,U_tau_gc1,
     & U_tau2,U_tau_gc2,
     & kappa,
     & wall_location_index,
     & ilower0,iupper0,
     & ilower1,iupper1,
     & ilower2,iupper2,
     & dx)

c
      implicit none
c
c     input variables
      INTEGER U_tau_gc0, U_tau_gc1, U_tau_gc2, tau_w_gc
      INTEGER ilower0, iupper0, ilower1, iupper1, ilower2, iupper2
      REAL kappa
      INTEGER wall_location_index
      REAL dx(0:NDIM-1)
c
      REAL tau_w(NODE3d(ilower, iupper, tau_w_gc))
      REAL U_tau0(SIDE3d0VECG(ilower,iupper,U_tau_gc))
      REAL U_tau1(SIDE3d1VECG(ilower,iupper,U_tau_gc))
      REAL U_tau2(SIDE3d2VECG(ilower,iupper,U_tau_gc))

c
c     output variables
      INTEGER P_k_gc
      REAL P_k(CELL3d(ilower,iupper,P_k_gc))


c     local variables
c      INTEGER i0, i1
c      REAL distance
c
c      do i1 = ilower1,iupper1
c        do i0 = ilower0,iupper0
c
c           if ((wall_location_index .eq. 0) .or.
c     &            (wall_location_index .eq. 1)) then
c            tau_w_cc = 0.5 * (tau_w(i0, i1) + tau_w(i0, i1+1))
c            U_tau_cc = 0.5 * (U_tau1(i0, i1) + U_tau1(i0, i1+1))
c            distance = dx(0) / 2.d0
c           else
c            tau_w_cc = 0.5 * (tau_w(i0, i1) + tau_w(i0+1, i1))
c            U_tau_cc = 0.5 * (U_tau0(i0, i1) + U_tau0(i0+1, i1))
c            distance = dx(1) / 2.d0
c           endif
c           P_k(i0, i1) = tau_w_cc * U_tau_cc / (kappa * distance)
c          print*, 'production from wall law'
c         print*, i0, i1, tau_w_cc, U_tau_cc, P_k(i0, i1)
c        enddo
c      enddo
      return
      end

c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute mu_w = tau_w*d/U_mag. Not implemented.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
      subroutine wall_viscosity_3d(
     & mu, mu_gc,
     & tau_w,tau_w_gc,
     & U0,U_gc0,
     & U1,U_gc1,
     & U2,U_gc2,
     & wall_location_index,
     & number_of_indices,
     & ilower0,iupper0,
     & ilower1,iupper1,
     & ilower2,iupper2,
     & dx)
c
      implicit none
c
c     input variables
      INTEGER U_gc0, U_gc1, U_gc2, tau_w_gc
      INTEGER ilower0, iupper0, ilower1, iupper1, ilower2, iupper2
      INTEGER wall_location_index

      REAL dx(0:NDIM-1)
      REAL number_of_indices(0:NDIM-1)
c
      REAL tau_w(EDGE3d0(ilower, iupper, tau_w_gc))
      REAL U0(SIDE3d0VECG(ilower, iupper, U_gc))
      REAL U1(SIDE3d1VECG(ilower, iupper, U_gc))
      REAL U2(SIDE3d2VECG(ilower, iupper, U_gc))
c
c     output variables
      INTEGER mu_gc
      REAL mu(EDGE3d0(ilower, iupper, mu_gc))

c

c
      return
      end

