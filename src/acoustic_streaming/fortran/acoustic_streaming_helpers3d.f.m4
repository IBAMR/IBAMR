c ---------------------------------------------------------------------
c
c Copyright (c) 2008 - 2023 by the IBAMR developers
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
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Copy staggered velocity field onto another staggered velocity field.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine copy_staggered_data_3d(U0, U1, U2,
     &     U_gcw, V0, V1, V2, V_gcw, ilower0, iupper0,
     &     ilower1, iupper1, ilower2, iupper2)
c
      implicit none
c
c     Input.
c
      INTEGER U_gcw, V_gcw
      INTEGER ilower0, iupper0
      INTEGER ilower1, iupper1
      INTEGER ilower2, iupper2

      REAL U0(SIDE3d0(ilower,iupper,U_gcw))
      REAL U1(SIDE3d1(ilower,iupper,U_gcw))
      REAL U2(SIDE3d2(ilower,iupper,U_gcw))
c
c     Input/Output.
c
      REAL V0(SIDE3d0(ilower,iupper,V_gcw))
      REAL V1(SIDE3d1(ilower,iupper,V_gcw))
      REAL V2(SIDE3d2(ilower,iupper,V_gcw))
c
c     Local variables.
c
      INTEGER i0,i1

c
c     Copy the components of U to V.
c

c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Copy cell-centered field onto another cell-centered field.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine copy_cell_data_3d(U,
     &     U_gcw, V, V_gcw, ilower0, iupper0,
     &     ilower1, iupper1, ilower2, iupper2)
c
      implicit none
c
c     Input.
c
      INTEGER U_gcw, V_gcw
      INTEGER ilower0, iupper0
      INTEGER ilower1, iupper1
      INTEGER ilower2, iupper2

      REAL U(CELL3d(ilower,iupper,U_gcw))

c
c     Input/Output.
c
      REAL V(CELL3d(ilower,iupper,V_gcw))
c
c     Local variables.
c
      INTEGER i0,i1

c
c     Copy the components of U to V.
c


c
      return
      end
c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the forcing term for the RHS of second order momentum equation due   
c     to the first order system f1 = - < div . (rho (U1 x U1))>
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      


      subroutine acoustic_momentum_coupling_3d(U0r, U1r, U2r, 
     &     Ur_gcw, U0i, U1i, U2i, Ui_gcw, rho0, rho1, rho2,
     &     rho_gcw, f0, f1, f2, f_gcw, ilower0, iupper0, 
     &     ilower1, iupper1, ilower2, iupper2, dx)

c
      implicit none
c
c     Input.
c
      INTEGER Ur_gcw, Ui_gcw, rho_gcw, f_gcw
      INTEGER ilower0, iupper0
      INTEGER ilower1, iupper1
      INTEGER ilower2, iupper2

      REAL U0r(SIDE3d0(ilower,iupper,Ur_gcw))
      REAL U1r(SIDE3d1(ilower,iupper,Ur_gcw))
      REAL U2r(SIDE3d2(ilower,iupper,Ur_gcw))
      REAL U0i(SIDE3d0(ilower,iupper,Ui_gcw))
      REAL U1i(SIDE3d1(ilower,iupper,Ui_gcw))
      REAL U2i(SIDE3d2(ilower,iupper,Ui_gcw))
      REAL rho0(SIDE3d0(ilower,iupper,rho_gcw))
      REAL rho1(SIDE3d1(ilower,iupper,rho_gcw))
      REAL rho2(SIDE3d2(ilower,iupper,rho_gcw))

      REAL dx(0:NDIM-1)

c
c     Input/Output.
c
      REAL f0(SIDE3d0(ilower,iupper,f_gcw))
      REAL f1(SIDE3d1(ilower,iupper,f_gcw))
      REAL f2(SIDE3d2(ilower,iupper,f_gcw))

c     Local variables
      
      REAL fac0, fac1, g1, g2
      INTEGER i0, i1

      fac0 = 1.0d0/dx(0)
      fac1 = 1.0d0/dx(1)
      
c
c     Compute the coupling source term.
c

           
c     
      return
      end
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the forcing term for the RHS of second order continuity equation due   
c     to the first order system f2 =  < div . (rho1 U1)>
c     Note that we are removing the -ve sign from f2 because the second order system 
c     considers -div (rho0 U2) as the continuity equation
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    


      subroutine acoustic_mass_coupling_3d(U0r, U1r, U2r,
     &     Ur_gcw, U0i, U1i, U2i, Ui_gcw, pr,
     &     pr_gcw, pi,  pi_gcw, c0,
     &     m, m_gcw, ilower0, iupper0, ilower1, iupper1, 
     &     ilower2, iupper2, dx)

c
      implicit none
c
c     Input.
c
      INTEGER Ur_gcw, Ui_gcw, pr_gcw, pi_gcw, m_gcw
      INTEGER ilower0, iupper0
      INTEGER ilower1, iupper1
      INTEGER ilower2, iupper2


      REAL U0r(SIDE3d0(ilower,iupper,Ur_gcw))
      REAL U1r(SIDE3d1(ilower,iupper,Ur_gcw))
      REAL U2r(SIDE3d2(ilower,iupper,Ur_gcw))
      REAL U0i(SIDE3d0(ilower,iupper,Ui_gcw))
      REAL U1i(SIDE3d1(ilower,iupper,Ui_gcw))
      REAL U2i(SIDE3d2(ilower,iupper,Ui_gcw))
      REAL pr(CELL3d(ilower,iupper,pr_gcw))
      REAL pi(CELL3d(ilower,iupper,pi_gcw))
      REAL c0
      REAL dx(0:NDIM-1)

c
c     Input/Output.
c
      REAL m(CELL3d(ilower,iupper,m_gcw))

c     Local variables

      REAL fac0, fac1, fac2, g1, g2
      INTEGER i0,i1

      fac0 = 1.d0/dx(0)
      fac1 = 1.d0/dx(1) 
      fac2 = 0.5d0/c0**2

c
c     Compute the coupling source term.
c

c     
      return
      end

