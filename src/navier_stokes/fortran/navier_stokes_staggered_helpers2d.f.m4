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
define(NDIM,2)dnl
define(REAL,`double precision')dnl
define(INTEGER,`integer')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim2d.i)dnl
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Interpolate the components of a staggered velocity field onto the
c     faces of zones centered about the components of the velocity
c     field.
c
c     NOTES:
c
c     U0 and U1 are standard side-centered staggered grid velocities for
c     the patch [(ifirst0,ilast0),(ifirst1,ilast1)].
c
c     V00 and V01 are face-centered staggered grid velocities defined at
c     the faces of the control volumes centered about the x components
c     of the side-centered velocity, i.e., face-centered staggered grid
c     velocities for the patch [(ifirst0,ilast0+1),(ifirst1,ilast1)].
c
c     V10 and V11 are face-centered staggered grid velocities defined at
c     the faces of the control volumes centered about the y components
c     of the side-centered velocity, i.e., face-centered staggered grid
c     velocities for the patch [(ifirst0,ilast0),(ifirst1,ilast1+1)].
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine navier_stokes_interp_comps2d(
     &     patch_ifirst0,patch_ilast0,
     &     patch_ifirst1,patch_ilast1,
     &     n_U_gc0,n_U_gc1,
     &     U0,U1,
     &     side0_ifirst0,side0_ilast0,
     &     side0_ifirst1,side0_ilast1,
     &     n_V0_gc0,n_V0_gc1,
     &     V00,V01,
     &     side1_ifirst0,side1_ilast0,
     &     side1_ifirst1,side1_ilast1,
     &     n_V1_gc0,n_V1_gc1,
     &     V10,V11)
c
      implicit none
c
c     Input.
c
      INTEGER patch_ifirst0,patch_ilast0
      INTEGER patch_ifirst1,patch_ilast1

      INTEGER n_U_gc0,n_U_gc1

      INTEGER side0_ifirst0,side0_ilast0
      INTEGER side0_ifirst1,side0_ilast1

      INTEGER n_V0_gc0,n_V0_gc1

      INTEGER side1_ifirst0,side1_ilast0
      INTEGER side1_ifirst1,side1_ilast1

      INTEGER n_V1_gc0,n_V1_gc1

      REAL U0(
     &     SIDE2d0VECG(patch_ifirst,patch_ilast,n_U_gc)
     &     )
      REAL U1(
     &     SIDE2d1VECG(patch_ifirst,patch_ilast,n_U_gc)
     &     )
c
c     Input/Output.
c
      REAL V00(
     &     FACE2d0VECG(side0_ifirst,side0_ilast,n_V0_gc)
     &     )
      REAL V01(
     &     FACE2d1VECG(side0_ifirst,side0_ilast,n_V0_gc)
     &     )
      REAL V10(
     &     FACE2d0VECG(side1_ifirst,side1_ilast,n_V1_gc)
     &     )
      REAL V11(
     &     FACE2d1VECG(side1_ifirst,side1_ilast,n_V1_gc)
     &     )
c
c     Local variables.
c
      INTEGER i0,i1
      INTEGER gc0,gc1
c
c     Interpolate the components of the velocity at each zone face.
c
      gc0 = min(n_U_gc0-1,n_V0_gc0)
      gc1 = min(n_U_gc1  ,n_V0_gc1)

      do    i1 = side0_ifirst1-gc1,side0_ilast1+gc1
         do i0 = side0_ifirst0-gc0,side0_ilast0+gc0
            V00(i0,i1) = 0.5d0*(U0(i0-1,i1)+U0(i0,i1))
            V01(i1,i0) = 0.5d0*(U1(i0-1,i1)+U1(i0,i1))
         enddo
      enddo

      gc0 = min(n_U_gc0  ,n_V1_gc0)
      gc1 = min(n_U_gc1-1,n_V1_gc1)

      do    i0 = side1_ifirst0-gc0,side1_ilast0+gc0
         do i1 = side1_ifirst1-gc1,side1_ilast1+gc1
            V10(i0,i1) = 0.5d0*(U0(i0,i1-1)+U0(i0,i1))
            V11(i1,i0) = 0.5d0*(U1(i0,i1-1)+U1(i0,i1))
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Reset the face-centered advection velocity about the control
c     volumes for each component of the velocity.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine navier_stokes_reset_adv_velocity2d(
     &     side0_ifirst0,side0_ilast0,
     &     side0_ifirst1,side0_ilast1,
     &     n_U_adv0_gc0,n_U_adv0_gc1,
     &     U_adv00,U_adv01,
     &     n_U_half0_gc0,n_U_half0_gc1,
     &     U_half00,U_half01,
     &     side1_ifirst0,side1_ilast0,
     &     side1_ifirst1,side1_ilast1,
     &     n_U_adv1_gc0,n_U_adv1_gc1,
     &     U_adv10,U_adv11,
     &     n_U_half1_gc0,n_U_half1_gc1,
     &     U_half10,U_half11)
c
      implicit none
c
c     Input.
c
      INTEGER side0_ifirst0,side0_ilast0
      INTEGER side0_ifirst1,side0_ilast1

      INTEGER n_U_adv0_gc0,n_U_adv0_gc1
      INTEGER n_U_half0_gc0,n_U_half0_gc1

      INTEGER side1_ifirst0,side1_ilast0
      INTEGER side1_ifirst1,side1_ilast1

      INTEGER n_U_adv1_gc0,n_U_adv1_gc1
      INTEGER n_U_half1_gc0,n_U_half1_gc1

      REAL U_half00(
     &     FACE2d0VECG(side0_ifirst,side0_ilast,n_U_half0_gc)
     &     )
      REAL U_half01(
     &     FACE2d1VECG(side0_ifirst,side0_ilast,n_U_half0_gc)
     &     )
      REAL U_half10(
     &     FACE2d0VECG(side1_ifirst,side1_ilast,n_U_half1_gc)
     &     )
      REAL U_half11(
     &     FACE2d1VECG(side1_ifirst,side1_ilast,n_U_half1_gc)
     &     )
c
c     Input/Output.
c
      REAL U_adv00(
     &     FACE2d0VECG(side0_ifirst,side0_ilast,n_U_adv0_gc)
     &     )
      REAL U_adv01(
     &     FACE2d1VECG(side0_ifirst,side0_ilast,n_U_adv0_gc)
     &     )
      REAL U_adv10(
     &     FACE2d0VECG(side1_ifirst,side1_ilast,n_U_adv1_gc)
     &     )
      REAL U_adv11(
     &     FACE2d1VECG(side1_ifirst,side1_ilast,n_U_adv1_gc)
     &     )
c
c     Local variables.
c
      INTEGER i0,i1
      INTEGER gc0,gc1
c
c     Reset the advection velocity.
c
      gc0 = min(n_U_adv0_gc0,n_U_adv1_gc0,n_U_half0_gc0,n_U_half1_gc0)
      gc1 = min(n_U_adv0_gc1,n_U_adv1_gc1,n_U_half0_gc1,n_U_half1_gc1)

      do    i1 = side0_ifirst1-gc1,side0_ilast1+gc1
         do i0 = side0_ifirst0-gc0,side0_ilast0+gc0
            U_adv00(i0,i1) = U_half00(i0,i1)
            U_adv01(i1,i0) = U_half10(i0,i1)
         enddo
      enddo

      do    i1 = side1_ifirst1-gc1,side1_ilast1+gc1
         do i0 = side1_ifirst0-gc0,side1_ilast0+gc0
            U_adv10(i0,i1) = U_half01(i1,i0)
            U_adv11(i1,i0) = U_half11(i1,i0)
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Convert a side-centered vector field into a face-centered vector
c     field.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine navier_stokes_side_to_face2d(
     &     ifirst0,ilast0,
     &     ifirst1,ilast1,
     &     u_sc0,u_sc1,u_sc_gcw,
     &     u_fc0,u_fc1,u_fc_gcw)
c
      implicit none
c
c     Input.
c
      INTEGER ifirst0,ilast0
      INTEGER ifirst1,ilast1

      INTEGER u_sc_gcw,u_fc_gcw

      REAL u_sc0(SIDE2d0(ifirst,ilast,u_sc_gcw))
      REAL u_sc1(SIDE2d1(ifirst,ilast,u_sc_gcw))
c
c     Input/Output.
c
      REAL u_fc0(FACE2d0(ifirst,ilast,u_fc_gcw))
      REAL u_fc1(FACE2d1(ifirst,ilast,u_fc_gcw))
c
c     Local variables.
c
      INTEGER i0,i1
      INTEGER gcw
c
c     Swap the velocity values.
c
      gcw = min(u_sc_gcw,u_fc_gcw)

      do i1 = ifirst1-gcw,ilast1+gcw
         do i0 = ifirst0-gcw,ilast0+gcw+1
            u_fc0(i0,i1) = u_sc0(i0,i1)
         enddo
      enddo

      do i1 = ifirst1-gcw,ilast1+gcw+1
         do i0 = ifirst0-gcw,ilast0+gcw
            u_fc1(i1,i0) = u_sc1(i0,i1)
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Convert a face-centered vector field into a side-centered vector
c     field.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine navier_stokes_face_to_side2d(
     &     ifirst0,ilast0,
     &     ifirst1,ilast1,
     &     u_fc0,u_fc1,u_fc_gcw,
     &     u_sc0,u_sc1,u_sc_gcw)
c
      implicit none
c
c     Input.
c
      INTEGER ifirst0,ilast0
      INTEGER ifirst1,ilast1

      INTEGER u_fc_gcw,u_sc_gcw

      REAL u_fc0(FACE2d0(ifirst,ilast,u_fc_gcw))
      REAL u_fc1(FACE2d1(ifirst,ilast,u_fc_gcw))
c
c     Input/Output.
c
      REAL u_sc0(SIDE2d0(ifirst,ilast,u_sc_gcw))
      REAL u_sc1(SIDE2d1(ifirst,ilast,u_sc_gcw))
c
c     Local variables.
c
      INTEGER i0,i1
      INTEGER gcw
c
c     Swap the velocity values.
c
      gcw = min(u_sc_gcw,u_fc_gcw)

      do i1 = ifirst1-gcw,ilast1+gcw
         do i0 = ifirst0-gcw,ilast0+gcw+1
            u_sc0(i0,i1) = u_fc0(i0,i1)
         enddo
      enddo

      do i1 = ifirst1-gcw,ilast1+gcw+1
         do i0 = ifirst0-gcw,ilast0+gcw
            u_sc1(i0,i1) = u_fc1(i1,i0)
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Copy SAMRAI velocity and pressure patch data to PETSc Vec  
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine copy_to_patchlevel_vec_mac2d(
     &     ifirst0,ilast0,
     &     ifirst1,ilast1,
     &     first_local,last_local,
     &     p_cc,p_cc_gcw,
     &     u_sc0,u_sc1,u_sc_gcw,
     &     p_dof_cc,p_dof_cc_gcw,
     &     u_dof_sc0, u_dof_sc1,u_dof_sc_gcw,
     &     arr)
c
      implicit none
c
c     Input.
c
      INTEGER ifirst0,ilast0
      INTEGER ifirst1,ilast1

      INTEGER first_local,last_local
      INTEGER p_cc_gcw,u_sc_gcw
      INTEGER p_dof_cc_gcw,u_dof_sc_gcw
c
c     Input/Output.
c
      REAL p_cc(CELL2d(ifirst,ilast,p_cc_gcw))
      REAL u_sc0(SIDE2d0(ifirst,ilast,u_sc_gcw))
      REAL u_sc1(SIDE2d1(ifirst,ilast,u_sc_gcw))

      INTEGER p_dof_cc(CELL2d(ifirst,ilast,p_dof_cc_gcw))
      INTEGER u_dof_sc0(SIDE2d0(ifirst,ilast,u_dof_sc_gcw))
      INTEGER u_dof_sc1(SIDE2d1(ifirst,ilast,u_dof_sc_gcw))

      REAL arr(first_local:last_local-1)
c
c     Local variables.
c
      INTEGER i0,i1
      INTEGER dof_idx
c
c     Copy velocity to array.
c
    
      do i1 = ifirst1,ilast1
         do i0 = ifirst0,ilast0+1
	    dof_idx = u_dof_sc0(i0,i1)
	    if ((dof_idx .ge. first_local) .and.
     &   	(dof_idx .lt. last_local)) then
	       arr(dof_idx) = u_sc0(i0,i1)
	    end if	  
         enddo
      enddo

      do i1 = ifirst1,ilast1+1
         do i0 = ifirst0,ilast0
	    dof_idx = u_dof_sc1(i0,i1)
	    if ((dof_idx .ge. first_local) .and.
     &      (dof_idx .lt. last_local)) then
     	       arr(dof_idx) = u_sc1(i0,i1)
	    end if	    
         enddo
      enddo
c

c
c     Copy pressure to array
c
      
      do i1 = ifirst1,ilast1
         do i0 = ifirst0,ilast0
	    dof_idx = p_dof_cc(i0,i1)
	    if ((dof_idx .ge. first_local) .and.
     &      (dof_idx .lt. last_local)) then
	       arr(dof_idx) = p_cc(i0,i1)
	    end if	  
         enddo
      enddo     
 
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Copy PETSc Vec velocity and pressure to SAMRAI patch data  
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine copy_from_patchlevel_vec_mac2d(
     &     ifirst0,ilast0,
     &     ifirst1,ilast1,
     &     first_local,last_local,
     &     p_cc,p_cc_gcw,
     &     u_sc0,u_sc1,u_sc_gcw,
     &     p_dof_cc,p_dof_cc_gcw,
     &     u_dof_sc0, u_dof_sc1,u_dof_sc_gcw,
     &     arr)
c
      implicit none
c
c     Input.
c
      INTEGER ifirst0,ilast0
      INTEGER ifirst1,ilast1

      INTEGER first_local,last_local
      INTEGER p_cc_gcw,u_sc_gcw
      INTEGER p_dof_cc_gcw,u_dof_sc_gcw
c
c     Input/Output.
c
      REAL p_cc(CELL2d(ifirst,ilast,p_cc_gcw))
      REAL u_sc0(SIDE2d0(ifirst,ilast,u_sc_gcw))
      REAL u_sc1(SIDE2d1(ifirst,ilast,u_sc_gcw))

      INTEGER p_dof_cc(CELL2d(ifirst,ilast,p_dof_cc_gcw))
      INTEGER u_dof_sc0(SIDE2d0(ifirst,ilast,u_dof_sc_gcw))
      INTEGER u_dof_sc1(SIDE2d1(ifirst,ilast,u_dof_sc_gcw))

      REAL arr(first_local:last_local-1)
c
c     Local variables.
c
      INTEGER i0,i1
      INTEGER dof_idx
c
c     Copy array velocity  to SAMRAI patch data.
c
    
      do i1 = ifirst1,ilast1
         do i0 = ifirst0,ilast0+1
	    dof_idx = u_dof_sc0(i0,i1)
	    if ((dof_idx .ge. first_local) .and.
     &      (dof_idx .lt. last_local)) then
	       u_sc0(i0,i1) = arr(dof_idx)
	    end if	  
         enddo
      enddo

      do i1 = ifirst1,ilast1+1
         do i0 = ifirst0,ilast0
	    dof_idx = u_dof_sc1(i0,i1)
	    if ((dof_idx .ge. first_local) .and.
     &      (dof_idx .lt. last_local)) then
	       u_sc1(i0,i1) =  arr(dof_idx)
	    end if	    
         enddo
      enddo
c

c
c     Copy array pressure to SAMRAI patch data.
c
      
      do i1 = ifirst1,ilast1
         do i0 = ifirst0,ilast0
	    dof_idx = p_dof_cc(i0,i1)
	    if ((dof_idx .ge. first_local) .and.
     &      (dof_idx .lt. last_local)) then
	       p_cc(i0,i1) = arr(dof_idx)
	    end if	  
         enddo
      enddo     
 
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

