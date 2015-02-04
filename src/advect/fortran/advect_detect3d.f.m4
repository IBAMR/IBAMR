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
dnl Process this file with m4 to produce FORTRAN source code
define(NDIM,3)dnl
define(REAL,`double precision')dnl
define(INTEGER,`integer')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim3d.i)dnl
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Detect sharp spatial gradients.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine advect_detectgrad3d(
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  vghost0,tagghost0,ttagghost0,
     &  vghost1,tagghost1,ttagghost1,
     &  vghost2,tagghost2,ttagghost2,
     &  dx,
     &  gradtol,
     &  dotag,
     &  var,
     &  tags,temptags)
c
      implicit none
c
      INTEGER
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  dotag,
     &  vghost0,vghost1,vghost2,
     &  tagghost0,tagghost1,tagghost2,
     &  ttagghost0,ttagghost1,ttagghost2
      REAL
     &  dx(0:NDIM-1),
     &  gradtol
      REAL
     &  var(CELL3dVECG(ifirst,ilast,vghost))
      INTEGER
     &  tags(CELL3dVECG(ifirst,ilast,tagghost)),
     &  temptags(CELL3dVECG(ifirst,ilast,ttagghost))
c
      REAL tol
      REAL facejump, loctol
      REAL presm1,presp1
      REAL diag(0:NDIM-1),diag012
      logical tagcell
      INTEGER ic0,ic1,ic2
c
      tol = gradtol
      diag(0) = sqrt(dx(2)**2+dx(1)**2)
      diag(1) = sqrt(dx(0)**2+dx(2)**2)
      diag(2) = sqrt(dx(0)**2+dx(1)**2)
      diag012 = sqrt(dx(0)**2+dx(1)**2+dx(2)**2)

      do ic2=ifirst2,ilast2
        do ic1=ifirst1,ilast1
          do ic0=ifirst0,ilast0

            if (tags(ic0,ic1,ic2) .ne. 0) then
              loctol = 0.125d0*tol
            else
              loctol = tol
            endif

            tagcell = .false.
c
c     One-dimensional diagonals.
c
            presm1 = var(ic0-1,ic1,ic2)
            presp1 = var(ic0+1,ic1,ic2)
            facejump = abs(var(ic0,ic1,ic2)-presm1)
            facejump = max(facejump,abs(var(ic0,ic1,ic2)-presp1))
            tagcell = ((facejump).gt.(loctol*dx(0)))
            if (.not.tagcell) then
              presm1 = var(ic0,ic1-1,ic2)
              presp1 = var(ic0,ic1+1,ic2)
              facejump = abs(var(ic0,ic1,ic2)-presm1)
              facejump = max(facejump,abs(var(ic0,ic1,ic2)-presp1))
              tagcell = ((facejump).gt.(loctol*dx(1)))
            endif
            if (.not.tagcell) then
              presm1 = var(ic0,ic1,ic2-1)
              presp1 = var(ic0,ic1,ic2+1)
              facejump = abs(var(ic0,ic1,ic2)-presm1)
              facejump = max(facejump,abs(var(ic0,ic1,ic2)-presp1))
              tagcell = ((facejump).gt.(loctol*dx(2)))
            endif
c
c     Two-dimensional diagonals.
c
            if (.not.tagcell) then
              presm1 = var(ic0,ic1-1,ic2-1)
              presp1 = var(ic0,ic1+1,ic2+1)
              facejump = abs(var(ic0,ic1,ic2)-presm1)
              facejump = max(facejump,abs(var(ic0,ic1,ic2)-presp1))
              tagcell = ((facejump).gt.(loctol*diag(0)))
            endif
            if (.not.tagcell) then
              presm1 = var(ic0,ic1+1,ic2-1)
              presp1 = var(ic0,ic1-1,ic2+1)
              facejump = abs(var(ic0,ic1,ic2)-presm1)
              facejump = max(facejump,abs(var(ic0,ic1,ic2)-presp1))
              tagcell = ((facejump).gt.(loctol*diag(0)))
            endif

            if (.not.tagcell) then
              presm1 = var(ic0-1,ic1,ic2-1)
              presp1 = var(ic0+1,ic1,ic2+1)
              facejump = abs(var(ic0,ic1,ic2)-presm1)
              facejump = max(facejump,abs(var(ic0,ic1,ic2)-presp1))
              tagcell = ((facejump).gt.(loctol*diag(1)))
            endif
            if (.not.tagcell) then
              presm1 = var(ic0-1,ic1,ic2+1)
              presp1 = var(ic0+1,ic1,ic2-1)
              facejump = abs(var(ic0,ic1,ic2)-presm1)
              facejump = max(facejump,abs(var(ic0,ic1,ic2)-presp1))
              tagcell = ((facejump).gt.(loctol*diag(1)))
            endif

            if (.not.tagcell) then
              presm1 = var(ic0-1,ic1-1,ic2)
              presp1 = var(ic0+1,ic1+1,ic2)
              facejump = abs(var(ic0,ic1,ic2)-presm1)
              facejump = max(facejump,abs(var(ic0,ic1,ic2)-presp1))
              tagcell = ((facejump).gt.(loctol*diag(2)))
            endif
            if (.not.tagcell) then
              presm1 = var(ic0-1,ic1+1,ic2)
              presp1 = var(ic0+1,ic1-1,ic2)
              facejump = abs(var(ic0,ic1,ic2)-presm1)
              facejump = max(facejump,abs(var(ic0,ic1,ic2)-presp1))
              tagcell = ((facejump).gt.(loctol*diag(2)))
            endif
c
c     Three-dimensional diagonals.
c
            if (.not.tagcell) then
              presm1 = var(ic0-1,ic1-1,ic2-1)
              presp1 = var(ic0+1,ic1+1,ic2+1)
              facejump = abs(var(ic0,ic1,ic2)-presm1)
              facejump = max(facejump,abs(var(ic0,ic1,ic2)-presp1))
              tagcell = ((facejump).gt.(loctol*diag012))
            endif
            if (.not.tagcell) then
              presm1 = var(ic0-1,ic1-1,ic2+1)
              presp1 = var(ic0+1,ic1+1,ic2-1)
              facejump = abs(var(ic0,ic1,ic2)-presm1)
              facejump = max(facejump,abs(var(ic0,ic1,ic2)-presp1))
              tagcell = ((facejump).gt.(loctol*diag012))
            endif
            if (.not.tagcell) then
              presm1 = var(ic0-1,ic1+1,ic2-1)
              presp1 = var(ic0+1,ic1-1,ic2+1)
              facejump = abs(var(ic0,ic1,ic2)-presm1)
              facejump = max(facejump,abs(var(ic0,ic1,ic2)-presp1))
              tagcell = ((facejump).gt.(loctol*diag012))
            endif
            if (.not.tagcell) then
              presm1 = var(ic0-1,ic1+1,ic2+1)
              presp1 = var(ic0+1,ic1-1,ic2-1)
              facejump = abs(var(ic0,ic1,ic2)-presm1)
              facejump = max(facejump,abs(var(ic0,ic1,ic2)-presp1))
              tagcell = ((facejump).gt.(loctol*diag012))
            endif

            if ( tagcell ) then
              temptags(ic0,ic1,ic2) = dotag
            endif
          enddo
        enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
