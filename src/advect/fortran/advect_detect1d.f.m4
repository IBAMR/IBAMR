c
c     Copyright (c) 2002-2010, Boyce Griffith
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
c        * Neither the name of New York University nor the names of its
c          contributors may be used to endorse or promote products
c          derived from this software without specific prior written
c          permission.
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
define(NDIM,1)dnl
define(REAL,`double precision')dnl
define(INTEGER,`integer')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim1d.i)dnl
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Detect sharp spatial gradients.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine advect_detectgrad1d(
     &  ifirst0,ilast0,
     &  vghost0,tagghost0,ttagghost0,
     &  dx,
     &  gradtol,
     &  dotag,donttag,
     &  var,
     &  tags,temptags)
c
      implicit none
include(TOP_SRCDIR/src/fortran/const.i)dnl
c
      INTEGER
     &  ifirst0,ilast0,
     &  dotag,donttag,
     &  vghost0,
     &  tagghost0,
     &  ttagghost0
      REAL
     &  dx(0:NDIM-1),
     &  gradtol
      REAL
     &  var(CELL1dVECG(ifirst,ilast,vghost))
      INTEGER
     &  tags(CELL1dVECG(ifirst,ilast,tagghost)),
     &  temptags(CELL1dVECG(ifirst,ilast,ttagghost))
c
      REAL tol
      REAL facejump, loctol
      REAL presm1,presp1
      logical tagcell
      INTEGER ic0
c
      tol = gradtol

      do ic0=ifirst0,ilast0

         if (tags(ic0) .ne. 0) then
            loctol = 0.125*tol
         else
            loctol = tol
         endif

         tagcell = .false.

         presm1 = var(ic0-1)
         presp1 = var(ic0+1)
         facejump = abs(var(ic0)-presm1)
         facejump = max(facejump,abs(var(ic0)-presp1))
         tagcell = ((facejump).gt.(loctol*dx(0)))

         if ( tagcell ) then
            temptags(ic0) = dotag
         endif
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
