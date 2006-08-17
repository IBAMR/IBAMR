dnl Process this file with m4 to produce FORTRAN source code
define(REAL,`double precision')dnl
define(INTEGER,`integer')dnl
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     The MUSCL limited 4th order centered difference approximation to
c     dQ/dx.
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
      REAL function muscldiff(Q)
c     
      implicit none
include(TOP_SRCDIR/src/fortran/const.i)dnl     
c     
c     Input.
c     
      REAL Q(-2:2)
c     
c     Local variables.
c     
      REAL dQ_lim,dQf_left,dQf_rght
c     
c     The MUSCL limited 4th order centered difference approximation to
c     dQ/dx.
c     
      if ( ((Q(1)-Q(0))*(Q(0)-Q(-1))).gt.0.d0 ) then
         if ( ((Q(0)-Q(-1))*(Q(-1)-Q(-2))).gt.0.d0 ) then
            dQ_lim = 2.d0*dmin1(dabs(Q(0)-Q(-1)),dabs(Q(-1)-Q(-2)))
            dQf_left = dmin1(0.5d0*dabs(Q(0)-Q(-2)),dQ_lim)
            
            if ( (Q(0)-Q(-2)).lt.0.d0 ) dQf_left = -dQf_left
         else
            dQf_left = 0.d0
         endif
         
         if ( ((Q(2)-Q(1))*(Q(1)-Q(0))).gt.0.d0 ) then
            dQ_lim = 2.d0*dmin1(dabs(Q(2)-Q(1)),dabs(Q(1)-Q(0)))
            dQf_rght = dmin1(0.5d0*dabs(Q(2)-Q(0)),dQ_lim)
            
            if ( (Q(2)-Q(0)).lt.0.d0 ) dQf_rght = -dQf_rght
         else
            dQf_rght = 0.d0
         endif
         
         dQ_lim = 2.d0*dmin1(dabs(Q(1)-Q(0)),dabs(Q(0)-Q(-1)))
         
         muscldiff = dmin1(twothird*dabs( Q( 1)-0.25d0*dQf_rght
     &                                   -Q(-1)-0.25d0*dQf_left),dQ_lim)
         
         if ( (Q(1)-Q(-1)).lt.0.d0 ) then
            muscldiff = -muscldiff
         endif
      else
         muscldiff = 0.d0
      endif
c     
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     The minmod function of three arguments.
c
c     The minmod function is a function of two or more arguments that
c     takes the value of the argument with the smallest modulus if all
c     arguments have the same sign.  Otherwise it takes the value zero.
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
      REAL function minmod3(a,b,c)
c     
      implicit none
include(TOP_SRCDIR/src/fortran/const.i)dnl     
c     
c     Input.
c     
      REAL a,b,c
c     
c     minmod(a,b,c)
c     
      if     ( (a.ge.0.0d0).and.(b.ge.0.0d0).and.(c.ge.0.d0) ) then
         minmod3 = dmin1(a,b,c)
      elseif ( (a.le.0.0d0).and.(b.le.0.0d0).and.(c.le.0.d0) ) then
         minmod3 = dmax1(a,b,c)
      else
         minmod3 = 0.d0
      endif
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
