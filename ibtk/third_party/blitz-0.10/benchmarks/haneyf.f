CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 
C  AUTHOR
C     Scott Haney, LLNL, swhaney@llnl.gov, (510) 423-6308
C 
C  COPYRIGHT
C     Copyright (c) 1995.  The Regents of the University of California.
C                          All rights reserved.
C 
C     Permission to use, copy, modify, and distribute this software for any
C     purpose without fee is hereby granted, provided that this entire
C     notice is included in all copies of any software which is or includes
C     a copy or modification of this software and in all copies of the
C     a copy or modification of this software and in all copies of the
C     supporting documentation for such software.
C 
C     This work was produced at the University of California, Lawrence
C     Livermore National Laboratory under contract no. W-7405-ENG-48 between
C     the U.S. Department of Energy and The Regents of the University of
C     California for the operation of UC LLNL.
C 
C  DISCLAIMER
C     This software was prepared as an account of work sponsored by an
C     agency of the United States Government.  Neither the United States
C     Government nor the University of California nor any of their
C     employees, makes any warranty, express or implied, or assumes any
C     liability or responsibility for the accuracy, completeness, or
C     usefulness of any information, apparatus, product, or process
C     disclosed, or represents that its use would not infringe
C     privately-owned rights.  Reference herein to any specific commercial
C     products, process, or service by trade name, trademark, manufacturer,
C     privately-owned rights.  Reference herein to any specific commercial
C     products, process, or service by trade name, trademark, manufacturer,
C     or otherwise, does not necessarily constitute or imply its
C     endorsement, recommendation, or favoring by the United States
C     Government or the University of California.  The views and opinions of
C     authors expressed herein do not necessarily state or reflect those of
C     the United States Government or the University of California, and
C     shall not be used for advertising or product endorsement purposes.
C 
C 
C  CREATED
C     July 22, 1994
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine vecopsf(li, R, a, n, iters)

      integer i
      real li(n),R(n),a(n)
      real xmu0, xm24
      parameter(xmu0 = 4.0 * 3.1415926535 * 1.0e-7)
      parameter(xm24 = 1.0 / 24.0)

      do 2000 j = 1, iters
      do 1000 i = 1, n
         li(i) = xmu0 * R(i) * 
     .        (0.5 * (1.0 + xm24 * (a(i) / R(i))**2) *
     .        log(32.0 * (R(i) / a(i))**2) + 
     .        0.84 * (a(i) / R(i))**2 - 0.2)

 1000 continue
 2000 continue

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Overhead timing
C     Added January 1997, Todd Veldhuizen
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine vecopsfo(li, R, a, n, iters)

      integer i
      real li(n),R(n),a(n)
      real xmu0, xm24
      parameter(xmu0 = 4.0 * 3.1415926535 * 1.0e-7)
      parameter(xm24 = 1.0 / 24.0)

      do 2000 j = 1, iters
 2000 continue

      return
      end
