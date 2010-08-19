c
c     Copyright (c) 2002-2010 Boyce Griffith
c
c     Permission is hereby granted, free of charge, to any person
c     obtaining a copy of this software and associated documentation
c     files (the "Software"), to deal in the Software without
c     restriction, including without limitation the rights to use, copy,
c     modify, merge, publish, distribute, sublicense, and/or sell copies
c     of the Software, and to permit persons to whom the Software is
c     furnished to do so, subject to the following conditions:
c
c     The above copyright notice and this permission notice shall be
c     included in all copies or substantial portions of the Software.
c
c     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
c     EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
c     MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
c     NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
c     HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
c     WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
c     OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
c     DEALINGS IN THE SOFTWARE.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the principal square root X of an N-by-N real matrix A,
c     specialized to the case N=3.
c
c     The algorithm implemented herein is described in: A. Bjorck and
c     S. Hammarling, A Schur method for the square root of a matrix, Lin
c     Algebra Appl, 52-53:127--140 (1983).
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine dsqrtm(X,A)
c
      implicit none
c
c     External functions.
c
      EXTERNAL DUMMY_SELECT
c
c     Parameters.
c
      INTEGER N,LWORK
      PARAMETER(N=3,LWORK=16*N)
      COMPLEX*16 ONE,ZERO
      PARAMETER(ONE=(1.d0,0.d0),ZERO=(0.d0,0.d0))
c
c     Input.
c
      REAL*8 A(N,N)
c
c     Ouput.
c
      REAL*8 X(N,N)
c
c     Local variables.
c
      INTEGER i,j,k,sdiag
      INTEGER SDIM,INFO
      COMPLEX*16 Q(N,N),S(N,N),U(N,N),W(N),WORK(LWORK)
      COMPLEX*16 QUQH(N,N),QU(N,N)
      REAL*8 RWORK(N)
      LOGICAL BWORK(N)
c
c     Compute the Schur form Q*S*(Q**H)
c
      do j = 1,N
         do i = 1,N
            S(i,j) = dcmplx(A(i,j))
         enddo
      enddo
      call ZGEES('V', 'N', DUMMY_SELECT, N, S, N, SDIM, W, Q, N, WORK,
     &     LWORK, RWORK, BWORK, INFO)
c
c     Compute the square root U of S one super-diagonal at a time.
c
      do j = 1,N-1              ! set the lower triangle to zero
         do i = j+1,N
            U(i,j) = ZERO
         enddo
      enddo

      do i = 1,N                ! set the diagonal elements
         U(i,i) = zsqrt(S(i,i))
      enddo

      do sdiag = 1,N-1          ! loop over the N-1 super-diagonals
         do i = 1,N-sdiag
            j = i+sdiag
            U(i,j) = S(i,j)
            do k = i+1,j-1
               U(i,j) = U(i,j) - U(i,k)*U(k,j)
            enddo
            U(i,j) = U(i,j) / (U(i,i) + U(j,j))
         enddo
      enddo
c
c     Compute X = Q*U*(Q**H).
c
      call ZGEMM('N', 'N', N, N, N, ONE, Q, N, U, N, ZERO, QU, N)
      call ZGEMM('N', 'C', N, N, N, ONE, QU, N, Q, N, ZERO, QUQH, N)
      do j = 1,N
         do i = 1,N
            X(i,j) = dble(QUQH(i,j))
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Dummy select function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      LOGICAL FUNCTION DUMMY_SELECT(ARG)
c
      IMPLICIT NONE
c
c     Input.
c
      COMPLEX*16 ARG
C
C     The return value is always .TRUE. unless ARG is nan.  It is
C     computed in the peculiar way below to avoid compiler messages
C     about unused function arguments.
c
      PRINT 1000
      DUMMY_SELECT = (ARG .EQ. ARG)
 1000 FORMAT (///1X, '**** ERROR:  ',
     $     'DUMMY_SELECT FUNCTION CALLED BUT NOT AVAILABLE. ****')
c
      stop
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
