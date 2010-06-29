ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the principal square root XX of a 3x3 matrix AA.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine sqrtm(XX,AA)
c
      implicit none
c
      integer*8 N, LDA, LDVS, LWORK
      parameter (N = 3)
      parameter (LDA = N)
      parameter (LDVS = 2 * N)
      parameter (LWORK = 3 * N)
c
c     Input.
c
      real*8 AA(N,N)
c
c     Ouput.
c
      real*8 XX(N,N)
c
c     Local variables.
c
      complex*16 A(LDA,N), VS(LDVS,N)
      complex*16 WORK(LWORK), W(N)

      double precision RWORK(N)
      integer*8 SDIM,INFO

      logical*8 BWORK(1:N)
      logical DUMMY_SELECT
      external DUMMY_SELECT

      integer*8   i,j,k
      complex*16  Q(N,N), QR(N,N), QH(N,N), R(N,N), X(N,N)
      complex*16  summ
c
c     Copy AA into the complex matrix A.
c
      do j = 1,N
         do i = 1,N
            A(i,j) = dcmplx(AA(i,j))
         enddo
      enddo
c
c     Compute the eigenvalues and Schur vectors of A.
c
      call ZGEES('V','N',DUMMY_SELECT, N, A, LDA, SDIM, W, VS, LDVS,
     $     WORK, LWORK, RWORK, BWORK, INFO)
c
c     Compute the square root R of T one column at a time.
c
      do j = 1,N
         do i = 1,N
            R(i,j)=0.0
         enddo
      enddo
      summ=0.0
      do j = 1,N
         R(j,j) = sqrt(A(j,j))
         do i = j-1,1,-1
            do k = i+1,j-1
               summ = summ + R(i,k)*R(k,j)
            enddo
            R(i,j) = (A(i,j) - summ)/(R(i,i) + R(j,j))
         enddo
      enddo
c
c     Copy VS into Q.
c
      do j = 1,N
         do i = 1,N
            Q(i,j) = VS(i,j)
         enddo
      enddo
c
c     Set X to be the square root of A : X = real(Q*R*Q**H).
c
      do j = 1,N
         do i = 1,N
            QH(i,j) = CONJG(Q(j,i))
         enddo
      enddo
      QR = MATMUL(Q,R)
      X = MATMUL(QR,QH)
c
c     Set XX to be the real part of X.
c
      do j = 1,N
         do i = 1,N
            XX(i,j) = dble(X(i,j))
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
      complex*16 ARG
C
C     The value computed is always .TRUE.  It is computed in the
C     peculiar way below to avoid compiler messages about unused
C     function arguments.
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
