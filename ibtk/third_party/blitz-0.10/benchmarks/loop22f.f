
      SUBROUTINE loop22_F77(N, x, a, b, c, d)
      INTEGER i, N
      REAL*8 x(N), a(N), b(N), c(N), d(N)

      DO i=1,N
          x(i) = x(i) + a(i)*b(i) + c(i)*d(i);
      END DO
      RETURN
      END


      SUBROUTINE loop22_F77Overhead(N, x, a, b, c, d)
      INTEGER i, N
      REAL*8 x(N), a(N), b(N), c(N), d(N)
      RETURN
      END
