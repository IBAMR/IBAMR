
      SUBROUTINE loop24_F77(N, x, a, b, c, d, y)
      INTEGER i, N
      REAL*8 x(N), a(N), b(N), c(N), d(N), y(N)

      DO i=1,N
          x(i) = a(i)*c(i) - b(i)*c(i); y(i) = a(i)*d(i) + b(i)+c(i);
      END DO
      RETURN
      END


      SUBROUTINE loop24_F77Overhead(N, x, a, b, c, d, y)
      INTEGER i, N
      REAL*8 x(N), a(N), b(N), c(N), d(N), y(N)
      RETURN
      END
