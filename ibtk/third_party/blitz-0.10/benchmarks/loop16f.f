
      SUBROUTINE loop16_F77(N, y, x, a, b, c, u)
      INTEGER i, N
      REAL*8 y(N), x(N), a(N), b(N), c(N), u

      DO i=1,N
          x(i) = a(i)+b(i)+c(i); y(i) = x(i)+c(i)+u;
      END DO
      RETURN
      END


      SUBROUTINE loop16_F77Overhead(N, y, x, a, b, c, u)
      INTEGER i, N
      REAL*8 y(N), x(N), a(N), b(N), c(N), u
      RETURN
      END
