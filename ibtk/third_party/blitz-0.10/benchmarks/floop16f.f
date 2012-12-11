
      SUBROUTINE floop16_F77(N, y, x, a, b, c, u)
      INTEGER i, N
      REAL*4 y(N), x(N), a(N), b(N), c(N), u

      DO i=1,N
          x(i) = a(i)+b(i)+c(i); y(i) = x(i)+c(i)+u;
      END DO
      RETURN
      END


      SUBROUTINE floop16_F77Overhead(N, y, x, a, b, c, u)
      INTEGER i, N
      REAL*4 y(N), x(N), a(N), b(N), c(N), u
      RETURN
      END
