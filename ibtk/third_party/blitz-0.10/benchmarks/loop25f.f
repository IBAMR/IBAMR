
      SUBROUTINE loop25_F77(N, x, a, b, c, y, u, v, w)
      INTEGER i, N
      REAL*8 x(N), a(N), b(N), c(N), y(N), u, v, w

      DO i=1,N
          x(i) = u*b(i); y(i) = v*b(i) + w*a(i) + u*c(i);
      END DO
      RETURN
      END


      SUBROUTINE loop25_F77Overhead(N, x, a, b, c, y, u, v, w)
      INTEGER i, N
      REAL*8 x(N), a(N), b(N), c(N), y(N), u, v, w
      RETURN
      END
