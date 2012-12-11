
      SUBROUTINE loop19_F77(N, y, x, a, b, u, v)
      INTEGER i, N
      REAL*8 y(N), x(N), a(N), b(N), u, v

      DO i=1,N
          x(i) = u*a(i); y(i) = v*b(i);
      END DO
      RETURN
      END


      SUBROUTINE loop19_F77Overhead(N, y, x, a, b, u, v)
      INTEGER i, N
      REAL*8 y(N), x(N), a(N), b(N), u, v
      RETURN
      END
