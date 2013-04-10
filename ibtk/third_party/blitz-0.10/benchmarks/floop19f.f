
      SUBROUTINE floop19_F77(N, y, x, a, b, u, v)
      INTEGER i, N
      REAL*4 y(N), x(N), a(N), b(N), u, v

      DO i=1,N
          x(i) = u*a(i); y(i) = v*b(i);
      END DO
      RETURN
      END


      SUBROUTINE floop19_F77Overhead(N, y, x, a, b, u, v)
      INTEGER i, N
      REAL*4 y(N), x(N), a(N), b(N), u, v
      RETURN
      END
