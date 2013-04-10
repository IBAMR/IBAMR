
      SUBROUTINE floop18_F77(N, x, a, b, u, v)
      INTEGER i, N
      REAL*4 x(N), a(N), b(N), u, v

      DO i=1,N
          x(i) = (u+a(i))*(v+b(i));
      END DO
      RETURN
      END


      SUBROUTINE floop18_F77Overhead(N, x, a, b, u, v)
      INTEGER i, N
      REAL*4 x(N), a(N), b(N), u, v
      RETURN
      END
