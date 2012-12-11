
      SUBROUTINE loop14_F77(N, y, x, a, b)
      INTEGER i, N
      REAL*8 y(N), x(N), a(N), b(N)

      DO i=1,N
          x(i) = a(i)+b(i); y(i) = a(i)-b(i);
      END DO
      RETURN
      END


      SUBROUTINE loop14_F77Overhead(N, y, x, a, b)
      INTEGER i, N
      REAL*8 y(N), x(N), a(N), b(N)
      RETURN
      END
