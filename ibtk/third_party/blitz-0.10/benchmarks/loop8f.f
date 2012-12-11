
      SUBROUTINE loop8_F77(N, a, x, u)
      INTEGER i, N
      REAL*8 a(N), x(N), u

      DO i=1,N
          x(i) = u/a(i);
      END DO
      RETURN
      END


      SUBROUTINE loop8_F77Overhead(N, a, x, u)
      INTEGER i, N
      REAL*8 a(N), x(N), u
      RETURN
      END
