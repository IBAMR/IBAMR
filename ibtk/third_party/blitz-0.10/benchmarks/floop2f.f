
      SUBROUTINE floop2_F77(N, x, y, u)
      INTEGER i, N
      REAL*4 x(N), y(N), u

      DO i=1,N
          x(i) = y(i)/u;
      END DO
      RETURN
      END


      SUBROUTINE floop2_F77Overhead(N, x, y, u)
      INTEGER i, N
      REAL*4 x(N), y(N), u
      RETURN
      END
