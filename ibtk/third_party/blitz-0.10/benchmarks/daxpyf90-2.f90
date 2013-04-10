! Fortran 90 DAXPY using arrays
SUBROUTINE f90daxpy(a, x, y, n, iters)
  IMPLICIT NONE
  INTEGER, INTENT( IN ) :: n, iters
  DOUBLE PRECISION, DIMENSION (n) :: x, y
  DOUBLE PRECISION, INTENT( IN ) :: a
  DOUBLE PRECISION :: b

  b = - a

  DO i=1,iters
    y = y + a * y
    y = y + b * y
  END DO
END SUBROUTINE
