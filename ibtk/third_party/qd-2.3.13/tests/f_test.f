! program fortran_test
subroutine f_main
! A simple test of the fortran wrappers

  use qdmodule
  implicit none
  integer*4 old_cw
  integer i
  type (qd_real) x, y, z

  call f_fpu_fix_start (old_cw)

  ! Test for read/write
  z = "3.14159265358979323846264338327950288419716939937510582097494459230"
  call write_scalar(6, z)

  ! Test for atan/write
  do i=1,3
    x = qdreal(dble(i))
    call write_scalar(6, x)
    y = atan(x)
    call write_scalar(6, y)
  end do

  call write_scalar(6, nan(x))
  call write_scalar(6, qdcomplex(x, y))

  call f_fpu_fix_end (old_cw)
end

