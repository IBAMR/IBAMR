! program to time the various routines
subroutine f_main

  use ddmodule
  implicit none
  integer*4 old_cw
  double precision t
  double precision second
  double precision time_thresh
  parameter (time_thresh = 0.5d0)
  type (dd_real) a, b, c, d
  integer n, i, k

  call f_fpu_fix_start (old_cw)

  write (6, *) 'Timing addition / subtraction ...'
  n = 512
  do k = 1, 25
    n = n * 2
    a = ddpi()
    b = sqrt(a)
    c = sqrt(b)
    d = sqrt(c)
    t = second()
    do i = 1, n
      a = b + c
      b = a - d
      a = b + c
      b = a - d
    enddo
    t = second() - t
    if (t .ge. time_thresh) exit
  enddo
  n = n * 4
  write (6, *) n, ' operations in ', t, ' seconds.'
  write (6, *) t/n*1.0d6, ' usec'
  call ddwrite(6, a)

  write (6, *) 'Timing multiplication ...'
  n = 512
  do k = 1, 25
    n = n * 2
    a = 1.0d0 + ddpi() * 1.0d-7
    b = a + 1.0d-7
    c = b + 1.0d-8
    d = c + 1.0d-9
    t = second()
    do i = 1, n
      a = b * c
      b = a * d
      a = b * c
      b = a * d
    enddo
    t = second() - t
    if (t .ge. time_thresh) exit
  enddo
  n = n * 4
  write (6, *) n, ' operations in ', t, ' seconds.'
  write (6, *) t/n*1.0d6, ' usec'
  call ddwrite(6, a)

  write (6, *) 'Timing division ...'
  n = 512
  do k = 1, 25
    n = n * 2
    a = 1.0d0 + ddpi()
    b = 2.0d0 + ddpi()
    c = 1.0d0 + 1.0d-8
    d = 1.0d0 + 1.0d-9
    t = second()
    do i = 1, n
      a = b / c
      b = a / d
      a = b / c
      b = a / d
    enddo
    t = second() - t
    if (t .ge. time_thresh) exit
  enddo
  n = n * 4
  write (6, *) n, ' operations in ', t, ' seconds.'
  write (6, *) t/n*1.0d6, ' usec'
  call ddwrite(6, a)

  write (6, *) 'Timing square root ...'
  n = 512
  do k = 1, 25
    n = n * 2
    a = 0.0d0
    b = 2.0d0 + ddpi()
    t = second()
    do i = 1, n
      a = sqrt(a + b)
    enddo
    t = second() - t
    if (t .ge. time_thresh) exit
  enddo
  write (6, *) n, ' operations in ', t, ' seconds.'
  write (6, *) t/n*1.0d6, ' usec'
  call ddwrite(6, a)

  write (6, *) 'Timing sin ...'
  n = 512
  do k = 1, 25
    n = n * 2
    a = 0.0d0
    c = 1.7d0 * ddreal(1.0d0) / dble(n)
    d = 2.45d0 * ddpi() / dble(n + 3)
    t = second()
    do i = 1, n
      a = a + sin(c)
      c = c + d
    enddo
    t = second() - t
    if (t .ge. time_thresh) exit
  enddo
  write (6, *) n, ' operations in ', t, ' seconds.'
  write (6, *) t/n*1.0d6, ' usec'
  call ddwrite(6, a)

  write (6, *) 'Timing log ...'
  n = 512
  do k = 1, 25
    n = n * 2
    a = 0.0d0
    c = exp(ddreal(-50.1d0));
    d = exp(ddreal(100.2d0) / dble(n))
    t = second()
    do i = 1, n
      a = a + log(c)
      c = c * d
    enddo
    t = second() - t
    if (t .ge. time_thresh) exit
  enddo
  write (6, *) n, ' operations in ', t, ' seconds.'
  write (6, *) t/n*1.0d6, ' usec'
  call ddwrite(6, a)

end

