      subroutine fidaxpy(a, x, y, n, iters)

      integer i
      double precision x(n), y(n), a
      double precision b
      b = - a

      do 2000 j = 1, iters
      do 1000 i = 1, n
         y(i) = y(i) + a * x(i)
 1000 continue

      do 3000 i = 1, n
         y(i) = y(i) + b * x(i)
 3000 continue
 2000 continue

      return
      end

      subroutine fidaxpyo(a, x, y, n, iters)

      double precision x(n), y(n), a

      return
      end

