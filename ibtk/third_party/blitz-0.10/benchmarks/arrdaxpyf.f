       subroutine arrdaxpyf(A, B, N, c1)
       integer N, iters
       double precision A(N,N,N), B(N,N,N)
       double precision c1
       double precision c2
       c2 = - c1

         do k=2,N-1
           do j=2,N-1
             do i=2,N-1
               A(i,j,k) = A(i,j,k) + c1 * B(i,j,k);
             enddo
           enddo
         enddo

         do k=2,N-1
           do j=2,N-1
             do i=2,N-1
               A(i,j,k) = A(i,j,k) + c2 * B(i,j,k);
             enddo
           enddo
         enddo

       return
       end

