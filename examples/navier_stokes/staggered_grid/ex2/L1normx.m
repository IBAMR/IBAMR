function n = L1normx(u)
   N1 = size(u,1);
   N2 = size(u,2);
   h = 1/N2;
   n = 0.5*sum(abs(u(1,:)))*h^2 + sum(sum(abs(u(2:N1-1,:))))*h^2 + 0.5*sum(abs(u(N1,:)))*h^2;
