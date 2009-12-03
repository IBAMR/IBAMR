function n = L2normx(u)
   N1 = size(u,1);
   N2 = size(u,2);
   h = 1/N2;
   n = sqrt(0.5*sum(u(1,:).^2)*h^2 + sum(sum(u(2:N1-1,:).^2))*h^2 + 0.5*sum(u(N1,:).^2)*h^2);
