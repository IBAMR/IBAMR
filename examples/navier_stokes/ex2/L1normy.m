function n = L1normy(v)
   N1 = size(v,1);
   N2 = size(v,2);
   h = 1/N1;
   n = 0.5*sum(abs(v(:,1)))*h^2 + sum(sum(abs(v(:,2:N2-1))))*h^2 + 0.5*sum(abs(v(:,N2)))*h^2;
