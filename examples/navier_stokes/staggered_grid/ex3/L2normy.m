function n = L2normy(v)
   N1 = size(v,1);
   N2 = size(v,2);
   h = 1/N1;
   n = sqrt(0.5*sum(v(:,1).^2)*h^2 + sum(sum(v(:,2:N2-1).^2))*h^2 + 0.5*sum(v(:,N2).^2)*h^2);
