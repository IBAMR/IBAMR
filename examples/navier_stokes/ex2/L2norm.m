function n = L2norm(p)
   N1 = size(p,1);
   N2 = size(p,2);
   if (N1 ~= N2)
      error('invalid grid spacing')
   end
   h = 1/N1;
   n = sqrt(sum(sum(p.^2))*h^2);
